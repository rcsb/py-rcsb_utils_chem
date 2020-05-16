##
# File:  OeDepictAlign.py
# Date:  27-Oct-2019  J. Westbrook
#
# Updates:
#
##
"""
Classes to depict aligned 2D chemical diagrams in a variety of media formats.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import logging
import os
import os.path

from openeye import oechem
from openeye import oedepict
from rcsb.utils.chem.OeDepict import OeDepictBase
from rcsb.utils.io.decorators import TimeoutException
from rcsb.utils.io.decorators import timeout

logger = logging.getLogger(__name__)


class OeDepictAlignBase(OeDepictBase):
    """  Base class for aligned 2D molecular renderings containing molecular
         object data and common display preferences.
    """

    def __init__(self):
        super(OeDepictAlignBase, self).__init__()
        #
        #
        self._refId = None
        self._refPath = None
        self._refMol = None
        self._refTitle = None
        self._refImagePath = None
        #
        self._fitId = None
        self._fitPath = None
        self._fitMol = None
        self._fitTitle = None
        self._fitImagePath = None
        #
        self._pairMolList = []
        #
        self._searchType = "default"
        self._minAtomMatchFraction = 0.50
        self._mcss = None
        self._miter = None
        #
        self._ss = None

    def setRefMol(self, oeMol, ccId, title=None, imagePath=None):
        try:
            self._refId = ccId
            self._refMol = oeMol
            #
            if title is not None:
                self._refMol.SetTitle(title)
                self._refTitle = title
            else:
                self._refMol.SetTitle(self._refId)
                self._refTitle = None
            #
            self._refImagePath = imagePath if imagePath is not None else self._refId + ".svg"
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return False

    def setFitMol(self, oeMol, ccId, title=None, imagePath=None):
        try:
            self._fitId = ccId
            self._fitMol = oeMol
            #
            if title is not None:
                self._fitMol.SetTitle(title)
                self._fitTitle = title
            else:
                self._fitMol.SetTitle(self._fitId)
                self._fitTitle = None
            #
            self._fitImagePath = imagePath if imagePath is not None else self._fitId + ".svg"
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def addFitMolList(self, oeMolList, suppressHydrogens=False, imageDirPath=".", imageFilePrefix=None):
        """Set the list of molecules to be compared with reference molecule by MCSS.

           From the input list build the internal pair list (self._pairMolList) of
           tuples  [(refId,refMol,refTitle,refImagePath,fitId,fitMol,fitTitle,fitImagePath),(),...]
        """
        self._pairMolList = []
        try:
            imageFilePrefix = imageFilePrefix + "-" if imageFilePrefix else ""
            for oeMol in oeMolList:
                refId = self._refId
                refImagePath = self._refImagePath
                refMol = oechem.OESuppressHydrogens(oechem.OEGraphMol(self._refMol)) if suppressHydrogens else self._refMol
                #
                fitId = str(oeMol.GetTitle()).upper()
                fitMol = oechem.OESuppressHydrogens(oechem.OEGraphMol(oeMol)) if suppressHydrogens else oeMol
                fitTitle = fitId + "/" + refId
                refTitle = refId + "/" + fitId
                fitImagePath = os.path.join(imageDirPath, imageFilePrefix + fitId + "-with-ref-" + refId + ".svg")
                self._pairMolList.append((refId, refMol, refTitle, refImagePath, fitId, fitMol, fitTitle, fitImagePath))
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return False

    def setPairMolList(self, pairOeMolList, suppressHydrogens=False, imageDirPath=".", imageFilePrefix=None):
        """Set the pairs of molecules to be compared by MCSS.

           From the input list build the internal pair list (self._pairMolList) of
           tuples  [(refId,refMol,refTitle,refImagePath,fitId,fitMol,fitTitle,fitImagePath),(),...]
        """
        self._pairMolList = []
        try:
            imageFilePrefix = imageFilePrefix + "-" if imageFilePrefix else ""
            for refId, refOeMol, fitId, fitOeMol in pairOeMolList:
                refImagePath = None
                refMol = oechem.OESuppressHydrogens(oechem.OEGraphMol(refOeMol)) if suppressHydrogens else refOeMol
                #
                fitMol = oechem.OESuppressHydrogens(oechem.OEGraphMol(fitOeMol)) if suppressHydrogens else fitOeMol
                fitTitle = fitId + "/" + refId
                refTitle = refId + "/" + fitId
                fitImagePath = os.path.join(imageDirPath, imageFilePrefix + "pair-fit-" + fitId + "-with-ref-" + refId + ".svg")
                refImagePath = os.path.join(imageDirPath, imageFilePrefix + "pair-fit-" + refId + "-with-ref-" + fitId + ".svg")
                self._pairMolList.append((refId, refMol, refTitle, refImagePath, fitId, fitMol, fitTitle, fitImagePath))
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return False

    def setSearchType(self, sType="default", minAtomMatchFraction=0.50):
        self._searchType = sType
        self._minAtomMatchFraction = minAtomMatchFraction
        return self._searchType

    def _setupMCSS(self, refmol):
        """ Internal initialization for the MCSS comparison.
        """
        # self._mcss = oechem.OEMCSSearch(oechem.OEMCSType_Approximate)
        self._mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
        #
        if self._searchType in ["default", "graph-strict"]:
            atomexpr = oechem.OEExprOpts_DefaultAtoms
            bondexpr = oechem.OEExprOpts_DefaultBonds
        elif self._searchType in ["relaxed", "graph-relaxed"]:
            # atomexpr = oechem.OEExprOpts_AtomicNumber
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge
            bondexpr = oechem.OEExprOpts_BondOrder
            bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_EqSingleDouble
            # bondexpr = 0
        elif self._searchType == "exact":
            atomexpr = oechem.OEExprOpts_ExactAtoms
            bondexpr = oechem.OEExprOpts_ExactBonds
            # OEAddExplicitHydrogens(refmol)
        else:
            atomexpr = oechem.OEExprOpts_DefaultAtoms
            bondexpr = oechem.OEExprOpts_DefaultBonds
        #
        # atomexpr = oechem.OEExprOpts_AtomicNumber|oechem.OEExprOpts_EqAromatic
        # bondexpr = 0
        #
        # atomexpr = oechem.OEExprOpts_AtomicNumber|oechem.OEExprOpts_Aromaticity
        # bondexpr = oechem.OEExprOpts_BondOrder|oechem.OEExprOpts_EqNotAromatic
        #
        self._mcss.Init(refmol, atomexpr, bondexpr)
        #
        # self._mcss.SetMCSFunc(OEMCSMaxBondsCompleteCycles())
        # self._mcss.SetMCSFunc(OEMCSMaxAtoms())
        #
        # Half of the reference molecule --
        #
        # nAtomsRef=refmol.NumAtoms()
        # self._mcss.SetMinAtoms(nAtomsRef/2)

    def _setupSubStructure(self, refmol):
        """ Internal initialization for a substructure comparison.
        """

        #
        if self._searchType == "default":
            atomexpr = oechem.OEExprOpts_DefaultAtoms
            bondexpr = oechem.OEExprOpts_DefaultBonds
        elif self._searchType == "relaxed":
            atomexpr = oechem.OEExprOpts_AtomicNumber
            bondexpr = 0
        elif self._searchType == "exact":
            atomexpr = oechem.OEExprOpts_ExactAtoms
            bondexpr = oechem.OEExprOpts_ExactBonds
            # OEAddExplicitHydrogens(refmol)
        else:
            atomexpr = oechem.OEExprOpts_DefaultAtoms
            bondexpr = oechem.OEExprOpts_DefaultBonds
        #
        self._ss = oechem.OESubSearch(refmol, atomexpr, bondexpr)
        #

    def _setHighlightStyleRef(self, matchObj, refMol):
        refdisp = oedepict.OE2DMolDisplay(refMol, self._opts)
        #
        optInverseFit = False
        if self._params["highlightStyleRef"] == "ballAndStick":
            hstyle = oedepict.OEHighlightStyle_BallAndStick
        elif self._params["highlightStyleRef"] == "stick":
            hstyle = oedepict.OEHighlightStyle_Stick
        elif self._params["highlightStyleRef"] == "stickInverse":
            hstyle = oedepict.OEHighlightStyle_Stick
            optInverseFit = True
        elif self._params["highlightStyleRef"] == "ballAndStickInverse":
            hstyle = oedepict.OEHighlightStyle_BallAndStick
            optInverseFit = True
        else:
            hstyle = oedepict.OEHighlightStyle_BallAndStick

        if self._params["highLightMatchColorRef"] == "blue":
            myHighLightMatchColor = oechem.OEBlueTint
        elif self._params["highLightMatchColorRef"] == "green":
            myHighLightMatchColor = oechem.OEGreenTint
        elif self._params["highLightMatchColorRef"] == "pink":
            myHighLightMatchColor = oechem.OEPinkTint
        else:
            myHighLightMatchColor = oechem.OEBlueTint

        if self._params["highLightNotMatchColorRef"] == "blue":
            myHighLightNotMatchColor = oechem.OEBlueTint
        elif self._params["highLightNotMatchColorRef"] == "green":
            myHighLightNotMatchColor = oechem.OEGreenTint
        elif self._params["highLightNotMatchColorRef"] == "pink":
            myHighLightNotMatchColor = oechem.OEPinkTint
        else:
            myHighLightNotMatchColor = oechem.OEBlueTint

        #
        matchedatoms = oechem.OEIsAtomMember(matchObj.GetPatternAtoms())
        matchedbonds = oechem.OEIsBondMember(matchObj.GetPatternBonds())
        if optInverseFit:
            oedepict.OEAddHighlighting(refdisp, myHighLightNotMatchColor, hstyle, oechem.OENotAtom(matchedatoms), oechem.OENotBond(matchedbonds))
        else:
            oedepict.OEAddHighlighting(refdisp, myHighLightMatchColor, hstyle, matchedatoms, matchedbonds)

        return refdisp

    def _setHighlightStyleFit(self, matchObj, fitMol):
        fitdisp = oedepict.OE2DMolDisplay(fitMol, self._opts)

        optInverseFit = False
        if self._params["highlightStyleFit"] == "ballAndStick":
            hstyle = oedepict.OEHighlightStyle_BallAndStick
        elif self._params["highlightStyleFit"] == "stick":
            hstyle = oedepict.OEHighlightStyle_Stick
        elif self._params["highlightStyleFit"] == "stickInverse":
            hstyle = oedepict.OEHighlightStyle_Stick
            optInverseFit = True
        elif self._params["highlightStyleFit"] == "ballAndStickInverse":
            hstyle = oedepict.OEHighlightStyle_BallAndStick
            optInverseFit = True
        else:
            hstyle = oedepict.OEHighlightStyle_BallAndStick

        if self._params["highLightMatchColorFit"] == "blue":
            myHighLightMatchColor = oechem.OEBlueTint
        elif self._params["highLightMatchColorFit"] == "green":
            myHighLightMatchColor = oechem.OEGreenTint
        elif self._params["highLightMatchColorFit"] == "pink":
            myHighLightMatchColor = oechem.OEPinkTint
        else:
            myHighLightMatchColor = oechem.OEBlueTint

        if self._params["highLightNotMatchColorFit"] == "blue":
            myHighLightNotMatchColor = oechem.OEBlueTint
        elif self._params["highLightNotMatchColorFit"] == "green":
            myHighLightNotMatchColor = oechem.OEGreenTint
        elif self._params["highLightNotMatchColorFit"] == "pink":
            myHighLightNotMatchColor = oechem.OEPinkTint
        else:
            myHighLightNotMatchColor = oechem.OEBlueTint

        #
        # matchedatoms = oechem.OEIsAtomMember(matchObj.GetPatternAtoms())
        # matchedbonds = oechem.OEIsBondMember(matchObj.GetPatternBonds())
        matchedatoms = oechem.OEIsAtomMember(matchObj.GetTargetAtoms())
        matchedbonds = oechem.OEIsBondMember(matchObj.GetTargetBonds())

        if optInverseFit:
            oedepict.OEAddHighlighting(fitdisp, myHighLightNotMatchColor, hstyle, oechem.OENotAtom(matchedatoms), oechem.OENotBond(matchedbonds))
        else:
            oedepict.OEAddHighlighting(fitdisp, myHighLightMatchColor, hstyle, matchedatoms, matchedbonds)

        return fitdisp


class OeDepictMCSAlignMultiPage(OeDepictAlignBase):
    """ Create 2D depictions of MCSS alignments. Inputs can be in the the form of pairs,
        lists, and pair lists of molecule object instances.

        Output images are rendered in a grid layout that can span multiple pages.
    """

    def __init__(self):
        super(OeDepictMCSAlignMultiPage, self).__init__()
        self.__grid = None
        self.__gridRows = None
        self.__gridCols = None
        self.__multi = None
        self.__image = None
        self.__citer = None
        #

    def alignPairListMulti(self, imagePath="multi.pdf"):
        aM = []
        self._params["gridCols"] = 2
        try:
            aM = self.__alignListMultiWorker(imagePath=imagePath, layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return aM

    def alignOneWithListMulti(self, imagePath="multi.pdf"):
        aM = []
        try:
            aM = self.__alignListMultiWorker(imagePath=imagePath, layout="list")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return aM

    def __setupImageMulti(self):
        """ Internal method to configure a multipage image.
        """
        #
        self.__gridRows = self._params["gridRows"]
        self.__gridCols = self._params["gridCols"]
        #
        if self._params["pageOrientation"] == "landscape":
            self.__multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Landscape, oedepict.OEPageSize_US_Letter)
        else:
            self.__multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Portrait, oedepict.OEPageSize_US_Letter)

        self.__newPage()

    def __newPage(self):
        """ Internal method to advance to a new page in a multipage configuration.
        """
        rows = self.__gridRows
        cols = self.__gridCols
        self.__image = self.__multi.NewPage()
        self.__grid = oedepict.OEImageGrid(self.__image, rows, cols)
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        logger.debug("Num columns %d", self.__grid.NumCols())
        logger.debug("Num rows    %d", self.__grid.NumRows())
        self._opts = oedepict.OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), oedepict.OEScale_AutoScale)
        self._assignDisplayOptions()
        self.__citer = self.__grid.GetCells()

    @timeout(500)
    def __alignListMultiWorker(self, imagePath="multi.pdf", layout="pairs"):
        """ Working method comparing a reference molecule with a list of fit molecules.

            pairMolList = (refId,refMol,refTitle,refImagePath,fitId,fitMol,fitTitle,fitImagePath)

            Map of corresponding atoms is returned.

            Image Output is in multipage layout.
        """
        #
        self.__setupImageMulti()
        #
        atomMap = []
        firstOne = True
        iCount = 0
        for (refId, refMol, _, _, fitId, fitMol, fitTitle, _) in self._pairMolList:
            iCount += 1
            #
            oedepict.OEPrepareDepiction(refMol)
            self._setupMCSS(refMol)
            #
            #
            fitMol.SetTitle(fitTitle)
            #
            oedepict.OEPrepareDepiction(fitMol)
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            mcssMinAtoms = int(minAtoms * self._minAtomMatchFraction)
            self._mcss.SetMinAtoms(mcssMinAtoms)

            # scaling
            refscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            fitscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)
            logger.debug("mcss match completed for refId %s fitId %s", refId, fitId)
            if self._miter.IsValid():
                match = self._miter.Target()
                oedepict.OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and layout in ["list"]) or layout in ["pairs"]:
                    firstOne = False
                    if layout in ["pairs"]:
                        refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    else:
                        refdisp = oedepict.OE2DMolDisplay(refMol, self._opts)
                    if not self.__citer.IsValid():
                        self.__newPage()
                    cell = self.__citer.Target()
                    self.__citer.Next()
                    oedepict.OERenderMolecule(cell, refdisp)
                    if self._params["cellBorders"]:
                        oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                if not self.__citer.IsValid():
                    self.__newPage()
                cell = self.__citer.Target()
                self.__citer.Next()
                oedepict.OERenderMolecule(cell, fitdisp)

                if self._params["cellBorders"]:
                    oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

                logger.debug("mcss match completed for refId %s fitId %s total map length %d", refId, fitId, len(atomMap))

        logger.debug("writing image %s", imagePath)
        oedepict.OEWriteMultiPageImage(imagePath, self.__multi)
        logger.debug("completed with map lenth %d", len(atomMap))
        #
        return atomMap


class OeDepictMCSAlignPage(OeDepictAlignBase):
    """ Create 2D depictions of MCSS alignments. Inputs can be in the the form of pairs,
        lists, and pair lists of molecule object instances.

        Output images are rendered to a single page image with grid layout.
    """

    def __init__(self):
        super(OeDepictMCSAlignPage, self).__init__()
        self.__grid = None
        # self.__gridRows = None
        # self.__gridCols = None
        # self.__multi = None
        self.__image = None
        self.__citer = None
        #

    def __setupImage(self):
        """ Internal method to configure a single pair alignment image.
        """
        #
        self.__image = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        rows = self._params["gridRows"]
        cols = self._params["gridCols"]
        self.__grid = oedepict.OEImageGrid(self.__image, rows, cols)
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        logger.debug("Num columns %d", self.__grid.NumCols())
        logger.debug("Num rows    %d", self.__grid.NumRows())
        self._opts = oedepict.OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), oedepict.OEScale_AutoScale)
        self._assignDisplayOptions()
        self.__citer = self.__grid.GetCells()

    def alignPair(self, imagePath="single-pair.png"):
        """  Depict a single aligned ref/fit molecule pair or the first ref/fit molecule pair on the
             current _pairMolList.  Display options set for a single grid row with two columns.
        """
        self._params["gridCols"] = 2
        self._params["gridRows"] = 1
        self.__setupImage()
        #
        aM = []
        self._pairMolList = []
        self._pairMolList.append((self._refId, self._refMol, self._refTitle, None, self._fitId, self._fitMol, self._fitTitle, None))
        try:
            aM = self.__alignListWorker(imagePath=imagePath, layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignPairList(self, imagePath="single-pair-list.png"):
        self._params["gridCols"] = 2
        try:
            aM = self.__alignListWorker(imagePath=imagePath, layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignOneWithList(self, imagePath="single-list.png"):
        try:
            aM = self.__alignListWorker(imagePath=imagePath, layout="list")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    # @timeout(15)
    def __alignListWorker(self, imagePath="single.pdf", layout="pairs"):
        """ Working method comparing a reference molecule with a list of fit molecules.

            pairMolList = (refId,refMol,refTitle,refImgPath,fitId,fitMol,fitTitle,fitImgPath)

            Map of corresponding atoms is returned.

            Output image is a single-page with grid layout.
        """
        #
        self.__setupImage()
        #
        atomMap = []

        firstOne = True

        for (refId, refMol, _, _, fitId, fitMol, fitTitle, _) in self._pairMolList:
            #
            oedepict.OEPrepareDepiction(refMol)
            self._setupMCSS(refMol)
            #
            title = fitTitle if fitTitle else fitId
            fitMol.SetTitle(title)

            oedepict.OEPrepareDepiction(fitMol)
            #
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            # scaling
            refscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            fitscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)

            if self._miter.IsValid():
                match = self._miter.Target()
                oedepict.OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and layout in ["list"]) or layout in ["pairs"]:
                    firstOne = False
                    if layout in ["pairs"]:
                        refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    else:
                        refdisp = oedepict.OE2DMolDisplay(refMol, self._opts)
                    if not self.__citer.IsValid():
                        break
                    cell = self.__citer.Target()
                    self.__citer.Next()
                    oedepict.OERenderMolecule(cell, refdisp)
                    if self._params["cellBorders"]:
                        oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                if not self.__citer.IsValid():
                    break
                cell = self.__citer.Target()
                self.__citer.Next()
                oedepict.OERenderMolecule(cell, fitdisp)

                if self._params["cellBorders"]:
                    oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        oedepict.OEWriteImage(imagePath, self.__image)
        return atomMap


class OeDepictMCSAlign(OeDepictAlignBase):
    """ Create 2D depictions of MCSS alignments. Inputs can be in the the form of pairs,
        lists, and pair lists of molecule object instances.

        Outputs are separate image files with a single diagram per file.
    """

    def __init__(self):
        super(OeDepictMCSAlign, self).__init__()
        #
        self.__imageRef = None
        self.__imageFit = None
        #

    def __setupImage(self):
        """ Internal method to configure a single page image.
        """
        #
        self.__imageRef = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self.__imageFit = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self._opts = oedepict.OE2DMolDisplayOptions(self.__imageRef.GetWidth(), self.__imageRef.GetHeight(), oedepict.OEScale_AutoScale)
        self._assignDisplayOptions()

    def alignPair(self):
        """  Depict a single aligned ref/fit molecule pair or the first ref/fit molecule pair on the
             current _pairMolList.  Display options set for a single grid row with two columns.
        """
        self._pairMolList = []
        self._pairMolList.append((self._refId, self._refMol, self._refTitle, self._refImagePath, self._fitId, self._fitMol, self._fitTitle, self._fitImagePath))
        try:
            aM = self.__alignListWorker(layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignPairList(self):
        try:
            aM = self.__alignListWorker(layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignOneWithList(self):
        try:
            aM = self.__alignListWorker(layout="list")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    @timeout(15)
    def __alignListWorker(self, layout="pairs"):
        """ Working method comparing a reference molecule with a list of fit molecules.

            pairMolList = (refId,refMol,refTitle,refImagePath, fitId,fitMol,fitTitle, fitImagePath)

            Map of corresponding atoms is returned.

            Writes separate output images for the reference and fit molecules for each comparison pair.
        """
        #
        atomMap = []
        firstOne = True

        for (refId, refMol, _, refImagePath, fitId, fitMol, fitTitle, fitImagePath) in self._pairMolList:
            #
            self.__setupImage()
            #
            oedepict.OEPrepareDepiction(refMol)
            #
            self._setupMCSS(refMol)

            fitMol.SetTitle(fitTitle)
            oedepict.OEPrepareDepiction(fitMol)
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            # scaling
            refscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            fitscale = oedepict.OEGetMoleculeScale(fitMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)
            if self._miter.IsValid():
                match = self._miter.Target()
                oedepict.OEPrepareAlignedDepiction(fitMol, self._mcss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and layout in ["list"]) or layout in ["pairs"]:
                    firstOne = False
                    if layout in ["pairs"]:
                        refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._mcss.GetPattern())
                    else:
                        refdisp = oedepict.OE2DMolDisplay(refMol, self._opts)
                    oedepict.OERenderMolecule(self.__imageRef, refdisp)
                    oedepict.OEWriteImage(refImagePath, self.__imageRef)

                # Depict fit molecule with MCS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                oedepict.OERenderMolecule(self.__imageFit, fitdisp)
                oedepict.OEWriteImage(fitImagePath, self.__imageFit)

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        return atomMap


class OeMCSAlignUtil(OeDepictAlignBase):
    """ Perform MCSS alignments.  Inputs can be in the the form of pairs,
        lists, and pair lists of molecule object instances.
    """

    def __init__(self, maxMatches=2048):
        super(OeMCSAlignUtil, self).__init__()
        #
        self.__maxMatches = maxMatches

    def doAlign(self):
        """ Test the MCSS comparison between the current reference and fit molecules -

            Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        minAtoms = min(nAtomsRef, nAtomsFit)
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = True
        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))

        return atomMap

    def __getAtomIndex(self, oeMol):
        atD = {}
        try:
            for ii, atom in enumerate(oeMol.GetAtoms()):
                atD[atom.GetName().strip()] = ii
        except Exception:
            pass
        return atD

    def __getChargeIndex(self, oeMol):
        atD = {}
        try:
            for atom in oeMol.GetAtoms():
                atD[atom.GetName().strip()] = atom.GetFormalCharge()
        except Exception:
            pass
        return atD

    def __getNeighbors(self, oeMol, atNameList):
        nL = []
        try:
            for atom in oeMol.GetAtoms():
                aN = atom.GetName().strip()
                if aN in atNameList:
                    tL = []
                    for nbr in atom.GetAtoms():
                        tL.append(nbr.GetName())
                        nL.append((aN, tL))
        except Exception:
            pass
        return nL

    def doTestAlign(self):
        """ Test the MCSS comparison between the current reference and fit molecules -

            Return list of corresponding atoms on success or an empty list otherwise.
        """
        self._minAtomMatchFraction = 0.9
        unique = True
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        #
        #
        minAtoms = min(nAtomsRef, nAtomsFit)
        self._mcss.SetMCSFunc(oechem.OEMCSMaxAtoms())
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = False

        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))

        logger.debug("nAtomsRef %d nAtomsFit %d match length %d", nAtomsRef, nAtomsFit, len(atomMap))

        return atomMap

    def doAlignWithAnal(self):
        """ Test the MCSS comparison between the current reference and fit molecules -

            Return list of corresponding atoms on success or an empty list otherwise.
        """
        self._minAtomMatchFraction = 0.9
        unique = True
        atomMap = []
        self._setupMCSS(self._refMol)
        #
        nAtomsRef = self._refMol.NumAtoms()
        nAtomsFit = self._fitMol.NumAtoms()
        #
        atomRefD = self.__getAtomIndex(self._refMol)
        atomFitD = self.__getAtomIndex(self._fitMol)
        #
        chgRefD = self.__getChargeIndex(self._refMol)
        chgFitD = self.__getChargeIndex(self._fitMol)
        #
        minAtoms = min(nAtomsRef, nAtomsFit)
        #
        self._mcss.SetMaxMatches(self.__maxMatches)
        nlim = oechem.OEGetMCSExhaustiveSearchTruncationLimit()
        logger.debug("search limit %d max matches %d", nlim, self.__maxMatches)
        self._mcss.SetMCSFunc(oechem.OEMCSMaxAtoms())
        self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))
        unique = False

        atomRefMatchD = {}
        atomFitMatchD = {}
        self._miter = self._mcss.Match(self._fitMol, unique)
        if self._miter.IsValid():
            match = self._miter.Target()
            for mAt in match.GetAtoms():
                atomMap.append((self._refId, mAt.pattern.GetName(), self._fitId, mAt.target.GetName()))
                atomRefMatchD[mAt.pattern.GetName()] = mAt.target.GetName().strip()
                atomFitMatchD[mAt.target.GetName()] = mAt.pattern.GetName().strip()
        logger.debug("nAtomsRef %d nAtomsFit %d match length %d", nAtomsRef, nAtomsFit, len(atomMap))

        #
        # Get unmapped reference and fit atoms -
        #
        uRefAtomList = []
        uRefNList = []
        for atN in atomRefD:
            if atN not in atomRefMatchD:
                uRefAtomList.append(atN)
        #
        if len(uRefAtomList) > 0:
            uRefNList = self.__getNeighbors(self._refMol, uRefAtomList)
        #
        uFitAtomList = []
        uFitNList = []
        for atN in atomFitD:
            if atN not in atomFitMatchD:
                uFitAtomList.append(atN)
        #
        #
        if len(uFitAtomList) > 0:
            uFitNList = self.__getNeighbors(self._fitMol, uFitAtomList)
        #
        chgDifRefD = {}
        for atN0 in chgRefD:
            if atN0 in atomRefMatchD:
                atN = atomRefMatchD[atN0]
                # if ((atN not in chgFitD) or (chgRefD[atN0] != chgFitD[atN])):
                if chgRefD[atN0] != chgFitD[atN]:
                    chgDifRefD[atN] = chgRefD[atN0]
        #
        chgDifFitD = {}
        for atN0 in chgFitD:
            if atN0 in atomFitMatchD:
                atN = atomFitMatchD[atN0]
                # if ((atN not in chgRefD) or (chgFitD[atN0] != chgRefD[atN])):
                if chgFitD[atN0] != chgRefD[atN]:
                    chgDifFitD[atN] = chgFitD[atN0]

        return atomMap, uRefAtomList, uRefNList, chgDifRefD, uFitAtomList, uFitNList, chgDifFitD

    def doAlignList(self):
        """ Test MCSS comparison between the current reference molecule with a list of fit molecules.

            pairMolList = (refId,refMol,refTitle,fitId,fitMol,fitTitle)

            Map of corresponding atoms is returned.

        """
        atomMap = []
        for (refId, refMol, _, _, fitId, fitMol, _, _) in self._pairMolList:
            #
            self._setupMCSS(refMol)
            #
            #
            nAtomsRef = refMol.NumAtoms()
            nAtomsFit = fitMol.NumAtoms()
            minAtoms = min(nAtomsRef, nAtomsFit)
            self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            unique = True
            self._miter = self._mcss.Match(fitMol, unique)

            if self._miter.IsValid():
                match = self._miter.Target()
                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        return atomMap


class OeDepictSubStructureAlign(OeDepictAlignBase):
    """ Create 2D depictions of substructure alignments. Inputs can be in the the form of pairs,
        lists, and pair lists of molecule object instances.

        Outputs are separate image files with a single diagram per file.
    """

    def __init__(self):
        super(OeDepictSubStructureAlign, self).__init__()
        #
        self.__imageRef = None
        self.__imageFit = None
        #

    def __setupImage(self):
        """ Internal method to configure a single page image.
        """
        #
        self.__imageRef = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self.__imageFit = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self._opts = oedepict.OE2DMolDisplayOptions(self.__imageRef.GetWidth(), self.__imageRef.GetHeight(), oedepict.OEScale_AutoScale)
        self._assignDisplayOptions()

    def alignPair(self):
        """  Depict a single aligned ref/fit molecule pair or the first ref/fit molecule pair on the
             current _pairMolList.  Display options set for a single grid row with two columns.
        """
        self._pairMolList = []
        self._pairMolList.append((self._refId, self._refMol, self._refTitle, self._refImagePath, self._fitId, self._fitMol, self._fitTitle, self._fitImagePath))
        try:
            aM = self.__alignListWorker(layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignPairList(self):
        try:
            aM = self.__alignListWorker(layout="pairs")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    def alignOneWithList(self):
        try:
            aM = self.__alignListWorker(layout="list")
        except TimeoutException:
            logger.info("Timeout exception")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aM

    @timeout(15)
    def __alignListWorker(self, layout="pairs"):
        """ Working method comparing a reference molecule with a list of fit molecules.

            pairMolList = (refId,refMol,refTitle,refImagePath, fitId,fitMol,fitTitle, fitImagePath)

            Map of corresponding atoms is returned.

            Writes separate output images for the reference and fit molecules for each comparison pair.
        """
        #
        atomMap = []
        firstOne = True

        for (refId, refMol, _, refImagePath, fitId, fitMol, fitTitle, fitImagePath) in self._pairMolList:
            #
            self.__setupImage()
            #
            oedepict.OEPrepareDepiction(refMol)
            #
            self._setupSubStructure(refMol)
            fitMol.SetTitle(fitTitle)
            oedepict.OEPrepareDepiction(fitMol)
            #
            # nAtomsRef = refMol.NumAtoms()
            # nAtomsFit = fitMol.NumAtoms()
            # minAtoms = min(nAtomsRef, nAtomsFit)
            # self._mcss.SetMinAtoms(int(minAtoms * self._minAtomMatchFraction))

            # scaling
            refscale = oedepict.OEGetMoleculeScale(refMol, self._opts)
            fitscale = oedepict.OEGetMoleculeScale(fitMol, self._opts)
            self._opts.SetScale(min(refscale, fitscale))
            #
            unique = True
            self._miter = self._ss.Match(fitMol, unique)
            if self._miter.IsValid():
                match = self._miter.Target()
                oedepict.OEPrepareAlignedDepiction(fitMol, self._ss.GetPattern(), match)

                # Depict reference molecule with MCS highlighting
                if (firstOne and layout in ["list"]) or layout in ["pairs"]:
                    firstOne = False
                    if layout in ["pairs"]:
                        refdisp = self._setHighlightStyleRef(matchObj=match, refMol=self._ss.GetPattern())
                    else:
                        refdisp = oedepict.OE2DMolDisplay(refMol, self._opts)
                    logger.info("writing %r", refImagePath)
                    oedepict.OERenderMolecule(self.__imageRef, refdisp)
                    oedepict.OEWriteImage(refImagePath, self.__imageRef)

                # Depict fit molecule with SS highlighting
                fitdisp = self._setHighlightStyleFit(matchObj=match, fitMol=fitMol)
                oedepict.OERenderMolecule(self.__imageFit, fitdisp)
                oedepict.OEWriteImage(fitImagePath, self.__imageFit)

                for mAt in match.GetAtoms():
                    atomMap.append((refId, mAt.pattern.GetName(), fitId, mAt.target.GetName()))

        return atomMap
