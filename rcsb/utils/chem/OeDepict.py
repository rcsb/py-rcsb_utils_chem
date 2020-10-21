##
# File:  OeDepict.py
# Date:  25-Oct-2019  J. Westbrook  adapted from earlier code module.
#
# Updates:
#  2-April-2020 jdw enlarge the base image size to avoid dropping labels.
#                   Tweak stereo style to include W/H but avoid undefined stereo marking
#
##
"""
Utilities to depict individual OE molecules.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os

from openeye import oechem
from openeye import oedepict

logger = logging.getLogger(__name__)


class LabelAtoms(oedepict.OEDisplayAtomPropBase):
    def __init__(self):
        oedepict.OEDisplayAtomPropBase.__init__(self)

    def __call__(self, atom):
        return atom.GetName()

    def CreateCopy(self):
        # __disown__ is required to allow C++ to take
        # ownsership of this object and its memory
        copy = LabelAtoms()
        return copy.__disown__()


class OeDepictBase(object):
    """Base class for 2D depictions in single and multi-page format."""

    def __init__(self):
        super(OeDepictBase, self).__init__()
        self._molTitleList = []
        self._opts = None
        #
        # internal dictionary of display parameters -
        #
        self._params = {
            "imageSizeX": 2500,
            "imageSizeY": 2500,
            "cellBorders": True,
            "suppressHydrogens": False,
            "labelAtomName": False,
            "labelAtomCIPStereo": False,
            "labelAtomIndex": False,
            "labelBondCIPStereo": False,
            "labelBondIndex": False,
            "bondDisplayWidth": None,
            "gridRows": 2,
            "gridCols": 2,
            "cellGap": 5,
            "cellMargin": 10,
            "pageOrientation": "landscape",
            "highlightStyleRef": "none",
            "highlightStyleFit": "ballAndStick",
            "highLightMatchColorRef": "blue",
            "highLightMatchColorFit": "blue",
            "highLightNotMatchColorRef": "pink",
            "highLightNotMatchColorFit": "pink",
        }

    def setDisplayOptions(self, **kwargs):
        self._params.update(kwargs)

    def setGridOptions(self, rows=1, cols=1, cellBorders=True):
        self._params["gridRows"] = rows
        self._params["gridCols"] = cols
        self._params["cellBorders"] = cellBorders

    def setMolTitleList(self, oeMolTitleList):
        """Set the list of OE Mols to be depicted as a list of tuples containing
        [(ccId,oeMol,titleString),(ccId,oeMol,titleString),...]
        """
        self._molTitleList = oeMolTitleList

    def _assignDisplayOptions(self):
        if self._params["labelAtomCIPStereo"]:
            #
            self._opts.SetAtomStereoStyle(oedepict.OEAtomStereoStyle_Display_All)
            # self._opts.SetAtomStereoStyle(oedepict.OEAtomStereoStyle_Display_CIPAtomStereo)

        if self._params["labelBondCIPStereo"]:
            # will include bowties for undefined stereo
            # self._opts.SetBondStereoStyle(oedepict.OEBondStereoStyle_Display_All)
            self._opts.SetBondStereoStyle(oedepict.OEBondStereoStyle_Display_CIPBondStereo)
            self._opts.SetAtomPropLabelFontScale(0.650)

        if self._params["labelAtomIndex"]:
            self._opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
            self._opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))
            self._opts.SetAtomPropLabelFontScale(0.650)

        if self._params["labelBondIndex"]:
            self._opts.SetBondPropertyFunctor(oedepict.OEDisplayBondIdx())
            self._opts.SetBondPropLabelFont(oedepict.OEFont(oechem.OEDarkBlue))

        if self._params["labelAtomName"]:
            atomlabel = LabelAtoms()
            self._opts.SetAtomPropertyFunctor(atomlabel)
            self._opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))
            self._opts.SetAtomPropLabelFontScale(0.650)

        if self._params["bondDisplayWidth"] is not None:
            pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_On, self._params["bondDisplayWidth"])
            self._opts.SetDefaultBondPen(pen)
            # remove for the moment.  not supported on all platforms
            # self._opts.SetBondWidthScaling(False)
        #
        # 5.0 is minimum size -
        self._opts.SetTitleHeight(5.0)
        #


class OeDepictMultiPage(OeDepictBase):

    """Create 2D depictions in multipage format from a list of OE molecules and title strings"""

    def __init__(self, useTitle=True):
        super(OeDepictMultiPage, self).__init__()
        self.__useTitle = useTitle
        self.__multi = None
        self.__image = None
        #

    def __setupImage(self):
        if self._params["pageOrientation"] == "landscape":
            self.__multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Landscape, oedepict.OEPageSize_US_Letter)
        else:
            self.__multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Portrait, oedepict.OEPageSize_US_Letter)
        self.__image = self.__multi.NewPage()
        self._opts = oedepict.OE2DMolDisplayOptions()

    def prepare(self):
        self.__setupImage()
        rows = self._params["gridRows"]
        cols = self._params["gridCols"]
        grid = oedepict.OEImageGrid(self.__image, rows, cols)

        citer = grid.GetCells()

        for ccId, oeMol, title in self._molTitleList:
            logger.debug("Preparing %s %r", ccId, title)
            if not citer.IsValid():
                # go to next page
                self.__image = self.__multi.NewPage()
                grid = oedepict.OEImageGrid(self.__image, rows, cols)
                grid.SetCellGap(self._params["cellGap"])
                grid.SetMargins(self._params["cellMargin"])
                citer = grid.GetCells()

            cell = citer.Target()
            #
            if self._params["suppressHydrogens"]:
                # mol = oeMol.getGraphMolSuppressH()
                #  OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
                mol = oechem.OESuppressHydrogens(oechem.OEGraphMol(oeMol))
            else:
                mol = oeMol

            if self.__useTitle and title:
                mol.SetTitle(title)
                self._opts.SetTitleHeight(5.0)
            else:
                mol.SetTitle("")
            #
            #
            oedepict.OEPrepareDepiction(mol)
            self._opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
            self._assignDisplayOptions()

            disp = oedepict.OE2DMolDisplay(mol, self._opts)
            oedepict.OERenderMolecule(cell, disp)
            oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

            citer.Next()

    def write(self, imagePath):
        try:
            dirPath, _ = os.path.split(imagePath)
            if not os.access(dirPath, os.W_OK):
                os.makedirs(dirPath, mode=0o755)
        except Exception:
            pass
        oedepict.OEWriteMultiPageImage(imagePath, self.__multi)


class OeDepict(OeDepictBase):

    """Create 2D depictions in single-page format from a list of OE molecules & title strings"""

    def __init__(self, useTitle=True):
        super(OeDepict, self).__init__()
        self.__useTitle = useTitle
        self.__image = None
        self.__grid = None

    def __setupImage(self):
        """Internal method to configure a single page image."""
        #
        self.__image = oedepict.OEImage(self._params["imageSizeX"], self._params["imageSizeY"])
        self.__grid = oedepict.OEImageGrid(self.__image, self._params["gridRows"], self._params["gridCols"])
        self.__grid.SetCellGap(self._params["cellGap"])
        self.__grid.SetMargins(self._params["cellMargin"])
        self._opts = oedepict.OE2DMolDisplayOptions(self.__grid.GetCellWidth(), self.__grid.GetCellHeight(), oedepict.OEScale_AutoScale)
        #
        logger.debug("Num columns %d", self.__grid.NumCols())
        logger.debug("Num rows    %d", self.__grid.NumRows())

    def prepare(self):
        """[summary]
         # OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
        oechem.OESuppressHydrogens(self.__oeMol)
        """
        self.__setupImage()
        for idx, cell in enumerate(self.__grid.GetCells()):
            ccId, oeMol, title = self._molTitleList[idx]
            logger.debug("Preparing %s %r", ccId, title)
            #
            if self._params["suppressHydrogens"]:
                # mol = oeMol.getGraphMolSuppressH()
                #  OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
                mol = oechem.OESuppressHydrogens(oechem.OEGraphMol(oeMol))
            else:
                mol = oeMol
            #
            if self.__useTitle and title:
                mol.SetTitle(title)
                self._opts.SetTitleHeight(5.0)
            else:
                mol.SetTitle("")
            #
            #
            oedepict.OEPrepareDepiction(mol)
            self._opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)

            self._assignDisplayOptions()

            disp = oedepict.OE2DMolDisplay(mol, self._opts)
            oedepict.OERenderMolecule(cell, disp)
            if self._params["cellBorders"]:
                oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OEBlackPen))

    def write(self, imagePath):
        try:
            dirPath, _ = os.path.split(imagePath)
            if not os.access(dirPath, os.W_OK):
                os.makedirs(dirPath, mode=0o755)
        except Exception:
            pass
        oedepict.OEWriteImage(imagePath, self.__image)
