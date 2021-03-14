##
# File:  OeAlignUtils.py
# Date:  17-Dec-2020  J. Westbrook
#
# Updates:
#  22-Jan-2021 jdw rename to OeAlignUtils() and add substructure mode
#   2-Feb-2021 jdw suppress OE warnings unless verbose mode is set
#   9-Mar-2021 jdw add setRefObj/setFitObj methods for setting molecule objects directly
##
"""
Utilities for performing substructure and maximum common substructure molecular alignments.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
from collections import namedtuple

from openeye import oechem
from wrapt_timeout_decorator import timeout

from rcsb.utils.chem.OeCommonUtils import OeCommonUtils
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)

ComponentAtomDetails = namedtuple("ComponentAtomDetails", "atIdx atNo atName atType x y z atFormalCharge")
AlignAtomMap = namedtuple("AlignAtomMap", "refId refAtIdx refAtNo refAtName fitId fitAtIdx fitAtNo fitAtName")
AlignAtomUnMapped = namedtuple("AlignAtomUnMapped", "fitId fitAtIdx fitAtNo fitAtType fitAtName fitAtFormalCharge x y z fitNeighbors")


class OeAlignUtils(object):
    """Perform substructure and maximum common substructure molecular alignments.  Targets can be chemical component
    identifiers or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    """

    def __init__(self, workPath=None, verbose=True, timeOut=None):
        #
        self.__verbose = verbose
        #
        self.__refId = None
        self.__refmol = None
        self.__refTitle = None
        #
        self.__fitId = None
        self.__fitmol = None
        self.__fitTitle = None
        #
        self.__pairTupleList = []
        #
        self.__searchType = "relaxed"
        #
        self.__refFD = {}
        self.__fitFD = {}
        #
        self.__mcss = None
        self.__ss = None
        self.__refPath = None
        self.__fitPath = None
        self.__workPath = workPath if workPath else "."
        self.timeOut = timeOut

    def setSearchType(self, sType="relaxed"):
        self.__searchType = sType
        return self.__searchType

    def setRefId(self, ccId, title=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the query molecule for MCSS comparison using the input chemical component ID.
        It is assumed that the definition for this ID can be obtained from the chemical component
        repository.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__refId = ccId
        ccIdU = ccId.upper()
        self.__refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")

        _, self.__refmol, self.__refFD = self.getCCDefFile(self.__refPath, suppressHydrogens=suppressHydrogens)
        #
        if title is not None:
            self.__refmol.SetTitle(title)
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None

        return self.__refmol.NumAtoms() if self.__refmol else 0

    def setRefPath(self, ccPath, title=None, molBuildType="model-xyz", suppressHydrogens=False, fType="CC", importType="2D"):
        """Set the query molecule for MCSS comparison using the input file path.

        The file type is either ['CC'] for a chemical component definition or another file type
        supported by OE toolkit assumed to have a conventional file extension for this type.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__refPath = ccPath
        if fType in ["CC"]:
            (self.__refId, self.__refmol, self.__refFD) = self.getCCDefFile(ccPath, molBuildType=molBuildType, suppressHydrogens=suppressHydrogens)
        else:
            (self.__refId, self.__refmol, self.__refFD) = self.__getMiscFile(ccPath, suppressHydrogens=suppressHydrogens, importType=importType, title=title)

        if self.__verbose:
            logger.debug("Derived ref ID     = %s", self.__refId)
            logger.debug("SMILES (stereo)  = %s", self.__refFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None

        return self.__refmol.NumAtoms() if self.__refmol else 0

    def setRefObj(self, ccObj, title=None, molBuildType="model-xyz", suppressHydrogens=False):
        """Set the query molecule for MCSS comparison using the input chemical component definition object.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        (self.__refId, self.__refmol, self.__refFD) = self.getCCDefObj(ccObj, molBuildType=molBuildType, suppressHydrogens=suppressHydrogens)
        if self.__verbose:
            logger.debug("Derived ref ID     = %s", self.__refId)
            logger.debug("SMILES (stereo)  = %s", self.__refFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None

        return self.__refmol.NumAtoms() if self.__refmol else 0

    def setFitId(self, ccId, title=None, suppressHydrogens=False, cachePath="/data/components/ligand-dict-v3"):
        """Set the ID of the target/library molecule for MCSS comparison."""
        self.__fitId = ccId
        ccIdU = ccId.upper()
        self.__fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
        self.__fitId, self.__fitmol, self.__fitFD = self.getCCDefFile(self.__fitPath, suppressHydrogens=suppressHydrogens)
        if self.__verbose:
            logger.info("Fit ID             = %s", self.__fitId)
            logger.info("SMILES (isomeric)  = %s", self.__fitFD["SMILES_STEREO"])
        if title is not None:
            self.__fitmol.SetTitle(title)
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            self.__fitTitle = None
        return self.__fitmol.NumAtoms() if self.__fitmol else 0

    def setFitPath(self, ccPath, title=None, molBuildType="model-xyz", suppressHydrogens=False, fType="CC", importType="2D", largestPart=False):
        """Set the path to the target/library molecule for MCSS comparison using the input file path.

        The file type is either 'CC' for a chemical component definition or another file type
        supported by OE toolkit assumed to have a conventional file extension for this type.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__fitPath = ccPath
        if fType in ["CC"]:
            (self.__fitId, self.__fitmol, self.__fitFD) = self.getCCDefFile(ccPath, molBuildType=molBuildType, suppressHydrogens=suppressHydrogens)
        else:
            (self.__fitId, self.__fitmol, self.__fitFD) = self.__getMiscFile(
                ccPath, suppressHydrogens=suppressHydrogens, importType=importType, title=title, largestPart=largestPart
            )

        if self.__verbose:
            logger.debug("Derived fit ID     = %s", self.__fitId)
            logger.debug("SMILES (stereo)  = %s", self.__fitFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__fitmol.SetTitle(title)
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__refId)
            self.__fitTitle = None
        return self.__fitmol.NumAtoms() if self.__fitmol else 0

    def setFitObj(self, ccObj, title=None, molBuildType="model-xyz", suppressHydrogens=False):
        """Set the object for the target/library molecule for MCSS comparison using the input file path.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        (self.__fitId, self.__fitmol, self.__fitFD) = self.getCCDefObj(ccObj, molBuildType=molBuildType, suppressHydrogens=suppressHydrogens)
        if self.__verbose:
            logger.debug("Derived fit ID     = %s", self.__fitId)
            logger.debug("SMILES (stereo)  = %s", self.__fitFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__fitmol.SetTitle(title)
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__refId)
            self.__fitTitle = None
        return self.__fitmol.NumAtoms() if self.__fitmol else 0

    def __setFitIdList(self, ccIdList, cachePath="/data/components/ligand-dict-v3"):
        """Set the list of IDs to be compared with reference molecule by MCSS.

        From the input ID list build the internal pair list of
        tuples  [(refId,refPath,refTitle,fitId,fitPath,fitTitle),(),...]
        """
        self.__pairTupleList = []
        for ccId in ccIdList:
            refId = self.__refId
            refPath = self.__refPath
            #
            fitId = ccId
            ccIdU = ccId.upper()
            fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            refTitle = refId + "/" + fitId
            fitTitle = fitId + "/" + refId
            self.__pairTupleList.append((refId, refPath, refTitle, fitId, fitPath, fitTitle))

    def __setPairIdList(self, pairIdList, cachePath="/data/components/ligand-dict-v3"):
        """Set the list of ID pais to be aligned by MCSS.

        From the input ID list build the internal pair list of
        tuples  [(refId,refPath,refTitle,fitId,fitPath,fitTitle),(),...]
        """

        self.__pairTupleList = []
        for refId, fitId in pairIdList:
            ccIdU = refId.upper()
            refPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            ccIdU = fitId.upper()
            fitPath = os.path.join(cachePath, ccIdU[0], ccIdU, ccIdU + ".cif")
            #
            refTitle = refId + "/" + fitId
            fitTitle = fitId + "/" + refId
            self.__pairTupleList.append((refId, refPath, refTitle, fitId, fitPath, fitTitle))

    def getCCDefFile(self, ccFilePath, molBuildType="model-xyz", suppressHydrogens=False):
        """Fetch the molecule definition (ccPath) and build OE molecules
        for comparison.

        """
        #
        mU = MarshalUtil(workPath=self.__workPath)
        rdCcObjL = mU.doImport(ccFilePath, fmt="mmcif")
        oemf = OeMoleculeFactory()
        if not self.__verbose:
            oemf.setQuiet()
        ccId = oemf.setChemCompDef(rdCcObjL[0])
        oemf.build(molBuildType=molBuildType)

        if self.__verbose:
            logger.info("  CCId               = %s", ccId)
            logger.info("  Title              = %s", oemf.getTitle())
            logger.info("  SMILES             = %s", oemf.getCanSMILES())
            logger.info("  SMILES (stereo)    = %s", oemf.getIsoSMILES())
            logger.info("  Formula (Hill)     = %s", oemf.getFormula())
            logger.info("  InChI key          = %s", oemf.getInChIKey())
            logger.info("  InChI              = %s", oemf.getInChI())

        fD = {}
        fD = {"Formula": oemf.getFormula(), "SMILES": oemf.getCanSMILES(), "SMILES_STEREO": oemf.getIsoSMILES(), "InChI": oemf.getInChI(), "InChIKey": oemf.getInChIKey()}

        if suppressHydrogens:
            tMol = oemf.getGraphMolSuppressH()
        else:
            tMol = oemf.getMol()

        fD["OEMOL"] = tMol
        fD["xyz"] = oemf.getAtomDetails(xyzType="model")

        return (ccId, tMol, fD)

    def getCCDefObj(self, dataContainer, molBuildType="model-xyz", suppressHydrogens=False):
        """Build OE molecule from the input chemical component definition object."""
        #
        oemf = OeMoleculeFactory()
        if not self.__verbose:
            oemf.setQuiet()
        ccId = oemf.setChemCompDef(dataContainer)
        oemf.build(molBuildType=molBuildType)

        if self.__verbose:
            logger.info("  CCId               = %s", ccId)
            logger.info("  Title              = %s", oemf.getTitle())
            logger.info("  SMILES             = %s", oemf.getCanSMILES())
            logger.info("  SMILES (stereo)    = %s", oemf.getIsoSMILES())
            logger.info("  Formula (Hill)     = %s", oemf.getFormula())
            logger.info("  InChI key          = %s", oemf.getInChIKey())
            logger.info("  InChI              = %s", oemf.getInChI())

        fD = {}
        fD = {"Formula": oemf.getFormula(), "SMILES": oemf.getCanSMILES(), "SMILES_STEREO": oemf.getIsoSMILES(), "InChI": oemf.getInChI(), "InChIKey": oemf.getInChIKey()}

        if suppressHydrogens:
            tMol = oemf.getGraphMolSuppressH()
        else:
            tMol = oemf.getMol()

        fD["OEMOL"] = tMol
        fD["xyz"] = oemf.getAtomDetails(xyzType="model")

        return (ccId, tMol, fD)

    def __getMiscFile(self, filePath, suppressHydrogens=False, importType="2D", title=None, largestPart=False):
        """Fetch a miscellaneous chemical file (ccPath) and build OE molecules
        for comparison.

        """
        try:
            oeioU = OeIoUtils()
            oeMolL = oeioU.fileToMols(filePath, use3D=importType == "3D", largestPart=largestPart)
            logger.info("Read (%d) from %s ", len(oeMolL), filePath)
            oeMol = oeMolL[0]

            ccId = title if title else oeMol.GetTitle()
            if title:
                oeMol.SetTitle(ccId)
            #
            oemf = OeMoleculeFactory()
            if not self.__verbose:
                oemf.setQuiet()
            oemf.setOeMol(oeMol, ccId)
            #
            fD = oemf.getOeMoleculeFeatures()
            if self.__verbose:
                logger.info("  Title              = %s", title)
                logger.info("  Title OEMF         = %s", oemf.getTitle())
                logger.info("  SMILES             = %s", oemf.getCanSMILES())
                logger.info("  SMILES (stereo)    = %s", oemf.getIsoSMILES())
                logger.info("  Formula (Hill)     = %s", oemf.getFormula())
                logger.info("  InChI key          = %s", oemf.getInChIKey())
                logger.info("  InChI              = %s", oemf.getInChI())
            #
            ccId = oemf.getTitle()
            if suppressHydrogens:
                tMol = oemf.getGraphMolSuppressH()
            else:
                tMol = oemf.getMol()

            molXyzL = []
            if importType == "3D":
                for atm in tMol.GetAtoms():
                    xyzL = oechem.OEFloatArray(3)
                    tMol.GetCoords(atm, xyzL)
                    molXyzL.append(
                        ComponentAtomDetails(
                            atIdx=atm.GetIdx(),
                            atNo=atm.GetAtomicNum(),
                            atName=atm.GetName(),
                            atType=atm.GetType(),
                            x=xyzL[0],
                            y=xyzL[1],
                            z=xyzL[2],
                            atFormalCharge=atm.GetFormalCharge(),
                        )
                    )
            fD = {}
            fD = {
                "Formula": oemf.getFormula(),
                "SMILES": oemf.getCanSMILES(),
                "SMILES_STEREO": oemf.getIsoSMILES(),
                "InChI": oemf.getInChI(),
                "InChIKey": oemf.getInChIKey(),
                "xyz": molXyzL,
            }
            for atm in tMol.GetAtoms():
                xyzL = oechem.OEFloatArray(3)
                tMol.GetCoords(atm, xyzL)
                if self.__verbose:
                    logger.debug("atom  %s %s %s %s %r", atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL)

            fD["OEMOL"] = tMol
            return (ccId, tMol, fD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return None, None, None

    def __setupMCSS(self, refmol, useExhaustive=True):
        """Internal initialization for the MCSS comparison."""
        #
        mode = oechem.OEMCSType_Exhaustive if useExhaustive else oechem.OEMCSType_Approximate
        self.__mcss = oechem.OEMCSSearch(mode)
        atomexpr, bondexpr = OeCommonUtils.getAtomBondExprOpts(self.__searchType)
        self.__mcss.Init(refmol, atomexpr, bondexpr)
        if self.__verbose:
            logger.info("Initialize MCSS (%r)", self.__searchType)

    def __setupSubStructure(self, refmol):
        """Internal initialization for a substructure comparison."""
        #
        atomexpr, bondexpr = OeCommonUtils.getAtomBondExprOpts(self.__searchType)
        self.__ss = oechem.OESubSearch(refmol, atomexpr, bondexpr)
        if self.__verbose:
            logger.info("Initialize SS (%r)", self.__searchType)

    @timeout("instance.timeOut", use_signals=False, dec_allow_eval=True)
    def doAlignSs(self, unique=True, maxMatches=20):
        """Test the SS comparison between current reference and fit molecules -
        Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMapL = []
        fitAtomUnMappedL = []
        #
        nAtomsRef = self.__refmol.NumAtoms()
        nAtomsFit = self.__fitmol.NumAtoms()
        # -------
        oechem.OEAddExplicitHydrogens(self.__refmol)
        oechem.OEAddExplicitHydrogens(self.__fitmol)
        fitAtD = {}
        #
        for at in self.__fitmol.GetAtoms():
            nAtL = at.GetAtoms()
            neighbors = [nAt.GetName() for nAt in nAtL]
            atType = oechem.OEGetAtomicSymbol(at.GetAtomicNum())
            xyzL = oechem.OEFloatArray(3)
            self.__fitmol.GetCoords(at, xyzL)
            fitAtD[at.GetIdx()] = AlignAtomUnMapped(
                fitId=self.__fitId,
                fitAtIdx=at.GetIdx(),
                fitAtName=at.GetName(),
                fitAtType=atType,
                fitAtNo=at.GetAtomicNum(),
                fitAtFormalCharge=at.GetFormalCharge(),
                x=xyzL[0],
                y=xyzL[1],
                z=xyzL[2],
                fitNeighbors=neighbors,
            )

        #
        logger.debug("nAtomsRef %d nAtomsFit %d", nAtomsRef, nAtomsFit)
        #
        # --------
        self.__setupSubStructure(self.__refmol)
        self.__ss.SetMaxMatches(maxMatches)
        miter = self.__ss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            for mAt in match.GetAtoms():
                atomMapL.append(
                    AlignAtomMap(
                        refId=self.__refId,
                        refAtIdx=mAt.pattern.GetIdx(),
                        refAtNo=mAt.pattern.GetAtomicNum(),
                        refAtName=mAt.pattern.GetName(),
                        fitId=self.__fitId,
                        fitAtIdx=mAt.target.GetIdx(),
                        fitAtNo=mAt.target.GetAtomicNum(),
                        fitAtName=mAt.target.GetName(),
                    )
                )
                fitAtD.pop(mAt.target.GetIdx())
            logger.debug("fitAtD %r", fitAtD)
        fitAtomUnMappedL = list(fitAtD.values())
        return (nAtomsRef, self.__refFD, nAtomsFit, self.__fitFD, atomMapL, fitAtomUnMappedL)

    @timeout("instance.timeOut", use_signals=False, dec_allow_eval=True)
    def doAlignMcss(self, unique=True, minFrac=1.0, useExhaustive=True):
        """Test the MCSS comparison between current reference and fit molecules -
        Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMapL = []
        fitAtomUnMappedL = []
        #
        nAtomsRef = self.__refmol.NumAtoms()
        nAtomsFit = self.__fitmol.NumAtoms()
        fitAtD = {}
        for at in self.__fitmol.GetAtoms():
            nAtL = at.GetAtoms()
            neighbors = [nAt.GetName() for nAt in nAtL]
            atType = oechem.OEGetAtomicSymbol(at.GetAtomicNum())
            xyzL = oechem.OEFloatArray(3)
            self.__fitmol.GetCoords(at, xyzL)
            fitAtD[at.GetIdx()] = AlignAtomUnMapped(
                fitId=self.__fitId,
                fitAtIdx=at.GetIdx(),
                fitAtName=at.GetName(),
                fitAtType=atType,
                fitAtNo=at.GetAtomicNum(),
                fitAtFormalCharge=at.GetFormalCharge(),
                x=xyzL[0],
                y=xyzL[1],
                z=xyzL[2],
                fitNeighbors=neighbors,
            )
        minAtoms = int(min(nAtomsRef, nAtomsFit) * minFrac)
        # -------
        self.__setupMCSS(self.__refmol, useExhaustive=useExhaustive)
        self.__mcss.SetMCSFunc(oechem.OEMCSMaxAtoms())
        self.__mcss.SetMinAtoms(minAtoms)
        oechem.OEAddExplicitHydrogens(self.__refmol)
        oechem.OEAddExplicitHydrogens(self.__fitmol)
        #
        # --------
        miter = self.__mcss.Match(self.__fitmol, unique)
        if miter.IsValid():
            match = miter.Target()
            for mAt in match.GetAtoms():

                atomMapL.append(
                    AlignAtomMap(
                        refId=self.__refId,
                        refAtIdx=mAt.pattern.GetIdx(),
                        refAtNo=mAt.pattern.GetAtomicNum(),
                        refAtName=mAt.pattern.GetName(),
                        fitId=self.__fitId,
                        fitAtIdx=mAt.target.GetIdx(),
                        fitAtNo=mAt.target.GetAtomicNum(),
                        fitAtName=mAt.target.GetName(),
                    )
                )
                fitAtD.pop(mAt.target.GetIdx())
        fitAtomUnMappedL = list(fitAtD.values())
        return (nAtomsRef, self.__refFD, nAtomsFit, self.__fitFD, atomMapL, fitAtomUnMappedL)
