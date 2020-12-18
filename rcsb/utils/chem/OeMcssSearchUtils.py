##
# File:  OeMcssSearchUtils.py
# Date:  17-Dec-2020  J. Westbrook
#
# Updates:
##
"""
Utilities for performing MCSS comparisons.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os

from openeye import oechem

from rcsb.utils.chem.OeCommonUtils import OeCommonUtils
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.decorators import timeout
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class OeMcssSearchUtils(object):
    """Perform MCSS comparisons.  Targets can be chemical component identifiers
    or paths to chemical component definition files.  Inputs can be in the the form of pairs,
    lists, and pair lists of chemical component definitions.

    """

    def __init__(self, verbose=True):
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
        self.__minAtomMatchFraction = 0.50
        #
        self.__searchType = "relaxed"
        #
        self.__refFD = {}
        self.__fitFD = {}
        #
        self.__mcss = None
        self.__refPath = None
        self.__fitPath = None

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
        self.__setupMCSS(self.__refmol)

    def setRefPath(self, ccPath, title=None, suppressHydrogens=False, fType="CC", importType="2D"):
        """Set the query molecule for MCSS comparison using the input file path.

        The file type is either ['CC'] for a chemical component definition or another file type
        supported by OE toolkit assumed to have a conventional file extension for this type.

        Once the reference molecule is built, the MCSS calculation is initialized.

        A title is optionally provided otherwise the component Id will be used.

        The hydrogen flag can be used to perform the MCSS using only heavy atoms.
        """
        self.__refPath = ccPath
        if fType in ["CC"]:
            (self.__refId, self.__refmol, self.__refFD) = self.getCCDefFile(ccPath, suppressHydrogens=suppressHydrogens)
        else:
            (self.__refId, self.__refmol, self.__refFD) = self.__getMiscFile(ccPath, suppressHydrogens=suppressHydrogens, importType=importType)

        if self.__verbose:
            logger.info("Derived ref ID     = %s", self.__refId)
            logger.info("SMILES (stereo)  = %s", self.__refFD["SMILES_STEREO"])
        #
        # Insert title here -
        if title is not None:
            self.__refmol.SetTitle(title)
            self.__refTitle = title
        else:
            self.__refmol.SetTitle(self.__refId)
            self.__refTitle = None

        self.__setupMCSS(self.__refmol)

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

    def setFitPath(self, ccPath, title=None, suppressHydrogens=False):
        """Set the path of the target/library molecule for MCSS comparison."""
        (self.__fitId, self.__fitmol, self.__fitFD) = self.getCCDefFile(ccPath, suppressHydrogens=suppressHydrogens)
        if title is not None:
            self.__fitmol.SetTitle(title)
            self.__fitTitle = title
        else:
            self.__fitmol.SetTitle(self.__fitId)
            self.__fitTitle = None

    def setFitIdList(self, ccIdList, cachePath="/data/components/ligand-dict-v3"):
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

    def setPairIdList(self, pairIdList, cachePath="/data/components/ligand-dict-v3"):
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
        mU = MarshalUtil()
        rdCcObjL = mU.doImport(ccFilePath, fmt="mmcif")
        oemf = OeMoleculeFactory()
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

        return (ccId, tMol, fD)

    def __getMiscFile(self, filePath, suppressHydrogens=False, importType="2D"):
        """Fetch a miscellaneous chemical file (ccPath) and build OE molecules
        for comparison.

        """
        try:
            oeioU = OeIoUtils()
            oeMolL = oeioU.fileToMols(filePath, use3D=importType == "3D")
            oeMol = oeMolL[0]
            ccId = oeMol.GetTitle()
            #
            oemf = OeMoleculeFactory()
            oemf.setOeMol(oeMol, ccId)
            #
            fD = oemf.getOeMoleculeFeatures()
            if self.__verbose:
                logger.info("+OeMCSS.__getMiscFile()")
                logger.info("  Title              = %s", oemf.getTitle())
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
                for ii, atm in enumerate(tMol.GetAtoms()):
                    xyzL = oechem.OEFloatArray(3)
                    tMol.GetCoords(atm, xyzL)
                    molXyzL.append((ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), "%.3f" % xyzL[0], "%.3f" % xyzL[1], "%.3f" % xyzL[2]))
            fD = {}
            fD = {
                "Formula": oemf.getFormula(),
                "SMILES": oemf.getCanSMILES(),
                "SMILES_STEREO": oemf.getIsoSMILES(),
                "InChI": oemf.getInChI(),
                "InChIKey": oemf.getInChIKey(),
                "xyz": molXyzL,
            }

            for ii, atm in enumerate(tMol.GetAtoms()):
                xyzL = oechem.OEFloatArray(3)
                tMol.GetCoords(atm, xyzL)
                if self.__verbose:
                    logger.debug("atom  %d %s %s %s %s %r", ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL)

            fD["OEMOL"] = tMol
            return (ccId, tMol, fD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return None, None, None

    def __setupMCSS(self, refmol):
        """Internal initialization for the MCSS comparison."""
        #
        self.__mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
        # self.__mcss = oechem.OEMCSSearch(oechem.OEMCSType_Approximate)
        atomexpr, bondexpr = OeCommonUtils.getAtomBondExprOpts(self.__searchType)
        self.__mcss.Init(refmol, atomexpr, bondexpr)
        logger.info("Initialize MCSS (%r)", self.__searchType)

    @timeout(100)
    def doAlign(self, unique=True, minFrac=1.0):
        """Test the MCSS comparison between current reference and fit molecules -
        Return list of corresponding atoms on success or an empty list otherwise.
        """
        atomMap = []
        #
        nAtomsRef = self.__refmol.NumAtoms()
        nAtomsFit = self.__fitmol.NumAtoms()
        minAtoms = int(min(nAtomsRef, nAtomsFit) * minFrac)
        # -------
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
                atomMap.append(
                    (
                        self.__refId,
                        mAt.pattern.GetIdx(),
                        mAt.pattern.GetAtomicNum(),
                        mAt.pattern.GetName(),
                        self.__fitId,
                        mAt.target.GetIdx(),
                        mAt.target.GetAtomicNum(),
                        mAt.target.GetName(),
                    )
                )
        return (nAtomsRef, self.__refFD, nAtomsFit, self.__fitFD, atomMap)
