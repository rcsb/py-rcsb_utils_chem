##
# File:    OeIoUtils.py
# Author:  jdw
# Date:    21-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to manage OE specific IO and format conversion operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import time

from openeye import oechem
from openeye import oegraphsim

from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory


logger = logging.getLogger(__name__)


class OeIoUtils(object):
    """ Utility methods to manage OE specific IO and format conversion operations.
    """

    def __init__(self, verbose=False):
        self.__verbose = verbose
        #

    def smilesToMol(self, smiles):
        """Parse the input SMILES string and return a molecule object (OeGraphMol).

        Args:
            smiles (str): SMILES string

        Returns:
            object: OeGraphMol() object or None for failure
        """
        try:
            mol = oechem.OEGraphMol()
            # convert the SMILES string into a molecule
            smiles.strip()
            if oechem.OESmilesToMol(mol, smiles):
                return mol
            else:
                logger.error("Parsing failed for SMILES string %s", smiles)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def inchiToMol(self, inchi):
        """Parse the input InChI string and return a molecule object (OeGraphMol).

        Args:
            inchi (str): InChI string

        Returns:
            object: OeGraphMol() object or None for failure
        """
        try:
            mol = oechem.OEGraphMol()
            inchi = inchi.strip()
            if oechem.OEInChIToMol(mol, inchi):
                return mol
            else:
                logger.error("Parsing failed for InChI string %s", inchi)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def smartsToQmol(self, smarts):
        """Parse the input SMARTS query string and return a query molecule object (OeQMol).

        Args:
            smarts (str): SMARTS query string

        Returns:
            object : OeQMol() object or None for failure
        """
        try:
            qmol = oechem.OEQMol()
            if oechem.OEParseSmarts(qmol, smarts):
                return qmol
            else:
                logger.error("Parsing failed for SMARTS string %s", smarts)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def fileToMols(self, filePath):
        """Parse the input path returning a list of molecule objects (OeGraphMol).

        Args:
            filePath (str): file path must have strandard recognized extension ('mol', 'sdf', 'smi', 'oeb').

        Returns:
            list : list of OeGraphMol() objects
        """
        mL = []
        try:
            ifs = oechem.oemolistream(filePath)
            for mol in ifs.GetOEGraphMols():
                mL.append(self.__standardPerceptions(oechem.OEGraphMol(mol)))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return mL

    def stringToMols(self, txt, sType="mol2"):
        """Parse the input string as input format type (sType) returning a list of
        molecule objects (OeGraphMol)

        Args:
            txt (str): string text of molecule data
            sType (str, optional): string data format (mol2, sdf, smiles) . Defaults to "mol2".

        Returns:
            list: list of OeGraphMol() objects
        """
        #
        mL = []
        try:
            if sType not in ["mol2", "sdf", "smiles"]:
                logger.error("Unsupported string data format")
                return None
            fD = {"mol2": oechem.OEFormat_MOL2, "sdf": oechem.OEFormat_SDF, "smiles": oechem.OEFormat_SMI}
            ifs = oechem.oemolistream()
            ifs.SetFormat(fD["sType"])
            if not ifs.openstring(txt):
                logger.error("Unable open string data for molecule reader")
                return None
            for mol in ifs.GetOEGraphMols():
                mL.append(self.__standardPerceptions(oechem.OEGraphMol(mol)))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return mL

    def readOeBinaryMolCache(self, filePath):
        """Return a list of OeGraphMol() objects read from the cached binary file.

        Args:
            filePath (str): file path for the binary OeMol cache

        Returns:
            dict: dictionary of OeGraphMol()'s {<ccId>: OeGraphMol(), ... }
        """
        retD = {}
        startTime = time.time()
        try:
            ifs = oechem.oemolistream()
            if ifs.open(filePath):
                for oeMol in ifs.GetOEGraphMols():
                    tMol = oechem.OEGraphMol(oeMol)
                    retD[tMol.GetTitle()] = tMol
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return retD

    def createOeFingerPrintDatabase(self, oeMolDbFilePath, oeFpDbFilePath, fpType="TREE"):
        """ Create fast search fingerprint database from the input molecular database.

        Args:
            oeMolDbFilePath (str): path to the input molecular database
            oeFpDbFilePath (str): path to the output fingerprint database
            fpType (str):  finger print type

        Returns:
            bool: True for success or False otherwise

        Supports:
            OEFPType_Circular
            OEFPType_Path
            OEFPType_Tree

        Not currently supported by OE fp search -
            OEFPType_MACCS166
            OEFPType_Lingo
        """
        try:
            _ = fpType
            startTime = time.time()
            fpD = {"TREE": oegraphsim.OEFPType_Tree, "CIRCULAR": oegraphsim.OEFPType_Circular, "PATH": oegraphsim.OEFPType_Path}
            myFpType = fpD[fpType] if fpType in fpD else oegraphsim.OEFPType_Tree
            opts = oegraphsim.OECreateFastFPDatabaseOptions(oegraphsim.OEGetFPType(myFpType))
            ok = oegraphsim.OECreateFastFPDatabaseFile(oeFpDbFilePath, oeMolDbFilePath, opts)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return ok

    def loadOeFingerPrintDatabase(self, oeFpDbFilePath, inMemory=True, fpType="TREE"):
        #
        _ = fpType
        startTime = time.time()
        if inMemory:
            memType = oegraphsim.OEFastFPDatabaseMemoryType_InMemory
        else:
            memType = oegraphsim.OEFastFPDatabaseMemoryType_MemoryMapped
        fpDb = oegraphsim.OEFastFPDatabase(oeFpDbFilePath, memType)
        if not fpDb.IsValid():
            logger.error("Cannot open fingerprint database %r", oeFpDbFilePath)
        #
        lenFp = fpDb.NumFingerPrints()
        memTypeStr = fpDb.GetMemoryTypeString()
        endTime = time.time()
        logger.info("Loaded fingerprint data length %d loaded %s (%.4f seconds)", lenFp, memTypeStr, endTime - startTime)
        return fpDb

    def loadOeBinaryDatabaseAndIndex(self, oeMolDbFilePath):
        molDb = None
        try:
            moldb = oechem.OEMolDatabase()
            if not moldb.Open(oeMolDbFilePath):
                logger.error("Unable to open %r", oeMolDbFilePath)
        except Exception as e:
            logger.exception("Loading %r failing with %s", oeMolDbFilePath, str(e))
        return molDb

    def createOeBinaryDatabaseAndIndex(self, oebMolFilePath, oeMolDbFilePath):
        """Create OE binary database file and associated index from the input serial
        binary data file.

        Args:
            oebMolFilePath (str): input OeMol stream binary file path
            oeMolDbFilePath (str): output OeMolDatabase file path

        Returns:
           int:  number of molecules processed in the database.
        """
        molCount = 0
        try:
            startTime = time.time()
            moldb = oechem.OEMolDatabase()
            if not moldb.Open(oebMolFilePath):
                logger.error("Read fails for %r", oebMolFilePath)
                return molCount
            #
            logger.info("Opened database in format %r num mols %d max index %d", moldb.GetFormat(), moldb.NumMols(), moldb.GetMaxMolIdx())
            moldb.Save(oeMolDbFilePath)
            tL = list(moldb.GetTitles())
            logger.info("First and last titles: %r %r", tL[0], tL[-1])
            molCount = moldb.NumMols()
            endTime = time.time()
            logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return molCount

    def buildOeBinaryMolCache(self, filePath, ccObjL, coordType="model", quietFlag=False):
        """Build cache of OEGraphMol() objects from the input chemical component definition list.

        Args:
            ccObjL (list): Chemical component dataContainer object list

        Returns:
          (int, int) : chemical component count processed successes and failures
        """
        try:
            ccCount = 0
            errCount = 0
            startTime = time.time()
            ofs = oechem.oemolostream()
            ofs.SetFormat(oechem.OEFormat_OEB)
            if ofs.open(filePath):
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                for ccObj in ccObjL:
                    ccId = oemf.set(ccObj)
                    if ccId:
                        if coordType:
                            ok = oemf.build3D(coordType=coordType)
                        else:
                            ok = oemf.build2D()

                        if ok:
                            oeMol = oemf.getMol()
                            oechem.OEWriteMolecule(ofs, oeMol)
                            ccCount += 1
                    else:
                        # incomplete component (e.g. missing atoms or bonds)
                        errCount += 1
            else:
                logger.error("Unable to open cache database %s", filePath)
                errCount += 1
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return ccCount, errCount

    def __standardPerceptions(self, oeMol, dType="2D"):
        try:
            if dType == "2D":
                # run standard perceptions --
                oechem.OEFindRingAtomsAndBonds(oeMol)
                oechem.OEPerceiveChiral(oeMol)

                for oeAt in oeMol.GetAtoms():
                    st = oeAt.GetStringData("StereoInfo")
                    if st == "R":
                        oechem.OESetCIPStereo(oeMol, oeAt, oechem.OECIPAtomStereo_R)
                    elif st == "S":
                        oechem.OESetCIPStereo(oeMol, oeAt, oechem.OECIPAtomStereo_S)

                for oeBnd in oeMol.GetBonds():
                    st = oeBnd.GetStringData("StereoInfo")
                    if st == "E":
                        oechem.OESetCIPStereo(oeMol, oeBnd, oechem.OECIPBondStereo_E)
                    elif st == "Z":
                        oechem.OESetCIPStereo(oeMol, oeBnd, oechem.OECIPBondStereo_Z)
            elif dType == "3D":
                # run standard perceptions --
                #
                oeMol.SetDimension(3)
                oechem.OE3DToInternalStereo(oeMol)
                oechem.OEFindRingAtomsAndBonds(oeMol)
                # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
                oechem.OEAssignAromaticFlags(oeMol, oechem.OEAroModelOpenEye)
                for atom in oeMol.GetAtoms():
                    oechem.OEPerceiveCIPStereo(oeMol, atom)

                for bond in oeMol.GetBonds():
                    if bond.GetOrder() == 2:
                        oechem.OEPerceiveCIPStereo(oeMol, bond)
                oechem.OEAddExplicitHydrogens(oeMol)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oeMol

    def createOeSubSearchDatabase(self, oebMolFilePath, oeSubSearchFilePath, screenType="SMARTS", numProc=2):
        sort = True
        keepTitle = True
        myScreenType = None
        if screenType == "MOLECULE":
            myScreenType = oechem.OEGetSubSearchScreenType(oechem.OESubSearchScreenType_Molecule)
        elif screenType == "MDL":
            myScreenType = oechem.OEGetSubSearchScreenType(oechem.OESubSearchScreenType_MDL)
        elif screenType == "SMARTS":
            myScreenType = oechem.OEGetSubSearchScreenType(oechem.OESubSearchScreenType_SMARTS)

        opts = oechem.OECreateSubSearchDatabaseOptions(myScreenType)
        opts.SetSortByBitCounts(sort)
        opts.SetKeepTitle(keepTitle)
        opts.SetNumProcessors(numProc)

        screenStr = myScreenType.GetName()
        logger.info("Using %d processor(s) to generate database with %s", numProc, screenStr)

        tracer = oechem.OEConsoleProgressTracer()
        ok = oechem.OECreateSubSearchDatabaseFile(oeSubSearchFilePath, oebMolFilePath, opts, tracer)
        return ok

    def loadOeSubSearchDatabase(self, oeSubSearchFilePath, screenType=None, numProc=1):
        ssDb = None
        try:
            _ = screenType
            ssDb = oechem.OESubSearchDatabase(oechem.OESubSearchDatabaseType_Default, numProc)
            tracer = oechem.OEConsoleProgressTracer()
            if not ssDb.Open(oeSubSearchFilePath, tracer):
                logger.error("Unable to open %r", oeSubSearchFilePath)
        except Exception as e:
            logger.exception("Loading %r failing with %s", oeSubSearchFilePath, str(e))
        return ssDb
