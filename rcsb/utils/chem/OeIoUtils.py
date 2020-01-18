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
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class OeIoUtils(object):
    """ Utility methods to manage OE specific IO and format conversion operations.
    """

    def __init__(self, **kwargs):
        self.__dirPath = kwargs.get("dirPath", ".")
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        #

    def getComponentDefinitions(self, ccdFilePath):
        try:
            rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
            logger.info("Read %s with %d definitions", ccdFilePath, len(rdCcObjL))
        except Exception as e:
            logger.exception("Loading %s failing with %s", ccdFilePath, str(e))
        return rdCcObjL

    def chemCompToMol(self, ccdFilePath, coordType="model", quietFlag=False):
        retMolL = []
        try:
            rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
            logger.info("Read %s with %d definitions", ccdFilePath, len(rdCcObjL))
            oemf = OeMoleculeFactory()
            if quietFlag:
                oemf.setQuiet()
            for ccObj in rdCcObjL:
                ccId = oemf.set(ccObj)
                if ccId:
                    if coordType:
                        ok = oemf.build3D(coordType=coordType)
                    else:
                        ok = oemf.build2D()
                    if ok:
                        oeMol = oemf.getMol()
                        retMolL.append(oeMol)
        except Exception as e:
            logger.exception("Loading %s failing with %s", ccdFilePath, str(e))
        return retMolL

    def smilesToMol(self, smiles, limitPerceptions=True):
        """Parse the input SMILES string and return a molecule object (OeGraphMol).

        Args:
            smiles (str): SMILES string
            limitPerceptions (bool): flag to limit the perceptions/transformations of input SMILES

        Returns:
            object: OeGraphMol() object or None for failure
        """
        try:
            mol = oechem.OEGraphMol()
            smiles.strip()
            if limitPerceptions:
                # convert the SMILES string into a molecule
                if oechem.OEParseSmiles(mol, smiles, False, True):
                    return mol
                else:
                    logger.error("Parsing failed for SMILES string %s", smiles)
            else:
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

    def fileToMols(self, filePath, use3D=False):
        """Parse the input path returning a list of molecule objects (OeGraphMol).

        Args:
            filePath (str): file path must have strandard recognized extension ('mol', 'sdf', 'smi', 'oeb').

        Returns:
            list : list of OeGraphMol() objects
        """
        mL = []
        oemf = OeMoleculeFactory()
        try:
            ifs = oechem.oemolistream(filePath)
            for mol in ifs.GetOEGraphMols():
                if use3D:
                    mL.append(oemf.updateOePerceptions3D(mol, aromaticModel=oechem.OEAroModelOpenEye))
                else:
                    mL.append(oemf.updateOePerceptions2D(mol, aromaticModel=oechem.OEAroModelOpenEye))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return mL

    def stringToMols(self, txt, sType="mol2", use3D=False):
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
        oemf = OeMoleculeFactory()
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
                if use3D:
                    mL.append(oemf.updateOePerceptions3D(mol, aromaticModel=oechem.OEAroModelOpenEye))
                else:
                    mL.append(oemf.updateOePerceptions2D(mol, aromaticModel=oechem.OEAroModelOpenEye))

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
            molDb = oechem.OEMolDatabase()
            if not molDb.Open(oeMolDbFilePath):
                logger.error("Unable to open %r", oeMolDbFilePath)
            molCount = molDb.NumMols()
            logger.info("Loaded OE database file containing %d molecules", molCount)
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

    def write(self, filePath, oeMol, constantMol=False):
        """Write an oeMol with format type inferred from the filePath extension (e.g. .mol)

        Args:
            filePath (str): filepath with a chemical type extension
            constantMol (bool, optional): copies molecule before performing format specific perceptions

        Returns:
            bool: True for success or False otherwise
        """
        try:
            ofs = oechem.oemolostream()
            ofs.open(filePath)
            logger.info("Writing %s title %s\n", filePath, oeMol.GetTitle())
            if constantMol:
                oechem.OEWriteConstMolecule(ofs, oeMol)
            else:
                oechem.OEWriteMolecule(ofs, oeMol)
            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", filePath, str(e))
        return False

    def serializeOe(self, oeMol):
        """ Create a string representing the content of the current OE molecule.   This
            serialization uses the OE internal binary format.
        """
        try:
            oms = oechem.oemolostream()
            oms.SetFormat(oechem.OEFormat_OEB)
            oms.openstring()
            oechem.OEWriteMolecule(oms, oeMol)
            logger.debug("SMILES %s", oechem.OECreateCanSmiString(oeMol))
            logger.debug("Atoms = %d", oeMol.NumAtoms())
            return oms.GetString()
        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def deserializeOe(self, oeS):
        """ Reconstruct an OE molecule from the input string serialization (OE binary).

            The deserialized molecule is used to initialize the internal OE molecule
            within this object.

            Returns:
                list:  OE GraphMol list
        """
        molList = []
        try:
            ims = oechem.oemolistream()
            ims.SetFormat(oechem.OEFormat_OEB)
            ims.openstring(oeS)
            for mol in ims.GetOEGraphMols():
                logger.debug("SMILES %s", oechem.OECreateCanSmiString(mol))
                logger.debug("title  %s", mol.GetTitle())
                logger.debug("atoms  %d", mol.NumAtoms())
                molList.append(oechem.OEGraphMol(mol))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return molList
