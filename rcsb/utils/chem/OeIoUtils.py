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
import os
import time

from openeye import oechem
from openeye import oegraphsim
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class OeIoUtils(object):
    """Utility methods to manage OE specific IO and format conversion operations."""

    def __init__(self, **kwargs):
        self.__dirPath = kwargs.get("dirPath", ".")
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oeErrorLevel = oechem.OEErrorLevel_Info
        if kwargs.get("quietFlag", False):
            self.setQuiet()
        #

    def setQuiet(self):
        """Suppress OE warnings and processing errors"""
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Quiet)
        self.__oeErrorLevel = oechem.OEErrorLevel_Quiet

    def getComponentDefinitions(self, ccdFilePath):
        rdCcObjL = []
        try:
            rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
            logger.info("Read %s with %d definitions", ccdFilePath, len(rdCcObjL))
        except Exception as e:
            logger.exception("Loading %s failing with %s", ccdFilePath, str(e))
        return rdCcObjL

    def suppressHydrogens(self, oeMol):
        tMol = oechem.OEMol(oeMol) if oeMol else None
        if tMol:
            oechem.OESuppressHydrogens(tMol)
        return tMol

    def chemCompToMol(self, ccdFilePath, molBuildType="model-xyz", quietFlag=False):
        retMolL = []
        try:
            rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
            logger.info("Read %s with %d definitions", ccdFilePath, len(rdCcObjL))
            oemf = OeMoleculeFactory()
            if quietFlag:
                oemf.setQuiet()
            for ccObj in rdCcObjL:
                ccId = oemf.setChemCompDef(ccObj)
                if ccId:
                    ok = oemf.build(molBuildType=molBuildType)
                    if ok:
                        oeMol = oemf.getMol()
                        retMolL.append(oeMol)
        except Exception as e:
            logger.exception("Loading %s failing with %s", ccdFilePath, str(e))
        return retMolL

    def descriptorToSmiles(self, descr, descrType, limitPerceptions=False, messageTag=None):
        """Parse the input descriptor string and return an OE smiles.

        Args:
            descr (str): descriptor
            descrType (str): descriptor type
            limitPerceptions (bool): flag to limit the perceptions/transformations of input descriptor
            messageTag (srt, optional): prefix string for error messages. Defaults to None.

        Returns:
            str: SMILES string
        """
        try:
            if "SMILES" in descrType.upper() and "ISO" in descrType.upper():
                oeMol = self.smilesToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    return oechem.OECreateIsoSmiString(oeMol)
                else:
                    return None
            if "SMILES" in descrType.upper():
                oeMol = self.smilesToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    return oechem.OECreateCanSmiString(oeMol)
                else:
                    return None
            elif "INCHI" in descrType.upper():
                oeMol = self.inchiToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    return oechem.OECreateIsoSmiString(oeMol)
            else:
                return None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def descriptorToMol(self, descr, descrType, limitPerceptions=False, messageTag=None):
        """Parse the input descriptor string and return a molecule object (OeGraphMol/OeQMol).

        Args:
            descr (str): descriptor
            descrType (str): descriptor type
            limitPerceptions (bool): flag to limit the perceptions/transformations of input descriptor
            messageTag (srt, optional): prefix string for error messages. Defaults to None.

        Returns:
            object: OeGraphMol()/OeQmol() object or None for failure

            ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noq
        """
        try:
            if "SMILES" in descrType.upper() and "ISO" in descrType.upper():
                oeMol = self.smilesToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    isoSmiles = oechem.OECreateIsoSmiString(oeMol)
                    return self.smilesToMol(isoSmiles, limitPerceptions=limitPerceptions, messageTag=messageTag)
                else:
                    return None
            if "SMILES" in descrType.upper():
                oeMol = self.smilesToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    smiles = oechem.OECreateCanSmiString(oeMol)
                    return self.smilesToMol(smiles, limitPerceptions=limitPerceptions, messageTag=messageTag)
                else:
                    return None
            elif "INCHI" in descrType.upper():
                oeMol = self.inchiToMol(descr, limitPerceptions=limitPerceptions, messageTag=messageTag)
                if oeMol:
                    isoSmiles = oechem.OECreateIsoSmiString(oeMol)
                    return self.smilesToMol(isoSmiles, limitPerceptions=limitPerceptions, messageTag=messageTag)
            elif "SMARTS" in descrType.upper():
                return self.smartsToQmol(descr, messageTag=messageTag)
            else:
                return None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def smilesToMol(self, smiles, limitPerceptions=False, messageTag=None):
        """Parse the input SMILES string and return a molecule object (OeGraphMol).

        Args:
            smiles (str): SMILES string
            limitPerceptions (bool): flag to limit the perceptions/transformations of input SMILES

        Returns:
            object: OeGraphMol() object or None for failure
        """
        try:
            label = messageTag if messageTag else ""
            mol = oechem.OEGraphMol()
            smiles.strip()
            if limitPerceptions:
                # convert the SMILES string into a molecule
                if oechem.OEParseSmiles(mol, smiles, False, False):
                    return mol
                else:
                    logger.debug("%s parsing failed for input SMILES string %s", label, smiles)
                    logger.error("%s parsing failed for input SMILES string", label)
            else:
                if oechem.OESmilesToMol(mol, smiles):
                    return mol
                else:
                    logger.debug("%s converting failed for input SMILES string %s", label, smiles)
                    logger.error("%s converting failed for input SMILES string", label)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def inchiToMol(self, inchi, limitPerceptions=False, messageTag=None):
        """Parse the input InChI string and return a molecule object (OeGraphMol).

        Args:
            inchi (str): InChI string

        Returns:
            object: OeGraphMol() object or None for failure

        """
        try:
            label = messageTag if messageTag else ""
            mol = oechem.OEGraphMol()
            inchi = inchi.strip()
            if limitPerceptions:
                if oechem.OEParseInChI(mol, inchi):
                    return mol
                else:
                    logger.debug("%s parsing failed for InChI string %r", label, inchi)
                    logger.error("%s parsing failed for InChI string", label)
            else:
                if oechem.OEInChIToMol(mol, inchi):
                    return mol
                else:
                    logger.debug("%s converting failed for InChI string %r", label, inchi)
                    logger.error("%s converting failed for InChI string", label)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def descriptorToQMol(self, descr, descrType, limitPerceptions=False, messageTag=None):
        """Parse the input descriptor string and return a query molecule object (OeQMol).

        Args:
            descr (str): descriptor
            descrType (str): descriptor type
            limitPerceptions (bool): flag to limit the perceptions/transformations of input descriptor
            messageTag (srt, optional): prefix string for error messages. Defaults to None.

        Returns:
            object: OeQmol() object or None for failure

        """
        oeQMol = label = None
        try:
            label = messageTag if messageTag else ""
            tMol = self.descriptorToMol(descr, descrType, limitPerceptions=limitPerceptions, messageTag=messageTag)
            if tMol:
                oeQMol = oechem.OEQMol(tMol)

        except Exception as e:
            logger.error("%s Failing for with %s", label, str(e))
        return oeQMol if oeQMol else None

    def smartsToQmol(self, smarts, messageTag=None):
        """Parse the input SMARTS query string and return a query molecule object (OeQMol).

        Args:
            smarts (str): SMARTS query string

        Returns:
            object : OeQMol() object or None for failure
        """
        try:
            label = messageTag if messageTag else ""
            qmol = oechem.OEQMol()
            if oechem.OEParseSmarts(qmol, smarts):
                return qmol
            else:
                logger.debug("%s parsing failed for SMARTS string %s", label, smarts)
                logger.error("%s parsing failed for SMARTS string", label)
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
            ifs = oechem.oemolistream()
            if ifs.open(filePath):
                for tMol in ifs.GetOEGraphMols():
                    oeMol = oechem.OEGraphMol(tMol)
                    # if oechem.OEReadMolecule(ifs, oeMol):
                    if use3D:
                        mL.append(oemf.updateOePerceptions3D(oeMol, aromaticModel=oechem.OEAroModelOpenEye))
                    else:
                        mL.append(oemf.updateOePerceptions2D(oeMol, aromaticModel=oechem.OEAroModelOpenEye))
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
            for tMol in ifs.GetOEGraphMols():
                oeMol = oechem.OEGraphMol(tMol)
                if use3D:
                    mL.append(oemf.updateOePerceptions3D(oeMol, aromaticModel=oechem.OEAroModelOpenEye))
                else:
                    mL.append(oemf.updateOePerceptions2D(oeMol, aromaticModel=oechem.OEAroModelOpenEye))

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

    def createOeFingerPrintDatabase(self, oeMolDbFilePath, oeFpDbFilePath, fpType="TREE", dbType="FAST"):
        if dbType == "FAST":
            return self.__createOeFastFingerPrintDatabase(oeMolDbFilePath, oeFpDbFilePath, fpType=fpType)
        else:
            return True

    def __createOeFastFingerPrintDatabase(self, oeMolDbFilePath, oeFpDbFilePath, fpType="TREE"):
        """Create fast search fingerprint database from the input molecular database.

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
        startTime = time.time()
        ok = False
        try:
            _ = fpType
            fpD = {"TREE": oegraphsim.OEFPType_Tree, "CIRCULAR": oegraphsim.OEFPType_Circular, "PATH": oegraphsim.OEFPType_Path}
            myFpType = fpD[fpType] if fpType in fpD else oegraphsim.OEFPType_Tree
            opts = oegraphsim.OECreateFastFPDatabaseOptions(oegraphsim.OEGetFPType(myFpType))
            ok = oegraphsim.OECreateFastFPDatabaseFile(oeFpDbFilePath, oeMolDbFilePath, opts)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return ok

    def loadOeFingerPrintDatabase(self, oeMolDbFilePath, oeFpDbFilePath, inMemory=False, fpType="TREE", fpDbType="FAST"):
        if fpDbType == "FAST":
            return self.__loadOeFastFingerPrintDatabase(oeFpDbFilePath, inMemory=inMemory, fpType=fpType)
        else:
            return self.__loadOeFingerPrintDatabase(oeMolDbFilePath, fpType=fpType)

    def __loadOeFingerPrintDatabase(self, oeMolDbFilePath, fpType="TREE"):
        """Create conventional search fingerprint database from the input molecular database.

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
            OEFPType_MACCS166
            OEFPType_Lingo
        """
        fpDb = None
        ok = False
        startTime = time.time()
        try:
            fpD = {
                "TREE": oegraphsim.OEFPType_Tree,
                "CIRCULAR": oegraphsim.OEFPType_Circular,
                "PATH": oegraphsim.OEFPType_Path,
                "MACCS": oegraphsim.OEFPType_MACCS166,
                "LINGO": oegraphsim.OEFPType_Lingo,
            }
            fpType = fpType if fpType and fpType in fpD else "TREE"
            tag = "FP_" + fpType
            oeFpType = fpD[fpType] if fpType in fpD else oegraphsim.OEFPType_Tree
            oeMolDb = self.loadOeBinaryDatabaseAndIndex(oeMolDbFilePath)
            #
            fpDb = oegraphsim.OEFPDatabase(oeFpType)
            numMols = oeMolDb.GetMaxMolIdx()
            logger.debug("fpType %r tag %r oeFpType %r", fpType, tag, oeFpType)
            oeMol = oechem.OEGraphMol()
            for idx in range(0, numMols):
                if oeMolDb.GetMolecule(oeMol, idx):
                    if oeMol.HasData(tag):
                        tfp = oeMol.GetData(tag)
                        fpDb.AddFP(tfp)
                    else:
                        fpDb.AddFP(oeMol)
                else:
                    logger.info("Missing molecule at index %r", idx)

            numFp = fpDb.NumFingerPrints()
            ok = numMols == numFp
            logger.info("Loaded molecules  %d %s fingerprints %d (%.4f seconds)", numMols, fpType, numFp, time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            fpDb = None
        endTime = time.time()
        logger.debug("Completed with status %r operation at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return fpDb

    def __loadOeFastFingerPrintDatabase(self, oeFpDbFilePath, inMemory=False, fpType="TREE"):
        #
        _ = fpType
        startTime = time.time()
        if inMemory:
            memType = oegraphsim.OEFastFPDatabaseMemoryType_InMemory
        else:
            memType = oegraphsim.OEFastFPDatabaseMemoryType_MemoryMapped
        if not self.__mU.exists(oeFpDbFilePath):
            logger.error("Missing fingerprint database file %r", oeFpDbFilePath)
        fpDb = oegraphsim.OEFastFPDatabase(oeFpDbFilePath, memType)
        if not fpDb.IsValid():
            logger.error("Cannot open fingerprint database %r", oeFpDbFilePath)
        #
        lenFp = fpDb.NumFingerPrints()
        memTypeStr = fpDb.GetMemoryTypeString()
        endTime = time.time()
        logger.info("Read fingerprint database length %d loaded %s (%.4f seconds)", lenFp, memTypeStr, endTime - startTime)
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

    def buildOeBinaryMolCache(self, filePath, ccObjD, molBuildType="model-xyz", quietFlag=False, fpTypeList=None, limitPerceptions=False, suppressHydrogens=False):
        """Build cache of OEMol() objects from the input chemical component definition list.

        Args:
            filePath (str): output cache file path
            ccObjD (dict):  chemical component object dictionary
            molBuildType (str, optional): [description]. Defaults to "model-xyz".
            quietFlag (bool, optional): [description]. Defaults to False.
            fpTypeList (list, optional): fingerprint type list. Defaults to None.
            limitPerceptions (bool, optional): suppress automatic chemical perceptions. Defaults to False.
            suppressHydrogens (bool, optional): suppress explicit hydrogen count. Defaults to False.

        Returns:
            (int, int, list): chem comp success count, error count, chem comp identifier failure list

        """
        ok = False
        startTime = time.time()
        failIdList = []
        ccCount = 0
        errCount = 0
        try:
            ofs = oechem.oemolostream()
            ofs.SetFormat(oechem.OEFormat_OEB)
            if ofs.open(filePath):
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                for ccId, ccObj in ccObjD.items():
                    tId = oemf.setChemCompDef(ccObj)
                    if tId and tId == ccId:
                        ok = oemf.build(molBuildType=molBuildType, limitPerceptions=limitPerceptions)
                        if ok and fpTypeList:
                            fpOk = oemf.addFingerPrints(fpTypeList)
                            if not fpOk:
                                logger.info("Fingerprint generation fails for %r", ccId)
                        if ok:
                            oeMol = oemf.getMol(suppressHydrogens=suppressHydrogens)
                            oechem.OEWriteMolecule(ofs, oeMol)
                            ccCount += 1
                    if not ok or not tId:
                        # build failed incomplete component (e.g. missing atoms or bonds)
                        errCount += 1
                        failIdList.append(ccId)
            else:
                logger.error("Unable to open cache database %s", filePath)
                errCount += 1
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return ccCount, errCount, failIdList

    #
    def buildOeBinaryMolCacheFromIndex(self, filePath, ccIdxD, quietFlag=False, fpTypeList=None, limitPerceptions=False, suppressHydrogens=False):
        """Build cache of OEGraphMol() objects from the input chemical component search index.

        Args:
            filePath (str): output cache file path
            ccIdxD (dict): search index dictionary
            quietFlag (bool, optional): suppress OE output. Defaults to False.
            fpTypeList (list, optional): list of fingerprint types. Defaults to None.
            limitPerceptions (bool, optional): suppress automatic chemical perceptions. Defaults to False.
            suppressHydrogens (bool, optional): suppress explicit hydrogen count. Defaults to False.

        Returns:
            (int, int, list): chem comp success count, error count, chem comp identifier failure list
        """
        failIdList = []
        ccCount = 0
        errCount = 0
        startTime = time.time()
        try:
            ofs = oechem.oemolostream()
            ofs.SetFormat(oechem.OEFormat_OEB)
            if ofs.open(filePath):
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                for searchCcId, ccIdx in ccIdxD.items():
                    oemf.setDescriptor(ccIdx["smiles"], "oe-iso-smiles", searchCcId)
                    ok = oemf.build(molBuildType="oe-iso-smiles", limitPerceptions=limitPerceptions)
                    if ok and fpTypeList:
                        fpOk = oemf.addFingerPrints(fpTypeList)
                        if not fpOk:
                            logger.info("Fingerprint generation fails for %r", searchCcId)
                    if ok:
                        if not suppressHydrogens:
                            oemf.addExplicitHydrogens()
                            oemf.setSimpleAtomNames()
                        oeMol = oemf.getMol(suppressHydrogens=suppressHydrogens)
                        oechem.OEWriteMolecule(ofs, oeMol)
                        ccCount += 1
                    if not ok:
                        # build failed incomplete component (e.g. missing atoms or bonds)
                        errCount += 1
                        failIdList.append(searchCcId)
            else:
                logger.error("Unable to open cache database %s", filePath)
                errCount += 1
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        endTime = time.time()
        logger.info("Completed operation at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        return ccCount, errCount, failIdList

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
            logger.info("Opened %r with %r molecules", oeSubSearchFilePath, ssDb.NumMolecules())
        except Exception as e:
            logger.exception("Loading %r failing with %s", oeSubSearchFilePath, str(e))
        return ssDb

    def write(self, filePath, oeMol, constantMol=False, addSdTags=True):
        """Write an oeMol with format type inferred from the filePath extension (e.g. .mol)

        Args:
            filePath (str): filepath with a chemical type extension
            constantMol (bool, optional): copies molecule before performing format specific perceptions

        Returns:
            bool: True for success or False otherwise
        """
        try:
            molId = os.path.splitext(os.path.basename(filePath))[0]
            fmt = os.path.splitext(os.path.basename(filePath))[1][1:].lower()
            #
            if addSdTags:
                oemf = OeMoleculeFactory()
                oemf.setOeMol(oeMol, molId)
                oemf.addSdTags()
                oeMol = oemf.getMol()
            #
            self.__mU.mkdir(os.path.dirname(filePath))
            ofs = oechem.oemolostream()
            ofs.open(filePath)
            logger.debug("Writing (fmt=%s) molId %s path %s title %s", fmt, molId, filePath, oeMol.GetTitle())
            #
            if constantMol:
                oechem.OEWriteConstMolecule(ofs, oeMol)
            else:
                oechem.OEWriteMolecule(ofs, oeMol)
            #
            # If this is a mol2 file, we need to replace the resname
            if fmt.startswith("mol2"):
                # If this is a mol2/mol2h substitute the default substructure id
                with open(filePath, "r") as ifh:
                    lines = ifh.readlines()
                lines = [line.replace("<0>", molId) for line in lines]
                with open(filePath, "w") as ofh:
                    ofh.writelines(lines)
            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", filePath, str(e))
        return False

    def serializeOe(self, oeMol):
        """Create a string representing the content of the current OE molecule.   This
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
        """Reconstruct an OE molecule from the input string serialization (OE binary).

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
