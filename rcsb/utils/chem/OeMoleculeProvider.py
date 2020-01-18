##
# File:    OeMoleculeProvider.py
# Author:  J. Westbrook
# Date:    24-Oct-2019
#
# Updates:
#  28-Oct-2019 jdw incorporate all of the public Bird cc definitions
##
"""
Utilities to read and process the dictionary of PDB chemical component definitions
and deliver OE molecule data.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import time
from collections import defaultdict

from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompDescriptorIt
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompIt
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompAtomIt
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class OeMoleculeProvider(object):
    """Utilities to read and process the dictionary of PDB chemical component definitions
    and deliver OE molecule data.
    """

    def __init__(self, **kwargs):
        self.__ccUrlTarget = kwargs.get("ccUrlTarget", "http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz")
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", "http://ftp.wwpdb.org/pub/pdb/data/bird/prd/prdcc-all.cif.gz")
        #
        self.__dirPath = kwargs.get("dirPath", ".")
        useCache = kwargs.get("useCache", True)
        numProc = kwargs.get("numProc", 4)
        self.__molLimit = kwargs.get("molLimit", None)
        oeMolFileName = kwargs.get("oeMolFileName", "oe-mol-components.oeb")
        oeMolDbFileName = kwargs.get("oeMolDbFileName", "oe-mol-db-components.oeb")
        #
        fpTypeList = kwargs.get("fpTypeList", ["TREE"])
        screenTypeList = kwargs.get("screenTypeList", ["SMARTS"])
        # screenTypeList = kwargs.get("screenTypeList", [])
        #
        ccIdxFileName = kwargs.get("ccIdxFileName", "cc-idx-components.pic")
        coordType = kwargs.get("coordType", None)
        #
        self.__oeMolDbFilePath = os.path.join(self.__dirPath, oeMolDbFileName)
        self.__ccIdxFilePath = os.path.join(self.__dirPath, ccIdxFileName)
        self.__fpDb = None
        self.__ssDb = None
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oeMolD = self.__reload(
            self.__ccUrlTarget,
            self.__birdUrlTarget,
            self.__dirPath,
            oeMolFileName,
            oeMolDbFileName,
            ccIdxFileName,
            coordType,
            fpTypeList,
            screenTypeList,
            useCache=useCache,
            numProc=numProc,
            molLimit=self.__molLimit,
        )
        #
        self.__oeMolDb = None
        self.__oeMolDbTitleD = None
        self.__ccIdxD = None

    def testCache(self, minCount=29000):
        num = self.__molLimit if self.__molLimit else minCount
        return self.__oeMolD and len(self.__oeMolD) >= num

    def getSubSearchDb(self, screenType="SMARTS", numProc=1):
        _ = screenType
        oeIo = OeIoUtils()
        if not self.__ssDb:
            fp = os.path.join(self.__dirPath, self.__getSubSearchFileName(screenType))
            self.__ssDb = oeIo.loadOeSubSearchDatabase(fp, screenType, numProc=numProc)
        return self.__ssDb

    def getFingerPrintDb(self, fpType="TREE"):
        oeIo = OeIoUtils()
        if not self.__fpDb:
            fp = os.path.join(self.__dirPath, self.__getFpDbFileName(fpType))
            self.__fpDb = oeIo.loadOeFingerPrintDatabase(fp, inMemory=True, fpType=fpType)
        return self.__fpDb

    def __getOeMolDbTitleIndex(self):
        oeMolDbTitleD = {}
        try:
            for idx in range(self.__oeMolDb.GetMaxMolIdx()):
                oeMolDbTitleD[self.__oeMolDb.GetTitle(idx)] = idx
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oeMolDbTitleD

    def getOeMolDatabase(self):
        if not self.__oeMolDb:
            oeIo = OeIoUtils()
            self.__oeMolDb = oeIo.loadOeBinaryDatabaseAndIndex(self.__oeMolDbFilePath)
            self.__oeMolDbTitleD = self.__getOeMolDbTitleIndex()
        #

        return self.__oeMolDb, self.__oeMolDbTitleD

    def getOeMolD(self):
        return self.__oeMolD

    def getMol(self, ccId):
        try:
            return self.__oeMolD[ccId]
        except Exception as e:
            logger.exception("Get molecule %r failing with %s", ccId, str(e))
        return None

    def getChemCompIdx(self):
        if not self.__ccIdxD:
            startTime = time.time()
            self.__ccIdxD = self.__mU.doImport(self.__ccIdxFilePath, fmt="pickle")
            endTime = time.time()
            logger.info("Loading %s with %d indexed records (%.4f seconds)", self.__ccIdxFilePath, len(self.__ccIdxD), endTime - startTime)
        #
        return self.__ccIdxD

    def __getFpDbFileName(self, fpType):
        return "oe-fp-database-%s.fpbin" % fpType

    def __getSubSearchFileName(self, screenType):
        return "oe-ss-database-%s.oeb" % screenType

    def __reload(
        self, ccUrlTarget, birdUrlTarget, dirPath, oeMolFileName, oeMolDbFileName, ccIdxFileName, coordType, fpTypeList, screenTypeList, useCache=True, numProc=1, molLimit=None
    ):
        """Reload input dictionary of chemical components and generate ...

        Args:
            ccUrlTarget (str): target url for chemical component dictionary resource file
            birdUrlTarget (str): target url for bird dictionary resource file (cc format)
            dirPath (str): path to the directory containing cache files
            #
            oeMolFileName (str): OE binary data file
            #
            eomolDbFileName (str): oeMol binary database file name
            ccIdxFileName (str): index file name
            coordType (str):  coordinates to use in building OE molecules from CIF components (model, ideal or None)
            #
            fpTypeList (list): fingerprint type (TREE, CIRCULAR, PATH)
            screenTypeList (list): fast sub search screen type (MOLECULE, SMARTS, MDL )
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict): something
        """
        reuseRemote = True
        retD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(ccUrlTarget)
        ccdFilePath = os.path.join(dirPath, fn)
        fn = fU.getFileName(birdUrlTarget)
        birdFilePath = os.path.join(dirPath, fn)
        #
        oeMolFilePath = os.path.join(dirPath, oeMolFileName)
        oeMolDbFilePath = os.path.join(dirPath, oeMolDbFileName)
        #
        fpPathD = {}
        for fpType in fpTypeList:
            fpPathD[fpType] = os.path.join(dirPath, self.__getFpDbFileName(fpType))

        screenPathD = {}
        for screenType in screenTypeList:
            screenPathD[screenType] = os.path.join(dirPath, self.__getSubSearchFileName(screenType))
        #
        ccIdxFilePath = os.path.join(dirPath, ccIdxFileName)
        #
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            pL = [oeMolFilePath, oeMolDbFilePath, ccIdxFilePath] + list(fpPathD.values())
            if not reuseRemote:
                pL.extend([ccdFilePath, birdFilePath])
            #
            logger.info("Clearing cache files %r", pL)
            for fp in pL:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(oeMolFilePath):
            oeIo = OeIoUtils()
            retD = oeIo.readOeBinaryMolCache(oeMolFilePath)
            logger.info("Loading OE binary molecule cache length %d", len(retD))
        else:
            ok = ok1 = ok2 = True
            if not (reuseRemote and fU.exists(ccdFilePath)):
                startTime = time.time()
                ok1 = fU.get(ccUrlTarget, ccdFilePath)
                endTime = time.time()
                if ok1:
                    logger.info("Fetched url %s for resource file %s (status = %r) (%.4f seconds)", ccUrlTarget, ccdFilePath, ok, endTime - startTime)
                else:
                    logger.error("Failing fetch of url %s for resource file %s (status = %r) (%.4f seconds)", ccUrlTarget, ccdFilePath, ok, endTime - startTime)
            if not (reuseRemote and fU.exists(birdFilePath)):
                startTime = time.time()
                ok2 = fU.get(birdUrlTarget, birdFilePath)
                endTime = time.time()
                if ok2:
                    logger.info("Fetched %s for resource file %s (status = %r) (%.4f seconds)", birdUrlTarget, birdFilePath, ok, endTime - startTime)
                else:
                    logger.error("Failing fetch for %s for resource file %s (status = %r) (%.4f seconds)", birdUrlTarget, birdFilePath, ok, endTime - startTime)

            if not (ok1 and ok2):
                logger.error("Cache rebuild failing")
            else:
                rdCcObjL = self.__getComponentDefinitions(ccdFilePath, birdFilePath)
                # Apply molecule limit
                ccObjL = rdCcObjL[:molLimit] if molLimit else rdCcObjL
                if molLimit:
                    logger.info("Using limited molecule subset count %d", molLimit)
                # -------
                startTime = time.time()
                oeIo = OeIoUtils()
                ccCount, errCount = oeIo.buildOeBinaryMolCache(oeMolFilePath, ccObjL, coordType=coordType, quietFlag=False)
                logger.info("Stored %d OeMols created with coordType %r (unconverted %d)", ccCount, coordType, errCount)
                retD = oeIo.readOeBinaryMolCache(oeMolFilePath)
                endTime = time.time()
                logger.info("Constructed %d cached OeMols (%.4f seconds)", len(retD), endTime - startTime)
                # --------
                startTime = time.time()
                molCount = oeIo.createOeBinaryDatabaseAndIndex(oeMolFilePath, oeMolDbFilePath)
                endTime = time.time()
                logger.info("Created and stored %d indexed OeMols in OE database format (%.4f seconds)", molCount, endTime - startTime)
                # --------
                for fpType in fpTypeList:
                    startTime = time.time()
                    ok = oeIo.createOeFingerPrintDatabase(oeMolDbFilePath, fpPathD[fpType], fpType=fpType)
                    endTime = time.time()
                    logger.info("Created and stored %s fingerprint database (%.4f seconds)", fpType, endTime - startTime)
                # ---------
                for screenType in screenTypeList:
                    startTime = time.time()
                    ok = oeIo.createOeSubSearchDatabase(oeMolFilePath, screenPathD[screenType], screenType=screenType, numProc=numProc)
                    endTime = time.time()
                    logger.info("Constructed substructure search database status %r with screenType %s (%.4f seconds)", ok, screenType, endTime - startTime)

                # ---------
                startTime = time.time()
                rD = self.__buildChemCompIndex(ccObjL)
                ok = self.__mU.doExport(ccIdxFilePath, rD, fmt="pickle")
                endTime = time.time()
                logger.info("Storing %s with %d raw indexed definitions (%.4f seconds)", ccIdxFilePath, len(rD), endTime - startTime)
        #
        logger.info("Reload OE binary molecule cache length %d", len(retD))
        return retD

    def getComponentDefinitions(self):
        fU = FileUtil()
        fn = fU.getFileName(self.__ccUrlTarget)
        ccdFilePath = os.path.join(self.__dirPath, fn)
        fn = fU.getFileName(self.__birdUrlTarget)
        birdFilePath = os.path.join(self.__dirPath, fn)
        return self.__getComponentDefinitions(ccdFilePath, birdFilePath)

    def __getComponentDefinitions(self, ccdFilePath, birdFilePath):
        startTime = time.time()
        logger.info("Reading %s", ccdFilePath)
        rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
        endTime = time.time()
        logger.info("Read %s with %d CCD definitions (%.4f seconds)", ccdFilePath, len(rdCcObjL), endTime - startTime)
        # -------
        startTime = time.time()
        logger.info("Reading %s", birdFilePath)
        birdCcObjL = self.__mU.doImport(birdFilePath, fmt="mmcif")
        endTime = time.time()
        logger.info("Read %s with %d BIRD definitions (%.4f seconds)", birdFilePath, len(birdCcObjL), endTime - startTime)
        rdCcObjL.extend(birdCcObjL)
        return rdCcObjL

    def __buildChemCompIndex(self, cL):
        """Internal method return a dictionary of extracted chemical component descriptors and formula.
        """
        rD = {}
        try:
            for dataContainer in cL:
                ccIt = PdbxChemCompIt(dataContainer)
                for cc in ccIt:
                    ccId = cc.getId()
                    formula = str(cc.getFormula()).replace(" ", "")
                    ambiguousFlag = cc.getAmbiguousFlag().upper() in ["Y", "YES"]
                    tch = cc.getFormalCharge()
                    fcharge = int(tch) if tch and tch not in [".", "?"] else 0
                #
                logger.debug("ccId %r formula %r ambiguous %r fcharge %r", ccId, formula, ambiguousFlag, fcharge)
                if fcharge:
                    sign = "+" if fcharge > 0 else "-"
                    mag = str(abs(fcharge)) if abs(fcharge) > 1 else ""
                    formula = formula + sign + mag
                #
                desIt = PdbxChemCompDescriptorIt(dataContainer)
                for des in desIt:
                    desType = des.getType().upper()
                    desProg = des.getProgram().upper()
                    if "OPEN" in desProg and desType == "SMILES_CANONICAL":
                        isoSmiles = des.getDescriptor()
                    elif "OPEN" in desProg and desType == "SMILES":
                        smiles = des.getDescriptor()
                    elif desType == "INCHI":
                        inchi = des.getDescriptor()
                    elif desType == "INCHIKEY":
                        inchiKey = des.getDescriptor()
                #
                atIt = PdbxChemCompAtomIt(dataContainer)
                typeCounts = defaultdict(int)
                for at in atIt:
                    aType = at.getType().upper()
                    typeCounts[aType] += 1
                #
                rD[ccId] = {
                    "FORMULA": formula,
                    "TYPE_COUNTS": typeCounts,
                    "AMBIGUOUS": ambiguousFlag,
                    "INCHI": inchi,
                    "INCHI_KEY": inchiKey,
                    "OE_ISO_SMILES": isoSmiles,
                    "OE_SMILES": smiles,
                }
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return rD
