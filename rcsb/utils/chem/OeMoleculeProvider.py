##
# File:    OeMoleculeProvider.py
# Author:  J. Westbrook
# Date:    24-Oct-2019
#
# Updates:
#  28-Oct-2019 jdw incorporate all of the public Bird cc definitions
##
"""
Utilities deliver OE molecule data for PDB chemical component definitions
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import time

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class OeMoleculeProvider(object):
    """Utilities build and deliver OE molecule databases from PDB chemical component definition data"""

    def __init__(self, **kwargs):
        """Utilities build and deliver OE molecule databases from PDB chemical component definition data
        Args:
            cachePath (str, optional): path to the directory containing cache files (default: '.')
            molBuildType (str,optional): data source for building OE molecules (default: "model-xyz")
            oeFileNamePrefix (str, optional) file name prefix for all generated databases (default: "oe")

        """
        # Database file names with be prefixed with base prefix plus the molecular build type and perception options
        oeFileNamePrefixBase = kwargs.get("oeFileNamePrefix", "oe")
        limitPerceptions = kwargs.get("limitPerceptions", False)
        molBuildType = kwargs.get("molBuildType", "model-xyz")
        if limitPerceptions and molBuildType in ["oe-smiles", "oe-iso-smiles", "inchi"]:
            self.__oeFileNamePrefix = oeFileNamePrefixBase + "-" + molBuildType + "-limit"
        else:
            self.__oeFileNamePrefix = oeFileNamePrefixBase + "-" + molBuildType
        #
        cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(cachePath, "oe_mol")
        #
        self.__fpDbD = {}
        self.__ssDb = None
        self.__oeMolD = {}
        self.__oeMolDb = None
        self.__oeMolDbTitleD = None
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__molCount = self.__reload(**kwargs)

    def testCache(self):
        return self.__mU.exists(os.path.join(self.__dirPath, self.__getOeMolFileName())) and self.__mU.exists(os.path.join(self.__dirPath, self.__getOeMolDbFileName()))

    def getSubSearchDb(self, screenType="SMARTS", numProc=1, forceRefresh=False):
        if not self.__ssDb or forceRefresh:
            oeIo = OeIoUtils()
            fp = os.path.join(self.__dirPath, self.__getSubSearchFileName(screenType))
            logger.info("Opening screened substructure search database %r", fp)
            self.__ssDb = oeIo.loadOeSubSearchDatabase(fp, screenType, numProc=numProc)
        return self.__ssDb

    def getFingerPrintDb(self, fpType, fpDbType="STANDARD", rebuild=False):
        if fpType not in self.__fpDbD or rebuild:
            oeIo = OeIoUtils()
            fastFpDbPath = os.path.join(self.__dirPath, self.__getFastFpDbFileName(fpType))
            oeMolDbFilePath = os.path.join(self.__dirPath, self.__getOeMolDbFileName())
            fpDb = oeIo.loadOeFingerPrintDatabase(oeMolDbFilePath, fastFpDbPath, inMemory=True, fpType=fpType, fpDbType=fpDbType)
            if fpDb:
                self.__fpDbD[fpType] = fpDb
        #
        return self.__fpDbD[fpType]

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
            self.__oeMolDb = oeIo.loadOeBinaryDatabaseAndIndex(os.path.join(self.__dirPath, self.__getOeMolDbFileName()))
            self.__oeMolDbTitleD = self.__getOeMolDbTitleIndex()
        return self.__oeMolDb, self.__oeMolDbTitleD

    def getOeMolD(self):
        try:
            if not self.__oeMolD:
                oeIo = OeIoUtils()
                self.__oeMolD = oeIo.readOeBinaryMolCache(os.path.join(self.__dirPath, self.__getOeMolFileName()))
                logger.info("Loading OE binary molecule cache length %d", len(self.__oeMolD))
            return self.__oeMolD
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def getMol(self, ccId):
        try:
            if not self.__oeMolD:
                oeIo = OeIoUtils()
                self.__oeMolD = oeIo.readOeBinaryMolCache(os.path.join(self.__dirPath, self.__getOeMolFileName()))
                logger.info("Loading OE binary molecule cache length %d", len(self.__oeMolD))
            return self.__oeMolD[ccId]
        except Exception as e:
            logger.exception("Get molecule %r failing with %s", ccId, str(e))
        return None

    def __getFastFpDbFileName(self, fpType):
        return "%s-fast-fp-database-%s.fpbin" % (self.__oeFileNamePrefix, fpType)

    def __getSubSearchFileName(self, screenType):
        return "%s-ss-database-%s.oeb" % (self.__oeFileNamePrefix, screenType)

    def __getOeMolDbFileName(self):
        return "%s-mol-db-components.oeb" % self.__oeFileNamePrefix

    def __getOeMolFileName(self):
        return "%s-mol-components.oeb" % self.__oeFileNamePrefix

    def __reload(self, **kwargs):
        """Reload the dictionary of OE molecules and related data artifacts for chemical component definitions.

        Args:
            molBuildType (str):  coordinates to use in building OE molecules from CIF components (model, ideal or None)
            limitPerceptions(bool): process input descriptors in essentially verbatim mode (default: True)
            fpTypeList (list): fingerprint type (TREE,PATH,MACCS,CIRCULAR,LINGO)
            screenTypeList (list): fast sub search screen type (MOLECULE, SMARTS, MDL, ... )
            useCache (bool, optional): flag to use cached files. Defaults to True.
            cachePath (str): path to the top cache directory. Defaults to '.'.
            numProc (int): number processors to engage in screen substructure search database generation.
            molLimit (int, optional): limiting number of molecules in data store (default: 0 no limit)
            suppressHydrogens (bool, optional): flag to suppress explicit hydrogens in the OE data store.

        Returns:
            (dict): dictionary of constructed OE molecules

        """
        useCache = kwargs.get("useCache", True)
        cachePath = kwargs.get("cachePath", ".")
        numProc = kwargs.get("numProc", 2)
        molLimit = kwargs.get("molLimit", 0)
        fpTypeList = kwargs.get("fpTypeList", ["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"])
        # screenTypeList = kwargs.get("screenTypeList", ["SMARTS"])
        screenTypeList = kwargs.get("screenTypeList", [])
        molBuildType = kwargs.get("molBuildType", "model-xyz")
        limitPerceptions = kwargs.get("limitPerceptions", False)
        quietFlag = kwargs.get("quietFlag", True)
        suppressHydrogens = kwargs.get("suppressHydrogens", False)
        logSizes = kwargs.get("logSizes", False)
        fpDbType = "STANDARD"
        #
        ccCount = 0
        oeCount = 0
        errCount = 0
        failIdList = []
        oeIo = OeIoUtils(quietFlag=quietFlag)
        # --------
        oeMolFilePath = os.path.join(self.__dirPath, self.__getOeMolFileName())
        if not useCache or (useCache and not self.__mU.exists(oeMolFilePath)):
            cmpKwargs = {k: v for k, v in kwargs.items() if k not in ["cachePath", "useCache", "molLimit"]}
            ccmP = ChemCompMoleculeProvider(cachePath=cachePath, useCache=True, molLimit=molLimit, **cmpKwargs)
            ok = ccmP.testCache(minCount=molLimit, logSizes=logSizes)
            ccObjD = ccmP.getMolD() if ok else {}
            ccCount = len(ccObjD)
            # -------
            startTime = time.time()
            oeCount, errCount, failIdList = oeIo.buildOeBinaryMolCache(
                oeMolFilePath, ccObjD, molBuildType=molBuildType, quietFlag=quietFlag, fpTypeList=fpTypeList, limitPerceptions=limitPerceptions, suppressHydrogens=suppressHydrogens
            )
            logger.info("Stored %d/%d OeMols (suppressH = %r) created with molBuildType %r (unconverted %d)", oeCount, ccCount, suppressHydrogens, molBuildType, errCount)
            if failIdList:
                logger.info("%r failures %r", molBuildType, failIdList)
            endTime = time.time()
            logger.info("Constructed %d/%d cached oeMols (%.4f seconds)", oeCount, ccCount, endTime - startTime)
        # --------
        oeMolDbFilePath = os.path.join(self.__dirPath, self.__getOeMolDbFileName())
        if not useCache or (useCache and not self.__mU.exists(oeMolDbFilePath)):
            startTime = time.time()
            molCount = oeIo.createOeBinaryDatabaseAndIndex(oeMolFilePath, oeMolDbFilePath)
            endTime = time.time()
            logger.info("Created and stored %d indexed OeMols in OE database format (%.4f seconds)", molCount, endTime - startTime)

        # --------
        if fpDbType == "FAST":
            for fpType in fpTypeList:
                startTime = time.time()
                #  Fast FP search database file names
                fpPath = os.path.join(self.__dirPath, self.__getFastFpDbFileName(fpType))
                if not useCache or (useCache and not self.__mU.exists(fpPath)):
                    ok = oeIo.createOeFingerPrintDatabase(oeMolDbFilePath, fpPath, fpType=fpType)
                    endTime = time.time()
                    logger.info("Created and stored %s fingerprint database (%.4f seconds)", fpType, endTime - startTime)
        # --------
        if molBuildType in ["oe-iso-smiles"]:
            for screenType in screenTypeList:
                startTime = time.time()
                fp = os.path.join(self.__dirPath, self.__getSubSearchFileName(screenType))
                if not useCache or (useCache and not self.__mU.exists(fp)):
                    ok = oeIo.createOeSubSearchDatabase(oeMolFilePath, fp, screenType=screenType, numProc=numProc)
                    endTime = time.time()
                    logger.info("Constructed screened substructure database (status %r) with screenType %s (%.4f seconds)", ok, screenType, endTime - startTime)
                    # ---------
                    ssDb = oeIo.loadOeSubSearchDatabase(fp, screenType=screenType, numProc=numProc)
                    ok = ssDb.NumMolecules() == oeCount
                    # ----------
        return oeCount
