##
# File:    OeSearchIndexUtilsTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for search modes using source molecular definitions coming from a search index.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.chem import __version__
from rcsb.utils.chem.ChemCompSearchWrapper import ChemCompSearchWrapper
from rcsb.utils.io.MarshalUtil import MarshalUtil


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSearchIndexUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__testFlagFull = False
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        # Run the bootstrap configuration
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__testBootstrapConfig()
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __testBootstrapConfig(self):
        """Test read/write search configuration.
        """
        try:
            cfgD = {}
            if self.__testFlagFull:
                os.environ["CHEM_SEARCH_CACHE_PATH"] = os.path.join(self.__cachePath)
                os.environ["CHEM_SEARCH_CC_PREFIX"] = "cc-full"
                ccFileNamePrefix = "cc-full"
                configDirPath = os.path.join(self.__cachePath, "config")
                configFilePath = os.path.join(configDirPath, ccFileNamePrefix + "-config.json")
                oeFileNamePrefix = cfgD.get("oeFileNamePrefix", "oe-full")
                ccFileNamePrefix = cfgD.get("oeFileNamePrefix", ccFileNamePrefix)
                ccUrlTarget = cfgD.get("ccUrlTarget", None)
                birdUrlTarget = cfgD.get("birdUrlTarget", None)
            else:
                os.environ["CHEM_SEARCH_CACHE_PATH"] = os.path.join(self.__cachePath)
                os.environ["CHEM_SEARCH_CC_PREFIX"] = "cc-abbrev"
                ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
                birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
                ccFileNamePrefix = "cc-abbrev"
                configDirPath = os.path.join(self.__cachePath, "config")
                configFilePath = os.path.join(configDirPath, ccFileNamePrefix + "-config.json")
                oeFileNamePrefix = cfgD.get("oeFileNamePrefix", "oe-abbrev")
            #
            molLimit = cfgD.get("molLimit", None)
            useCache = cfgD.get("useCache", False)
            logSizes = cfgD.get("logSizes", False)
            #
            numProc = cfgD.get("numProc", 12)
            maxProc = os.cpu_count()
            numProc = min(numProc, maxProc)
            maxChunkSize = cfgD.get("maxChunkSize", 10)
            logger.info("+++ >>> Using MAXPROC %d", numProc)
            #
            limitPerceptions = cfgD.get("limitPerceptions", False)
            quietFlag = cfgD.get("quietFlag", True)
            #
            fpTypeCuttoffD = {"TREE": 0.6, "MACCS": 0.9, "PATH": 0.6, "CIRCULAR": 0.6, "LINGO": 0.9}
            buildTypeList = ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
            #
            oesmpKwargs = {
                "ccUrlTarget": ccUrlTarget,
                "birdUrlTarget": birdUrlTarget,
                "cachePath": self.__cachePath,
                "useCache": useCache,
                "ccFileNamePrefix": ccFileNamePrefix,
                "oeFileNamePrefix": oeFileNamePrefix,
                "limitPerceptions": limitPerceptions,
                "minCount": None,
                "maxFpResults": 50,
                "fpTypeCuttoffD": fpTypeCuttoffD,
                "buildTypeList": buildTypeList,
                "screenTypeList": None,
                "quietFlag": quietFlag,
                "numProc": numProc,
                "maxChunkSize": maxChunkSize,
                "molLimit": molLimit,
                "logSizes": logSizes,
            }
            ccsiKwargs = {
                "ccUrlTarget": ccUrlTarget,
                "birdUrlTarget": birdUrlTarget,
                "cachePath": self.__cachePath,
                "useCache": useCache,
                "ccFileNamePrefix": ccFileNamePrefix,
                "oeFileNamePrefix": oeFileNamePrefix,
                "limitPerceptions": limitPerceptions,
                "minCount": None,
                "numProc": numProc,
                "quietFlag": quietFlag,
                "maxChunkSize": maxChunkSize,
                "molLimit": None,
                "logSizes": False,
            }
            configD = {"versionNumber": 0.20, "ccsiKwargs": ccsiKwargs, "oesmpKwargs": oesmpKwargs}
            self.__mU.mkdir(configDirPath)
            self.__mU.doExport(configFilePath, configD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAReadConfig(self):
        """Test read/access configuration
        """
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBUpdateChemCompIndex(self):
        """Test update chemical component/Bird basic index.
        """
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateChemCompIndex()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCUpdateSearchIndex(self):
        """Test update search index.
        """
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateSearchIndex()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testDUpdateSearchMoleculeProvider(self):
        """Test update of the search molecule provider.
        """
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateSearchMoleculeProvider()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testEReloadSearchDatabase(self):
        """Test reload search databases.
        """
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.reloadSearchDatabase()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testZoomMatchDescriptor(self):
        """Test descriptor matching
        """
        try:
            numMolsTest = 500
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateChemCompIndex(useCache=True)
            self.assertTrue(ok)
            ccIdx = ccsw.getChemCompIndex()
            ok = ccsw.reloadSearchDatabase()
            self.assertTrue(ok)
            #
            logger.debug("ccIdx (%d) keys %r entry %r", len(ccIdx), list(ccIdx.keys())[:10], ccIdx["000"])
            #
            logger.info("Dependencies loaded - Starting search test scan of (limit=%r)", numMolsTest)
            for ii, (ccId, ccD) in enumerate(ccIdx.items(), 1):
                if numMolsTest and ii > numMolsTest:
                    break
                for buildType in self.__buildTypeList:
                    if buildType in ccD:
                        startTime = time.time()
                        retStatus, ssL, fpL = ccsw.matchByDescriptor(ccD[buildType], buildType)
                        mOk = self.__resultContains(ccId, ssL)
                        self.assertTrue(mOk)
                        #
                        ssCcIdList = list(set([t.ccId.split("|")[0] for t in ssL]))
                        fpCcIdList = list(set([t.ccId.split("|")[0] for t in fpL]))
                        logger.info(
                            "%s (%d) for buildType %s (%d/%d) (%r) %r (%.4f secs)",
                            ccId,
                            ii,
                            buildType,
                            len(ssCcIdList),
                            len(fpCcIdList),
                            mOk and retStatus == 0,
                            ssCcIdList,
                            time.time() - startTime,
                        )
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testZoomMatchFormula(self):
        """Test formula matching
        """
        try:
            numMolsTest = 500
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateChemCompIndex(useCache=True)
            self.assertTrue(ok)
            ccIdx = ccsw.getChemCompIndex()
            #
            logger.debug("ccIdx (%d) keys %r entry %r", len(ccIdx), list(ccIdx.keys())[:10], ccIdx["000"])
            #
            logger.info("Dependencies loaded - Starting formula test scan of (limit=%r)", numMolsTest)
            for ii, (ccId, idxD) in enumerate(ccIdx.items(), 1):
                if numMolsTest and ii > numMolsTest:
                    break
                #
                startTime = time.time()
                elementRangeD = {el: {"min": eCount, "max": eCount} for el, eCount in idxD["type-counts"].items()}
                retStatus, rL = ccsw.matchByFormulaRange(elementRangeD, ccId)
                mOk = self.__resultContains(ccId, rL)
                self.assertTrue(mOk)
                logger.info(
                    "%s (%d) (%d) (%r) (%.4f secs)", ccId, ii, len(rL), mOk and retStatus == 0, time.time() - startTime,
                )
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def __resultContains(self, ccId, matchResultList):
        for matchResult in matchResultList:
            if ccId in matchResult.ccId:
                return True
        return False


def fullSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchIndexUtilsTests("testMatchSmiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fullSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
