##
# File:    ChemCompSearchWrapperTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2020
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
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class ChemCompSearchWrapperTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        #
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        # Set the external environment for the wrapper class-
        self.__testFlagFull = False
        self.__testStash = False
        if self.__testFlagFull:
            os.environ["CHEM_SEARCH_CACHE_PATH"] = os.path.join(self.__cachePath)
            os.environ["CHEM_SEARCH_CC_PREFIX"] = "cc-full"
        else:
            os.environ["CHEM_SEARCH_CACHE_PATH"] = os.path.join(self.__cachePath)
            os.environ["CHEM_SEARCH_CC_PREFIX"] = "cc-abbrev"
        #
        self.__numMolsTest = 20
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAABuildConfiguration(self):
        """Test case - build configuration -"""
        try:
            ccsw = ChemCompSearchWrapper()
            ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif") if not self.__testFlagFull else None
            birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif") if not self.__testFlagFull else None
            ok = ccsw.setConfig(ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAABuildDependencies(self):
        """Test case - build all dependencies (convenience method)"""
        try:
            ccsw = ChemCompSearchWrapper()
            ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif") if not self.__testFlagFull else None
            birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif") if not self.__testFlagFull else None
            ok = ccsw.buildDependenices(ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skip("Private test")
    def testAABuildDependenciesAndStash(self):
        """Test case - build, stash and restore dependencies -"""
        try:
            ccsw = ChemCompSearchWrapper()
            ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif") if not self.__testFlagFull else None
            birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif") if not self.__testFlagFull else None
            ok = ccsw.buildDependenices(ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget)
            self.assertTrue(ok)
            #
            if self.__testStash:
                url = "sftp://bl-east.rcsb.org"
                userName = ""
                pw = ""
                dirPath = "4-coastal"
                ok = ccsw.stashDependencies(url, dirPath, userName=userName, pw=pw)
                self.assertTrue(ok)
                #
                fileU = FileUtil()
                fileU.remove(self.__cachePath)
                #
                url = "http://bl-east.rcsb.org"
                ok = ccsw.restoreDependencies(url, dirPath)
                #
                fileU.remove(self.__cachePath)
                #
                url = "sftp://bl-east.rcsb.org"
                ok = ccsw.restoreDependencies(url, dirPath, userName=userName, pw=pw)
                self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAReadConfig(self):
        """Test read/access configuration"""
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBUpdateChemCompIndex(self):
        """Test update chemical component/Bird basic index."""
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
        """Test update search index."""
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
        """Test update of the search molecule provider."""
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
        """Test reload search databases."""
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
        """Test descriptor matching"""
        try:
            numMolsTest = self.__numMolsTest
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
                        retStatus, ssL, fpL = ccsw.searchByDescriptor(ccD[buildType], buildType, matchOpts="graph-relaxed")
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

    def testZoomSubStructSearch(self):
        """Test substructure search"""
        try:
            numMolsTest = self.__numMolsTest
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateChemCompIndex(useCache=True)
            self.assertTrue(ok)
            ccIdx = ccsw.getChemCompIndex()
            ok = ccsw.updateSearchIndex(useCache=True)
            self.assertTrue(ok)
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
                        retStatus, ssL, _ = ccsw.searchByDescriptor(ccD[buildType], buildType, matchOpts="sub-struct-graph-relaxed")
                        if retStatus == -100:
                            logger.warning("Descriptor error continuing...")
                            continue
                        mOk = self.__resultContains(ccId, ssL)
                        self.assertTrue(mOk)
                        #
                        ssCcIdList = list(set([t.ccId.split("|")[0] for t in ssL]))
                        logger.info(
                            "%s (%d) for buildType %s (%d) (%r) %r (%.4f secs)",
                            ccId,
                            ii,
                            buildType,
                            len(ssCcIdList),
                            mOk and retStatus == 0,
                            ssCcIdList,
                            time.time() - startTime,
                        )
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testZoomMatchFormula(self):
        """Test formula matching"""
        try:
            numMolsTest = self.__numMolsTest
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
                    "%s (%d) (%d) (%r) (%.4f secs)",
                    ccId,
                    ii,
                    len(rL),
                    mOk and retStatus == 0,
                    time.time() - startTime,
                )
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skip("Private test")
    def testZoomFingerprintRepeat(self):
        """Test substructure search"""
        try:
            ccsw = ChemCompSearchWrapper()
            ok = ccsw.readConfig()
            self.assertTrue(ok)
            ok = ccsw.updateChemCompIndex(useCache=True)
            self.assertTrue(ok)
            ccIdx = ccsw.getChemCompIndex()
            ok = ccsw.updateSearchIndex(useCache=True)
            self.assertTrue(ok)
            ok = ccsw.reloadSearchDatabase()
            self.assertTrue(ok)
            #
            logger.debug("ccIdx (%d)", len(ccIdx))
            #
            fpL = []
            descr = "InChI=1S/C9H15N5O3/c1-3(15)6(16)4-2-11-7-5(12-4)8(17)14-9(10)13-7/h3-4,6,12,15-16H,2H2,1H3,(H4,10,11,13,14,17)/t3-,4-,6-/m1/s1"
            for ii in range(100):
                startTime = time.time()
                retStatus, _, fpL = ccsw.searchByDescriptor(descr, "InChI", matchOpts="fingerprint-similarity")
                if retStatus == -100:
                    logger.warning("Descriptor error continuing...")
                    continue
                rD = {}
                for mr in fpL:
                    ccId = mr.ccId.split("|")[0]
                    rD[ccId] = max(rD[ccId], mr.fpScore) if ccId in rD else mr.fpScore
                rTupL = sorted(rD.items(), key=lambda kv: kv[1], reverse=True)
                rL = [rTup[0] for rTup in rTupL]
                scoreL = [rTup[1] for rTup in rTupL]
                logger.info("%4d (%3d) %r (%.4f secs)", ii, len(rL), retStatus == 0, time.time() - startTime)
            logger.info("rL %r", rL)
            logger.info("scoreL %r", scoreL)
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
    suiteSelect.addTest(ChemCompSearchWrapperTests("testMatchSmiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fullSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
