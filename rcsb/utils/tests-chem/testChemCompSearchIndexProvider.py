##
# File:    ChemCompSearchIndexProviderTests.py
# Author:  J. Westbrook
# Date:    3-Mar-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to generate a search index for reasonable tautomer and protomer
configurations for core PDB chemical component definitions.

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
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChemCompSearchIndexProviderTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testChemCompSearchIndexCacheFilesAbbrev(self):
        """Test search index constructure for an abbreviated chemical component resource file."""
        self.__testBuildSearchIndexCacheFiles(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, ccFileNamePrefix="cc-abbrev")

    def testChemCompSearchIndexCacheFilesAbbrevMp(self):
        """Test search index constructure for an abbreviated chemical component resource file (multi).
        Creates index length - 2837  411s on numproc 4 on macbookpro .383 GB resident mem
        """
        self.__testBuildSearchIndexCacheFiles(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, ccFileNamePrefix="cc-abbrev", numProc=4)

    @unittest.skipIf(skipFlag, "Long test")
    def testChemCompSearchIndexCacheFilesFullMp(self):
        """Test search index construction of full chemical component resource files.

        461 s (numproc=12 json format) 7.1 GB resident mem linux len=121079
        """
        self.__testBuildSearchIndexCacheFiles(ccFileNamePrefix="cc-full", numProc=4)

    @unittest.skipIf(skipFlag, "Long test")
    def testChemCompSearchIndexCacheFilesFiltered(self):
        """Test search index construction of a filtered subset of chemical component definitions."""
        self.__testBuildSearchIndexCacheFiles(ccFileNamePrefix="cc-filtered")

    def __testBuildSearchIndexCacheFiles(self, **kwargs):
        """Test build search index chemical component cache files from the input component dictionaries"""
        molLimit = kwargs.get("molLimit", None)
        useCache = kwargs.get("useCache", False)
        logSizes = kwargs.get("logSizes", False)
        limitPerceptions = kwargs.get("limitPerceptions", False)
        numProc = kwargs.get("numProc", 1)
        maxChunkSize = kwargs.get("maxChunkSize", 5)
        molLimit = kwargs.get("molLimit", None)
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        quietFlag = kwargs.get("quietFlag", True)
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        #
        ccsiP = ChemCompSearchIndexProvider(
            ccUrlTarget=ccUrlTarget,
            birdUrlTarget=birdUrlTarget,
            cachePath=self.__cachePath,
            useCache=useCache,
            molLimit=molLimit,
            ccFileNamePrefix=ccFileNamePrefix,
            limitPerceptions=limitPerceptions,
            numProc=numProc,
            maxChunkSize=maxChunkSize,
            quietFlag=quietFlag,
        )
        ok = ccsiP.testCache(minCount=molLimit, logSizes=logSizes)
        self.assertTrue(ok)
        logger.info(" ******* Completed operation ******** ")
        #
        return ccsiP

    def testFormulaMatch(self):
        """Test formula match   ..."""
        ccsidxP = self.__testBuildSearchIndexCacheFiles(
            ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, logSizes=False, useCache=True, ccFileNamePrefix="cc-abbrev"
        )
        ccidxD = ccsidxP.getIndex()
        logger.info("Matching formula for %d definitions", len(ccidxD))
        for ccId, idxD in ccidxD.items():
            startTime = time.time()
            fQueryD = {el: {"min": eCount, "max": eCount} for el, eCount in idxD["type-counts"].items()}
            if fQueryD:
                rL = ccsidxP.matchMolecularFormulaRange(fQueryD, matchSubset=False)
                logger.debug("%s formula matches %r (%.4f seconds)", ccId, rL, time.time() - startTime)
                ok = self.__resultContains(ccId, rL)
                if not ok:
                    logger.info("%s formula not matched %r %r  (%.4f seconds)", ccId, rL, fQueryD, time.time() - startTime)
                self.assertTrue(ok)

    def testFormulaSubsetMatchAbbrev(self):
        """Test formula range match on abbreviated data set  ..."""
        ccsidxP = self.__testBuildSearchIndexCacheFiles(
            ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, logSizes=False, useCache=True, ccFileNamePrefix="cc-abbrev"
        )
        startTime = time.time()
        fQueryD = {"C": {"min": 50, "max": 65}}
        rL = ccsidxP.matchMolecularFormulaRange(fQueryD, matchSubset=True)
        logger.info("C formula subset matches (%d) (%.4f seconds)", len(rL), time.time() - startTime)
        self.assertGreaterEqual(len(rL), 10)

    @unittest.skipIf(skipFlag, "Long test")
    def testFormulaSubsetMatchFull(self):
        """Test formula range match on the full index   ...
        Full subset search
        """
        ccsidxP = self.__testBuildSearchIndexCacheFiles(useCache=True, numProc=4, ccFileNamePrefix="cc-full")
        startTime = time.time()
        fQueryD = {"C": {"min": 10, "max": 65}, "O": {"min": 1, "max": 20}, "N": {"min": 2, "max": 12}}
        rL = ccsidxP.matchMolecularFormulaRange(fQueryD, matchSubset=True)
        logger.info("Formula subset matches (%d) (%.4f seconds)", len(rL), time.time() - startTime)
        self.assertGreaterEqual(len(rL), 10)

    @unittest.skipIf(skipFlag, "Long test")
    def testMinumumFormulaSubsetFull(self):
        """Test formula filter in minumum formula composition on the full index   ..."""
        ccsidxP = self.__testBuildSearchIndexCacheFiles(useCache=True, numProc=4, ccFileNamePrefix="cc-full")
        startTime = time.time()
        fQueryD = {"C": {"min": 10}, "O": {"min": 1}, "N": {"min": 2}}
        rL = ccsidxP.filterMinimumMolecularFormula(fQueryD)
        logger.info("Formula minimum filters (%d) (%.4f seconds)", len(rL), time.time() - startTime)
        self.assertGreaterEqual(len(rL), 10)

    def __resultContains(self, ccId, matchResultList):
        for matchResult in matchResultList:
            if ccId in matchResult.ccId:
                return True
        return False


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFiltered"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesAbbrev"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFullMp"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
