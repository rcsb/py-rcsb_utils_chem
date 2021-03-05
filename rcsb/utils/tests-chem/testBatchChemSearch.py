##
# File:    BatchChemSearchTests.py
# Author:  J. Westbrook
# Date:    3-Mar-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for batch search mode.

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
from rcsb.utils.chem.BatchChemSearch import BatchChemSearch


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class BatchChemSearchTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__smilesCod = os.path.join(self.__dataPath, "allcod.smi")
        #
        self.__queryResultsPath = os.path.join(self.__cachePath, "batch-descriptor-results.json")
        self.__testFlagFull = False
        #
        os.environ["CHEM_SEARCH_CACHE_PATH"] = os.path.join(self.__cachePath)
        os.environ["CHEM_SEARCH_CC_PREFIX"] = "cc-full" if self.__testFlagFull else "cc-abbrev"
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif") if not self.__testFlagFull else None
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif") if not self.__testFlagFull else None
        #
        self.__numMolsTest = 4000
        self.__numProc = 4
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBatchChemSearchSetup(self):
        """Test case - COD SMILES search """
        try:
            bsw = BatchChemSearch(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, cachePath=self.__cachePath, numProc=self.__numProc)
            ok = bsw.testCache()
            self.assertTrue(ok)
            smiL = bsw.fetchDescriptorList(self.__smilesCod, swap=True)
            logger.info("Query length (%d)", len(smiL))
            #
            smiL = bsw.splitSmiles(smiL)
            retL = bsw.doQuery(smiL[: self.__numMolsTest], "SMILES", matchOpts="graph-exact")
            logger.info("Result length (%d)", len(retL))
            self.assertGreater(len(retL), 2)
            for ii, ret in enumerate(retL, 1):
                logger.debug("%5d %8s %4s (%.3f) %s: %s", ii, ret.queryId, ret.ccId, ret.fpScore, ret.queryType, ret.query)
            #
            ok = bsw.storeMatchList(self.__queryResultsPath, retL)
            self.assertTrue(ok)
            mL = bsw.fetchMatchList(self.__queryResultsPath)
            self.assertEqual(len(mL), len(retL))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def batchChemSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(BatchChemSearchTests("testBatchChemSearchSetup"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = batchChemSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
