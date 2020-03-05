##
# File:    ChemCompIndexProviderTests.py
# Author:  J. Westbrook
# Date:    17-Feb-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to read and process the dictionary of PDB chemical component definitions.

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
from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChemCompIndexProviderTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testChemCompIndexCacheFilesAbbrev(self):
        """ Test construction of full chemical component resource files.
        """
        self.__testBuildMoleculeCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-abbrev")

    def testChemCompIndexCacheFilesFull(self):
        """ Test construction of full chemical component resource files.
        """
        self.__testBuildMoleculeCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-full")

    def testChemCompIndexCacheFilesFiltered(self):
        """ Test construction of a filtered subset of chemical component definitions.
        """
        self.__testBuildMoleculeCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-filtered")

    def __testBuildMoleculeCacheFiles(self, **kwargs):
        """ Test build chemical component cache files from the input component dictionaries
        """
        molLimit = kwargs.get("molLimit", None)
        useCache = kwargs.get("useCache", True)
        logSizes = kwargs.get("logSizes", False)
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        ccmP = ChemCompIndexProvider(cachePath=self.__cachePath, useCache=useCache, molLimit=molLimit, ccFileNamePrefix=ccFileNamePrefix)
        ok = ccmP.testCache(minCount=molLimit, logSizes=logSizes)
        self.assertTrue(ok)
        logger.info(" ******* Completed operation ******** ")
        #


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompIndexProviderTests("testChemCompIndexCacheFilesFull"))
    suiteSelect.addTest(ChemCompIndexProviderTests("testChemCompIndexCacheFilesFiltered"))
    suiteSelect.addTest(ChemCompIndexProviderTests("testChemCompIndexCacheFilesAbbrev"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
