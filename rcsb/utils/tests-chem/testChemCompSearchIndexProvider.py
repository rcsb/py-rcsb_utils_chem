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
    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f MB", rusageMax / 10 ** 6)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testChemCompSearchIndexCacheFilesAbbrev(self):
        """ Test search index constructure for an abbreviated chemical component resource file.
        """
        self.__testBuildSearchIndexCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-abbrev", molLimit=50)

    def testChemCompSearchIndexCacheFilesAbbrevMp(self):
        """ Test search index constructure for an abbreviated chemical component resource file (multi).
            Creates index length - 2837  411s on numproc 4 on macbookpro .383 GB resident mem
        """
        self.__testBuildSearchIndexCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-abbrev", molLimit=None, numProc=4)

    def testChemCompSearchIndexCacheFilesFull(self):
        """ Test search index construction of full chemical component resource files.
        """
        self.__testBuildSearchIndexCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-full")

    def testChemCompSearchIndexCacheFilesFiltered(self):
        """ Test search index construction of a filtered subset of chemical component definitions.
        """
        self.__testBuildSearchIndexCacheFiles(logSizes=True, useCache=False, ccFileNamePrefix="cc-filtered")

    def __testBuildSearchIndexCacheFiles(self, **kwargs):
        """ Test build search index chemical component cache files from the input component dictionaries
        """
        molLimit = kwargs.get("molLimit", None)
        useCache = kwargs.get("useCache", True)
        logSizes = kwargs.get("logSizes", False)
        limitPerceptions = kwargs.get("limitPerceptions", True)
        numProc = kwargs.get("numProc", 1)
        maxChunkSize = kwargs.get("maxChunkSize", 5)
        molLimit = kwargs.get("molLimit", None)
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        ccsiP = ChemCompSearchIndexProvider(
            cachePath=self.__cachePath,
            useCache=useCache,
            molLimit=molLimit,
            ccFileNamePrefix=ccFileNamePrefix,
            limitPerceptions=limitPerceptions,
            numProc=numProc,
            maxChunkSize=maxChunkSize,
        )
        ok = ccsiP.testCache(minCount=molLimit, logSizes=logSizes)
        self.assertTrue(ok)
        logger.info(" ******* Completed operation ******** ")
        #


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFull"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFiltered"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesAbbrev"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
