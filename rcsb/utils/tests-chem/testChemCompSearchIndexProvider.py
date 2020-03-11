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
        """ Test search index constructure for an abbreviated chemical component resource file.
        """
        self.__testBuildSearchIndexCacheFiles(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, ccFileNamePrefix="cc-abbrev")

    def testChemCompSearchIndexCacheFilesAbbrevMp(self):
        """ Test search index constructure for an abbreviated chemical component resource file (multi).
            Creates index length - 2837  411s on numproc 4 on macbookpro .383 GB resident mem
        """
        self.__testBuildSearchIndexCacheFiles(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, ccFileNamePrefix="cc-abbrev", numProc=4)

    @unittest.skipIf(skipFlag, "Long test")
    def testChemCompSearchIndexCacheFilesFullMp(self):
        """ Test search index construction of full chemical component resource files.

            461 s (numproc=12 json format) 7.1 GB resident mem linux len=121079
        """
        self.__testBuildSearchIndexCacheFiles(ccFileNamePrefix="cc-full", numProc=4)

    @unittest.skipIf(skipFlag, "Long test")
    def testChemCompSearchIndexCacheFilesFiltered(self):
        """ Test search index construction of a filtered subset of chemical component definitions.
        """
        self.__testBuildSearchIndexCacheFiles(ccFileNamePrefix="cc-filtered")

    def __testBuildSearchIndexCacheFiles(self, **kwargs):
        """ Test build search index chemical component cache files from the input component dictionaries
        """
        molLimit = kwargs.get("molLimit", None)
        useCache = kwargs.get("useCache", False)
        logSizes = kwargs.get("logSizes", False)
        limitPerceptions = kwargs.get("limitPerceptions", True)
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


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFiltered"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesAbbrev"))
    suiteSelect.addTest(ChemCompSearchIndexProviderTests("testChemCompSearchIndexCacheFilesFullMp"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
