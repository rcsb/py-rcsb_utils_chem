##
# File:    ChemCompMoleculeProviderTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#  5-Mar-2020 jdw all cases including filtered subsets working with naming prefixes.
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
import resource
import time
import unittest


from rcsb.utils.chem import __version__
from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChemCompMoleculeProviderTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        self.__missedIdsPath = os.path.join(self.__dataPath, "filteredIds.json")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f MB", rusageMax / 10 ** 6)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildMoleculeCacheFilesAbbrev(self):
        """Test construction of abbreviated of chemical component definitions from alternate source paths
        """
        self.__testBuildMoleculeCacheFiles(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, ccFileNamePrefix="cc-abbrev", useCache=False, molLimit=500)

    def testBuildMoleculeCacheFilesSubset(self):
        """ Test construction of a subset (1K) of chemical component definitions.
        """
        self.__testBuildMoleculeCacheFiles(molLimit=1000, useCache=False)
        #

    def testSubsetBuildMoleculeCacheFiltered(self):
        """ Test construction of a filtered selection of chemical component definitions.
        """
        mU = MarshalUtil()
        fD = mU.doImport(self.__missedIdsPath, fmt="json")
        filterIdD = {ccId: True for ccId in fD["filteredIdList"]}
        self.__testBuildMoleculeCacheFiles(useCache=False, filterIdD=filterIdD, ccFileNamePrefix="cc-filtered")
        #

    def testBuildMoleculeCacheFilesFull(self):
        """ Test construction and reload of full chemical component resource files.
            Running 138 sec sp macbookpro 6.6G resident mem
        """
        self.__testBuildMoleculeCacheFiles(useCache=False, logSizes=True)
        logger.info("Reloading cache file")
        self.__testBuildMoleculeCacheFiles(useCache=True, logSizes=False)

    def __testBuildMoleculeCacheFiles(self, **kwargs):
        """ Test build chemical component cache files from the input component dictionaries
        """
        try:
            ccUrlTarget = kwargs.get("ccUrlTarget", None)
            birdUrlTarget = kwargs.get("birdUrlTarget", None)
            molLimit = kwargs.get("molLimit", None)
            minCount = kwargs.get("minCount", 500)
            useCache = kwargs.get("useCache", True)
            logSizes = kwargs.get("logSizes", False)
            filterIdD = kwargs.get("filterIdD", None)
            ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
            #
            if ccUrlTarget and birdUrlTarget:
                ccmP = ChemCompMoleculeProvider(
                    ccUrlTarget=ccUrlTarget,
                    birdUrlTarget=birdUrlTarget,
                    ccFileNamePrefix=ccFileNamePrefix,
                    cachePath=self.__cachePath,
                    useCache=useCache,
                    molLimit=molLimit,
                    filterIdD=filterIdD,
                )
                ok = ccmP.testCache(minCount=molLimit, logSizes=logSizes)
                self.assertTrue(ok)
            else:
                ccmP = ChemCompMoleculeProvider(cachePath=self.__cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix, molLimit=molLimit, filterIdD=filterIdD)
                ok = ccmP.testCache(minCount=minCount, logSizes=logSizes)
                self.assertTrue(ok)
        except Exception as e:
            logger.info("Failing with %s", str(e))
            self.fail()


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompMoleculeProviderTests("testBuildMoleculeCacheFilesAbbrev"))
    suiteSelect.addTest(ChemCompMoleculeProviderTests("testBuildMoleculeCacheFilesSubset"))
    suiteSelect.addTest(ChemCompMoleculeProviderTests("testBuildMoleculeCacheFilesFiltered"))
    suiteSelect.addTest(ChemCompMoleculeProviderTests("testBuildMoleculeCacheFilesFull"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
