##
# File:    ChemAxonDescriptorProviderTests.py
# Author:  J. Westbrook
# Date:    17-Aug-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Utilities to deliver ChemAxon rendered chemical descriptors for chemical component definitions.
"""

__docformat__ = "google en"
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
from rcsb.utils.chem.ChemAxonDescriptorProvider import ChemAxonDescriptorProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChemAxonDescriptorProviderTests(unittest.TestCase):
    buildFull = False

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

    def testChemAxonDescriptorAbbrev(self):
        """Test construction of abbreviated ChemAxon mapping file."""
        self.__buildMoleculeCacheFiles(buildFlag=True, ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, molLimit=10, useCache=False, ccFileNamePrefix="cc-abbrev")

    def testChemAxonDescriptorCacheFull(self):
        """Test get cached full ChemAxon mapping file."""
        self.__buildMoleculeCacheFiles(buildFlag=False, useCache=True, ccFileNamePrefix="cc-full")

    @unittest.skipUnless(buildFull, "Troubleshooting test")
    def testChemAxonDescriptorFull(self):
        """Test construction of full ChemAxon mapping file."""
        self.__buildMoleculeCacheFiles(buildFlag=True, useCache=False, ccFileNamePrefix="cc-full")

    def __buildMoleculeCacheFiles(self, buildFlag=False, **kwargs):
        """Test build mapping and chemical component cache files from the input component dictionaries"""
        molLimit = kwargs.get("molLimit", None)
        useCache = kwargs.get("useCache", True)
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc-full")
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        #
        caxP = ChemAxonDescriptorProvider(
            ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget, cachePath=self.__cachePath, useCache=useCache, molLimit=molLimit, ccFileNamePrefix=ccFileNamePrefix
        )
        if buildFlag:
            caxP.buildDescriptors()
            ok = caxP.testCache()
            #
            ok = caxP.updateDescriptors(useCache=True)
            self.assertTrue(ok)
        #
        caxP = ChemAxonDescriptorProvider(
            ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget, cachePath=self.__cachePath, useCache=True, molLimit=molLimit, ccFileNamePrefix=ccFileNamePrefix
        )
        ok = caxP.testCache()
        self.assertTrue(ok)

        return ok


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemAxonDescriptorProviderTests("testChemAxonDescriptorAbbrev"))
    suiteSelect.addTest(ChemAxonDescriptorProviderTests("testChemAxonDescriptorCacheFull"))
    suiteSelect.addTest(ChemAxonDescriptorProviderTests("testChemAxonDescriptorFull"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
