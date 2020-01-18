##
# File:    OeMoleculeProviderTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
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
import time
import unittest


from rcsb.utils.chem import __version__
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeMoleculeProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__startTime = time.time()
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        # self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildMoleculeCacheFiles(self):
        """ Test build OE cache files from full component dictionary
        """
        minCount = 500
        oemp = OeMoleculeProvider(
            ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType="model", useCache=False
        )
        ok = oemp.testCache(minCount=minCount)
        self.assertTrue(ok, minCount)
        logger.info(" ******* Completed initial load ******** ")
        #
        oemp = OeMoleculeProvider(
            ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType="model", useCache=True
        )
        ok = oemp.testCache(minCount=minCount)
        self.assertTrue(ok)
        ccId = "004"
        oeMol = oemp.getMol(ccId)
        logger.info("%s atom count %d", ccId, len(list(oeMol.GetAtoms())))
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 20)
        #
        oeDb, oeDbIdx = oemp.getOeMolDatabase()
        logger.info("Type db %r idx %r", type(oeDb), type(oeDbIdx))
        #
        #
        #
        ssDb = oemp.getSubSearchDb()
        logger.info("Type %r", type(ssDb))
        fpDb = oemp.getFingerPrintDb()
        logger.info("Type %r", type(fpDb))


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFiles"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
