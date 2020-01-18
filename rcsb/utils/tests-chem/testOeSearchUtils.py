##
# File:    OeSearchUtilsTests.py
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
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSearchUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        self.__minCount = 500
        self.__startTime = time.time()
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSubStructureSearch(self):
        oemp = OeMoleculeProvider(ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, dirPath=os.path.join(self.__cachePath, "chem_comp"), useCache=True)
        ok = oemp.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        oesU = OeSearchUtils(oemp)
        for ccId in ["001", "004"]:
            oeMol = oemp.getMol(ccId)
            matchMolL = oesU.searchSubStructure(oeMol)
            logger.info("%s Match length %d", ccId, len(matchMolL))


def fingerprintSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchUtilsTests("testSubStructureSearch"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = fingerprintSearch()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
