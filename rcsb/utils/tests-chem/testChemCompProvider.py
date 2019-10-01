##
# File:    ChemCompProviderTests.py
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
import unittest

from rcsb.utils.chem.ChemCompProvider import ChemCompProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChemCompProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReadChemCompDictionary(self):
        ccm = ChemCompProvider(dirPath=os.path.join(self.__workPath, "chem_comp"), useCache=False)
        ccm.testCache()
        rD = {}
        # rD = ccm.getMapping()
        #
        # logger.info("Chemical component length %d", len(rD))
        for ccId in rD:
            if len(rD[ccId]) > 1:
                logger.info("Match length for %r %d", ccId, len(rD[ccId]))


def readChemCompDictionary():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompProviderTests("testReadChemCompDictionary"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompDictionary()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
