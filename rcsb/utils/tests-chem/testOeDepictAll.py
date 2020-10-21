##
#
# File:    OeDepictAllTests.py
# Author:  jdw
# Date:    5-May-2013
# Version: 0.001
#
# Updates:
#  4-May-2014 jdw add example for depiction from SMILES input
#  6-Jun-2016 jdw general cleanup
##
"""
A collection of tests for the OEDepict and related classes.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import logging
import os
import os.path
import time
import unittest

from rcsb.utils.chem.OeDepict import OeDepict
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeDepictAllTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        #
        self.__oeMolD = self.__getCache()

    def tearDown(self):
        pass

    def __getCache(self):
        oemp = OeMoleculeProvider(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            ccFileNamePrefix="cc-abbrev",
            cachePath=self.__cachePath,
            molBuildType="model-xyz",
            useCache=True,
            oeFileNamePrefix="oe-abbrev",
        )
        ok = oemp.testCache()
        self.assertTrue(ok)
        return oemp.getOeMolD()

    def testDepictAll(self):
        """Test case -  single OE molecule depiction."""
        try:
            for ccId, oeMol in self.__oeMolD.items():
                imagePath = os.path.join(self.__workPath, "image", ccId[0], ccId + ".svg")
                oed = OeDepict()
                title = ""
                oed.setMolTitleList([(ccId, oeMol, title)])
                oed.setDisplayOptions(
                    labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, labelBondCIPStereo=True, cellBorders=False, bondDisplayWidth=0.5
                )
                oed.setGridOptions(rows=1, cols=1)
                oed.prepare()
                oed.write(imagePath)
            for ccId, oeMol in self.__oeMolD.items():
                imagePath = os.path.join(self.__workPath, "image_labeled", ccId[0], ccId + ".svg")
                oed = OeDepict()
                title = ""
                oed.setMolTitleList([(ccId, oeMol, title)])
                oed.setDisplayOptions(
                    labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, abelBondCIPStereo=True, cellBorders=False, bondDisplayWidth=0.5
                )
                oed.setGridOptions(rows=1, cols=1)
                oed.prepare()
                oed.write(imagePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteDepict():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictAllTests("testDepictAll"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteDepict()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
