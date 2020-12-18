##
# File:    testOeMcssSearchUtils.py
# Author:  jdw
# Date:    17-Dec-2020
# Version: 0.001
#
# Updates:
##
"""
A collection of tests for MCSS comparison operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.chem import __version__
from rcsb.utils.chem.OeMcssSearchUtils import OeMcssSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeMcssSearchUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        #
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testMCSSAlignCCPair(self):
        """Test case -  Simple pairwise MCSS alignment   (cif) """
        try:
            refPath = os.path.join(self.__dataPath, "001.cif")
            fitPath = os.path.join(self.__dataPath, "001.cif")
            oed = OeMcssSearchUtils(verbose=self.__verbose)
            oed.setSearchType(sType="relaxed")
            oed.setRefPath(refPath)
            logger.info("Ref path %s", refPath)
            oed.setFitPath(fitPath)
            logger.info("Fit path %s", fitPath)
            (nAtomsRef, refFD, nAtomsFit, fitFD, atomMapL) = oed.doAlign(minFrac=0.8)
            self.assertEqual(nAtomsRef, nAtomsFit)
            self.assertEqual(refFD["InChIKey"], fitFD["InChIKey"])
            if len(atomMapL) > 0:
                for (rCC, rIdx, rAtNum, rAt, fCC, fIdx, fAtNum, fAt) in atomMapL:
                    logger.debug("%5s %5s %5s %-5s %5s %5s %5s %-5s", rCC, rIdx, rAtNum, rAt, fCC, fIdx, fAtNum, fAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSSAlignMixedPair(self):
        """Test case -  Simple pairwise MCSS alignment   (sdf/cif) """
        try:
            refPath = os.path.join(self.__dataPath, "001.mol2")
            fitPath = os.path.join(self.__dataPath, "001.cif")

            oed = OeMcssSearchUtils(verbose=self.__verbose)
            oed.setSearchType(sType="relaxed")
            logger.info("Ref path %s", refPath)
            oed.setRefPath(refPath, title=None, suppressHydrogens=False, fType="mol2", importType="3D")
            #
            logger.info("Fit path %s", fitPath)
            oed.setFitPath(fitPath)

            (nAtomsRef, refFD, nAtomsFit, fitFD, atomMapL) = oed.doAlign(minFrac=0.8)
            self.assertEqual(nAtomsRef, nAtomsFit)
            self.assertEqual(refFD["SMILES"], fitFD["SMILES"])
            if len(atomMapL) > 0:
                for (rCC, rIdx, rAtNum, rAt, fCC, fIdx, fAtNum, fAt) in atomMapL:
                    logger.debug("%5s %5s %5s %-5s %5s %5s %5s %-5s", rCC, rIdx, rAtNum, rAt, fCC, fIdx, fAtNum, fAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteAlignPair():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMcssSearchUtilsTests("testMCSSAlignPair"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteAlignPair()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
