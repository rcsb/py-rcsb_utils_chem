##
# File:    OeChemCompTests.py
# Author:  jdw
# Date:    28-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests of OeChemCompUtils.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeChemCompUtils import OeChemCompUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeChemCompTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccCifPath = os.path.join(self.__cachePath, "cc-cif")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        self.__molLimit = None

    def tearDown(self):
        pass

    def __getChemCompDefs(self, molLimit=None):
        ccMolD = {}
        try:
            useCache = True
            ccFileNamePrefix = "cc-abbrev"
            ccmP = ChemCompMoleculeProvider(
                ccUrlTarget=self.__ccUrlTarget,
                birdUrlTarget=self.__birdUrlTarget,
                cachePath=self.__cachePath,
                useCache=useCache,
                ccFileNamePrefix=ccFileNamePrefix,
                molLimit=molLimit,
            )
            ok = ccmP.testCache(minCount=molLimit)
            self.assertTrue(ok)
            ccMolD = ccmP.getMolD()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ccMolD

    def testBuildCifFromOE(self):
        """Build chemical component definitions from OE Mol object"""
        try:
            ccMolD = self.__getChemCompDefs()
            oemf = OeMoleculeFactory()
            #
            for ccId, ccObj in list(ccMolD.items())[:10]:
                # ----
                tId = oemf.setChemCompDef(ccObj)
                self.assertEqual(tId, ccId)
                ok = oemf.build(molBuildType="model-xyz")
                self.assertTrue(ok)
                fp = os.path.join(self.__ccCifPath, ccId + "-gen.cif")
                oeMol = oemf.getMol()
                oeccU = OeChemCompUtils()
                ok = oeccU.addOeMol(ccId, oeMol, missingModelXyz=False, writeIdealXyz=False)
                self.assertTrue(ok)
                ok = oeccU.write(fp)

                # ----
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteBuildCifTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeChemCompTests("testBuildCifFromOE"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteBuildCifTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
