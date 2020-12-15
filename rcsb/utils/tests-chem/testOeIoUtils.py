##
# File:    OeOeIoUtilsTests.py
# Author:  jdw
# Date:    28-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests for OeIoUtils() operations.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.chem.PdbxChemComp import PdbxChemCompIt

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeIoUtilsTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__molfileDirPath = os.path.join(self.__cachePath, "molfiles")
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

    def testIpOps(self):
        """Test IO operation on generated related molecules"""
        try:
            oeIoU = OeIoUtils()
            mU = MarshalUtil()
            mU.mkdir(self.__molfileDirPath)
            ccMolD = self.__getChemCompDefs()
            oemf = OeMoleculeFactory()
            for ccId, ccObj in list(ccMolD.items())[:10]:
                # ----
                tId = oemf.setChemCompDef(ccObj)
                self.assertEqual(tId, ccId)
                relatedIdxD = oemf.buildRelated(limitPerceptions=False)
                logger.info("%s generated %d molecular forms", ccId, len(relatedIdxD))
                for sId, idxD in relatedIdxD.items():
                    logger.info("sId %r smiles %r", sId, idxD["smiles"])
                    mol2Path = os.path.join(self.__molfileDirPath, sId + ".mol2")
                    oeMol = oeIoU.descriptorToMol(idxD["smiles"], "oe-iso-smiles", limitPerceptions=False, messageTag=None)
                    oeIoU.write(mol2Path, oeMol, constantMol=True, addSdTags=True)
                    sdfPath = os.path.join(self.__molfileDirPath, sId + ".mol")
                    oeMol = oeIoU.descriptorToMol(idxD["smiles"], "oe-iso-smiles", limitPerceptions=False, messageTag=None)
                    oeIoU.write(sdfPath, oeMol, constantMol=True, addSdTags=True)

                # ----
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteIoOpsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeIoUtilsTests("testIoOps"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteIoOpsTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
