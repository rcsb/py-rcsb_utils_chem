##
# File:    CactvsMolecularFactoryTests.py
# Author:  jdw
# Date:    25-Jan-2020
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests of CactvsMolecularFactory to compare assigned and computed features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.IoUtils import IoUtils
from rcsb.utils.chem.CactvsMoleculeFactory import CactvsMoleculeFactory


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class CactvsMolecularFactoryTests(unittest.TestCase):
    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__cactvsPythonPath = "/apps/lib/cspy"
        self.__molLimit = 50
        self.__exportFlag = True
        #

    def tearDown(self):
        pass

    def testFromSdf(self):
        ctvsmf = CactvsMoleculeFactory(self.__cactvsPythonPath, self.__workPath)
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL:
            ccId = rdCcObj.getName()
            sdfS, atomIdxD = ioU.makeSdf(rdCcObj)
            if self.__exportFlag:
                fp = os.path.join(self.__workPath, rdCcObj.getName() + ".sdf")
                with open(fp, "w") as ofh:
                    ofh.write("%s" % sdfS)
                ctvsmf.setFile(ccId, fp, atomIdxD=atomIdxD)
                ctvsmf.getMoleculeFeatures()
                #
                cactvsJsonPath = os.path.join(self.__workPath, rdCcObj.getName() + ".json")
                ctvsmf.annotate(cactvsJsonPath)
                #
                mD = ctvsmf.getMoleculeFeatures()
                #
            sdfL.append(sdfS)
        #
        # logger.info("\n%s", "\n".join(sdfL))


def suiteCactvsBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CactvsMolecularFactoryTests("testFromSdf"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCactvsBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
