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
from rcsb.utils.chem.MoleculeAnnotationsCompare import MoleculeAnnotationsCompare


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class CactvsMoleculeFactoryTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cactvsPythonPath = "/apps/lib/cspy"
        self.__molLimit = 50
        self.__exportFlag = True
        #

    def tearDown(self):
        pass

    @unittest.skipIf(skipFlag, "Special install CACTVS dependency")
    def testFromSdf(self):
        aroModel = "daylight"
        retD = {}
        maC = MoleculeAnnotationsCompare()
        ctvsmf = CactvsMoleculeFactory(self.__cactvsPythonPath, self.__workPath)
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL:
            ccId = rdCcObj.getName()
            _, sdfS, atomIdxD = ioU.makeSdf(rdCcObj)
            if self.__exportFlag:
                fp = os.path.join(self.__workPath, rdCcObj.getName() + ".sdf")
                with open(fp, "w", encoding="utf-8") as ofh:
                    ofh.write("%s" % sdfS)
                ctvsmf.setFile(ccId, fp, atomIdxD=atomIdxD)
                cactvsJsonPath = os.path.join(self.__workPath, rdCcObj.getName() + ".json")
                ctvsmf.annotate(cactvsJsonPath, aroModel=aroModel)
                #
                tstFD = ctvsmf.getMoleculeFeatures(aroModel=aroModel)
                refFD = maC.getChemCompFeatures(rdCcObj, descriptorProgram="CACTVS")
                ok, retCmp = maC.compare(refFD, tstFD, tstInfo="cactvs default")
                if not ok:
                    retD[ccId] = retCmp
                #
            sdfL.append(sdfS)
        #
        logger.info("Components processed %d differences %d", len(rdCcObjL), len(retD))
        # logger.info("\n%s", "\n".join(sdfL))


def suiteCactvsBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CactvsMoleculeFactoryTests("testFromSdf"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCactvsBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
