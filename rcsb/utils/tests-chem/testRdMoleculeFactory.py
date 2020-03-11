##
# File:    RdMolecularFactoryTests.py
# Author:  jdw
# Date:    05-Jan-2020
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests of RdMolecularFactory to compare assigned and computed features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.IoUtils import IoUtils

try:
    from rcsb.utils.chem.RdMoleculeFactory import RdMoleculeFactory
except Exception:
    pass


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class RdMolecularFactoryTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__molLimit = 50
        self.__exportFlag = True
        #

    def tearDown(self):
        pass

    @unittest.skipIf(skipFlag, "Special install RDKit dependency")
    def testFromSdf(self):
        rdmf = RdMoleculeFactory()
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL:
            ccId = rdCcObj.getName()
            _, sdfS, atomIdxD = ioU.makeSdf(rdCcObj)
            if self.__exportFlag:
                fp = os.path.join(self.__workPath, rdCcObj.getName() + ".sdf")
                with open(fp, "w") as ofh:
                    ofh.write("%s" % sdfS)
                rdmf.setFile(ccId, fp, atomIdxD=atomIdxD)
                rdmf.updateProperties()
                rdmf.updateStereoChemAssignments()
                rdmf.updateAromaticAssignments()
                _ = rdmf.getMoleculeFeatures()
                #
            sdfL.append(sdfS)
        #
        logger.info("\n%s", "\n".join(sdfL))


def suiteRdBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(RdMolecularFactoryTests("testFromSdf"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteRdBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
