##
# File:    ObMolecularFactoryTests.py
# Author:  jdw
# Date:    11-Jan-2020
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests of ObMolecularFactory to compare assigned and computed features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.IoUtils import IoUtils
from rcsb.utils.chem.ObMoleculeFactory import ObMoleculeFactory


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ObMolecularFactoryTests(unittest.TestCase):
    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__molLimit = 50
        self.__exportFlag = True
        #

    def tearDown(self):
        pass

    def testFromSdfFile(self):
        obmf = ObMoleculeFactory()
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
                obmf.setFile(ccId, fp, molFormat="mol", atomIdxD=atomIdxD)
                # obmf.updateProperties()
                # obmf.updateStereoChemAssignments()
                # obmf.updateAromaticAssignments()
                # _ = obmf.getRdMoleculeFeatures()
                #
            sdfL.append(sdfS)
        #
        logger.debug("\n%s", "\n".join(sdfL))

    def testFromSdfString(self):
        obmf = ObMoleculeFactory()
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL:
            ccId = rdCcObj.getName()
            _, sdfS, atomIdxD = ioU.makeSdf(rdCcObj)
            obmf.setString(ccId, sdfS, molFormat="mol", atomIdxD=atomIdxD)
            obD = obmf.getMoleculeFeatures()
            logger.info("dictionary %r", obD)
            #
            sdfL.append(sdfS)
        #
        logger.debug("\n%s", "\n".join(sdfL))


def suiteObBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ObMolecularFactoryTests("testFromSdfFile"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteObBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
