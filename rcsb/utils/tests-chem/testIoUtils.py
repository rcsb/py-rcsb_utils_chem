##
# File:    IoUtilsTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to read and write chemical component definitions.

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
from rcsb.utils.chem.IoUtils import IoUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class IoUtilsTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__startTime = time.time()
        self.__exportFlag = True
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSdfGenerator(self):
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL:
            ok, sdfS, _ = ioU.makeSdf(rdCcObj, molBuildType="ideal-xyz")
            if not ok:
                continue
            if self.__exportFlag:
                fp = os.path.join(self.__workPath, rdCcObj.getName() + ".sdf")
                with open(fp, "w") as ofh:
                    ofh.write("%s" % sdfS)
            sdfL.append(sdfS)
        #

    @unittest.skipIf(skipFlag, "Long test")
    def testSdfBulkGenerator(self):
        ioU = IoUtils()
        sdfL = []
        rdCcObjL = ioU.getComponentDefinitions(os.path.join(self.__dataPath, "components.cif.gz"))
        self.assertGreaterEqual(len(rdCcObjL), 4)
        for rdCcObj in rdCcObjL[:500]:
            ok, sdfS, _ = ioU.makeSdf(rdCcObj, molBuildType="ideal-xyz")
            if not ok:
                continue
            sdfL.append(sdfS)
        #
        fp = os.path.join(self.__workPath, "components.sdf")
        with open(fp, "w") as ofh:
            ofh.write("%s" % "\n".join(sdfL))


def sdfWriterSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(IoUtilsTests("testSdfGenerator"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = sdfWriterSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
