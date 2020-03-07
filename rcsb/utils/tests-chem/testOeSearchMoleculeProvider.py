##
# File:    OeSearchMoleculeProviderTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to read and process search OE molecule instances and related databases created
from chemical component search indices.

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
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeSearchMoleculeProviderTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        #
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSearchBuildMoleculeCacheFilesAbbrev(self):
        self.__testBuildSearchMoleculeCacheFiles(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            ccFileNamePrefix="cc-abbrev",
            oeFileNamePrefix="oe-abbrev",
            molLimit=None,
            fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"],
            screenTypeList=["SMARTS"],
            numProc=2,
        )

    def testSearchBuildMoleculeCacheFilesFull(self):
        self.__testBuildSearchMoleculeCacheFiles(
            ccFileNamePrefix="cc-full", oeFileNamePrefix="oe-full", molLimit=None, fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"], screenTypeList=[], numProc=2,
        )

    def __testBuildSearchMoleculeCacheFiles(self, **kwargs):
        """ Test build OE cache files from full component dictionary
        """
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        molLimit = kwargs.get("molLimit", 0)
        quietFlag = kwargs.get("quietFlag", True)
        fpTypeList = kwargs.get("fpTypeList", ["TREE"])
        screenTypeList = kwargs.get("screenTypeList", [])
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        oeFileNamePrefix = kwargs.get("oeFileNamePrefix", "oe")
        numProc = kwargs.get("numProc", 2)
        #
        startTime = time.time()
        if ccUrlTarget and birdUrlTarget and molLimit:
            # Using abbreviated reference source files
            oesmp = OeSearchMoleculeProvider(
                ccUrlTarget=ccUrlTarget,
                birdUrlTarget=birdUrlTarget,
                cachePath=self.__cachePath,
                ccFileNamePrefix=ccFileNamePrefix,
                oeFileNamePrefix=oeFileNamePrefix,
                useCache=False,
                quietFlag=quietFlag,
                fpTypeList=fpTypeList,
                screenTypeList=screenTypeList,
                numProc=numProc,
            )
        else:
            oesmp = OeSearchMoleculeProvider(
                cachePath=self.__cachePath,
                ccFileNamePrefix=ccFileNamePrefix,
                oeFileNamePrefix=oeFileNamePrefix,
                useCache=False,
                quietFlag=quietFlag,
                molLimit=molLimit,
                fpTypeList=fpTypeList,
                screenTypeList=screenTypeList,
                numProc=numProc,
            )
        ok = oesmp.testCache()
        self.assertTrue(ok)
        #
        endTime = time.time()
        logger.info(">> Completed load molLimit %r (%.4f seconds)", molLimit, endTime - startTime)
        #
        # ---
        oesmp = OeSearchMoleculeProvider(cachePath=self.__cachePath, ccFileNamePrefix=ccFileNamePrefix, oeFileNamePrefix=oeFileNamePrefix, useCache=True)
        #
        deltaMol = 2
        minMol = minNumFp = molLimit - deltaMol if molLimit else 500
        for fpType in fpTypeList:
            fpDb = oesmp.getFingerPrintDb(fpType="TREE")
            logger.debug("fpType %r length %d", fpType, fpDb.NumFingerPrints())
            self.assertGreaterEqual(fpDb.NumFingerPrints(), minNumFp)
        #
        ccId = "004"
        oeMol = oesmp.getMol(ccId)
        logger.debug("%s atom count %d", ccId, len(list(oeMol.GetAtoms())))
        #
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)
        #
        oeDb, oeDbIdx = oesmp.getOeMolDatabase()
        logger.debug("Type db %r length %d type idx %r length %d", type(oeDb), oeDb.NumMols(), type(oeDbIdx), len(oeDbIdx))
        self.assertGreaterEqual(oeDb.NumMols(), minMol)
        self.assertGreaterEqual(len(oeDbIdx), minMol)
        #
        if screenTypeList:
            ssDb = oesmp.getSubSearchDb()
            self.assertGreaterEqual(ssDb.NumMolecules(), minMol)


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchMoleculeProviderTests("testBuildSearchMoleculeCacheFilesAbbrev"))

    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
