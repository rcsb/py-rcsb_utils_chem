##
# File:    OeMoleculeProviderTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to read and process OE molecule instances and related databases created
from of PDB chemical component definitions.

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
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeMoleculeProviderTests(unittest.TestCase):
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

    def testBuildMoleculeCacheFilesAbbrev(self):
        for molBuildType in ["oe-iso-smiles", "model-xyz", "ideal-xyz"]:
            self.__testBuildMoleculeCacheFiles(
                ccUrlTarget=self.__ccUrlTarget,
                birdUrlTarget=self.__birdUrlTarget,
                ccFileNamePrefix="cc-abbrev",
                oeFileNamePrefix="oe-abbrev",
                molLimit=500,
                fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"],
                screenTypeList=["SMARTS"],
                molBuildType=molBuildType,
            )

    def testBuildMoleculeCacheFilesSubset(self):
        """ Test construction of OE resource files for a subset (1K) of full public chemical reference dictionary
        """
        for molBuildType in ["oe-iso-smiles", "model-xyz", "ideal-xyz"]:
            self.__testBuildMoleculeCacheFiles(
                molBuildType=molBuildType,
                molLimit=1000,
                ccFileNamePrefix="cc-1k",
                oeFileNamePrefix="oe-1k",
                screenTypeList=["SMARTS"],
                fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"],
            )
        #

    def testBuildMoleculeCacheFilesFiltered(self):
        """ Test construction of OE resource files for a filtered subset of full public chemical reference dictionary
        """
        for molBuildType in ["oe-iso-smiles", "model-xyz", "ideal-xyz"]:
            self.__testBuildMoleculeCacheFiles(
                molBuildType=molBuildType,
                molLimit=1000,
                ccFileNamePrefix="cc-filtered",
                oeFileNamePrefix="oe-filtered",
                screenTypeList=["SMARTS"],
                fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"],
            )
        #

    def testBuildMoleculeCacheFilesFull(self):
        """ Test construction of OE resource files for the full public chemical reference dictionary

            For 31K cc's the build time is ~3123 secs (4 procs)

        """
        for molBuildType in ["oe-iso-smiles", "model-xyz", "ideal-xyz"]:
            self.__testBuildMoleculeCacheFiles(
                molBuildType=molBuildType, ccFileNamePrefix="cc-full", oeFileNamePrefix="oe-full", fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"]
            )

    def testBuildMoleculeCacheFilesFullScreened(self):
        """ Test construction of OE resource files with screened substructure database for the full public chemical reference dictionary
        """
        for molBuildType in ["oe-iso-smiles"]:
            self.__testBuildMoleculeCacheFiles(
                molBuildType=molBuildType,
                ccFileNamePrefix="cc-full",
                oeFileNamePrefix="oe-full",
                screenTypeList=["SMARTS"],
                fpTypeList=["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"],
            )

    def __testBuildMoleculeCacheFiles(self, **kwargs):
        """ Test build OE cache files from full component dictionary
        """
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        molLimit = kwargs.get("molLimit", 0)
        quietFlag = kwargs.get("quietFlag", True)
        molBuildType = kwargs.get("molBuildType", "ideal-xyz")
        fpTypeList = kwargs.get("fpTypeList", ["TREE"])
        screenTypeList = kwargs.get("screenTypeList", [])
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        oeFileNamePrefix = kwargs.get("oeFileNamePrefix", "oe")
        #
        startTime = time.time()
        if ccUrlTarget and birdUrlTarget and molLimit:
            # Using abbreviated reference source files
            oemp = OeMoleculeProvider(
                ccUrlTarget=ccUrlTarget,
                birdUrlTarget=birdUrlTarget,
                cachePath=self.__cachePath,
                ccFileNamePrefix=ccFileNamePrefix,
                oeFileNamePrefix=oeFileNamePrefix,
                molBuildType=molBuildType,
                useCache=False,
                quietFlag=quietFlag,
                fpTypeList=fpTypeList,
                screenTypeList=screenTypeList,
            )
        else:
            oemp = OeMoleculeProvider(
                cachePath=self.__cachePath,
                ccFileNamePrefix=ccFileNamePrefix,
                oeFileNamePrefix=oeFileNamePrefix,
                molBuildType=molBuildType,
                useCache=False,
                quietFlag=quietFlag,
                molLimit=molLimit,
                fpTypeList=fpTypeList,
                screenTypeList=screenTypeList,
            )
        ok = oemp.testCache()
        self.assertTrue(ok)
        #
        endTime = time.time()
        logger.info(">> Completed load molBuildType %r molLimit %r (%.4f seconds)", molBuildType, molLimit, endTime - startTime)
        #
        # ---
        oemp = OeMoleculeProvider(cachePath=self.__cachePath, ccFileNamePrefix=ccFileNamePrefix, oeFileNamePrefix=oeFileNamePrefix, molBuildType=molBuildType, useCache=True)
        #
        deltaMol = 2
        minMol = minNumFp = molLimit - deltaMol if molLimit else 30000
        for fpType in fpTypeList:
            fpDb = oemp.getFingerPrintDb(fpType="TREE")
            logger.debug("fpType %r length %d", fpType, fpDb.NumFingerPrints())
            self.assertGreaterEqual(fpDb.NumFingerPrints(), minNumFp)
        #
        ccId = "004"
        oeMol = oemp.getMol(ccId)
        logger.debug("%s atom count %d", ccId, len(list(oeMol.GetAtoms())))
        #
        if molBuildType in ["oe-iso-smiles"]:
            self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)
        else:
            self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 20)
        #
        oeDb, oeDbIdx = oemp.getOeMolDatabase()
        logger.debug("Type db %r length %d type idx %r length %d", type(oeDb), oeDb.NumMols(), type(oeDbIdx), len(oeDbIdx))
        self.assertGreaterEqual(oeDb.NumMols(), minMol)
        self.assertGreaterEqual(len(oeDbIdx), minMol)
        #
        if molBuildType in ["oe-iso-smiles"] and screenTypeList:
            ssDb = oemp.getSubSearchDb()
            self.assertGreaterEqual(ssDb.NumMolecules(), minMol)


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFilesAbbrev"))
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFilesSubset"))
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFilesFiltered"))
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFilesFull"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
