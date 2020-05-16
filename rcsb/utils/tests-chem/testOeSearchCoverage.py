##
# File:    OeSearchCoverageTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for search modes using source molecular definitions coming from a search index. (COVERAGE Tests)

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.chem import __version__
from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSearchIndexUtilsTests(unittest.TestCase):
    skipFlag = False

    def setUp(self):
        self.__useCacheFlag = False
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        self.__fpTypeCuttoffList = [("TREE", 0.6), ("PATH", 0.6), ("MACCS", 0.9), ("CIRCULAR", 0.6), ("LINGO", 0.9)]
        self.__screenTypeList = ["SMARTS"]
        self.__numProc = 1
        self.__startTime = time.time()
        #
        # self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "acdlabs-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        self.__numMols = 31
        self.__myKwargs = {
            "ccUrlTarget": self.__ccUrlTarget,
            "birdUrlTarget": self.__birdUrlTarget,
            "cachePath": self.__cachePath,
            "useCache": self.__useCacheFlag,
            "ccFileNamePrefix": "cc-abbrev",
            "oeFileNamePrefix": "oe-abbrev",
            "limitPerceptions": False,
            "minCount": 30,
            "maxFpResults": 50,
            "fpTypeCuttoffList": self.__fpTypeCuttoffList,
            "buildTypeList": self.__buildTypeList,
            "screenTypeList": self.__screenTypeList,
        }
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __resultContains(self, ccId, matchResultList):
        for matchResult in matchResultList:
            if ccId in matchResult.ccId:
                return True
        return False

    def __getSearchDataProviders(self, **kwargs):
        oesmP = OeSearchMoleculeProvider(**kwargs)
        ok = oesmP.testCache()
        ccIdxP = ChemCompIndexProvider(**kwargs)
        ok = ccIdxP.testCache()
        self.assertTrue(ok)
        ccIdxD = ccIdxP.getIndex()
        return oesmP, ccIdxD

    def testSubStructureSearchExhaustiveAbbrev(self):
        """Exhaustive substructure search. (abbrev)
        """
        return self.__exhaustiveSubStructureSearch(self.__numMols, **self.__myKwargs)

    def __exhaustiveSubStructureSearch(self, numMols, **kwargs):
        """Exhaustive substructure search.
        """
        try:
            limitPerceptions = kwargs.get("limitPerceptions", False)
            buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles"])
            oesmP, ccIdxD = self.__getSearchDataProviders(**kwargs)
            oesU = OeSearchUtils(oesmP, fpTypeList=[])
            oeioU = OeIoUtils()
            #
            for ccId, ccD in list(ccIdxD.items())[:numMols]:
                matchCount = 0
                mtS = set()
                for buildType in buildTypeList:
                    if buildType in ccD:
                        oeMol = oeioU.descriptorToMol(ccD[buildType], buildType, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildType)
                        if not oeMol:
                            logger.error("%s %s build query molecule build fails (skipping)", ccId, buildType)
                            continue
                        # ----
                        startTime = time.time()
                        retStatus, mL = oesU.searchSubStructure(oeMol, matchOpts="graph-strict")
                        if not retStatus:
                            logger.info("%s match fails for build type %s", ccId, buildType)
                        elif not self.__resultContains(ccId, mL):
                            logger.info("%s failed match length %d build type %s in (%.4f seconds)", ccId, len(mL), buildType, time.time() - startTime)
                        elif self.__resultContains(ccId, mL):
                            mtS.update([m.ccId for m in mL])
                            matchCount += 1
                        self.assertTrue(retStatus)
                        self.assertTrue(self.__resultContains(ccId, mL))
                if matchCount:
                    logger.info("%s MATCHES %d: %r", ccId, matchCount, mtS)
                else:
                    logger.info("%s NO MATCHES", ccId)
                # ----
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()
        return False

    def testFingerprintSearchAbbrev(self):
        """Fingerprint search.
        """
        return self.__fingerPrintSearch(self.__numMols, **self.__myKwargs)

    def __fingerPrintSearch(self, numMols, **kwargs):
        maxFpResults = kwargs.get("maxFpResults", 50)
        limitPerceptions = kwargs.get("limitPerceptions", False)
        fpTypeCuttoffList = kwargs.get("fpTypeCuttoffList", [("TREE", 0.6)])
        buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles"])
        #
        oesmP, ccIdxD = self.__getSearchDataProviders(**kwargs)
        oesU = OeSearchUtils(oesmP, fpTypeList=[tup[0] for tup in fpTypeCuttoffList])
        oeioU = OeIoUtils()
        # This will reload the oe binary cache.
        oeMol = oesmP.getMol("004")
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)
        missedFpD = {}
        missedBuildD = {}
        numMols = min(len(ccIdxD), numMols) if numMols else len(ccIdxD)
        logger.info("Begin finger print search on %d molecules", numMols)
        # ----
        startTime = time.time()
        for ccId, ccD in list(ccIdxD.items())[:numMols]:
            for buildType in buildTypeList:
                if buildType in ccD:
                    oeMol = oeioU.descriptorToMol(ccD[buildType], buildType, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildType)
                    if not oeMol:
                        continue
                    selfHit = False
                    for fpType, minFpScore in fpTypeCuttoffList:
                        retStatus, mL = oesU.searchFingerPrints(oeMol, fpType=fpType, minFpScore=minFpScore, maxFpResults=maxFpResults)
                        self.assertTrue(retStatus)
                        #
                        matchedSelf = self.__resultContains(ccId, mL)
                        selfHit = selfHit or matchedSelf
                        if not matchedSelf:
                            missedFpD.setdefault(ccId, []).append((buildType, fpType, len(mL)))
                    #
                    if not selfHit:
                        missedBuildD.setdefault(ccId, []).append(buildType)
        # ------
        for ccId, bTL in missedBuildD.items():
            logger.info("%s missed all fptypes:  buildtype list %r", ccId, bTL)

        if ccId in missedFpD:
            logger.info("%s unmatched by fpTypes %r", ccId, missedFpD[ccId])

        # ----
        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(fpTypeCuttoffList), numMols, time.time() - startTime)
        # ----
        return True


def subStructureSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchExhaustiveAbbrev"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchWithFpAbbrev"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = subStructureSearch()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
