##
# File:    OeSearchIndexUtilsTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for search modes using source molecular definitions coming from a search index.

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
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSearchIndexUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        # self.__fpTypeList = ["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"]
        self.__fpTypeCuttoffList = [("TREE", 0.6), ("PATH", 0.6), ("MACCS", 0.9), ("CIRCULAR", 0.6), ("LINGO", 0.9)]
        self.__screenType = "SMARTS"
        self.__numProc = 1
        self.__minCount = 500
        self.__startTime = time.time()
        #
        # self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "acdlabs-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        self.__buildTypeList = ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]
        self.__numMols = 3000
        self.__myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "ccFileNamePrefix": "cc-abbrev",
            "oeFileNamePrefix": "oe-abbrev",
            "limitPerceptions": False,
            "minCount": 500,
            "maxFpResults": 50,
            "fpTypeCuttoffList": self.__fpTypeCuttoffList,
            "buildTypeList": self.__buildTypeList,
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
        minCount = kwargs.get("minCount", 500)
        oesmP = OeSearchMoleculeProvider(**kwargs)
        ok = oesmP.testCache()
        ccIdxP = ChemCompIndexProvider(**kwargs)
        ok = ccIdxP.testCache(minCount=minCount)
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
                for buildType in buildTypeList:
                    if buildType in ccD:
                        oeMol = oeioU.descriptorToMol(ccD[buildType], buildType, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildType)
                        if not oeMol:
                            continue
                        # ----
                        startTime = time.time()
                        retStatus, mL = oesU.searchSubStructure(oeMol, matchOpts="simple")
                        logger.info("%s match length %d build type %s in (%.4f seconds)", ccId, len(mL), buildType, time.time() - startTime)
                        self.assertTrue(retStatus)
                        self.assertTrue(self.__resultContains(ccId, mL))
                # ----
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

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

    def testFingerPrintScoresAbbrev(self):
        """Fingerprint scores. (abbreviated)
        """
        return self.__fingerPrintScores(self.__numMols, **self.__myKwargs)

    def testFingerPrintScoresFull(self):
        """Fingerprint scores. (full)
        """
        numMols = 300000
        myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "ccFileNamePrefix": "cc-full",
            "oeFileNamePrefix": "oe-full",
            "limitPerceptions": False,
            "minCount": 500,
            "maxFpResults": 50,
            "fpTypeCuttoffList": self.__fpTypeCuttoffList,
            "buildTypeList": self.__buildTypeList,
        }
        return self.__fingerPrintScores(numMols, **myKwargs)

    def __fingerPrintScores(self, numMols, **kwargs):
        maxFpResults = kwargs.get("maxResults", 50)
        limitPerceptions = kwargs.get("limitPerceptions", False)
        fpTypeCuttoffList = kwargs.get("fpTypeCuttoffList", [("TREE", 0.6)])
        buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles"])
        doDisplay = kwargs.get("doDisplay", False)
        #
        oesmP, ccIdxD = self.__getSearchDataProviders(**kwargs)
        oesU = OeSearchUtils(oesmP, fpTypeList=[tup[0] for tup in fpTypeCuttoffList])
        oeioU = OeIoUtils()
        # This will reload the oe binary cache.
        oeMol = oesmP.getMol("004")
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)
        #
        missedFpD = {}
        missedBuildD = {}
        numMols = min(len(ccIdxD), numMols) if numMols else len(ccIdxD)
        logger.info("Begin finger print score search on %d molecules", numMols)
        # ----
        startTime = time.time()
        for ii, ccId, in enumerate(list(ccIdxD.keys())[:numMols]):
            ccD = ccIdxD[ccId]
            for buildType in buildTypeList:
                if buildType in ccD:
                    oeMol = oeioU.descriptorToMol(ccD[buildType], buildType, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildType)
                    if not oeMol:
                        logger.debug("%s build failed for %s - skipping", ccId, buildType)
                        continue
                    maxHits = 0
                    minHits = maxFpResults
                    selfHit = False
                    #
                    startTime1 = time.time()
                    for fpType, minFpScore in fpTypeCuttoffList:
                        retStatus, mL = oesU.getFingerPrintScores(oeMol, fpType, minFpScore, maxFpResults)
                        self.assertTrue(retStatus)
                        logger.debug("%s fpType %r hits %d", ccId, fpType, len(mL))
                        maxHits = max(maxHits, len(mL))
                        minHits = min(minHits, len(mL))
                        matchedSelf = self.__resultContains(ccId, mL)
                        selfHit = selfHit or matchedSelf
                        if not matchedSelf:
                            missedFpD.setdefault(ccId, []).append((buildType, fpType, len(mL)))
                    #
                    if not selfHit:
                        missedBuildD.setdefault(ccId, []).append(buildType)
                    #
                    if maxHits < 1 or not selfHit:
                        logger.info("%s buildType %r min hits %d max hits %d (%.4f seconds)", ccId, buildType, minHits, maxHits, time.time() - startTime1)
                else:
                    logger.debug("%s missing descriptor %r", ccId, buildType)
            if ii % 100 == 0:
                logger.info("Completed %d of %d missed count %d", ii, numMols, len(missedBuildD))

        # ------
        for ccId, bTL in missedBuildD.items():
            logger.info("%s missed all fptypes:  buildtype list %r", ccId, bTL)

        if ccId in missedFpD:
            logger.info("%s unmatched by fpTypes %r", ccId, missedFpD[ccId])

        #
        if doDisplay:
            for ccId, bTL in missedBuildD.items():
                idxD = ccIdxD[ccId]
                if "oe-iso-smiles" in idxD:
                    for bT in bTL:
                        self.__displayAlignedDescriptorPair(ccId, idxD["oe-iso-smiles"], "oe-iso-smiles", idxD[bT], bT, title=None, limitPerceptions=True)

        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(fpTypeCuttoffList), numMols, time.time() - startTime)
        # ----                                  ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):

    def testSubStructureSearchWithFpAbbrev(self):
        """Substructure search with fingerprint prefilter. (abbreviated data)
        """
        return self.__sssWithFingerPrintFromDescriptor(self.__numMols, **self.__myKwargs)

    def testSubStructureSearchWithFpFull(self):
        """Substructure search with fingerprint prefilter. (full)
        """
        numMols = 300000
        myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "ccFileNamePrefix": "cc-full",
            "oeFileNamePrefix": "oe-full",
            "limitPerceptions": False,
            "minCount": 500,
            "maxFpResults": 50,
            "fpTypeCuttoffList": self.__fpTypeCuttoffList,
            "buildTypeList": self.__buildTypeList,
        }
        return self.__sssWithFingerPrintFromDescriptor(numMols, **myKwargs)

    def __sssWithFingerPrintFromDescriptor(self, numMols, **kwargs):
        maxFpResults = kwargs.get("maxResults", 50)
        limitPerceptions = kwargs.get("limitPerceptions", False)
        fpTypeCuttoffList = kwargs.get("fpTypeCuttoffList", [("TREE", 0.6)])
        buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles"])
        doDisplay = kwargs.get("doDisplay", False)
        #
        oesmP, ccIdxD = self.__getSearchDataProviders(**kwargs)
        oesU = OeSearchUtils(oesmP, fpTypeList=[tup[0] for tup in fpTypeCuttoffList])
        oeioU = OeIoUtils()
        # This will reload the oe binary cache.
        oeMol = oesmP.getMol("004")
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)

        matchOpts = "simple"
        missTupL = []
        missedD = {}
        missedFpD = {}
        numMols = min(len(ccIdxD), numMols) if numMols else len(ccIdxD)
        logger.info("Begin substructure search w/ finger print filter on %d molecules", numMols)
        # ----
        startTime = time.time()
        for ii, ccId, in enumerate(list(ccIdxD.keys())[:numMols]):
            ccD = ccIdxD[ccId]
            for buildType in buildTypeList:
                if buildType in ccD:
                    startTime1 = time.time()
                    oeMol = oeioU.descriptorToMol(ccD[buildType], buildType, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildType)
                    if not oeMol:
                        logger.debug("%s build failed for %s - skipping", ccId, buildType)
                        continue
                    maxHits = 0
                    minHits = maxFpResults
                    selfHit = False
                    for fpType, minFpScore in fpTypeCuttoffList:
                        retStatus, mL = oesU.searchSubStructureWithFingerPrint(oeMol, fpType, minFpScore, maxFpResults, matchOpts=matchOpts)
                        self.assertTrue(retStatus)
                        logger.debug("%s fpType %r hits %d", ccId, fpType, len(mL))
                        maxHits = max(maxHits, len(mL))
                        minHits = min(minHits, len(mL))
                        matchedSelf = self.__resultContains(ccId, mL)
                        selfHit = selfHit or matchedSelf
                        if not matchedSelf:
                            missedFpD.setdefault(ccId, []).append((buildType, fpType, len(mL)))
                    if not selfHit:
                        missedD.setdefault(ccId, []).append(buildType)

                    if maxHits < 1 or not selfHit:
                        logger.info("%s (%r) buildType %r min hits %d max hits %d (%.4f seconds)", ccId, selfHit, buildType, minHits, maxHits, time.time() - startTime1)
                else:
                    logger.debug("%s missing descriptor %r", ccId, buildType)
            if ii % 100 == 0:
                logger.info("Completed %d of %d missed count %d", ii, numMols, len(missedD))
        #
        for ccId, missL in missedD.items():
            logger.info("%s missed list %r", ccId, missL)
            if ccId in missedFpD:
                logger.info("%s unmatched for fpTypes %r", ccId, missedFpD[ccId])
        # ----
        if doDisplay:
            mD = {}
            for missTup in missTupL:
                mD.setdefault(missTup[0], []).append(missTup[1])

            for ccId, buildTypeL in mD.items():
                idxD = ccIdxD[ccId]
                if "oe-iso-smiles" in idxD:
                    for buildType in buildTypeL:
                        self.__displayAlignedDescriptorPair(ccId, idxD["oe-iso-smiles"], "oe-iso-smiles", idxD[buildType], buildType, title=None, limitPerceptions=True)

        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(fpTypeCuttoffList), numMols, time.time() - startTime)
        # ----                                  ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):

    def testSubStructureSearchScreenedAbbrev(self):
        """Screened substructure search.
        """
        numMols = 3000
        myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "ccFileNamePrefix": "cc-abbrev",
            "oeFileNamePrefix": "oe-abbrev",
            "limitPerceptions": False,
            "buildTypeList": ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles"],
            "screenTypeList": ["SMARTS"],
        }
        return self.__subStructureSearchScreened(numMols, **myKwargs)

    def testSubStructureSearchScreenedFiltered(self):
        """Screened substructure search.
        """
        numMols = 5000
        myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "ccFileNamePrefix": "cc-filtered",
            "oeFileNamePrefix": "oe-filtered",
            "limitPerceptions": False,
            "buildTypeList": ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles"],
            "screenTypeList": ["SMARTS"],
        }
        return self.__subStructureSearchScreened(numMols, **myKwargs)

    def __subStructureSearchScreened(self, numMols, **kwargs):
        #
        buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles"])
        screenTypeList = kwargs.get("screenTypeList", ["SMARTS"])
        oesmP, ccIdxD = self.__getSearchDataProviders(**kwargs)
        for screenType in screenTypeList:
            oesU = OeSearchUtils(oesmP, screenType=screenType, numProc=self.__numProc)
            oeioU = OeIoUtils()
            #
            missL = []
            numMols = min(len(ccIdxD), numMols) if numMols else len(ccIdxD)
            for ii, ccId, in enumerate(list(ccIdxD.keys())[:numMols]):
                ccD = ccIdxD[ccId]
                for buildType in buildTypeList:
                    if buildType in ccD:
                        if screenType == "SMARTS":
                            oeQMol = oeioU.descriptorToMol(ccD[buildType], "SMARTS", messageTag=ccId + ":" + buildType)
                        else:
                            oeQMol = oeioU.descriptorToQMol(ccD[buildType], "SMARTS", messageTag=ccId + ":" + buildType)
                        if not oeQMol:
                            logger.debug("%s build failed for %s - skipping", ccId, buildType)
                            continue
                        # ----
                        startTime = time.time()
                        retStatus, mL = oesU.searchSubStructureScreened(oeQMol, maxMatches=100)
                        if retStatus:
                            logger.debug("%s - %s - %s (status=%r) match length %d in (%.4f seconds)", ccId, buildType, screenType, retStatus, len(mL), time.time() - startTime)
                        if not self.__resultContains(ccId, mL):
                            missL.append((ccId, buildType, screenType))
                        # ----
                if ii % 100 == 0:
                    logger.info("Completed %d of %d missed count %d", ii, numMols, len(missL))
            logger.info("Screen %r missed searches (%d) %r", screenType, len(missL), missL)

    def __displayAlignedDescriptorPair(self, ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):
        oeioU = OeIoUtils()
        oeMolRef = oeioU.descriptorToMol(descrRef, buildTypeRef, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildTypeRef)
        oeMolFit = oeioU.descriptorToMol(descrFit, buildTypeFit, limitPerceptions=limitPerceptions, messageTag=ccId + ":" + buildTypeFit)
        #
        oed = OeDepictMCSAlignPage()
        oed.setSearchType(sType="relaxed", minAtomMatchFraction=0.50)
        oed.setDisplayOptions(
            labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
        )
        oed.setRefMol(oeMolRef, ccId)
        oed.setFitMol(oeMolFit, ccId)
        myTitle = title if title else buildTypeRef + "-" + buildTypeFit
        imgPath = os.path.join(self.__workPath, myTitle + "-" + ccId + ".svg")
        logger.info("Using image path %r", imgPath)
        aML = oed.alignPair(imagePath=imgPath)
        if aML:
            logger.info("%s aligned image path %r", ccId, imgPath)
            for (rCC, rAt, tCC, tAt) in aML:
                logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)


def subStructureSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchExhaustiveAbbrev"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchWithFpAbbrev"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchScreenedAbbrev"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchScreenedFiltered"))
    return suiteSelect


def fingerprintSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchIndexUtilsTests("testFingerPrintSearchAbbrev"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testFingerPrintScoresAbbrev"))
    return suiteSelect


def fullSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchIndexUtilsTests("testSubStructureSearchWithFpFull"))
    suiteSelect.addTest(OeSearchIndexUtilsTests("testFingerPrintScoresFull"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = fullSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
