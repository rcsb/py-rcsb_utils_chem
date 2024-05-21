##
# File:    OeSearchUtilsTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for various search modes for core PDB chemical component definitions.

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
from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSearchUtilsTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        self.__fpTypeList = ["TREE", "PATH", "MACCS", "CIRCULAR", "LINGO"]
        self.__fpTypeCuttoffList = [("TREE", 0.6), ("PATH", 0.6), ("MACCS", 0.9), ("CIRCULAR", 0.6), ("LINGO", 0.9)]
        self.__screenType = "SMARTS"
        self.__numProc = 1
        self.__minCount = 500
        self.__startTime = time.time()
        self.__myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "fpTypeList": self.__fpTypeList,
            "ccFileNamePrefix": "cc-abbrev",
            "oeFileNamePrefix": "oe-abbrev",
            "molBuildType": "model-xyz",
            "limitPerceptions": False,
        }
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __resultContains(self, ccId, matchResultList):
        for matchResult in matchResultList:
            if matchResult.ccId == ccId:
                return True
        return False

    @unittest.skipIf(skipFlag, "deprecated test")
    def testSubStructureSearch(self):
        oemp = OeMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        oesU = OeSearchUtils(oemp, fpTypeList=self.__fpTypeList)
        numMols = 10
        for ccId, _ in list(ccIdxD.items())[:numMols]:
            # ----
            startTime = time.time()
            oeMol = oemp.getMol(ccId)
            retStatus, mL = oesU.searchSubStructure(oeMol, matchOpts="relaxed")
            logger.info("%s match length %d in (%.4f seconds)", ccId, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            self.assertTrue(self.__resultContains(ccId, mL))
            # self.assertGreaterEqual(len(mL), 1)
            # ----

    @unittest.skipIf(skipFlag, "deprecated test")
    def testFingerPrintSearch(self):
        oemp = OeMoleculeProvider(**self.__myKwargs)
        # This will reload the oe binary cache.
        oeMol = oemp.getMol("004")
        self.assertGreaterEqual(len(list(oeMol.GetAtoms())), 12)
        #
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        minScore = 0.50
        maxResults = 50
        numMols = 50
        oesU = OeSearchUtils(oemp, fpTypeList=self.__fpTypeList)
        # ----
        startTime = time.time()
        for ccId, _ in list(ccIdxD.items())[:numMols]:
            for fpType in self.__fpTypeList:
                oeMol = oemp.getMol(ccId)
                retStatus, mL = oesU.searchFingerPrints(oeMol, fpType=fpType, minFpScore=minScore, maxFpResults=maxResults)
                self.assertTrue(retStatus)
                self.assertTrue(self.__resultContains(ccId, mL))
                # self.assertGreaterEqual(len(mL), 1)
        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(self.__fpTypeList), numMols, time.time() - startTime)
        # ----

    @unittest.skipIf(skipFlag, "deprecated test")
    def testFingerPrintScores(self):
        oemp = OeMoleculeProvider(**self.__myKwargs)
        #
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        limitPerceptions = False
        maxResults = 100
        numMols = 20
        oeioU = OeIoUtils()
        oesU = OeSearchUtils(oemp, fpTypeList=self.__fpTypeList)
        missedFpD = {}
        missedBuildD = {}
        # ----
        startTime = time.time()
        for ccId, ccD in list(ccIdxD.items())[:numMols]:
            for buildType in ["oe-iso-smiles", "oe-smiles", "acdlabs-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]:
                if buildType in ccD:
                    oeMol = None
                    logger.debug("Search %s %r", ccId, ccD[buildType])
                    if buildType in ["inchi"]:
                        # oeMol = oeioU.inchiToMol(ccD[buildType], limitPerceptions=limitPerceptions)
                        #
                        oemf = OeMoleculeFactory()
                        oemf.setDescriptor(ccD["inchi"], "inchi", ccId)
                        ok = oemf.build(molBuildType="inchi", limitPerceptions=limitPerceptions)
                        if not ok:
                            logger.info("%s build failed with InChI %r", ccId, ccD["inchi"])
                        else:
                            oeMol = oemf.getMol()
                            if oemf.getInChI() != ccD["inchi"]:
                                logger.info("%s regenerated InChI differs\n%r\n%s", ccId, ccD["inchi"], oemf.getInChI())
                        #
                    else:
                        oeMol = oeioU.smilesToMol(ccD[buildType], limitPerceptions=limitPerceptions)
                    if not oeMol:
                        continue
                    maxHits = 0
                    minHits = maxResults
                    selfHit = False
                    #
                    for fpType, minFpScore in self.__fpTypeCuttoffList:
                        retStatus, mL = oesU.getFingerPrintScores(oeMol, fpType, minFpScore, maxResults)
                        self.assertTrue(retStatus)
                        logger.info("%s fpType %r hits %d", ccId, fpType, len(mL))
                        maxHits = max(maxHits, len(mL))
                        minHits = min(minHits, len(mL))
                        matchedSelf = self.__resultContains(ccId, mL)
                        selfHit = selfHit or matchedSelf
                        if not matchedSelf:
                            missedFpD.setdefault(ccId, []).append((buildType, fpType, len(mL)))
                    #
                    if not selfHit:
                        missedBuildD.setdefault(ccId, []).append(buildType)

                    logger.info("%s buildType %r min hits %d max hits %d", ccId, buildType, minHits, maxHits)
                else:
                    logger.info("%s missing descriptor %r", ccId, buildType)

        #
        for ccId, bTL in missedBuildD.items():
            logger.info("%s missed build type list %r", ccId, bTL)

        if ccId in missedFpD:
            logger.info("%s unmatched fpTypes %r", ccId, missedFpD[ccId])

        #
        doDepict = False
        if doDepict:
            for ccId, bTL in missedBuildD.items():
                idxD = ccIdxD[ccId]
                if "oe-iso-smiles" in idxD:
                    for bT in bTL:
                        self.__displayAlignedDescriptorPair(ccId, idxD["oe-iso-smiles"], "oe-iso-smiles", idxD[bT], bT, title=None, limitPerceptions=True)

        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(self.__fpTypeList), numMols, time.time() - startTime)
        # ----                                  ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):

    @unittest.skipIf(skipFlag, "deprecated test")
    def testSssWithFingerPrintFromDescriptor(self):
        oemp = OeMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        limitPerceptions = False
        # minFpScore = 0.5
        maxFpResults = 50
        matchOpts = "graph-relaxed"
        numMols = 20
        oeioU = OeIoUtils()
        oesU = OeSearchUtils(oemp, fpTypeList=self.__fpTypeList)
        missTupL = []
        missedD = {}
        missedFpD = {}
        # ----
        startTime = time.time()
        for ccId, ccD in list(ccIdxD.items())[:numMols]:
            for buildType in ["oe-iso-smiles", "oe-smiles", "acdlabs-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"]:
                if buildType in ccD:
                    oeMol = None
                    logger.debug("Search %s %r", ccId, ccD[buildType])
                    if buildType in ["inchi"]:
                        oemf = OeMoleculeFactory()
                        oemf.setDescriptor(ccD["inchi"], "inchi", ccId)
                        ok = oemf.build(molBuildType="inchi", limitPerceptions=limitPerceptions)
                        if not ok:
                            logger.info("%s build failed with InChI %r", ccId, ccD["inchi"])
                        else:
                            oeMol = oemf.getMol()
                            if oemf.getInChI() != ccD["inchi"]:
                                logger.info("%s regenerated InChI differs\n%r\n%s", ccId, ccD["inchi"], oemf.getInChI())
                    else:
                        oeMol = oeioU.smilesToMol(ccD[buildType], limitPerceptions=limitPerceptions)
                    if not oeMol:
                        continue
                    maxHits = 0
                    minHits = maxFpResults
                    selfHit = False
                    for fpType, minFpScore in self.__fpTypeCuttoffList:
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

                    logger.info("%s (%r) buildType %r min hits %d max hits %d", ccId, selfHit, buildType, minHits, maxHits)
                else:
                    logger.info("%s missing descriptor %r", ccId, buildType)
        #
        for ccId, missL in missedD.items():
            logger.info("%s missed list %r", ccId, missL)
            if ccId in missedFpD:
                logger.info("%s unmatched for fpTypes %r", ccId, missedFpD[ccId])
        # ----
        doDepict = False
        if doDepict:
            mD = {}
            for missTup in missTupL:
                mD.setdefault(missTup[0], []).append(missTup[1])

            for ccId, buildTypeL in mD.items():
                idxD = ccIdxD[ccId]
                if "oe-iso-smiles" in idxD:
                    for buildType in buildTypeL:
                        self.__displayAlignedDescriptorPair(ccId, idxD["oe-iso-smiles"], "oe-iso-smiles", idxD[buildType], buildType, title=None, limitPerceptions=True)

        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(self.__fpTypeList), numMols, time.time() - startTime)
        # ----                                  ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):

    @unittest.skipIf(skipFlag, "deprecated test")
    def testSubStructureSearchWithFingerPrint(self):
        oemp = OeMoleculeProvider(**self.__myKwargs)
        #
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        minFpScore = 0.40
        maxFpResults = 50
        numMols = 20
        matchOpts = "graph-relaxed"
        oesU = OeSearchUtils(oemp, fpTypeList=self.__fpTypeList)
        # ----
        startTime = time.time()
        for ccId, _ in list(ccIdxD.items())[:numMols]:
            for fpType in self.__fpTypeList:
                oeMol = oemp.getMol(ccId)
                retStatus, mL = oesU.searchSubStructureWithFingerPrint(oeMol, fpType, minFpScore, maxFpResults, matchOpts=matchOpts)
                self.assertTrue(retStatus)
                self.assertTrue(self.__resultContains(ccId, mL))

        logger.info("%s fingerprints search on %d in (%.4f seconds)", len(self.__fpTypeList), numMols, time.time() - startTime)
        # ----

    @unittest.skipIf(skipFlag, "deprecated test")
    def testSubStructureSearchScreened(self):
        oeioU = OeIoUtils()
        oemp = OeMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**self.__myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        oesU = OeSearchUtils(oemp, screenType=self.__screenType, numProc=self.__numProc)
        numMols = 20
        missL = []
        for ccId, ccD in list(ccIdxD.items())[:numMols]:
            # ----
            startTime = time.time()
            if "oe-smiles" not in ccD:
                continue
            logger.info("Search %s %r", ccId, ccD["oe-smiles"])
            oeQMol = oeioU.smartsToQmol(ccD["oe-smiles"])
            retStatus, mL = oesU.searchSubStructureScreened(oeQMol, maxMatches=100)
            if retStatus:
                logger.info("%s (status=%r) match length %d in (%.4f seconds)", ccId, retStatus, len(mL), time.time() - startTime)
            if not self.__resultContains(ccId, mL):
                missL.append(ccId)
            #
            # self.assertGreaterEqual(len(mL), 1)
            # ----
        logger.info("Missed searches (%d) %r", len(missL), missL)

    @unittest.skipIf(skipFlag, "Long troubleshooting test")
    def testSubStructureSearchScreenedFiltered(self):
        myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            "fpTypeList": self.__fpTypeList,
            "ccFileNamePrefix": "cc-filtered",
            "oeFileNamePrefix": "oe-filtered",
            "molBuildType": "oe-iso-smiles",
            "limitPerceptions": False,
        }
        oeioU = OeIoUtils()
        oemp = OeMoleculeProvider(**myKwargs)
        ok = oemp.testCache()
        ccmP = ChemCompIndexProvider(**myKwargs)
        ccIdxD = ccmP.getIndex()
        ok = ccmP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        oesU = OeSearchUtils(oemp, screenType=self.__screenType, numProc=self.__numProc)
        numMols = 5000
        missL = []
        for ccId, ccD in list(ccIdxD.items())[:numMols]:
            # ----
            startTime = time.time()
            if "oe-smiles" not in ccD:
                continue
            logger.info("Search %s %r", ccId, ccD["oe-smiles"])
            oeQMol = oeioU.smartsToQmol(ccD["oe-smiles"])
            retStatus, mL = oesU.searchSubStructureScreened(oeQMol, maxMatches=100)
            logger.info("%s (status=%r)match length %d in (%.4f seconds)", ccId, retStatus, len(mL), time.time() - startTime)
            if not self.__resultContains(ccId, mL):
                missL.append(ccId)

            # self.assertGreaterEqual(len(mL), 1)
            # ----
        logger.info("Missed searches (%d) %r", len(missL), missL)

    def __displayAlignedDescriptorPair(self, ccId, descrRef, buildTypeRef, descrFit, buildTypeFit, title=None, limitPerceptions=True):
        oemfRef = OeMoleculeFactory()
        oemfRef.setDescriptor(descrRef, buildTypeRef, ccId)
        oemfRef.build(molBuildType=buildTypeRef, limitPerceptions=limitPerceptions)
        oeMolRef = oemfRef.getMol()
        #
        oemfFit = OeMoleculeFactory()
        oemfFit.setDescriptor(descrFit, buildTypeFit, ccId)
        oemfFit.build(molBuildType=buildTypeFit, limitPerceptions=limitPerceptions)
        oeMolFit = oemfFit.getMol()
        #
        oed = OeDepictMCSAlignPage()
        oed.setSearchType(sType="graph-relaxed", minAtomMatchFraction=0.50)
        oed.setDisplayOptions(labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5)
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


def rawSubStructureSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchUtilsTests("testSubStructureSearch"))
    suiteSelect.addTest(OeSearchUtilsTests("testSubStructureSearchWithFingerPrint"))
    suiteSelect.addTest(OeSearchUtilsTests("testSubStructureSearchScreened"))
    return suiteSelect


def fingerprintSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSearchUtilsTests("testFingerPrintSearch"))
    suiteSelect.addTest(OeSearchUtilsTests("testFingerPrintScores"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = rawSubStructureSearch()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
