##
# File:    OeSubStructSearchUtilsTests.py
# Author:  J. Westbrook
# Date:    2-Oct-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for substructure search modes on core and search index PDB chemical component definitions.

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
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.chem.OeSubStructSearchUtils import OeSubStructSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSubStructSearchUtilsTests(unittest.TestCase):
    useFull = False

    def setUp(self):

        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        self.__numProcPrep = 8
        self.__numProcSearch = 8
        self.__minCount = None
        self.__startTime = time.time()
        #
        if OeSubStructSearchUtilsTests.useFull:
            self.__myKwargs = {
                "cachePath": self.__cachePath,
                "useCache": True,
                "ccFileNamePrefix": "cc-full",
                "oeFileNamePrefix": "oe-full",
                "molBuildType": "model-xyz",
                "limitPerceptions": False,
                "screenTypeList": None,
                "numProc": self.__numProcPrep,
                "suppressHydrogens": True,
                "matchOpts": "sub-struct-graph-relaxed",
            }
        else:
            self.__myKwargs = {
                "ccUrlTarget": self.__ccUrlTarget,
                "birdUrlTarget": self.__birdUrlTarget,
                "cachePath": self.__cachePath,
                "useCache": True,
                "ccFileNamePrefix": "cc-abbrev",
                "oeFileNamePrefix": "oe-abbrev",
                "molBuildType": "model-xyz",
                "limitPerceptions": False,
                "screenTypeList": None,
                "numProc": self.__numProcPrep,
                "suppressHydrogens": True,
                "matchOpts": "sub-struct-graph-relaxed",
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

    @unittest.skipIf(not useFull, "Requires full data set")
    def testSubStructureSearchFromIndexSelected(self):
        matchOpts = self.__myKwargs.get("matchOpts", "sub-struct-graph-relaxed")
        numProc = self.__numProcSearch
        oemp = OeSearchMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        self.assertTrue(ok)
        oesU = OeSubStructSearchUtils(oemp)
        #
        ccIdxP = ChemCompSearchIndexProvider(**self.__myKwargs)
        ok = ccIdxP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        ccIdxD = ccIdxP.getIndex()
        ky = next(iter(ccIdxD))
        oeMol = oemp.getMol(ky)
        #
        for ccId in ["BNZ", "ALA"]:
            # ----
            startTime = time.time()
            oeMol = oemp.getMol(ccId)
            #
            ccIdL = oesU.prefilterIndex(oeMol, ccIdxP, matchOpts=matchOpts)
            logger.info("%s search length %d in (%.4f seconds)", ccId, len(ccIdL), time.time() - startTime)
            #
            retStatus, mL = oesU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
            logger.info("%s status %r result length %d in (%.4f seconds)", ccId, retStatus, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            self.assertTrue(self.__resultContains(ccId, mL))
            # ----

    def testSubStructureSearchFromIndexBase(self):
        matchOpts = self.__myKwargs.get("matchOpts", "sub-struct-graph-relaxed")
        numProc = self.__numProcSearch
        oemp = OeSearchMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        self.assertTrue(ok)
        oesU = OeSubStructSearchUtils(oemp)
        #
        ccIdxP = ChemCompSearchIndexProvider(**self.__myKwargs)
        ok = ccIdxP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        ccIdxD = ccIdxP.getIndex()
        ky = next(iter(ccIdxD))
        oeMol = oemp.getMol(ky)

        numMols = 10
        for ccId, _ in list(ccIdxD.items())[:numMols]:
            # ----
            startTime = time.time()
            oeMol = oemp.getMol(ccId)
            #
            ccIdL = oesU.prefilterIndex(oeMol, ccIdxP, matchOpts=matchOpts)
            logger.info("%s search length %d in (%.4f seconds)", ccId, len(ccIdL), time.time() - startTime)
            #
            retStatus, mL = oesU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
            logger.info("%s status %r result length %d in (%.4f seconds)", ccId, retStatus, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            self.assertTrue(self.__resultContains(ccId, mL))
            # ----

    @unittest.skipIf(not useFull, "Requires full data set")
    def testSubStructureSearchBaseSelected(self):
        matchOpts = self.__myKwargs.get("matchOpts", "sub-struct-graph-relaxed")
        numProc = self.__numProcSearch
        oemp = OeMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        self.assertTrue(ok)
        oesU = OeSubStructSearchUtils(oemp)
        #
        #
        ccIdxP = ChemCompIndexProvider(**self.__myKwargs)
        ok = ccIdxP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        ccIdxD = ccIdxP.getIndex()
        #
        ky = next(iter(ccIdxD))
        oeMol = oemp.getMol(ky)
        #
        for ccId in ["BNZ", "ALA"]:
            # ----
            startTime = time.time()
            oeMol = oemp.getMol(ccId)
            #
            ccIdL = oesU.prefilterIndex(oeMol, ccIdxP, matchOpts=matchOpts)

            logger.info("%s search length %d in (%.4f seconds)", ccId, len(ccIdL), time.time() - startTime)
            #
            retStatus, mL = oesU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
            logger.info("%s result length %d in (%.4f seconds)", ccId, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            self.assertTrue(self.__resultContains(ccId, mL))
            # ----

    def testSubStructureSearchBase(self):

        matchOpts = self.__myKwargs.get("matchOpts", "sub-struct-graph-relaxed")
        numProc = self.__numProcSearch
        oemp = OeMoleculeProvider(**self.__myKwargs)
        ok = oemp.testCache()
        self.assertTrue(ok)
        oesU = OeSubStructSearchUtils(oemp)
        #
        ccIdxP = ChemCompIndexProvider(**self.__myKwargs)
        ok = ccIdxP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        ccIdxD = ccIdxP.getIndex()
        #
        ky = next(iter(ccIdxD))
        oeMol = oemp.getMol(ky)
        #
        numMols = 10
        for ccId, _ in list(ccIdxD.items())[:numMols]:
            # ----
            startTime = time.time()
            oeMol = oemp.getMol(ccId)
            ccIdL = oesU.prefilterIndex(oeMol, ccIdxP, matchOpts=matchOpts)
            logger.info("%s search length %d in (%.4f seconds)", ccId, len(ccIdL), time.time() - startTime)
            #
            retStatus, mL = oesU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
            logger.info("%s result length %d in (%.4f seconds)", ccId, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            self.assertTrue(self.__resultContains(ccId, mL))
            # ----

    #  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

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


def exhaustiveSubStructureSearch():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSubStructSearchUtilsTests("testSubStructureSearchFromIndex"))
    suiteSelect.addTest(OeSubStructSearchUtilsTests("testSubStructureSearchBase"))
    suiteSelect.addTest(OeSubStructSearchUtilsTests("testSubStructureSearchFromIndexBase"))
    suiteSelect.addTest(OeSubStructSearchUtilsTests("testSubStructureSearchFromIndexSelected"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = exhaustiveSubStructureSearch()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
