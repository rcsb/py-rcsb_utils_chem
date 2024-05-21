##
# File:    OeSubStructSearchCompareTests.py
# Author:  J. Westbrook
# Date:    26-Feb-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Match and Substructure search comparison tests on the full data including display.
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
from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage, OeDepictMCSAlignMultiPage, OeDepictSubStructureAlignMultiPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils
from rcsb.utils.chem.OeSubStructSearchUtils import OeSubStructSearchUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class OeSubStructSearchCompareTests(unittest.TestCase):
    useFull = False

    def setUp(self):

        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        self.__doDisplay = True
        self.__numProcPrep = 6
        self.__numProcSearch = 6
        self.__minCount = None
        self.__startTime = time.time()
        #
        if OeSubStructSearchCompareTests.useFull:
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
                "fpTypeCuttoffD": {"TREE": 0.6, "MACCS": 0.9},
                "maxFpResults": 50,
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
                "fpTypeCuttoffD": {"TREE": 0.6, "MACCS": 0.9},
                "maxFpResults": 50,
            }
        #
        self.__oesmP = OeSearchMoleculeProvider(**self.__myKwargs)
        ok = self.__oesmP.testCache()
        self.assertTrue(ok)
        #
        self.__ccmP = ChemCompMoleculeProvider(**self.__myKwargs)
        self.__ccmP.testCache()
        #
        self.__ccsidxP = ChemCompSearchIndexProvider(**self.__myKwargs)
        ok = self.__ccsidxP.testCache(minCount=self.__minCount)
        self.assertTrue(ok)
        self.__oessU = OeSubStructSearchUtils(self.__oesmP)
        ok = self.__oessU.testCache()
        self.assertTrue(ok)
        #
        fpTypeCuttoffD = self.__myKwargs.get("fpTypeCuttoffD", {})
        fpTypeList = [k for k, v in fpTypeCuttoffD.items()]
        self.__oesU = OeSearchUtils(self.__oesmP, fpTypeList=fpTypeList)
        ok = self.__oesU.testCache()
        self.assertTrue(ok)
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(not useFull, "Requires full data set")
    def testSubStructSearchDescriptor(self):
        #
        query = "n1ccccc1"
        queryId = "query-smiles"
        queryType = "oe-iso-smiles"
        #
        limitPerceptions = self.__myKwargs.get("limitPerceptions", False)
        suppressHydrogens = self.__myKwargs.get("suppressHydrogens", True)
        numProc = self.__myKwargs.get("numProc", 4)
        # for matchOpts in ["sub-struct-graph-relaxed", "sub-struct-graph-relaxed-stereo", "sub-struct-graph-strict"]:
        for matchOpts in ["sub-struct-graph-strict"]:
            #
            oeMol = self.__getMol(query, queryType, queryId, limitPerceptions=limitPerceptions, suppressHydrogens=suppressHydrogens)
            startTime = time.time()
            retStatus, mL = self.__search(oeMol, matchOpts, numProc)
            logger.info("%s status (%r) matchOpts %s result %d in (%.4f seconds)", queryId, retStatus, matchOpts, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            if queryType == "CC":
                self.assertTrue(self.__resultContains(queryId, mL))
            #
            if self.__doDisplay:
                self.__display(mL, query, queryId, queryType, matchOpts)

    @unittest.skipIf(not useFull, "Requires full data set")
    def testSubStructSearchSelected(self):
        #
        query = queryId = "STI"
        queryType = "CC"
        #
        limitPerceptions = self.__myKwargs.get("limitPerceptions", False)
        suppressHydrogens = self.__myKwargs.get("suppressHydrogens", True)
        numProc = self.__myKwargs.get("numProc", 4)
        for matchOpts in ["sub-struct-graph-relaxed", "sub-struct-graph-relaxed-stereo", "sub-struct-graph-strict"]:
            #
            oeMol = self.__getMol(query, queryType, queryId, limitPerceptions=limitPerceptions, suppressHydrogens=suppressHydrogens)
            startTime = time.time()
            retStatus, mL = self.__search(oeMol, matchOpts, numProc)
            logger.info("%s status (%r) matchOpts %s result %d in (%.4f seconds)", queryId, retStatus, matchOpts, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            if queryType == "CC":
                self.assertTrue(self.__resultContains(queryId, mL))
            #
            if self.__doDisplay:
                self.__display(mL, query, queryId, queryType, matchOpts)

    @unittest.skipIf(not useFull, "Requires full data set")
    def testSubStructSearchAll(self):
        #
        ccD = self.__ccmP.getMolD()
        for ccId in ccD:
            query = queryId = ccId
            if ccId in ["UNX", "UNL", "UNK", "DUM"]:
                continue
            queryType = "CC"
            #
            limitPerceptions = self.__myKwargs.get("limitPerceptions", False)
            suppressHydrogens = self.__myKwargs.get("suppressHydrogens", True)
            numProc = self.__myKwargs.get("numProc", 5)
            for matchOpts in ["sub-struct-graph-relaxed", "sub-struct-graph-relaxed-stereo", "sub-struct-graph-strict"]:
                #
                oeMol = self.__getMol(query, queryType, queryId, limitPerceptions=limitPerceptions, suppressHydrogens=suppressHydrogens)
                if oeMol.NumAtoms() < 3:
                    continue
                #
                startTime = time.time()
                retStatus, mL = self.__search(oeMol, matchOpts, numProc)
                logger.info("%s status (%r) matchOpts %s result %d in (%.4f seconds)", queryId, retStatus, matchOpts, len(mL), time.time() - startTime)
                self.assertTrue(retStatus)
                if queryType == "CC":
                    self.assertTrue(self.__resultContains(queryId, mL))
                #
                if self.__doDisplay:
                    self.__display(mL, query, queryId, queryType, matchOpts)

    @unittest.skipIf(not useFull, "Requires full data set")
    def testMatchSearchSelected(self):
        #
        query = queryId = "STI"
        queryType = "CC"
        #
        limitPerceptions = self.__myKwargs.get("limitPerceptions", False)
        suppressHydrogens = self.__myKwargs.get("suppressHydrogens", True)
        numProc = self.__myKwargs.get("numProc", 4)
        for matchOpts in ["fingerprint-similarity", "graph-relaxed", "graph-relaxed-stereo", "graph-strict"]:
            #
            oeMol = self.__getMol(query, queryType, queryId, limitPerceptions=limitPerceptions, suppressHydrogens=suppressHydrogens)
            startTime = time.time()
            retStatus, mL = self.__search(oeMol, matchOpts, numProc)
            logger.info("%s status (%r) matchOpts %s result %d in (%.4f seconds)", queryId, retStatus, matchOpts, len(mL), time.time() - startTime)
            self.assertTrue(retStatus)
            if queryType == "CC":
                self.assertTrue(self.__resultContains(queryId, mL))
            #
            if self.__doDisplay:
                self.__display(mL, query, queryId, queryType, matchOpts)

    def __search(self, oeMol, matchOpts, numProc):
        fpL = mL = None
        if matchOpts.startswith("sub-struct-"):
            retStatus, mL = self.__subStructureSearch(oeMol, matchOpts=matchOpts, numProc=numProc)
        else:
            retStatus, mL, fpL = self.__matchSearch(oeMol, matchOpts=matchOpts)
        #
        rL = fpL if matchOpts in ["fingerprint-similarity"] else mL
        return retStatus, rL

    #
    def __subStructureSearch(self, oeMol, matchOpts, numProc):
        ##
        ccIdL = self.__oessU.prefilterIndex(oeMol, self.__ccsidxP, matchOpts=matchOpts, skipFeatures=False)
        retStatus, mL = self.__oessU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
        return retStatus, mL

    def __matchSearch(self, oeMol, matchOpts="graph-relaxed"):
        ssL = fpL = []
        try:
            fpTypeCuttoffD = self.__myKwargs.get("fpTypeCuttoffD", {})
            maxFpResults = self.__myKwargs.get("maxFpResults", 50)
            retStatus, ssL, fpL = self.__oesU.searchSubStructureAndFingerPrint(oeMol, list(fpTypeCuttoffD.items())[:2], maxFpResults, matchOpts=matchOpts)
            # logger.info("fpL %r", fpL)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            #
        return retStatus, ssL, fpL

    def __getMol(self, query, queryType, queryId, limitPerceptions=False, suppressHydrogens=True):
        oeioU = OeIoUtils()
        if queryType == "CC":
            oeMol = self.__oesmP.getMol(query)
        else:
            oeMol = oeioU.descriptorToMol(query, queryType, limitPerceptions=limitPerceptions, messageTag=queryId)
        #
        if suppressHydrogens:
            oeMol = oeioU.suppressHydrogens(oeMol)
        oeMol.SetTitle(queryId)
        return oeMol

    def __resultContains(self, ccId, matchResultList):
        for matchResult in matchResultList:
            if matchResult.ccId == ccId:
                return True
        return False

    #
    #  ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
    def __display(self, mL, query, queryId, queryType, matchOpts):
        smL = sorted(mL, key=lambda kv: kv.fpScore, reverse=True)
        # ----
        tD = {}
        for sm in smL:
            ccId = sm.ccId.split("|")[0]
            tD.setdefault(ccId, []).append(sm)
        dL = []
        for ccId, ttL in tD.items():
            if len(ttL) == 1:
                dL.append(ttL[0])
            else:
                parent = False
                for tt in ttL:
                    if tt.ccId == ccId:
                        dL.append(tt)
                        parent = True
                        break
                if not parent:
                    dL.append(ttL[0])
        # ----
        pdfImagePath = os.path.join(self.__workPath, queryId + "-" + matchOpts + ".pdf")
        self.__displayPaginatedAlignments(pdfImagePath, query, queryType, queryId, dL, matchOpts=matchOpts)

    def __displayPaginatedAlignments(self, pdfImagePath, query, queryType, queryId, matchResultList, matchOpts="relaxed-stereo", alignMode="SS"):
        refId = queryId
        oeMolRef = self.__getMol(query, queryType, queryId, limitPerceptions=False, suppressHydrogens=True)
        pairList = []
        for mr in sorted(matchResultList, key=lambda kv: kv.fpScore, reverse=True):
            fitId = mr.ccId.split("|")[0]
            if len(mr.ccId) > 4:
                fitId = fitId + " (tautomer/protomer)"
            oeMolFit = self.__oesmP.getMol(mr.ccId)
            pairList.append((refId, oeMolRef, fitId, oeMolFit))
        #
        self.__depictFitList(pdfImagePath, pairList, matchOpts=matchOpts, alignMode=alignMode)

    def __pairDepictPage(self, imagePath, refId, refTitle, refMol, fitId, fitTitle, fitMol, matchOpts="strict"):
        """Depict pairwise alignment of the input reference and fit molecules.

        Args:
            imagePath (str): path to image (format by path extension)
            refId (str): reference molecule identifier
            refTitle (str): reference molecule title
            refMol (obj): reference OE molecule object
            fitId (str): fit molecule identifier
            fitTitle (str): fit molecule title
            fitMol (obj): fit OE molecule object
            matchOpts (str, optional): alignment criteria (relaxed|relaxed-stereo|strict). Defaults to "strict".

        Returns:
            (list): atom mapping in all aligned figures
                    [(reference component Id, reference atom name, fit chemical component Id, fit atom name)
        """
        aML = []
        try:
            oed = OeDepictMCSAlignPage()
            oed.setSearchType(sType=matchOpts)

            oed.setRefMol(refMol, refId, title=refTitle)
            oed.setFitMol(fitMol, fitId, title=fitTitle)
            oed.setDisplayOptions(
                imageSizeX=2000,
                imageSizeY=1000,
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                bondDisplayWidth=0.5,
                highLightMatchColorRef="green",
                highLightNotMatchColorRef="pink",
            )
            aML = oed.alignPair(imagePath=imagePath)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aML

    def __depictFitList(self, pdfImagePath, pairList, matchOpts="exact", alignMode="SS"):
        """Depict pairwise alignments with multi-page layout in PDF format.

        Args:
            pdfImagePath (str): PDF image path
            pairList (list): [(refId, refOeMol, fitId, fitOeMol)]

        Returns:
            (list): atom mapping in all aligned figures
                    [(reference component Id, reference atom name, fit chemical component Id, fit atom name)
        """
        aML = []
        try:
            if alignMode == "MCSS":
                oed = OeDepictMCSAlignMultiPage()
            else:
                oed = OeDepictSubStructureAlignMultiPage()
            oed.setSearchType(sType=matchOpts)
            oed.setPairMolList(pairList)

            oed.setDisplayOptions(
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                pageOrientation="portrait",
                gridRows=4,
                bondDisplayWidth=0.5,
                highLightMatchColorRef="green",
                highLightNotMatchColorRef="pink",
            )
            aML = oed.alignPairListMulti(imagePath=pdfImagePath)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aML


def searchComparison():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeSubStructSearchCompareTests("testSubStructureSearchFromIndex"))
    suiteSelect.addTest(OeSubStructSearchCompareTests("testSubStructureSearchBase"))
    suiteSelect.addTest(OeSubStructSearchCompareTests("testSubStructureSearchFromIndexBase"))
    suiteSelect.addTest(OeSubStructSearchCompareTests("testSubStructureSearchFromIndexSelected"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = searchComparison()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
