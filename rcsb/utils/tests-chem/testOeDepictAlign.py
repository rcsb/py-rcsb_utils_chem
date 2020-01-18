##
# File:    OeDepictAlignTests.py
# Author:  jdw
# Date:    28-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests for the OEDepictAlign and related classes which perform
MCSS comparison and aligned depiction.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlign, OeDepictMCSAlignMultiPage, OeDepictMCSAlignPage, OeDepictSubStructureAlign, OeMCSAlignUtil
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeDepictAlignTests(unittest.TestCase):
    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        #
        self.__refId = "C"
        #
        self.__idList = ["cg1", "atp", "gtp", "A", "C", "G", "DG"]
        self.__pairIdList = [("c", "cg1"), ("c", "atp"), ("c", "gtp"), ("c", "A"), ("c", "C"), ("c", "G"), ("c", "DG")]
        self.__ssPairIdList = [("A", "A"), ("C", "C"), ("G", "G"), ("ALA", "TRP")]
        #
        self.__rnaPairIdList = [
            ("ZAD", "A"),
            ("1MA", "A"),
            ("MA6", "A"),
            ("SRA", "A"),
            ("12A", "A"),
            ("2MA", "A"),
            ("LCA", "A"),
            ("T6A", "A"),
            ("MAD", "A"),
            ("RIA", "A"),
            ("A3P", "A"),
            ("A2L", "A"),
            ("A23", "A"),
            ("8AN", "A"),
            ("PU", "A"),
            ("6IA", "A"),
            ("P5P", "A"),
            ("PPU", "A"),
            ("5FA", "A"),
            ("A2M", "A"),
            ("MTU", "A"),
            ("A44", "A"),
            ("AVC", "A"),
            ("MGQ", "A"),
            ("AP7", "A"),
            ("MIA", "A"),
            ("AET", "A"),
            ("A5O", "A"),
            ("N5M", "C"),
            ("OMC", "C"),
            ("10C", "C"),
            ("ZCY", "C"),
            ("ZBC", "C"),
            ("PMT", "C"),
            ("M5M", "C"),
            ("IC", "C"),
            ("S4C", "C"),
            ("M4C", "C"),
            ("CBV", "C"),
            ("CH", "C"),
            ("CCC", "C"),
            ("CSF", "C"),
            ("C43", "C"),
            ("C31", "C"),
            ("C2L", "C"),
            ("1SC", "C"),
            ("A5M", "C"),
            ("5IC", "C"),
            ("5MC", "C"),
            ("4OC", "C"),
            ("TPG", "G"),
            ("OMG", "G"),
            ("1MG", "G"),
            ("YYG", "G"),
            ("YG", "G"),
            ("XTS", "G"),
            ("PGP", "G"),
            ("23G", "G"),
            ("G2L", "G"),
            ("2MG", "G"),
            ("GRB", "G"),
            ("QUO", "G"),
            ("ZGU", "G"),
            ("M2G", "G"),
            ("CG1", "G"),
            ("KAG", "G"),
            ("MGV", "G"),
            ("IG", "G"),
            ("GTP", "G"),
            ("GOM", "G"),
            ("GH3", "G"),
            ("7MG", "G"),
            ("G7M", "G"),
            ("G46", "G"),
            ("G48", "G"),
            ("GDP", "G"),
            ("G25", "G"),
            ("O2G", "G"),
            ("GAO", "G"),
            ("N6G", "G"),
            ("T41", "T"),
            ("T38", "T"),
            ("T2S", "T"),
            ("T23", "T"),
            ("SMT", "T"),
            ("T39", "T"),
            ("ZTH", "T"),
            ("UAR", "U"),
            ("UR3", "U"),
            ("UMP", "U"),
            ("4SU", "U"),
            ("UD5", "U"),
            ("U8U", "U"),
            ("U37", "U"),
            ("U36", "U"),
            ("URD", "U"),
            ("US5", "U"),
            ("5MU", "U"),
            ("125", "U"),
            ("OMU", "U"),
            ("126", "U"),
            ("ZBU", "U"),
            ("127", "U"),
            ("ONE", "U"),
            ("5BU", "U"),
            ("70U", "U"),
            ("U34", "U"),
            ("U31", "U"),
            ("5FU", "U"),
            ("2OM", "U"),
            ("CNU", "U"),
            ("RSQ", "U"),
            ("RUS", "U"),
            ("SUR", "U"),
            ("SSU", "U"),
            ("2MU", "U"),
            ("3AU", "U"),
            ("PYO", "U"),
            ("2AU", "U"),
            ("U2P", "U"),
            ("U2L", "U"),
            ("IU", "U"),
            ("PSU", "U"),
            ("3MU", "U"),
            ("FHU", "U"),
            ("MNU", "U"),
            ("H2U", "U"),
            ("3TD", "U"),
        ]
        #
        self.__oeMolD = self.__getCache(coordType="model", useCache=True)

    def tearDown(self):
        pass

    def __getCache(self, coordType="model", useCache=True):
        oemp = OeMoleculeProvider(dirPath=os.path.join(self.__cachePath, "chem_comp"), screenTypeList=[], coordType=coordType, useCache=useCache)
        ok = oemp.testCache()
        self.assertTrue(ok)
        return oemp.getOeMolD()

    def __getMolDepictList(self, ccIdList):
        dL = []
        for ccId in ccIdList:
            ccId = ccId.upper()
            if ccId in self.__oeMolD:
                dL.append((self.__oeMolD[ccId], ccId, "Title: " + ccId))
            else:
                logger.info("Missing molecule %r", ccId)
        return dL

    def testMCSAlignPairDepict(self):
        """ Test case -  Simple pairwise MCSS alignment  -  Each aligned pair output to a separate image file
        """
        try:
            oed = OeDepictMCSAlignPage()
            oed.setDisplayOptions(
                labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
            )

            oed.setRefMol(self.__oeMolD[self.__refId], self.__refId)
            fitTupList = self.__getMolDepictList(self.__idList)

            for fitTup in fitTupList:
                oed.setFitMol(fitTup[0], fitTup[1])
                imgPath = os.path.join(self.__workPath, "ref-" + self.__refId + "-trg-" + fitTup[1] + ".svg")
                logger.info("Using image path %r", imgPath)
                aML = oed.alignPair(imagePath=imgPath)
                self.assertGreater(len(aML), 2)
                if aML:
                    for (rCC, rAt, tCC, tAt) in aML:
                        logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSRelaxAlignPairDepict(self):
        """Test case -  Relaxed pairwise MCSS alignment
        """
        try:
            oed = OeDepictMCSAlignPage()
            oed.setSearchType(sType="relaxed")
            oed.setDisplayOptions(
                labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
            )

            oed.setRefMol(self.__oeMolD[self.__refId], self.__refId)
            fitTupList = self.__getMolDepictList(self.__idList)
            for fitTup in fitTupList:
                oed.setFitMol(fitTup[0], fitTup[1])
                fName = os.path.join(self.__workPath, "relaxed-ref-" + self.__refId + "-trg-" + fitTup[1] + ".svg")

                oed.setDisplayOptions(
                    imageSizeX=2000,
                    imageSizeY=1000,
                    labelAtomName=True,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    highlightStyleFit="ballAndStickInverse",
                    bondDisplayWidth=1.0,
                )
                aML = oed.alignPair(imagePath=fName)
                self.assertGreater(len(aML), 2)
                if aML:
                    for (rCC, rAt, tCC, tAt) in aML:
                        logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignFitListDepictMulti(self):
        """Test case -  List view of pairwise MCS alignment - multipage output
        """
        try:
            oed = OeDepictMCSAlignMultiPage()
            oed.setRefMol(self.__oeMolD[self.__refId.upper()], self.__refId)
            fitMolList = [self.__oeMolD[ccId.upper()] for ccId in self.__idList]
            oed.addFitMolList(fitMolList, suppressHydrogens=False, imageDirPath=self.__workPath)
            oed.setDisplayOptions(
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                gridRows=3,
                gridCols=3,
                bondDisplayWidth=0.5,
            )
            imagePath = os.path.join(self.__workPath, "list-example-mcs-alignment.pdf")
            aML = oed.alignOneWithListMulti(imagePath=imagePath)
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignPairListDepictPortrait(self):
        """Test case -  List view MCS alignment using pair id list input
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__pairIdList]
            oed = OeDepictMCSAlignMultiPage()
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
            )

            imagePath = os.path.join(self.__workPath, "pair-list-example-mcs-alignment-portrait.pdf")
            aML = oed.alignPairListMulti(imagePath=imagePath)
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignRnaPairListDepict(self):
        """Test case -  Modified RNA nucleotide alignment with parent nucleotied using pair list input
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__rnaPairIdList]
            oed = OeDepictMCSAlignMultiPage()
            oed.setPairMolList(pairList)
            oed.setDisplayOptions(
                labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
            )

            imagePath = os.path.join(self.__workPath, "rna-modified-pair-alignment.pdf")
            aML = oed.alignPairListMulti(imagePath=imagePath)
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignAtomMap(self):
        """Test case -  match test with return of atom maps
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__rnaPairIdList]
            oed = OeMCSAlignUtil()
            for refId, refMol, fitId, fitMol in pairList:
                oed.setRefMol(refMol, refId)
                oed.setFitMol(fitMol, fitId)
                aML = oed.doAlign()
                self.assertGreater(len(aML), 2)
                if aML:
                    logger.debug("Match length %3d for: %s %s", len(aML), refId, fitId)
                else:
                    logger.info("Match failed for: %s %s", refId, fitId)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignListMultiDepict(self):
        """Test case -  List view of MCS alignment --- on multi-pages --
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__pairIdList]
            oed = OeDepictMCSAlignMultiPage()
            oed.setDisplayOptions(
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                gridRows=3,
                gridCols=3,
                bondDisplayWidth=0.5,
            )
            oed.setPairMolList(pairList)
            imageFile = os.path.join(self.__workPath, "mcs-align-with-list-multi.pdf")
            aML = oed.alignOneWithListMulti(imagePath=imageFile)
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignListDepict(self):
        """Test case -  List view of MCS alignment on single image -
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__pairIdList]
            oed = OeDepictMCSAlignPage()
            oed.setDisplayOptions(
                imageSizeX=1500,
                imageSizeY=1500,
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                gridRows=3,
                gridCols=3,
                bondDisplayWidth=1.0,
            )
            oed.setPairMolList(pairList)
            imageFile = os.path.join(self.__workPath, "mcs-align-with-list-single.svg")
            aML = oed.alignOneWithList(imagePath=imageFile)
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignListDepictSingle(self):
        """Test case -  View MCS alignment in multiple single image files.
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__pairIdList]
            oed = OeDepictMCSAlign()
            oed.setDisplayOptions(
                imageSizeX=500,
                imageSizeY=500,
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                bondDisplayWidth=1.0,
            )
            oed.setPairMolList(pairList, imageDirPath=self.__workPath, imageFilePrefix="single")
            aML = oed.alignOneWithList()
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMCSAlignPathListDepictSingle(self):
        """Test case -  View MCS alignments in single image files.

                        Input details specified chemical component and image file paths.

                        Image file paths must end with a recognized image format (e.g. svg, png, jpg)
        """
        try:
            # Use inverse  highlighting of matching/non-matching atoms/bonds -
            hOpt = "inverse"
            #
            oed = OeDepictMCSAlign()
            if hOpt == "inverse":
                oed.setDisplayOptions(
                    imageSizeX=500,
                    imageSizeY=500,
                    labelAtomName=True,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    highlightStyleFit="ballAndStickInverse",
                    bondDisplayWidth=1.0,
                )
            elif hOpt == "match":
                oed.setDisplayOptions(
                    imageSizeX=500,
                    imageSizeY=500,
                    labelAtomName=True,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    highlightStyleFit="ballAndStick",
                    bondDisplayWidth=1.0,
                )

            else:
                oed.setDisplayOptions(
                    imageSizeX=500,
                    imageSizeY=500,
                    labelAtomName=True,
                    labelAtomCIPStereo=True,
                    labelAtomIndex=False,
                    labelBondIndex=False,
                    highlightStyleFit="None",
                    bondDisplayWidth=1.0,
                )
            oed.setRefMol(self.__oeMolD[self.__refId], self.__refId)
            fitMolList = [self.__oeMolD[ccId.upper()] for ccId in self.__idList]
            oed.addFitMolList(fitMolList, suppressHydrogens=False, imageDirPath=self.__workPath, imageFilePrefix=hOpt)

            aML = oed.alignOneWithList()
            self.assertGreater(len(aML), 2)
            #
            # Write out atom correspondences --
            #
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSSAlignListDepictSingle(self):
        """Test case -  View substructure alignment in multiple single image files.
        """
        try:
            pairList = [(tup[0].upper(), self.__oeMolD[tup[0].upper()], tup[1].upper(), self.__oeMolD[tup[1].upper()]) for tup in self.__ssPairIdList]
            oed = OeDepictSubStructureAlign()
            oed.setDisplayOptions(
                imageSizeX=500,
                imageSizeY=500,
                labelAtomName=False,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStick",
                highLightMatchColorFit="green",
                bondDisplayWidth=0.5,
            )
            oed.setPairMolList(pairList, imageDirPath=self.__workPath, imageFilePrefix="single")
            aML = oed.alignOneWithList()
            self.assertGreater(len(aML), 2)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteAlignTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignAtomMap"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSRelaxAlignPairDepict"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignPairDepict"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignFitListDepictMulti"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignPairListDepictPortrait"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignRnaPairListDepict"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignListMultiDepict"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignListDepict"))
    suiteSelect.addTest(OeDepictAlignTests("testMCSAlignPathListDepictSingle"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteAlignTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
