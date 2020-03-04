##
# File:    OeMolecularFactoryTests.py
# Author:  jdw
# Date:    28-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests of OeMolecularFactory to compare assigned and computed features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import unittest

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.MoleculeAnnotationsCompare import MoleculeAnnotationsCompare
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompIt

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeMolecularFactoryTests(unittest.TestCase):
    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        # self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        self.__molLimit = 50

    def tearDown(self):
        pass

    def __getChemCompDefs(self, molLimit=500):
        ccMolD = {}
        try:
            useCache = True
            ccFileNamePrefix = "cc-abbrev"
            ccmP = ChemCompMoleculeProvider(
                ccUrlTarget=self.__ccUrlTarget,
                birdUrlTarget=self.__birdUrlTarget,
                cachePath=self.__cachePath,
                useCache=useCache,
                ccFileNamePrefix=ccFileNamePrefix,
                molLimit=molLimit,
            )
            ok = ccmP.testCache(minCount=molLimit)
            self.assertTrue(ok)
            ccMolD = ccmP.getMolD()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ccMolD

    def testBuildRelated(self):
        """Compare constructed molecule with underlying chemical definitions -
        """
        try:
            ccMolD = self.__getChemCompDefs()
            oemf = OeMoleculeFactory()
            for ccId, ccObj in list(ccMolD.items())[:10]:
                # ----
                tId = oemf.setChemCompDef(ccObj)
                self.assertEqual(tId, ccId)
                smiD = oemf.buildRelated(limitPerceptions=False)
                logger.info("%s related molecular forms %d", ccId, len(smiD))

                # ----
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSelfConsistency(self):
        """Compare constructed molecule with underlying chemical definitions -
        """
        try:
            failL = []
            ccMolD = self.__getChemCompDefs()
            #
            # molBuildTypeL = ["model-xyz", "ideal-xyz", None]
            # molBuildTypeL = [None]
            #
            oemf = OeMoleculeFactory()
            macmp = MoleculeAnnotationsCompare()

            limitPerceptions = False
            # buildTypeRef = "oe-iso-smiles"
            buildTypeRef = "model-xyz"
            filterHydrogens = False
            if buildTypeRef in ["oe-iso-smiles", "oe-smiles", "cactvs-smiles", "cactvs-iso-smiles", "acdlabs-smiles", "inchi"]:
                filterHydrogens = True
            #
            for ccId, ccObj in ccMolD.items():
                # ----
                tId = oemf.setChemCompDef(ccObj)
                self.assertEqual(tId, ccId)
                ok = oemf.build(molBuildType=buildTypeRef, limitPerceptions=limitPerceptions, normalize=False)
                if not ok:
                    logger.info("Build using %r failed for %s", buildTypeRef, ccId)
                    continue
                #
                doTautomers = False
                if doTautomers:
                    tautomerMolL = oemf.getTautomerList()
                    logger.info("%s number reasonable tautomers %d", ccId, len(tautomerMolL))
                #
                refFD = macmp.getChemCompFeatures(ccObj, descriptorProgram="OPENEYE", filterHydrogens=filterHydrogens)
                tstFD = oemf.getOeMoleculeFeatures(filterHydrogens=filterHydrogens)
                # logger.info("tstFD %r", tstFD)
                ok, retCmp = macmp.compare(refFD, tstFD, tstInfo="Openeye ISO SMILES")
                if not ok:
                    logger.info("Comparison failed build type %r and %r", buildTypeRef, ccId)
                    logger.debug(
                        "diff -> atomatic atoms %r stereo atoms %r bond types %r aromatic bonds %r",
                        retCmp.difAromaticAtoms,
                        retCmp.difStereoAtoms,
                        retCmp.difTypeBonds,
                        retCmp.difAromaticBonds,
                    )
                    failL.append(ccId)
                #
            logger.info("Failures (%d) %r: ", len(failL), failL)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testRoundTrip(self):
        """Round trip smiles comparisons -
        """
        try:
            ccMolD = self.__getChemCompDefs()
            # useCache = True
            # quietFlag = False
            # molBuildTypeL = ["model-xyz", "ideal-xyz", None]
            # molBuildTypeL = [None]
            buildTypeRef = "oe-iso-smiles"
            oemf1 = OeMoleculeFactory()
            oemf2 = OeMoleculeFactory()
            #
            for ccId, ccObj in ccMolD.items():
                # ----
                ccIt = iter(PdbxChemCompIt(ccObj))
                cc = next(ccIt)
                formula = cc.getFormulaWithCharge()
                # ccId = cc.getId()
                ccName = cc.getName()
                ifCharge = cc.getFormalChargeAsInt()
                isAmbiguous = cc.getAmbiguousFlag() in ["Y", "y"]
                isCurrent = cc.getReleaseStatus() in ["REL"]
                logger.debug("%s name %r formula %r charge %d", ccId, ccName, formula, ifCharge)
                # ----
                ccId = oemf1.setChemCompDef(ccObj)
                ok = oemf1.build(molBuildType=buildTypeRef, limitPerceptions=True)
                if not ok:
                    logger.info("Build using %r failed for %s (ambiguous flag %r current %r)", buildTypeRef, ccId, isAmbiguous, isCurrent)
                #
                isDiff = False
                #
                if isDiff:
                    genIsoSmi = oemf1.getCanSMILES()
                    oemf2 = OeMoleculeFactory()
                    oemf2.setDescriptor(genIsoSmi, "oe-iso-smiles", ccId)
                    oemf2.build(molBuildType="oe-iso-smiles", limitPerceptions=True)
                    regenIsoSmi = oemf2.getIsoSMILES()
                    if genIsoSmi != regenIsoSmi:
                        logger.info("%s  regenerated ISOSMILES differ \n -- INP: %s\n -- OUT: %s", ccId, genIsoSmi, regenIsoSmi)

                    oed = OeDepictMCSAlignPage()
                    oed.setDisplayOptions(
                        labelAtomName=True,
                        labelAtomCIPStereo=True,
                        labelAtomIndex=False,
                        labelBondIndex=False,
                        labelBondCIPStereo=True,
                        highlightStyleFit="ballAndStickInverse",
                        highLightNotMatchColorRef="pink",
                        bondDisplayWidth=0.5,
                    )
                    oed.setRefMol(oemf1.getGraphMol(), ccId)
                    oed.setFitMol(oemf2.getGraphMol(), ccId)
                    imgPath = os.path.join(self.__workPath, "compare-assigned-" + ccId + "-calc-" + ccId + ".svg")
                    logger.info("Using image path %r", imgPath)
                    aML = oed.alignPair(imagePath=imgPath)
                    if aML:
                        for (rCC, rAt, tCC, tAt) in aML:
                            logger.info("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
                else:
                    logger.info("%s matched all cases", ccId)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testBuilders(self):
        try:
            ccMolD = self.__getChemCompDefs()
            quietFlag = False
            molBuildTypeL = ["model-xyz", "ideal-xyz", None]
            molBuildTypeL = [None]
            for molBuildType in molBuildTypeL:
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                #
                eCount = 0
                for tId, ccObj in ccMolD.items():
                    ccId = oemf.setChemCompDef(ccObj)
                    self.assertEqual(tId, ccId)
                    logger.debug("Building %s using molBuildType %r", ccId, molBuildType)
                    if ccId:
                        ok = oemf.build(molBuildType=molBuildType)
                        logger.debug("Comparing built component %s using molBuildType %r", ccId, molBuildType)
                        # ok = oemf.compare()
                        ok = True
                        if not ok:
                            logger.info("Failing on %s molBuildType %r", ccId, molBuildType)
                            eCount += 1
                        # self.assertTrue(ok)
                    else:
                        logger.error("Cannot process %r", ccObj.getName())
                logger.info("Processed %d components molBuildType %r errors %d", len(ccMolD), molBuildType, eCount)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteCompareDescriptorsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMolecularFactoryTests("testCompareDescriptors"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCompareDescriptorsTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
