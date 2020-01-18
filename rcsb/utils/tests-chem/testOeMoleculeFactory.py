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
from collections import defaultdict

from openeye import oechem
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompDescriptorIt, PdbxChemCompIt

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
        #

    def tearDown(self):
        pass

    def testRoundTrip(self):
        """Round trip smiles comparisons -
        """
        try:
            useAbbrev = True
            # useCache = True
            quietFlag = False
            molLimit = 100
            # coordTypeL = ["model", "ideal", None]
            # coordTypeL = [None]
            oeIo = OeIoUtils()
            if useAbbrev:
                rdCcObjL = oeIo.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
            #
            # oemp = OeMoleculeProvider(dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType=coordType, useCache=useCache, molLimit=self.__molLimit)
            # ok = oemp.testCache()
            # self.assertTrue(ok)
            # rdCcObjL = oemp.getComponentDefinitions()
            #
            self.assertGreater(len(rdCcObjL), 4)
            ccObjL = rdCcObjL[:molLimit] if molLimit else rdCcObjL
            for ccObj in ccObjL:
                ccIt = PdbxChemCompIt(ccObj)
                for cc in ccIt:
                    formula = cc.getFormulaWithCharge()
                    ccId = cc.getId()
                    ccName = cc.getName()
                    ifCharge = cc.getFormalChargeAsInt()
                    isAmbiguous = cc.getAmbiguousFlag() in ["Y", "y"]
                    isCurrent = cc.getReleaseStatus() in ["REL"]
                logger.debug("%s name %r formula %rcharge %d", ccId, ccName, formula, ifCharge)
                desIt = PdbxChemCompDescriptorIt(ccObj)
                isoSmiles = smiles = inchi = inchiKey = None
                for des in desIt:
                    desType = des.getType().upper()
                    desProg = des.getProgram().upper()
                    desText = des.getDescriptor().strip()
                    if "OPEN" in desProg and desType == "SMILES_CANONICAL":
                        isoSmiles = desText
                    elif "CACTVS" in desProg and desType == "SMILES_CANONICAL":
                        isoSmilesCactvs = desText
                    elif "OPEN" in desProg and desType == "SMILES":
                        smiles = desText
                    elif desType == "INCHI":
                        inchi = desText
                    elif desType == "INCHIKEY":
                        inchiKey = desText
                    logger.debug("%s type %r prog %r text %r", ccId, desType, desProg, desText)
                #

                oemf1 = OeMoleculeFactory()
                if quietFlag:
                    oemf1.setQuiet()
                #
                if not isoSmiles:
                    logger.info("%s No OE ISOSMILES ambiguous %r current %r", ccId, isAmbiguous, isCurrent)
                    continue
                #
                if not isoSmilesCactvs:
                    logger.info("%s No CACTVS ISOSMILES ambiguous %r current %r", ccId, isAmbiguous, isCurrent)
                    continue
                #
                isDiff = False
                genMol = oeIo.smilesToMol(isoSmiles, limitPerceptions=True)
                oemf1.setOeMol(genMol, ccId)
                #
                genIsoSmi = oemf1.getIsoSMILES()
                if genIsoSmi != isoSmiles:
                    isDiff = True
                    logger.info("%s ISOSMILES differ \n -- INP: %s\n -- OUT: %s", ccId, isoSmiles, genIsoSmi)

                genSmi = oemf1.getCanSMILES()
                if genSmi != smiles:
                    logger.info("%s SMILES differ \n -- INP: %s\n -- OUT: %s", ccId, smiles, genSmi)

                genKey = oemf1.getInChIKey()
                if inchiKey != genKey:
                    logger.info("%s InChIKeys differ \n -- INP: %s\n -- OUT: %s", ccId, inchiKey, genKey)
                    # logger.info("%s InChiKeys differ", ccId)

                genInChI = oemf1.getInChI()
                if inchi != genInChI:
                    logger.info("%s InChI differ \n -- INP: %s\n -- OUT: %s", ccId, inchi, genInChI)
                #
                if isDiff:
                    regenMol = oeIo.smilesToMol(genIsoSmi, limitPerceptions=True)
                    oemf2 = OeMoleculeFactory()
                    oemf2.setOeMol(regenMol, ccId)
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
                    oed.setRefMol(genMol, ccId)
                    oed.setFitMol(regenMol, ccId)
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
            quietFlag = False
            molLimit = 40000
            coordTypeL = ["model", "ideal", None]
            coordTypeL = [None]
            oeIo = OeIoUtils()
            rdCcObjL = oeIo.getComponentDefinitions(os.path.join(self.__dataPath, "components-abbrev.cif"))
            #
            self.assertGreater(len(rdCcObjL), 4)
            ccObjL = rdCcObjL[:molLimit] if molLimit else rdCcObjL
            logger.info("Processing %d components", len(ccObjL))
            for coordType in coordTypeL:
                logger.debug("Processing %d of %d components with coordType %r", len(ccObjL), len(rdCcObjL), coordType)
                #
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                #
                eCount = 0
                for ii, ccObj in enumerate(ccObjL, 1):
                    ccId = oemf.set(ccObj)
                    logger.debug("Building %s using coordType %r", ccId, coordType)
                    if ccId:
                        if coordType:
                            ok = oemf.build3D(coordType=coordType)
                        else:
                            ok = oemf.build2D()
                        logger.debug("Comparing built component %s using coordType %r", ccId, coordType)
                        ok = oemf.compare()
                        if not ok:
                            logger.info("Failing on %s coordType %r component number %d", ccId, coordType, ii)
                            eCount += 1
                        # self.assertTrue(ok)
                    else:
                        logger.error("Cannot process %r", ccObj.getName())
                logger.info("Processing %d components coordType %r errors %d", len(ccObjL), coordType, eCount)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCompareDescriptors(self):

        molLimit = self.__molLimit
        # for coordType in ["model", "ideal", None]:
        for coordType in [None]:
            logger.info("Rebuilding cache using coordType %r (molLimit=%r)", coordType, molLimit)
            self.__testReproduceDescriptors(coordType, useCache=False, molLimit=molLimit)

    def __testReproduceDescriptors(self, coordType, useCache=True, molLimit=None):
        try:
            oed = OeDepictMCSAlignPage()
            oed.setDisplayOptions(
                labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
            )
            #
            oeIo = OeIoUtils()
            oemp = OeMoleculeProvider(
                ccUrlTarget=self.__ccUrlTarget,
                birdUrlTarget=self.__birdUrlTarget,
                dirPath=os.path.join(self.__cachePath, "chem_comp"),
                coordType=coordType,
                useCache=useCache,
                molLimit=molLimit,
            )
            ok = oemp.testCache()
            self.assertTrue(ok)
            oeMolD = oemp.getOeMolD()
            ccIdxD = oemp.getChemCompIdx()
            oemf = OeMoleculeFactory()
            countD = defaultdict(int)
            for ccId, oeMol in oeMolD.items():
                countD["total components"] += 1
                if ccId not in ccIdxD:
                    logger.info("Missing ccIndex entry for %s", ccId)
                    continue
                ccdD = ccIdxD[ccId]
                if ccdD["AMBIGUOUS"]:
                    countD["ambiguous component"] += 1
                    continue
                #
                countD["total molecules"] += 1
                oemf.setOeMol(oeMol, ccId)

                nativeCanIsoSmiles = oechem.OECreateIsoSmiString(oeMol)
                canIsoSmiles = oechem.OEMolToSmiles(oeMol)
                isoSmiles = oemf.getIsoSMILES()
                canSmiles = oemf.getCanSMILES()
                # check interal consistency
                if nativeCanIsoSmiles != isoSmiles:
                    logger.error("%s stored and calculated OE smiles differ %s %s", ccId, nativeCanIsoSmiles, isoSmiles)
                if canIsoSmiles != isoSmiles:
                    logger.error("%s calculated OE ISO and canonical smiles differ %s %s", ccId, isoSmiles, canIsoSmiles)

                # compare with archived values
                if isoSmiles != ccdD["OE_ISO_SMILES"]:
                    logger.debug("%s ISO SMILES differ \nccd: %r  \nOE:  %r", ccId, ccdD["OE_ISO_SMILES"], isoSmiles)
                    countD["iso_smiles_diff"] += 1
                    tMol = oeIo.smilesToMol(ccdD["OE_ISO_SMILES"])
                    #
                    oed.setRefMol(tMol, ccId)
                    oed.setFitMol(oeMol, ccId)
                    imgPath = os.path.join(self.__workPath, "compare-assigned-" + ccId + "-calc-" + ccId + ".svg")
                    logger.info("Using image path %r", imgPath)
                    aML = oed.alignPair(imagePath=imgPath)
                    if aML:
                        for (rCC, rAt, tCC, tAt) in aML:
                            logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)

                # ----------
                if canSmiles != ccdD["OE_SMILES"]:
                    logger.debug("%s CAN SMILES differ \nccd: %r  \nOE:  %r", ccId, ccdD["OE_SMILES"], canSmiles)
                    countD["can_smiles_diff"] += 1
                    tMol = oeIo.smilesToMol(ccdD["OE_SMILES"])

                formula = oemf.getFormula()
                if formula.upper() != ccdD["FORMULA"].upper():
                    logger.debug("%s formulas differ \nccd: %r  \nOE:  %r", ccId, ccdD["FORMULA"], formula)
                    countD["formula_diff"] += 1
                # ---------
                inchiKey = oemf.getInChIKey()
                if inchiKey != ccdD["INCHI_KEY"]:
                    logger.debug("%s InChI keys differ \nccd: %r  \nOE:  %r", ccId, ccdD["INCHI_KEY"], inchiKey)
                    countD["inchi_key_diff"] += 1
                #
                inchi = oemf.getInChI()
                if inchi != ccdD["INCHI"]:
                    logger.debug("%s InChIs differ \nccd: %r  \nOE:  %r", ccId, ccdD["INCHI"], inchi)
                    countD["inchi_diff"] += 1
            #
            #
            for ky, vl in countD.items():
                logger.info("%-12s %6d", ky, vl)
        except Exception as e:
            logger.exception("Failing with %s", str(e))


def suiteCompareDescriptorsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMolecularFactoryTests("testCompareDescriptors"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCompareDescriptorsTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
