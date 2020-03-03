##
# File:    OeDepictCompareTests.py
# Author:  jdw
# Date:    28-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
A collection of tests to compare assigned and computed features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import os
import time
import unittest
from collections import defaultdict

from openeye import oechem
from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider
from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeDepict import OeDepict
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeDepictCompareTests(unittest.TestCase):
    def setUp(self):
        #
        self.__startTime = time.time()
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-all.cif")
        self.__molLimit = 50
        self.__ccIdList = ["002", "004", "PRD_000921"]
        self.__myKwargs = {
            "cachePath": self.__cachePath,
            "useCache": True,
            # "ccFileNamePrefix": "cc-filtered",
            "ccFileNamePrefix": "cc-full",
            "molLimit": None,
        }
        #

    def tearDown(self):
        pass

    def __getChemCompDefs(self):
        ccMolD = {}
        ccIdxD = {}
        try:
            ccmP = ChemCompMoleculeProvider(**self.__myKwargs)
            ok = ccmP.testCache()
            ccMolD = ccmP.getMolD()
            ccmP = ChemCompIndexProvider(**self.__myKwargs)
            ccIdxD = ccmP.getIndex()
            ok = ccmP.testCache(minCount=500)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ccMolD, ccIdxD

    def testCompareByBuildType(self):
        """Compare depictions constructed molecules with various builds from chemical defintions -
        all build types 8769 (all)
        connect - smiles 6743
        model vs iso smiles 5937
        ideal va iso smiles  7047
        """
        ccResultD = {}
        genResultD = {}
        smilesByBuildTypeD = {}
        try:
            ccMolD, ccIdxD = self.__getChemCompDefs()
            #
            limitPerceptions = True
            # molBuildTypeL = ["model-xyz", "ideal-xyz", "connection-table", "oe-iso-smiles"]
            molBuildTypeL = ["ideal-xyz", "oe-iso-smiles"]
            #
            startTime = time.time()
            oefm = OeMoleculeFactory()
            oefm.setQuiet()
            for molBuildType in molBuildTypeL:
                for ccId, idxD in ccIdxD.items():
                    ccObj = ccMolD[ccId]
                    # ----
                    ccIsoSmiles = idxD["oe-iso-smiles"]
                    ccSmiles = idxD["oe-smiles"]
                    # ----
                    tId = oefm.setChemCompDef(ccObj)
                    if not tId:
                        logger.info("Skipping bad component %r", ccId)
                        continue
                    self.assertEqual(tId, ccId)
                    ok = oefm.build(molBuildType=molBuildType, limitPerceptions=limitPerceptions)
                    if not ok:
                        logger.info("Build using %r failed for %s", molBuildType, ccId)
                        continue
                    # ------
                    oeMol = oefm.getGraphMol()
                    oeIsoSmiles = oefm.getIsoSMILES()
                    oeSmiles = oefm.getCanSMILES()
                    ccEq = oeIsoSmiles == ccIsoSmiles and oeSmiles == ccSmiles
                    #
                    oefmR = OeMoleculeFactory()
                    oefmR.setQuiet()
                    ccIdGen = ccId + "_gen"
                    oefmR.setDescriptor(oeIsoSmiles, "oe-iso-smiles", ccIdGen)
                    ok = oefmR.build(molBuildType="oe-iso-smiles", limitPerceptions=limitPerceptions)
                    if not ok:
                        logger.info("Build using %r failed for %s", molBuildType, ccIdGen)
                        continue
                    # ------
                    #
                    # oeMolGen = oefmR.getGraphMol()
                    oeIsoSmilesGen = oefmR.getIsoSMILES()
                    oeSmilesGen = oefmR.getCanSMILES()
                    genEq = oeIsoSmiles == oeIsoSmilesGen and oeSmiles == oeSmilesGen
                    smilesByBuildTypeD.setdefault(ccId, {}).setdefault(molBuildType, []).append(oeIsoSmilesGen)
                    #
                    logger.debug("%s buildType %s ccEq %r genEq %r", ccId, molBuildType, ccEq, genEq)
                    if not ccEq:
                        ccResultD.setdefault(molBuildType, []).append(ccId)
                    if not genEq:
                        genResultD.setdefault(molBuildType, []).append(ccId)

                    if False:
                        pS = "-limited" if limitPerceptions else ""
                        imagePath = os.path.join(self.__workPath, ccId + "-%s%s.svg" % (molBuildType, pS))
                        oed = OeDepict()
                        title = ""
                        oed.setMolTitleList([(ccId, oeMol, title)])
                        oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, cellBorders=False, bondDisplayWidth=0.5)
                        oed.setGridOptions(rows=1, cols=1)
                        oed.prepare()
                        oed.write(imagePath)
            logger.info("Completed comparing %d molecules in %d builds (%.4f seconds)", len(ccIdxD), len(molBuildTypeL), time.time() - startTime)
            #
            #
            for molBuildType in molBuildTypeL:
                if molBuildType in genResultD:
                    logger.info("GEN %s (%d) %r", molBuildType, len(genResultD[molBuildType]), genResultD[molBuildType])

            numDiff = 0
            for ccId, btD in smilesByBuildTypeD.items():
                tS = set()
                for molBuildType, sL in btD.items():
                    tS.add(sL[0])
                if len(tS) > 1:
                    numDiff += 1
                    logger.debug("%s diff smiles (%d) %r", ccId, len(tS), tS)
            logger.info("Components with inconsistent SMILES %d", numDiff)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCompareDescriptors(self):
        molLimit = self.__molLimit
        # for molBuildType in ["model-xyz", "ideal-xyz", None]:
        for molBuildType in ["model-xyz"]:
            logger.info("Rebuild cache using molBuildType %r (molLimit=%r)", molBuildType, molLimit)
            self.__testReproduceDescriptors(molBuildType, limitPerceptions=True)

    def __displayAlignedPair(self, ccIdRef, oeMolRef, ccIdFit, oeMolFit, title=None):
        oed = OeDepictMCSAlignPage()
        oed.setDisplayOptions(
            labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
        )
        oed.setRefMol(oeMolRef, ccIdRef)
        oed.setFitMol(oeMolFit, ccIdFit)
        myTitle = title if title else "compare"
        imgPath = os.path.join(self.__workPath, myTitle + "-" + ccIdRef + "-" + ccIdFit + ".svg")
        # logger.info("Using image path %r", imgPath)
        aML = oed.alignPair(imagePath=imgPath)
        if aML:
            for (rCC, rAt, tCC, tAt) in aML:
                logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        oed = OeDepictMCSAlignPage()
        oed.setDisplayOptions(
            labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
        )

    def __testReproduceDescriptors(self, molBuildType, limitPerceptions=True):
        #
        ccMolD, ccIdxD = self.__getChemCompDefs()
        oemf = OeMoleculeFactory()
        countD = defaultdict(int)
        for ccId, ccDef in ccMolD.items():
            tId = oemf.setChemCompDef(ccDef)
            if ccId != tId:
                continue
            oemf.build(molBuildType=molBuildType, limitPerceptions=limitPerceptions)
            oeMol = oemf.getMol()
            #
            countD["total components"] += 1
            if ccId not in ccIdxD:
                logger.info("Missing ccIndex entry for %s", ccId)
                continue
            ccdD = ccIdxD[ccId]
            if ccdD["ambiguous"]:
                countD["ambiguous component"] += 1
                continue
            #
            countD["total molecules"] += 1

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
            if isoSmiles != ccdD["oe-iso-smiles"]:
                logger.info("%s ISO SMILES differ \nccd: %r  \nOE:  %r", ccId, ccdD["oe-iso-smiles"], isoSmiles)
                countD["iso_smiles_diff"] += 1
            # ----------
            if canSmiles != ccdD["oe-smiles"]:
                logger.info("%s CAN SMILES differ \nccd: %r  \nOE:  %r", ccId, ccdD["oe-smiles"], canSmiles)
                countD["smiles_diff"] += 1

            formula = oemf.getFormula()
            if formula.upper() != ccdD["formula"].upper():
                logger.debug("%s formulas differ \nccd: %r  \nOE:  %r", ccId, ccdD["formula"], formula)
                countD["formula_diff"] += 1
            # ---------
            inchiKey = oemf.getInChIKey()
            if inchiKey != ccdD["inchikey"]:
                logger.debug("%s InChI keys differ \nccd: %r  \nOE:  %r", ccId, ccdD["inchikey"], inchiKey)
                countD["inchikey_diff"] += 1
            #
            inchi = oemf.getInChI()
            if inchi != ccdD["inchi"]:
                logger.debug("%s InChIs differ \nccd: %r  \nOE:  %r", ccId, ccdD["inchi"], inchi)
                countD["inchi_diff"] += 1
        #
        #
        for ky, vl in countD.items():
            logger.info("%-12s %6d", ky, vl)


def suiteCompareDescriptorsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictCompareTests("testCompareDescriptors"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCompareDescriptorsTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
