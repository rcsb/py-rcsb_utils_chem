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

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeMolecularFactoryTests(unittest.TestCase):
    def setUp(self):
        #
        self.__workPath = os.path.join(HERE, "test-output")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__molLimit = 50
        #

    def tearDown(self):
        pass

    def testCompareDescriptors(self):
        molLimit = self.__molLimit
        # for coordType in ["model", "ideal", None]:
        for coordType in [None]:
            logger.info("Rebuild cache using coordType %r (molLimit=%r)", coordType, molLimit)
            self.__testReproduceDescriptors(coordType, useCache=False, molLimit=molLimit)

    def __testReproduceDescriptors(self, coordType, useCache=True, molLimit=None):
        oed = OeDepictMCSAlignPage()
        oed.setDisplayOptions(
            labelAtomName=True, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, highlightStyleFit="ballAndStickInverse", bondDisplayWidth=0.5
        )
        #
        oeIo = OeIoUtils()
        oemp = OeMoleculeProvider(dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType=coordType, useCache=useCache, molLimit=molLimit)
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
                if len(aML) > 0:
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


def suiteCompareDescriptorsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMolecularFactoryTests("testCompareDescriptors"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite1 = suiteCompareDescriptorsTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite1)
    #
