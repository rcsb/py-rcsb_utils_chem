##
# File:    OeMoleculeProviderTests.py
# Author:  J. Westbrook
# Date:    1-Oct-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities to read and process the dictionary of PDB chemical component definitions.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest
from collections import defaultdict

from openeye import oechem
from rcsb.utils.chem import __version__
from rcsb.utils.chem.OeMoleculeProvider import OeMoleculeProvider
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeMoleculeProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildMoleculeCacheFiles(self):
        oemp = OeMoleculeProvider(dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType="model", useCache=False)
        ok = oemp.testCache()
        self.assertTrue(ok)

    def testCompareMoleculeCacheFiles(self):
        for coordType in ["model", "ideal", None]:
            logger.info(">>>> Rebuild cache using coordType %r", coordType)
            self.__testReproduceDescriptors(coordType)

    def __testReproduceDescriptors(self, coordType, useCache=False):
        oemp = OeMoleculeProvider(dirPath=os.path.join(self.__cachePath, "chem_comp"), coordType=coordType, useCache=useCache)
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
            # ----------
            if canSmiles != ccdD["OE_SMILES"]:
                logger.debug("%s CAN SMILES differ \nccd: %r  \nOE:  %r", ccId, ccdD["OE_SMILES"], canSmiles)
                countD["can_smiles_diff"] += 1

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


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeMoleculeProviderTests("testBuildMoleculeCacheFiles"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
