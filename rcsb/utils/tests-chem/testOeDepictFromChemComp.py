##
#
# File:    OeDepictFromChemTests.py
# Author:  jdw
# Date:    5-May-2013
# Version: 0.001
#
# Updates:
#  4-May-2014 jdw add example for depiction from SMILES input
#  6-Jun-2016 jdw general cleanup
##
"""
A collection of tests for the OEDepict and related classes.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import logging
import os
import os.path
import time
import unittest

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeDepict import OeDepict
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class OeDepictFromChemCompTests(unittest.TestCase):
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
        #

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

    def testDepictByBuildType(self):
        """Compare depictions constructed molecules with various builds from chemical defintions -
        """
        try:
            ccIdList = self.__ccIdList
            ccMolD = self.__getChemCompDefs()
            #
            limitPerceptions = True
            molBuildTypeL = ["model-xyz", "ideal-xyz", "connection-table", "oe-iso-smiles"]
            #
            startTime = time.time()
            oefm = OeMoleculeFactory()
            for molBuildType in molBuildTypeL:
                for ccId in ccIdList:
                    ccObj = ccMolD[ccId]
                    # ----
                    tId = oefm.setChemCompDef(ccObj)
                    self.assertEqual(tId, ccId)
                    ok = oefm.build(molBuildType=molBuildType, limitPerceptions=limitPerceptions)
                    if not ok:
                        logger.info("Build using %r failed for %s", molBuildType, ccId)
                        continue
                    #
                    oeMol = oefm.getGraphMol()
                    pS = "-limited" if limitPerceptions else ""
                    imagePath = os.path.join(self.__workPath, ccId + "-%s%s.svg" % (molBuildType, pS))
                    oed = OeDepict()
                    title = ""
                    oed.setMolTitleList([(ccId, oeMol, title)])
                    oed.setDisplayOptions(labelAtomName=False, labelAtomCIPStereo=True, labelAtomIndex=False, labelBondIndex=False, cellBorders=False, bondDisplayWidth=0.5)
                    oed.setGridOptions(rows=1, cols=1)
                    oed.prepare()
                    oed.write(imagePath)
            logger.info("Completed depictions on %d molecules (%.4f seconds)", len(ccIdList) * len(molBuildTypeL), time.time() - startTime)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteDepictFromChemComp():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(OeDepictFromChemCompTests("testDepictByBuildType"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteDepictFromChemComp()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
