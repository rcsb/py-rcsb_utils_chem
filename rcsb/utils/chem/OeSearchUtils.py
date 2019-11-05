##
# File:    OeSearchUtils.py
# Author:  jdw
# Date:    22-Oct-2019
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to manage OE specific search operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import time

from openeye import oechem
from openeye import oegraphsim

# from rcsb.utils.chem.OeMolFactory import OeMolFactory


logger = logging.getLogger(__name__)


class OeSearchUtils(object):
    """ Utility methods to manage OE specific search operations.
    """

    def __init__(self, ccP, verbose=False):
        self.__verbose = verbose
        self.__fpDb = ccP.getFingerPrintDb()
        self.__oeMolDb, self.__oeMolDbTitleD = ccP.getOeMolDatabase()
        self.__ssDb = None

    def searchSubStructure(self, queryMol, idxList=None, reverseFlag=False):
        """Perform a substructure search for the input query molecule on the binary
        database of molecules.  The search optionally restricted to the input index
        list.   The sense of the search may be optionally reversed.

        Args:
            queryMol (object): query molecule OeGraphMol or OeQmol
            idxList ([type], optional): [description]. Defaults to None.
            reverseFlag (bool, optional): [description]. Defaults to False.

        Returns:
            [type]: [description]
        """
        mL = []
        try:
            ss = oechem.OESubSearch(queryMol)
            if not ss.IsValid():
                logger.error("Unable to initialize substructure search!")
                return mL
            #
            idxIt = idxList if idxList else range(self.__oeMolDb.GetMaxMolIdx())

            for idx in idxIt:
                mol = oechem.OEGraphMol()
                # title = self.__oeMolDb.GetTitle(idx)
                if not self.__oeMolDb.GetMolecule(mol, idx):
                    logger.error("Unable to read a molecule from index %r", idx)
                    continue
                oechem.OEPrepareSearch(mol, ss)
                if ss.SingleMatch(mol) != reverseFlag:
                    mL.append(mol)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return mL

    def searchFingerPrints(self, queryMol, numResults=10):
        try:
            opts = oegraphsim.OEFPDatabaseOptions(numResults, oegraphsim.OESimMeasure_Tanimoto)
            if not oegraphsim.OEAreCompatibleDatabases(self.__oeMolDb, self.__fpDb):
                logger.error("Molecule and fingerprint databases are not compatible!")

            fptype = self.__fpDb.GetFPTypeBase()
            lenFp = self.__fpDb.NumFingerPrints()
            memTypeStr = self.__fpDb.GetMemoryTypeString()
            logger.info("Using fingerprint type %s", fptype.GetFPTypeString())

            startTime = time.time()
            scores = self.__fpDb.GetSortedScores(queryMol, opts)
            endTime = time.time()
            logger.info("%5.2f sec to search %d fingerprints %s", endTime - startTime, lenFp, memTypeStr)

            startTime = time.time()
            hL = []
            hit = oechem.OEGraphMol()
            for si in scores:
                if self.__oeMolDb.GetMolecule(hit, si.GetIdx()):
                    oechem.OESetSDData(hit, "For %s index %r similarity score", "%.2f" % self.__oeMolDb.GetTitle(si.GetIdx()), si.GetIdx(), si.GetScore())
                    hL.append(hit)
            #
            endTime = time.time()
            logger.info("%5.2f sec to write %d hits", endTime - startTime, opts.GetLimit())
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return hL

    def searchSubStructureScreened(self, qmol, maxmatches=10):
        hL = []
        op = "search"
        # search database

        if op == "count":
            oechem.OEThrow.Info("Number of hits: %d" % self.__ssDb.NumMatches(qmol))
        else:
            query = oechem.OESubSearchQuery(qmol, maxmatches)
            result = oechem.OESubSearchResult()
            status = self.__ssDb.Search(result, query)

            logger.info("Search status = %r", oechem.OESubSearchStatusToName(status))
            logger.info("Number of targets  = %d", result.NumTargets())
            logger.info("Number of screened = %d", result.NumScreened())
            logger.info("Number of searched = %d", result.NumSearched())
            logger.info("Number of total matches = %d", result.NumTotalMatches())
            logger.info("Number of kept  matches =%d ", result.NumMatches())
            #
            for index in result.GetMatchIndices():
                hL.append(self.__ssDb.GetTitle(index))

        return hL
