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
from collections import namedtuple

from openeye import oechem
from openeye import oegraphsim


logger = logging.getLogger(__name__)

MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore", defaults=(None,) * 7)


class OeSearchUtils(object):
    """ Utility methods to manage OE specific search operations.
    """

    def __init__(self, oemP, fpTypeList=None, screenType=None, numProc=2, verbose=False):
        startTime = time.time()
        self.__verbose = verbose
        self.__fpDbD = {fpType: oemP.getFingerPrintDb(fpType) for fpType in fpTypeList} if fpTypeList else {}
        self.__oeMolDb, self.__oeMolDbTitleD = oemP.getOeMolDatabase()
        self.__idxTitleD = {v: k for k, v in self.__oeMolDbTitleD.items()}
        if screenType:
            self.__ssDb = oemP.getSubSearchDb(screenType=screenType, numProc=numProc, forceRefresh=True)
        endTime = time.time()
        logger.info("Loaded %d definitions (%.4f seconds)", self.__oeMolDb.NumMols(), endTime - startTime)
        logger.debug("self.__oeMolDbTitleD %s", list(self.__oeMolDbTitleD.items())[:5])

    def searchSubStructure(self, oeQueryMol, idxList=None, reverseFlag=False, matchOpts="default"):
        """Perform a substructure search for the input query molecule on the binary
        database of molecules.  The search optionally restricted to the input index
        list.   The sense of the search may be optionally reversed.

        Args:
            oeQueryMol (object): query molecule OeGraphMol or OeQmol
            idxList ([type], optional): [description]. Defaults to None.
            reverseFlag (bool, optional): [description]. Defaults to False.

        Returns:
            [type]: [description]
        """
        hL = []
        retStatus = True
        try:
            # logger.info("Query mol type %r", type(oeQueryMol))
            if matchOpts == "default":
                atomexpr = oechem.OEExprOpts_DefaultAtoms
                bondexpr = oechem.OEExprOpts_DefaultBonds
            elif matchOpts == "simple-stereo":
                atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_FormalCharge
                bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Chiral
            elif matchOpts == "simple":
                atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge
                bondexpr = oechem.OEExprOpts_BondOrder
            #
            ss = oechem.OESubSearch(oeQueryMol, atomexpr, bondexpr)
            if not ss.IsValid():
                retStatus = False
                logger.error("Unable to initialize substructure search!")
                return retStatus, hL
            #
            searchType = "exhaustive-substructure"
            if idxList:
                searchType = "prefilterd-substructure"
            idxIt = idxList if idxList else range(self.__oeMolDb.GetMaxMolIdx())

            for idx in idxIt:
                mol = oechem.OEGraphMol()
                ccId = self.__oeMolDb.GetTitle(idx)
                if not self.__oeMolDb.GetMolecule(mol, idx):
                    logger.error("Unable to read molecule %r at index %r", ccId, idx)
                    continue
                oechem.OEPrepareSearch(mol, ss)
                if ss.SingleMatch(mol) != reverseFlag:
                    hL.append(MatchResults(ccId=ccId, oeMol=mol, searchType=searchType, matchOpts=matchOpts))
            retStatus = True
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        #
        return retStatus, hL

    def searchFingerPrints(self, oeQueryMol, fpType, minFpScore=None, maxFpResults=100, annotateMols=False, verbose=False):
        hL = []
        retStatus = True
        try:
            fpDb = self.__fpDbD[fpType] if fpType in self.__fpDbD else {}
            if not fpDb:
                retStatus = False
                return retStatus, hL
            #
            opts = oegraphsim.OEFPDatabaseOptions(maxFpResults, oegraphsim.OESimMeasure_Tanimoto)
            if minFpScore:
                opts.SetCutoff(minFpScore)
            #
            if verbose:
                logger.info("Using %d fingerprint %s type %s", fpDb.NumFingerPrints(), fpType, fpDb.GetFPTypeBase().GetFPTypeString())
                startTime = time.time()
            #
            scores = fpDb.GetSortedScores(oeQueryMol, opts)
            oeMol = oechem.OEGraphMol()
            for si in scores:
                if self.__oeMolDb.GetMolecule(oeMol, si.GetIdx()):
                    ccId = self.__oeMolDb.GetTitle(si.GetIdx())
                    if annotateMols:
                        tS = "For %s index %r %r similarity score %.4f " % (ccId, si.GetIdx(), self.__idxTitleD[si.GetIdx()], si.GetScore())
                        oechem.OESetSDData(oeMol, fpType, tS)
                    hL.append(MatchResults(ccId=ccId, oeMol=oeMol, searchType="fp", fpType=fpType, fpScore=si.GetScore()))
            if verbose:
                endTime = time.time()
                logger.info("Fingerprint %s returning %d hits (%.4f sec)", fpType, len(hL), endTime - startTime)
        except Exception as e:
            retStatus = False
            logger.exception("Failing fpType %r with %s", fpType, str(e))

        return retStatus, hL

    def getFingerPrintScores(self, oeQueryMol, fpType, minFpScore, maxFpResults):
        hL = []
        retStatus = True
        try:
            fpDb = self.__fpDbD[fpType]
            opts = oegraphsim.OEFPDatabaseOptions(maxFpResults, oegraphsim.OESimMeasure_Tanimoto)
            if minFpScore:
                opts.SetCutoff(minFpScore)
            scores = fpDb.GetSortedScores(oeQueryMol, opts)
            # hL = [(self.__oeMolDb.GetTitle(si.GetIdx()), si.GetScore()) for si in scores]
            hL = [MatchResults(ccId=self.__oeMolDb.GetTitle(si.GetIdx()), searchType="fp", fpType=fpType, fpScore=si.GetScore()) for si in scores]
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        return retStatus, hL

    def searchSubStructureWithFingerPrint(self, oeQueryMol, fpType, minFpScore, maxFpResults, matchOpts="simple"):
        hL = []
        retStatus = True
        try:
            fpDb = self.__fpDbD[fpType]
            opts = oegraphsim.OEFPDatabaseOptions(maxFpResults, oegraphsim.OESimMeasure_Tanimoto)
            if minFpScore:
                opts.SetCutoff(minFpScore)
            scores = fpDb.GetSortedScores(oeQueryMol, opts)
            idxList = [si.GetIdx() for si in scores]
            retStatus, hL = self.searchSubStructure(oeQueryMol, idxList=idxList, reverseFlag=False, matchOpts=matchOpts)
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        return retStatus, hL

    def searchSubStructureScreened(self, qmol, maxMatches=10):
        """Perform screened substructure search on input query molecule (OEQMol) subject to maxMatches.

        Args:
            qmol (OEQmol): OE query molecule
            maxMatches (int, optional): maximum number matches returned . Defaults to 10.

        Returns:
            bool, list: search return status, list of namedtuple MatchResults.
        """
        retStatus = True
        hL = []
        try:
            query = oechem.OESubSearchQuery(qmol, maxMatches)
            result = oechem.OESubSearchResult()
            status = self.__ssDb.Search(result, query)
            statusText = oechem.OESubSearchStatusToName(status)
            if statusText.upper() != "FINISHED":
                retStatus = False
                logger.info("Search failing with status = %r", statusText)
                return retStatus, hL
            #
            if self.__verbose:
                logger.info("Number of targets  = %d", result.NumTargets())
                logger.info("Number of screened = %d", result.NumScreened())
                logger.info("Number of searched = %d", result.NumSearched())
                logger.info("Number of total matches = %d", result.NumTotalMatches())
                logger.info("Number of kept  matches = %d ", result.NumMatches())
            #
            for index in result.GetMatchIndices():
                hL.append(MatchResults(ccId=self.__ssDb.GetTitle(index), searchType="screened-substructure"))
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))

        return retStatus, hL
