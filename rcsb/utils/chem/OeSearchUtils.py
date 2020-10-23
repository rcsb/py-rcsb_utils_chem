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
Utilities to manage OE specific similarity search (match) operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import time
from collections import namedtuple
from collections import OrderedDict

from openeye import oechem
from openeye import oegraphsim

from rcsb.utils.chem.OeCommonUtils import OeCommonUtils

logger = logging.getLogger(__name__)

MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore oeIdx formula", defaults=(None,) * 9)


class OeSearchUtils(object):
    """Utilities to manage OE specific similarity search (match) operations."""

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

    def testCache(self):
        """Check for existence of data dependencies and non-zero content counts

        Returns:
            bool: True for success or False otherwise
        """
        okdb = self.__oeMolDb and self.__oeMolDb.NumMols() and self.__idxTitleD and len(self.__idxTitleD)
        if not okdb or len(self.__fpDbD) < 1:
            return False
        #
        ok = True
        for fpType, fpDb in self.__fpDbD.items():
            ok = ok and fpDb.NumFingerPrints() > 1
            logger.debug("fpType %r fp count %d status %r", fpType, fpDb.NumFingerPrints(), ok)
        logger.info("Return status %r", ok)
        return ok

    def searchSubStructure(self, oeQueryMol, idxList=None, reverseFlag=False, matchOpts="graph-relaxed"):
        """Perform a graph match for the input query molecule on the binary
        database of molecules.  The search optionally restricted to the input index
        list.   The sense of the search may be optionally reversed.

        Args:
            oeQueryMol (object): query molecule OeGraphMol or OeQmol
            idxList ([type], optional): [description]. Defaults to None.
            reverseFlag (bool, optional): [description]. Defaults to False.
            matchOpts (str, optional): graph match criteria type (graph-strict|graph-relaxed|graph-relaxed-stereo). Defaults to "graph-relaxed".

        Returns:
            [type]: [description]
        """
        hL = []
        retStatus = True
        try:
            # logger.info("Query mol type %r", type(oeQueryMol))
            atomexpr, bondexpr = OeCommonUtils.getAtomBondExprOpts(matchOpts)
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

    def searchFingerPrints(self, oeQueryMol, fpType, minFpScore=None, maxFpResults=50, annotateMols=False, verbose=False):
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
        """Return finger print search scores for the input OE molecule.

        Args:
            oeQueryMol (OEmol): OE graph molecule
            fpType (str): fingerprint type  [TREE,PATH,MACCS,CIRCULAR,LINGO]
            fpMinScore (float): min fingerprint match score (0.0-1.0)
            maxFpResults (int): maximum number of finger print results returned

        Returns:
            (bool, list): status, finger match lists of type (MatchResults)
        """
        hL = []
        retStatus = True
        try:
            fpDb = self.__fpDbD[fpType]
            opts = oegraphsim.OEFPDatabaseOptions(maxFpResults, oegraphsim.OESimMeasure_Tanimoto)
            if minFpScore:
                opts.SetCutoff(minFpScore)
            scores = fpDb.GetSortedScores(oeQueryMol, opts)
            hL = [MatchResults(ccId=self.__oeMolDb.GetTitle(si.GetIdx()), searchType="fp", fpType=fpType, fpScore=si.GetScore(), oeIdx=si.GetIdx()) for si in scores]
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        return retStatus, hL

    def searchSubStructureAndFingerPrint(self, oeQueryMol, fpTypeCutoffList, maxFpResults, matchOpts="graph-relaxed"):
        """Return graph match and finger print search results for the input OE molecule using finger print pre-filtering.

        Args:
            oeQueryMol (OEmol): OE graph molecule
            fpTypeCutoffList (list): [(finger print type, min score),...]
            maxFpResults (int): maximum number of finger print results returned
            matchOpts (str, optional): graph match criteria type (graph-strict|graph-relaxed|...). Defaults to "graph-relaxed".

        Returns:
            (bool, list, list): status, graph match and finger match lists of type (MatchResults)
        """
        fpL = []
        ssL = []
        retStatus = True
        try:
            for fpType, fpCutoff in fpTypeCutoffList:
                ok, tL = self.getFingerPrintScores(oeQueryMol, fpType, fpCutoff, maxFpResults)
                fpL.extend(tL)
            fpL = list(set(fpL))
            fpL = sorted(fpL, key=lambda nTup: nTup.fpScore, reverse=True)
            idxList = [nTup.oeIdx for nTup in fpL]
            idxList = list(OrderedDict.fromkeys(idxList))
            # -- only continue with a non-empty fingerprint result --
            if matchOpts not in ["fingerprint-similarity"] and idxList:
                # Save the maximum fp score
                fpScoreD = {}
                for fpTup in fpL:
                    fpScoreD[fpTup.ccId] = max(fpScoreD[fpTup.ccId], fpTup.fpScore) if fpTup.ccId in fpScoreD else fpTup.fpScore
                retStatus, rTupL = self.searchSubStructure(oeQueryMol, idxList=idxList, reverseFlag=False, matchOpts=matchOpts)
                for rTup in rTupL:
                    ssL.append(rTup._replace(fpScore=fpScoreD[rTup.ccId]))
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        return ok and retStatus, ssL, fpL

    def searchSubStructureWithFingerPrint(self, oeQueryMol, fpType, minFpScore, maxFpResults, matchOpts="graph-relaxed"):
        """Return graph match search results for the input OE molecule using finger print pre-filtering.

        Args:
            oeQueryMol (OEmol): OE graph molecule
            fpTypeCutoffList (list): [(finger print type, min score),...]
            maxFpResults (int): maximum number of finger print results returned
            matchOpts (str, optional): graph match criteria type (graph-strict|graph-relaxed). Defaults to "graph-relaxed".

        Returns:
            (bool, list, list): status, graph match and finger match lists of type (MatchResults)
        """
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
