##
# File:    OeSubStructSearchUtils.py
# Author:  jdw
# Date:    2-Oct-2020
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to manage OE specific substructure search operations (exhaustive and formula/feature prefiltered)
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import time
from collections import namedtuple

from openeye import oechem

from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore oeIdx formula", defaults=(None,) * 9)


class OeSubStructSearchWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for substructure search --
    """

    def __init__(self, oeQueryMol, oeMolDb, matchOpts="graph-relaxed", **kwargs):
        self.__oeMolDb = oeMolDb
        self.__oeQueryMol = oeQueryMol
        self.__matchOpts = matchOpts
        if matchOpts in ["default", "strict", "graph-strict", "graph-default"]:
            # atomexpr = oechem.OEExprOpts_DefaultAtoms
            # bondexpr = oechem.OEExprOpts_DefaultBonds
            self.__atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_Aromaticity
            self.__bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Aromaticity | oechem.OEExprOpts_Chiral
        elif matchOpts in ["relaxed-stereo", "graph-relaxed-stereo"]:
            self.__atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_FormalCharge
            self.__bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Chiral
        elif matchOpts in ["relaxed", "graph-relaxed", "simple"]:
            self.__atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge
            self.__bondexpr = oechem.OEExprOpts_BondOrder
        else:
            logger.error("Unanticipated match options %r", matchOpts)
        #
        #
        _ = kwargs

    def subStructureSearch(self, dataList, procName, optionsD, workingDir):
        """Search index"""
        _ = optionsD
        _ = workingDir
        successList = []
        diagList = []
        #
        try:
            retStatus, successList = self.__subStructureSearch(self.__oeQueryMol, dataList, reverseFlag=False)
            logger.debug("%s status %r found %d search candidates from %d definitions ", procName, retStatus, len(successList), len(dataList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, successList, diagList

    def __subStructureSearch(self, oeQueryMol, idxList, reverseFlag=False):
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
            ss = oechem.OESubSearch(oeQueryMol, self.__atomexpr, self.__bondexpr)
            if not ss.IsValid():
                retStatus = False
                logger.error("Unable to initialize substructure search!")
                return retStatus, hL
            #
            for idx in idxList:
                mol = oechem.OEGraphMol()
                if not self.__oeMolDb.GetMolecule(mol, idx):
                    ccId = self.__oeMolDb.GetTitle(idx)
                    logger.error("Unable to read molecule %r at index %r", ccId, idx)
                    continue
                oechem.OEPrepareSearch(mol, ss)
                if ss.SingleMatch(mol) != reverseFlag:
                    hL.append(idx)
            retStatus = True
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        #
        return retStatus, hL


class OeSubStructSearchUtils(object):
    """Utilities to manage OE specific substructure search operations (exhaustive and formula/feature prefiltered)"""

    def __init__(self, oemP, screenType=None, numProc=2, chunkSize=10, verbose=False):
        startTime = time.time()
        self.__verbose = verbose
        self.__oeMolDb, self.__oeMolDbTitleD = oemP.getOeMolDatabase()
        self.__idxTitleD = {v: k for k, v in self.__oeMolDbTitleD.items()}
        self.__numProc = numProc
        self.__chunkSize = chunkSize
        #
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
        ok = True
        okdb = self.__oeMolDb and self.__oeMolDb.NumMols() and self.__idxTitleD and len(self.__idxTitleD)
        if not okdb:
            ok = False
        #
        logger.info("Return status %r", ok)
        return ok

    def searchSubStructure(self, oeQueryMol, idxList=None, ccIdList=None, reverseFlag=False, matchOpts="graph-relaxed", numProc=1):
        if ccIdList:
            # logger.info("idxTitleD %r", list(self.__oeMolDbTitleD.keys())[:100])
            idxList = [self.__oeMolDbTitleD[ccId] for ccId in ccIdList if ccId in self.__oeMolDbTitleD]
        if numProc == 1:
            return self.__searchSubStructure(oeQueryMol, idxList=idxList, reverseFlag=reverseFlag, matchOpts=matchOpts)
        else:
            return self.__searchSubStructureMulti(oeQueryMol, idxList=idxList, matchOpts=matchOpts, numProc=numProc, maxChunkSize=self.__chunkSize)

    # ## JDW
    def __searchSubStructureMulti(self, oeQueryMol, idxList, matchOpts="graph-relaxed", numProc=2, maxChunkSize=10):
        #
        try:
            searchType = "exhaustive-substructure"
            if idxList:
                searchType = "prefilterd-substructure"
            idxList = idxList if idxList else list(range(self.__oeMolDb.GetMaxMolIdx()))
            #
            hL = []
            rWorker = OeSubStructSearchWorker(oeQueryMol, self.__oeMolDb, matchOpts=matchOpts)
            mpu = MultiProcUtil(verbose=True)
            optD = {"maxChunkSize": maxChunkSize}
            mpu.setOptions(optD)
            mpu.set(workerObj=rWorker, workerMethod="subStructureSearch")
            _, _, resultList, _ = mpu.runMulti(dataList=idxList, numProc=numProc, numResults=1, chunkSize=maxChunkSize)
            logger.info("Multi-proc result length %r", len(resultList[0]))
            for idx in resultList[0]:
                ccId = self.__oeMolDb.GetTitle(idx)
                hL.append(MatchResults(ccId=ccId, searchType=searchType, matchOpts=matchOpts))
                # hL.append(ccId)
            retStatus = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            retStatus = False
        return retStatus, hL

    # ## JDW
    def __searchSubStructure(self, oeQueryMol, idxList=None, reverseFlag=False, matchOpts="graph-relaxed"):
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
            if matchOpts in ["default", "strict", "graph-strict", "graph-default"]:
                # atomexpr = oechem.OEExprOpts_DefaultAtoms
                # bondexpr = oechem.OEExprOpts_DefaultBonds
                atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_Aromaticity
                bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Aromaticity | oechem.OEExprOpts_Chiral
            elif matchOpts in ["relaxed-stereo", "graph-relaxed-stereo"]:
                atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_FormalCharge
                bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Chiral
            elif matchOpts in ["relaxed", "graph-relaxed", "simple"]:
                atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge
                bondexpr = oechem.OEExprOpts_BondOrder
            else:
                logger.error("Unanticipated match options %r", matchOpts)
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
                    # hL.append(MatchResults(ccId=ccId, oeMol=mol, searchType=searchType, matchOpts=matchOpts))
                    hL.append(MatchResults(ccId=ccId, searchType=searchType, matchOpts=matchOpts))
            retStatus = True
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        #
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
