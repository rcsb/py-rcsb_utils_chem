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
Utilities to manage OE specific substructure search operations (w/ formula/feature prefiltering)
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
import time
from collections import namedtuple

from openeye import oechem

from rcsb.utils.chem.OeCommonUtils import OeCommonUtils
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
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
        self.__atomexpr, self.__bondexpr = OeCommonUtils.getAtomBondExprOpts(matchOpts)
        #
        _ = kwargs

    def subStructureSearch(self, dataList, procName, optionsD, workingDir):
        """Search index"""
        _ = optionsD
        _ = workingDir
        successList = []
        diagList = []
        sL = []
        # sucessList,resultList,diagList=workerFunc(runList=nextList,procName, optionsD, workingDir)
        try:
            retStatus, successList, sL = self.__subStructureSearch(self.__oeQueryMol, dataList, reverseFlag=False)
            logger.debug("%s status %r found %d search candidates from %d definitions ", procName, retStatus, len(successList), len(dataList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, successList, sL, diagList

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
        sL = []
        retStatus = True
        try:
            ss = oechem.OESubSearch(oeQueryMol, self.__atomexpr, self.__bondexpr)
            if not ss.IsValid():
                retStatus = False
                logger.error("Unable to initialize substructure search!")
                return retStatus, hL, sL
            #
            for idx in idxList:
                mol = oechem.OEGraphMol()
                ccId = self.__oeMolDb.GetTitle(idx)
                if not self.__oeMolDb.GetMolecule(mol, idx):
                    logger.error("Unable to read molecule %r at index %r", ccId, idx)
                    continue
                oechem.OEPrepareSearch(mol, ss)
                if ss.SingleMatch(mol) != reverseFlag:
                    score = float(oeQueryMol.NumAtoms()) / float(mol.NumAtoms())
                    hL.append(idx)
                    sL.append((idx, score))
                    # logger.info("%s queryAtoms %d molAtoms %d score %.4f", ccId, oeQueryMol.NumAtoms(), mol.NumAtoms(), score)
            retStatus = True
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        #
        return retStatus, hL, sL


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

    def prefilterIndex(self, oeQueryMol, idxP, matchOpts="relaxed", skipFeatures=False):
        """Filter the full search index base on minimum chemical formula an feature criteria.

        Args:
            oeQueryMol (object): search target moleculed (OEMol)
            idxP (object): instance ChemCompSearchIndexProvider()
            matchOpts (str, optional): search criteria options. Defaults to "default".
            skipFeatures (bool, optional): skip feature filters. Defaults to False.

        Returns:
            (list): list of chemical component identifiers in the filtered search space
        """
        startTime = time.time()
        oemf = OeMoleculeFactory()
        oemf.setOeMol(oeQueryMol, "queryTarget")
        typeCountD = oemf.getElementCounts(useSymbol=True)
        # ccIdL1 = idxP.filterMinimumMolecularFormula(typeCountD)
        #
        featureCountD = oemf.getFeatureCounts() if not skipFeatures else {}
        # Adjust filter according to search options
        if matchOpts in matchOpts in ["relaxed", "graph-relaxed", "simple", "sub-struct-graph-relaxed"]:
            for ky in ["rings_ar", "at_ar", "at_ch"]:
                featureCountD.pop(ky, None)
        elif matchOpts in ["relaxed-stereo", "graph-relaxed-stereo", "sub-struct-graph-relaxed-stereo", "graph-relaxed-stereo-sdeq", "sub-struct-graph-relaxed-stereo-sdeq"]:
            for ky in ["rings_ar", "at_ar"]:
                featureCountD.pop(ky, None)
        elif matchOpts in ["default", "strict", "graph-strict", "graph-default", "sub-struct-graph-strict"]:
            pass
        ccIdL = idxP.filterMinimumFormulaAndFeatures(typeCountD, featureCountD)
        logger.info("Pre-filtering results for formula+feature %d (%.4f seconds)", len(ccIdL), time.time() - startTime)
        return ccIdL

    def searchSubStructure(self, oeQueryMol, idxList=None, ccIdList=None, reverseFlag=False, matchOpts="graph-relaxed", numProc=1):
        if ccIdList:
            idxList = [self.__oeMolDbTitleD[ccId] for ccId in ccIdList if ccId in self.__oeMolDbTitleD]
        if numProc == 1:
            return self.__searchSubStructure(oeQueryMol, idxList=idxList, reverseFlag=reverseFlag, matchOpts=matchOpts)
        else:
            return self.__searchSubStructureMulti(oeQueryMol, idxList=idxList, matchOpts=matchOpts, numProc=numProc, maxChunkSize=self.__chunkSize)

    def __searchSubStructureMulti(self, oeQueryMol, idxList, matchOpts="graph-relaxed", numProc=2, maxChunkSize=10):
        #
        hL = []
        startTime = time.time()
        try:
            searchType = "exhaustive-substructure"
            if idxList:
                searchType = "prefilterd-substructure"
            idxList = idxList if idxList else list(range(self.__oeMolDb.GetMaxMolIdx()))
            #

            rWorker = OeSubStructSearchWorker(oeQueryMol, self.__oeMolDb, matchOpts=matchOpts)
            mpu = MultiProcUtil(verbose=True)
            optD = {"maxChunkSize": maxChunkSize}
            mpu.setOptions(optD)
            mpu.set(workerObj=rWorker, workerMethod="subStructureSearch")
            _, _, resultList, _ = mpu.runMulti(dataList=idxList, numProc=numProc, numResults=2, chunkSize=maxChunkSize)
            logger.debug("Multi-proc result length %d/%d", len(resultList[0]), len(resultList[1]))
            for idx, score in resultList[1]:
                ccId = self.__oeMolDb.GetTitle(idx)
                hL.append(MatchResults(ccId=ccId, searchType=searchType, matchOpts=matchOpts, fpScore=score))
            retStatus = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            retStatus = False
        logger.info("Substructure search returns %d (%.4f seconds)", len(hL), time.time() - startTime)
        return retStatus, hL

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
        startTime = time.time()
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
                    score = float(oeQueryMol.NumAtoms()) / float(oeQueryMol.NumAtoms())
                    hL.append(MatchResults(ccId=ccId, searchType=searchType, matchOpts=matchOpts, fpScore=score))
            retStatus = True
        except Exception as e:
            retStatus = False
            logger.exception("Failing with %s", str(e))
        #
        logger.info("Substructure search returns %d (%.4f seconds)", len(hL), time.time() - startTime)
        return retStatus, hL
