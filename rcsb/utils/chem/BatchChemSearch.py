##
# File:    BatchChemSearchWrapperTests.py
# Author:  jdw
# Date:    1-Mar-2021
# Version: 0.001
#
# Updates:
#
##
"""
Wrapper for batch chemical search operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import time

from collections import namedtuple

from rcsb.utils.chem.ChemCompSearchWrapper import ChemCompSearchWrapper
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

BatchResults = namedtuple("BatchResults", "queryId query queryType matchOpts ccId fpScore", defaults=(None,) * 6)


class BatchChemSearchWorker(object):
    def __init__(self, ccsW, **kwargs):
        self.__ccsW = ccsW
        _ = kwargs

    def searchDescriptorList(self, dataList, procName, optionsD, workingDir):
        """Multiprocessing method template to execute batch descriptor query"""
        _ = optionsD
        _ = workingDir
        queryType = optionsD.get("queryType", "SMILES")
        matchOpts = optionsD.get("matchOpts", "graph-relaxed")
        successList = []
        failList = []
        retList = []
        diagList = []
        #
        try:
            retList, failList = self.__queryDescriptorList(procName, dataList, queryType, matchOpts=matchOpts)
            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures", procName, len(failList))

            logger.debug("%s built %d search candidates from %d definitions with failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList

    def __queryDescriptorList(self, procName, queryPairList, queryType, matchOpts="graph-strict"):
        """Internal method to perform description search for the input descriptor/identifier list"""
        rL = []
        fL = []
        queryId = None
        try:
            for queryPair in queryPairList:
                queryId, descr = queryPair
                #
                retStatus, mL, _ = self.__ccsW.searchByDescriptor(descr, queryType, matchOpts=matchOpts)
                #
                if retStatus == 0 and len(mL) > 0:
                    rL.append((queryId, descr, queryType, mL))
                else:
                    fL.append(queryPair)
            #
        except Exception as e:
            logger.exception("%s %r failing with %s", procName, queryId, str(e))
        #
        return rL, fL


class BatchChemSearch(object):
    """Wrapper for batch chemical search operations."""

    def __init__(self, **kwargs):
        """Wrapper class for batch chemical search/depiction operations.

        Path and prefix data for wrapper class may be set as keyword arguments
        as environmental variables.

        Args:
            ccUrlTarget (str, optional): path to concatenated chemical component definition file. Defaults to public data file.
            birdUrlTarget (str, optional): path to the concatenated BIRD definition file.  Defaults to public data file.
            cachePath (str): path to top-level cache directory used to store search index file dependencies
                             (default environment variable CHEM_SEARCH_CACHE_PATH or ".")
            numProc (int): multi-process cores to reserve. Default to 6.
            chunkSize (int): multi-process batch size.  Defaults to 50.
        """
        self.__startTime = time.time()
        #
        self.__useCache = kwargs.get("useCache", True)
        self.__numProc = kwargs.get("numProc", 6)
        self.__chunkSize = kwargs.get("chunkSize", 50)
        #
        self.__ccUrlTarget = kwargs.get("ccUrlTarget", None)
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", None)
        self.__ccFileNamePrefix = kwargs.get("ccFileNamePrefix", os.environ.get("CHEM_SEARCH_CC_PREFIX", "cc-full"))
        #
        self.__cachePath = kwargs.get("cachePath", os.environ.get("CHEM_SEARCH_CACHE_PATH", "."))
        # ---
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        # ---
        self.__ccsw = self.__reload()
        #

    def __reload(self):
        ccsw = ChemCompSearchWrapper(cachePath=self.__cachePath)
        ok1 = ccsw.setConfig(
            self.__ccUrlTarget, self.__birdUrlTarget, ccFileNamePrefix=self.__ccFileNamePrefix, useCache=self.__useCache, numProc=self.__numProc, maxChunkSize=self.__chunkSize
        )
        ok2 = ccsw.readConfig()
        ok3 = ccsw.reloadSearchDatabase()
        return ccsw if ok1 and ok2 and ok3 else None

    def testCache(self):
        return self.__ccsw is not None

    def fetchDescriptorList(self, filePath, fmt="tdd", swap=False):
        """Fetch list of descriptor and identifier code pairs.

        Args:
            filePath (str): input descriptor/identifier file path
            fmt (str, optional):  tab delimited 'tdd' or comma-separated values 'csv'. Defaults to 'tdd'.
            swap (bool, optional): swap order of input fields. Defaults to False.

        Returns:
            (list): [description]
        """
        rL = []
        pairL = self.__mU.doImport(filePath, fmt=fmt, rowFormat="list")
        if swap:
            for pair in pairL:
                if len(pair) != 2:
                    continue
                rL.append((pair[1], pair[0]))
        return rL

    def splitSmiles(self, pairList):
        rL = []
        for pair in pairList:
            queryId, descr = pair
            if "." in descr:
                ff = descr.split(".")
                descr = max(ff, key=len)
                rL.append((queryId, descr))
            else:
                rL.append(pair)
        return rL

    def storeMatchList(self, filePath, matchList):
        mL = []
        # namedtuple("BatchResults", "queryId query queryType matchOpts ccId fpScore", defaults=(None,) * 6)
        for match in matchList:
            mL.append({"queryId": match.queryId, "query": match.query, "queryType": match.queryType, "matchOpts": match.matchOpts, "ccId": match.ccId, "fpScore": match.fpScore})
        return self.__mU.doExport(filePath, mL, fmt="json")

    def fetchMatchList(self, filePath):
        mDL = self.__mU.doImport(filePath, fmt="json")
        mL = []
        for mD in mDL:
            mL.append(BatchResults(queryId=mD["queryId"], query=mD["query"], queryType=mD["queryType"], matchOpts=mD["matchOpts"], ccId=mD["ccId"], fpScore=mD["fpScore"]))
        return mL

    def doQuery(self, queryPairList, queryType, matchOpts="graph-strict"):
        """Batch descriptor query.

        Args:
            queryPairList (list): [description]
            queryType ([type]): [description]
            matchOpts (str, optional): [description]. Defaults to "graph-strict".

        Returns:
            (list): [BatchResults(), BatchResults(), ...]
        """

        rL = []
        try:
            if queryType.upper() in ["SMILES", "INCHI"]:
                rL = self.__descriptorQueryMulti(queryPairList, queryType, matchOpts, numProc=self.__numProc, chunkSize=self.__chunkSize)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return rL

    def __descriptorQueryMulti(self, queryPairList, queryType, matchOpts, numProc=2, chunkSize=50):
        """Internal method to invoke descriptor query in multiprocess mode

        Args:
            queryPairList (list): [(identifer, descriptor), ... ]
            queryType (str): SMILES|InChI
            matchOpts (str): match criteria (e.g., graph-strict)
            numProc (int, optional): number of multiprocess cores. Defaults to 2.
            chunkSize (int, optional): multiprocess batch size. Defaults to 50.

        Returns:
            (list): [BatchResults(), ...]
        """
        logger.info("Input %r query length %d using %s numProc %d", queryType, len(queryPairList), matchOpts, numProc)
        rWorker = BatchChemSearchWorker(self.__ccsw)
        mpu = MultiProcUtil(verbose=True)
        optD = {"matchOpts": matchOpts, "queryType": queryType}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="searchDescriptorList")
        ok, failList, resultList, _ = mpu.runMulti(dataList=queryPairList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.debug("Search completed with failures (%d)", len(failList))
        logger.info("Multi-proc status %r failures %r result length %r", ok, len(failList), len(resultList[0]))
        #
        rL = []
        for tup in resultList[0]:
            for mr in tup[3]:
                rL.append(BatchResults(queryId=tup[0], query=tup[1], queryType=tup[2], matchOpts=matchOpts, ccId=mr.ccId, fpScore=mr.fpScore))
        #
        logger.info("Multi-proc status %r failures %r result length %r", ok, len(failList), len(rL))
        #
        return rL
