##
# File:    ChemCompSearchIndexProvider.py
# Author:  J. Westbrook
# Date:    3-Mar-2020
#
# Updates:
#  13-Mar-2020 jdw Add formula index search method.
##
"""
Utilities to read and process an index of PDB chemical component definitions.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import time
from collections import namedtuple

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.IoUtil import getObjSize
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.io.SingletonClass import SingletonClass
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

# from rcsb.utils.multiproc.MultiProcPoolUtil import MultiProcPoolUtil


logger = logging.getLogger(__name__)


MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore oeIdx formula", defaults=(None,) * 9)


class ChemCompSearchIndexWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for calculating search index candidates --
    """

    def __init__(self, ccObjD, **kwargs):
        self.__ccObjD = ccObjD
        _ = kwargs

    def buildRelatedList(self, dataList, procName, optionsD, workingDir):
        """Build search candidates for the input list of chemical component definitions
        and return index feature data.
        """
        _ = optionsD
        _ = workingDir
        # fetchLimit = optionsD.get("fetchLimit", None)
        limitPerceptions = optionsD.get("limitPerceptions", True)
        quietFlag = optionsD.get("quietFlag", False)
        successList = []
        failList = []
        retList = []
        diagList = []
        #
        try:
            retList, failList = self.__buildChemCompSearchIndex(procName, dataList, limitPerceptions=limitPerceptions, quietFlag=quietFlag)
            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built %d search candidates from %d definitions with failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList

    def __buildChemCompSearchIndex(self, procName, ccIdList, limitPerceptions=False, quietFlag=False):
        """Internal method return a dictionary of extracted chemical component descriptors and formula."""
        rL = []
        fL = []
        try:
            for ccId in ccIdList:
                if ccId not in self.__ccObjD:
                    logger.error("%s missing chemical definition for %s", procName, ccId)
                    fL.append(ccId)
                    continue
                dataContainer = self.__ccObjD[ccId]
                # ----
                oemf = OeMoleculeFactory()
                if quietFlag:
                    oemf.setQuiet()
                tId = oemf.setChemCompDef(dataContainer)
                if tId != ccId:
                    logger.error("%s %s chemical component definition import error", procName, ccId)
                    fL.append(ccId)
                    continue
                relD = oemf.buildRelated(limitPerceptions=limitPerceptions)
                logger.debug("%s %s related molecular forms %d", procName, ccId, len(relD))
                if relD:
                    rL.extend([relD[v] for v in relD])
                else:
                    fL.append(ccId)
        except Exception as e:
            logger.exception("%s failing with %s", procName, str(e))
        return rL, fL


class ChemCompSearchIndexProvider(object):
    """Utilities to read and process the index of chemical component definitions search targets"""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "chem_comp")
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        self.__searchIdx = self.__reload(**kwargs)

    def testCache(self, minCount=None, logSizes=False):
        if logSizes and self.__searchIdx:
            logger.info("searchIdxD (%.2f MB)", getObjSize(self.__searchIdx) / 1000000.0)
        ok = self.__searchIdx and len(self.__searchIdx) >= minCount if minCount else self.__searchIdx is not None
        return ok

    def getIndex(self):
        return self.__searchIdx

    def getIndexEntry(self, searchCcId):
        try:
            return self.__searchIdx[searchCcId]
        except Exception as e:
            logger.debug("Get index entry %r failing with %s", searchCcId, str(e))
        return None

    def getIndexFilePath(self):
        return os.path.join(self.__dirPath, "%s-search-idx-chemical-components.json" % self.__ccFileNamePrefix)

    def __reload(self, **kwargs):
        """Reload or created index of PDB chemical components.

        Args:
            cachePath (str): path to the directory containing cache files
            ccIdxFileName (str): serialized chemical component data index file name

         Returns:
            (list): chemical component data containers
        """
        #
        searchIdxD = {}
        useCache = kwargs.get("useCache", True)
        molLimit = kwargs.get("molLimit", 0)
        numProc = kwargs.get("numProc", 1)
        maxChunkSize = kwargs.get("maxChunkSize", 20)
        limitPerceptions = kwargs.get("limitPerceptions", True)
        quietFlag = kwargs.get("quietFlag", True)
        searchIdxFilePath = self.getIndexFilePath()
        #
        if useCache and self.__mU.exists(searchIdxFilePath):
            _, fExt = os.path.splitext(searchIdxFilePath)
            searchIdxFormat = "json" if fExt == ".json" else "pickle"
            rdCcIdxD = self.__mU.doImport(searchIdxFilePath, fmt=searchIdxFormat)
            searchIdxD = {k: rdCcIdxD[k] for k in sorted(rdCcIdxD.keys())[:molLimit]} if molLimit else rdCcIdxD
        else:
            cmpKwargs = {k: v for k, v in kwargs.items() if k not in ["cachePath", "useCache", "molLimit"]}
            ccmP = ChemCompMoleculeProvider(cachePath=self.__cachePath, useCache=True, molLimit=molLimit, **cmpKwargs)
            ok = ccmP.testCache(minCount=molLimit, logSizes=True)
            if ok:
                searchIdxD = self.__updateChemCompSearchIndex(ccmP.getMolD(), searchIdxFilePath, molLimit, limitPerceptions, numProc, maxChunkSize, quietFlag)
                logger.info("Storing %s with data for %d search candidates (status=%r) ", searchIdxFilePath, len(searchIdxD), ok)
        #
        #
        for idxD in searchIdxD.values():
            idxD["atom-types"] = set(idxD["type-counts"].keys()) if "type-counts" in idxD else set()

        return searchIdxD

    def __updateChemCompSearchIndex(self, ccObjD, filePath, molLimit, limitPerceptions, numProc, maxChunkSize, quietFlag):
        searchIdxD = {}
        try:
            # Serialized index of chemical component search targets
            startTime = time.time()
            _, fExt = os.path.splitext(filePath)
            fileFormat = "json" if fExt == ".json" else "pickle"
            if numProc <= 1:
                searchIdxD = self.__buildChemCompSearchIndex(ccObjD, limitPerceptions=limitPerceptions, molLimit=molLimit)
            else:
                searchIdxD = self.__buildChemCompSearchIndexMulti(
                    ccObjD, limitPerceptions=limitPerceptions, molLimit=molLimit, numProc=numProc, maxChunkSize=maxChunkSize, quietFlag=quietFlag
                )

            ok = self.__mU.doExport(filePath, searchIdxD, fmt=fileFormat)
            endTime = time.time()
            logger.info("Storing %s (%s) with %d search definitions (status=%r) (%.4f seconds)", filePath, fileFormat, len(searchIdxD), ok, endTime - startTime)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return searchIdxD

    def __buildChemCompSearchIndex(self, ccObjD, limitPerceptions=False, molLimit=None):
        """Internal method return a dictionary of extracted chemical component descriptors and formula."""
        rD = {}
        try:
            for ii, ccId in enumerate(ccObjD, 1):
                if molLimit and ii > molLimit:
                    break
                # ----
                oemf = OeMoleculeFactory()
                oemf.setQuiet()
                tId = oemf.setChemCompDef(ccObjD[ccId])
                if tId != ccId:
                    logger.error("%s chemical component definition import error", ccId)
                smiD = oemf.buildRelated(limitPerceptions=limitPerceptions)
                logger.debug("%s related molecular forms %d", ccId, len(smiD))
                rD.update(smiD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return rD

    def __buildChemCompSearchIndexMulti(self, ccObjD, limitPerceptions=False, molLimit=None, numProc=2, maxChunkSize=20, quietFlag=False):
        #
        ccIdList = sorted(ccObjD.keys())[:molLimit] if molLimit else sorted(ccObjD.keys())
        logger.info("Input definition length %d numProc %d limitPerceptions %r", len(ccIdList), numProc, limitPerceptions)
        #
        rWorker = ChemCompSearchIndexWorker(ccObjD)
        # mpu = MultiProcPoolUtil(verbose=True)
        mpu = MultiProcUtil(verbose=True)
        optD = {"maxChunkSize": maxChunkSize, "limitPerceptions": limitPerceptions, "quietFlag": quietFlag}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="buildRelatedList")
        ok, failList, resultList, _ = mpu.runMulti(dataList=ccIdList, numProc=numProc, numResults=1, chunkSize=maxChunkSize)
        if failList:
            logger.info("Index definitions with failures (%d): %r", len(failList), failList)
        logger.info("Multi-proc status %r failures %r result length %r", ok, len(failList), len(resultList[0]))
        # JDW
        rD = {vD["name"]: vD for vD in resultList[0]}
        return rD

    def matchMolecularFormulaRange(self, typeRangeD, matchSubset=False):
        """Find matching formula for the input atom type range query (evaluates min <= ff <= max).

        Args:
            typeRangeD (dict): dictionary of element ranges {'<element_name>: {'min': <int>, 'max': <int>}}
            matchSubset (bool, optional): test for formula subset (default: False)

        Returns:
            (list):  chemical component identifiers with matching formula (MatchResults)
        """
        rL = []
        try:
            if not typeRangeD:
                return rL
            myTypeRangeD = {k.upper(): v for k, v in typeRangeD.items()}
            queryTypeS = set(myTypeRangeD.keys())
            for ccId, idxD in self.__searchIdx.items():
                tD = idxD["type-counts"]
                # targetTypeS = set(tD.keys())
                if not matchSubset and idxD["atom-types"] != queryTypeS:
                    continue
                #
                if not queryTypeS.issubset(idxD["atom-types"]):
                    continue
                match = True
                for atomType, rangeD in myTypeRangeD.items():
                    try:
                        if ("min" in rangeD and rangeD["min"] > tD[atomType]) or ("max" in rangeD and rangeD["max"] < tD[atomType]):
                            match = False
                            break
                    except Exception:
                        match = False
                        break
                if match:
                    # logger.info("%s formula %r query %r", ccId, idxD["type-counts"], typeRangeD)
                    rL.append(MatchResults(ccId=ccId, searchType="formula", formula=idxD["formula"]))
        except Exception as e:
            logger.exception("Failing for %r with %s", typeRangeD, str(e))
        return rL

    def filterMinimumMolecularFormula(self, typeCountD):
        """Find molecules with the minumum formula composition for the input atom type query (evaluates min <= ff).

        Args:
            typeCountD (dict): dictionary of element minimum values {'<element_name>: #}

        Returns:
            (list):  chemical component identifiers
        """
        rL = []
        try:
            if not typeCountD:
                return list(self.__searchIdx.keys())

            queryTypeS = set(typeCountD.keys())
            for ccId, idxD in self.__searchIdx.items():
                tD = idxD["type-counts"]
                if not queryTypeS.issubset(tD):
                    continue
                match = True
                for atomType, minCount in typeCountD.items():
                    try:
                        if minCount > tD[atomType]:
                            match = False
                            break
                    except Exception:
                        match = False
                        break
                if match:
                    rL.append(ccId)
        except Exception as e:
            logger.exception("Failing for %r with %s", typeCountD, str(e))
        return rL

    def filterMinimumFormulaAndFeatures(self, typeCountD, featureCountD):
        """Find molecules with the minumum formula and feature composition.

        Args:
            typeCountD (dict): dictionary of element minimum values {'<element_name>: #}
            featureCountD (dict): dictionary of feature minimum values {'<element_name>: #}

        Returns:
            (list):  chemical component identifiers
        """
        rL = []
        try:
            if not typeCountD or not featureCountD:
                return list(self.__searchIdx.keys())
            # ----
            featureQueryS = set(featureCountD.keys())
            typeQueryS = set(typeCountD.keys())
            #
            for ccId, idxD in self.__searchIdx.items():
                tD = idxD["type-counts"]
                fD = idxD["feature-counts"]
                #
                if not typeQueryS.issubset(tD) or not featureQueryS.issubset(fD):
                    continue

                match = True
                for atomType, minCount in typeCountD.items():
                    try:
                        if minCount > tD[atomType]:
                            match = False
                            break
                    except Exception:
                        match = False
                        break

                if not match:
                    continue
                #
                for featureType, minCount in featureCountD.items():
                    try:
                        if minCount > fD[featureType]:
                            match = False
                            break
                    except Exception:
                        match = False
                        break
                #
                if match:
                    rL.append(ccId)
        except Exception as e:
            logger.exception("Failing for %r with %s", typeCountD, str(e))
        return rL
