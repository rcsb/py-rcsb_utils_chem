##
# File:    ChemCompSearchIndexProvider.py
# Author:  J. Westbrook
# Date:    3-Mar-2020
#
# Updates:
#
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

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.OeMoleculeFactory import OeMoleculeFactory
from rcsb.utils.io.IoUtil import getObjSize
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class ChemCompSearchIndexProvider(SingletonClass):
    """Utilities to read and process the index of chemical component definitions search targets
    """

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "chem_comp")
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__searchIdx = self.__reload(**kwargs)

    def testCache(self, minCount=29000, logSizes=False):
        if logSizes and self.__searchIdx:
            logger.info("ccIdxD (%.2f MB)", getObjSize(self.__searchIdx) / 1000000.0)
        ok = self.__searchIdx and len(self.__searchIdx) >= minCount if minCount else self.__searchIdx is not None
        return ok

    def getIndex(self):
        return self.__searchIdx

    def getMol(self, ccId):
        try:
            return self.__searchIdx[ccId]
        except Exception as e:
            logger.debug("Get molecule %r failing with %s", ccId, str(e))
        return None

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
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc")
        searchIdxFilePath = os.path.join(self.__dirPath, "%s-search-idx-components.json" % ccFileNamePrefix)
        #
        if useCache and self.__mU.exists(searchIdxFilePath):
            _, fExt = os.path.splitext(searchIdxFilePath)
            searchIdxFormat = "json" if fExt == "json" else "pickle"
            rdCcIdxD = self.__mU.doImport(searchIdxFilePath, fmt=searchIdxFormat)
            searchIdxD = {k: rdCcIdxD[k] for k in sorted(rdCcIdxD.keys())[:molLimit]} if molLimit else rdCcIdxD
        else:
            cmpKwargs = {k: v for k, v in kwargs.items() if k not in ["cachePath", "useCache", "molLimit"]}
            ccmP = ChemCompMoleculeProvider(cachePath=self.__cachePath, useCache=True, molLimit=molLimit, **cmpKwargs)
            ok = ccmP.testCache(minCount=molLimit, logSizes=True)
            if ok:
                searchIdxD = self.__updateChemCompSearchIndex(ccmP.getMolD(), searchIdxFilePath)
                logger.info("Storing %s with data for %d definitions (status=%r) ", searchIdxFilePath, len(searchIdxD), ok)
        #
        return searchIdxD

    def __updateChemCompSearchIndex(self, ccObjD, filePath):
        idxD = {}
        try:
            # Serialized index of chemical component search targets
            startTime = time.time()
            _, fExt = os.path.splitext(filePath)
            fileFormat = "json" if fExt == "json" else "pickle"
            idxD = self.__buildChemCompSearchIndex(ccObjD)
            ok = self.__mU.doExport(filePath, idxD, fmt=fileFormat)
            endTime = time.time()
            logger.info("Storing %s with %d raw indexed definitions (status=%r) (%.4f seconds)", filePath, len(idxD), ok, endTime - startTime)
        #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return idxD

    def __buildChemCompSearchIndex(self, cD):
        """Internal method return a dictionary of extracted chemical component descriptors and formula.
        """
        rD = {}
        try:
            for ccId, dataContainer in cD.items():
                # ----
                oemf = OeMoleculeFactory()
                tId = oemf.setChemCompDef(dataContainer)
                if tId != ccId:
                    logger.error("%s chemical component definition import error", ccId)
                smiD = oemf.buildRelated(limitPerceptions=False)
                logger.info("%s related molecular forms %d", ccId, len(smiD))
                rD.update(smiD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return rD
