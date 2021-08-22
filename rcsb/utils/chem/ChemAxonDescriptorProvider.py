##
# File:    ChemAxonDescriptorProvider.py
# Author:  J. Westbrook
# Date:    17-Aug-2021
#
# Updates:
#
##
"""
Utilities to deliver ChemAxon rendered chemical descriptors for chemical component definitions.
"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import datetime
import time

from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil
from rcsb.utils.io.StashableBase import StashableBase


logger = logging.getLogger(__name__)


class ChemAxonDescriptorProvider(StashableBase):
    """Utilities to deliver ChemAxon rendered chemical descriptors for chemical component definitions."""

    def __init__(self, **kwargs):
        #
        dirName = "chemaxon"
        if "cachePath" in kwargs:
            self.__cachePath = os.path.abspath(kwargs.get("cachePath", None))
            self.__dirPath = os.path.join(self.__cachePath, dirName)
        super(ChemAxonDescriptorProvider, self).__init__(self.__cachePath, [dirName])
        #
        self.__molLimit = kwargs.get("molLimit", 0)
        self.__ccUrlTarget = kwargs.get("ccUrlTarget", None)
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", None)
        useCache = kwargs.get("useCache", True)
        self.__chunkSize = kwargs.get("chunkSize", 100)
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__ccFileNamePrefix = kwargs.get("ccFileNamePrefix", "cc-full")
        self.__version = None
        self.__descrD = self.__reload(useCache)

    def testCache(self, minCount=None):
        ok = self.__descrD and len(self.__descrD) >= minCount if minCount else self.__descrD is not None
        logger.info("Loaded ChemAxon descriptors for (%d) components (success %r)", len(self.__descrD) if self.__descrD else 0, ok)
        return ok

    def getDescriptorIndex(self):
        return self.__descrD

    def getIndexFilePath(self):
        return os.path.join(self.__dirPath, "%s-chemaxon-descriptors.json" % self.__ccFileNamePrefix)

    def getVersion(self):
        return self.__version

    def __reload(self, useCache):
        """Reload or created Chemaxon descriptor mapping index.

        Args:
            cachePath (str): path to the directory containing cache files
            chunkSize (int, optional): number of SMILES per request. Defaults to 100.

         Returns:
            (dict): chemical component data containers for each indexed chemical component
        """
        #
        descrD = {}
        descrFilePath = self.getIndexFilePath()
        #
        if not (useCache and self.__mU.exists(descrFilePath)):
            url = "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/CHEMAXON/cc-full-chemaxon-descriptors.json"
            _ = self.__fetchUrl(url, self.__dirPath)
        #
        _, fExt = os.path.splitext(descrFilePath)
        descrFormat = "json" if fExt == ".json" else "pickle"
        if self.__mU.exists(descrFilePath):
            dD = self.__mU.doImport(descrFilePath, fmt=descrFormat)
            descrD = dD["smiles"]
            self.__version = dD["version"]
        #
        return descrD

    def __fetchUrl(self, urlTarget, dirPath, useCache=False):
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        if not (useCache and fU.exists(filePath)):
            startTime = time.time()
            ok2 = fU.get(urlTarget, filePath)
            endTime = time.time()
            if ok2:
                logger.info("Fetched %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok2, endTime - startTime)
            else:
                logger.error("Failing fetch for %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok2, endTime - startTime)
        #
        return filePath

    def buildDescriptors(self):
        descrFilePath = self.getIndexFilePath()
        ccidxP = ChemCompIndexProvider(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            cachePath=self.__cachePath,
            useCache=True,
            molLimit=self.__molLimit,
            ccFileNamePrefix=self.__ccFileNamePrefix,
        )
        ok = ccidxP.testCache()
        if ok:
            ccIdList = ccidxP.getIdList()
            self.__descrD = self.__fetchDescriptors(ccIdList, ccidxP, chunkSize=self.__chunkSize)
            tS = datetime.datetime.now().isoformat()
            vS = datetime.datetime.now().strftime("%Y-%m-%d")
            self.__version = vS
            dD = {"created": tS, "version": vS, "smiles": self.__descrD}
            ok = self.__mU.doExport(descrFilePath, dD, fmt="json", indent=3)
            logger.info("Stored %s descriptors for %d components (status=%r) ", descrFilePath, len(self.__descrD), ok)

    def updateDescriptors(self, useCache=True):

        ccidxP = ChemCompIndexProvider(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            cachePath=self.__cachePath,
            useCache=useCache,
            molLimit=None,
            ccFileNamePrefix=self.__ccFileNamePrefix,
        )
        ok = ccidxP.testCache()
        if ok:
            ccIdList = ccidxP.getIdList()
            curIdList = list(self.__descrD.keys())
            updIdList = list(set(ccIdList) - set(curIdList))
            if updIdList:
                logger.info("Updating Chemaxon descriptors for (%d) components", len(updIdList))
                uD = self.__fetchDescriptors(updIdList, ccidxP, chunkSize=self.__chunkSize)
                self.__descrD.update(uD)
                descrFilePath = self.getIndexFilePath()
                tS = datetime.datetime.now().isoformat()
                vS = datetime.datetime.now().strftime("%Y-%m-%d")
                self.__version = vS
                dD = {"created": tS, "version": vS, "smiles": self.__descrD}
                ok = self.__mU.doExport(descrFilePath, dD, fmt="json", indent=3)
        #
        return ok

    def __fetchDescriptors(self, ccIdList, ccidxP, chunkSize=100):
        """Fetch transformed SMILES descriptors from the ChemAxon webservice.

            Args:
                ccIdList (list, str): chemical component identifier list
                ccidxP (object): instance of the ChemCompIndexProvider()
                chunksize (int, optional): number of SMILES per request. Defaults to 100.

            Returns:
                (dict): dictionary {<ccId>: [<transformed SMILES>, ...], ...}

        Example API parameter data:
                            {
                            "errorHandlingMode": "FAIL_ON_ERROR",
                            "inputParams": "smiles",
                            "outputParams": "smiles",
                            "structures": [
                                "CC(C)[C@H](N)C=O",
                                "CC[C@H](C)[C@H](N)C=O",
                                "CC(C)C[C@H](N)C=O"
                            ]
                            }

        Example query:
        curl -X POST "https://jchem-microservices.chemaxon.com/jwsio/rest-v1/molconvert/batch" -H "accept: */*"
               -H "Content-Type: application/json" -d "{ \"errorHandlingMode\": \"FAIL_ON_ERROR\", \"inputParams\": \"smiles\",
               \"outputParams\": \"mrv\", \"structures\": [ \"CC(C)[C@H](N)C=O\", \"CC[C@H](C)[C@H](N)C=O\", \"CC(C)C[C@H](N)C=O\" ]}"
        """
        descrD = {}
        smilesCcIdD = {}
        smilesD = {}
        for ccId in ccIdList:
            smiL = list(set(ccidxP.getSMILES(ccId, smiTypeList=["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles"])))
            smilesCcIdD.setdefault(ccId, []).extend(smiL)
            for smi in smiL:
                smilesD.setdefault(smi, []).append(ccId)
        #
        logger.info("Translating (%d) SMILES for components (%d)", len(smilesD), len(smilesCcIdD))
        # ----
        smiLL = [list(smilesD.keys())[i : i + chunkSize] for i in range(0, len(smilesD), chunkSize)]
        # ---
        baseUrl = "https://jchem-microservices.chemaxon.com"
        endPoint = "jwsio/rest-v1/molconvert/batch"
        # hL = [("Accept", "application/json"), ("Content-Type", "application/json")]
        hD = {"Accept": "application/json", "Content-Type": "application/json"}
        try:
            pD = {"errorHandlingMode": "SKIP_ERROR", "inputParams": "smiles", "outputParams": "smiles"}
            #
            iCount = 0
            for smiL in smiLL:
                iCount += 1
                ureq = UrlRequestUtil()
                pD["structures"] = smiL
                logger.debug("pD %r", pD)
                rDL, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, headers=hD, sendContentType="application/json", returnContentType="application/json")
                logger.debug("API result (%r) %r", retCode, rDL)
                if rDL and len(rDL) == len(smiL):
                    for ii, rD in enumerate(rDL):
                        if "structure" in rD and "successful" in rD and rD["successful"]:
                            if smiL[ii] == rD["structure"]:
                                continue
                            for ccId in smilesD[smiL[ii]]:
                                if ccId in descrD and rD["structure"] in descrD[ccId]:
                                    continue
                                if rD["structure"] in smilesCcIdD[ccId]:
                                    continue
                                descrD.setdefault(ccId, []).append(rD["structure"])
                else:
                    logger.info("Chunk %d failed (%d)", iCount, len(rDL))
                if iCount % 10 == 0:
                    logger.info("Completed processing chunk (%d/%d)", iCount, len(smiLL))

            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return descrD
