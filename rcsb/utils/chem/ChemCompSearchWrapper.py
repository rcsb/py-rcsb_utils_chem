##
# File:    ChemCompSearchWrapper.py
# Author:  jdw
# Date:    9-Mar-2020
# Version: 0.001
#
# Updates:
#
##
"""
Wrapper for chemical component search operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import copy
import logging
import platform
import resource
import os
import time

from collections import namedtuple

from rcsb.utils.chem.ChemCompIndexProvider import ChemCompIndexProvider
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.MolecularFormula import MolecularFormula
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeSearchUtils import OeSearchUtils
from rcsb.utils.chem.OeSubStructSearchUtils import OeSubStructSearchUtils
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.SftpUtil import SftpUtil
from rcsb.utils.io.SingletonClass import SingletonClass

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logger = logging.getLogger(__name__)

MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore oeIdx formula", defaults=(None,) * 9)


class ChemCompSearchWrapper(SingletonClass):
    """Wrapper for chemical component search operations."""

    def __init__(self, **kwargs):
        """Wrapper class for chemical search/depiction operations.

        Path and prefix data for wrapper class may be set as keyword arguments
        as environmental variables.

        Args:
            cachePath (str): path to top-level cache directory used to store search index file dependencies
                             (default environment variable CHEM_SEARCH_CACHE_PATH or ".")
            ccFileNamePrefix (str): prefix code used to distinguish different subsets of chemical definitions
                                    (default environment variable CHEM_SEARCH_CC_PREFIX or "cc-full")

        """
        self.__startTime = time.time()
        #
        self.__cachePath = kwargs.get("cachePath", os.environ.get("CHEM_SEARCH_CACHE_PATH", "."))
        self.__ccFileNamePrefix = kwargs.get("ccFileNamePrefix", os.environ.get("CHEM_SEARCH_CC_PREFIX", "cc-full"))
        #
        self.__dependFileName = "ChemCompSearchWrapperData.tar.gz"
        self.__dependTarFilePath = os.path.join(self.__cachePath, self.__dependFileName)
        # ---
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        # ---
        self.__configD = {}
        self.__ccIdxP = None
        self.__siIdxP = None
        self.__siIdx = {}
        self.__oesmP = None
        self.__oesU = None
        self.__oesubsU = None
        # ---
        self.__statusDescriptorError = -100
        self.__searchError = -200
        self.__searchSuccess = 0

    def setConfig(self, ccUrlTarget, birdUrlTarget, **kwargs):
        """Provide the chemical definition source path details for rebuilding search
           index file dependencies.

        Args:
            ccUrlTarget (str): path to concatenated chemical component definition file
            birdUrlTarget (str): path to the concatenated BIRD definition file

            Other options are propagated to configurations of the wrapped classes in __bootstrapConfig()

        """
        kwargs["ccUrlTarget"] = ccUrlTarget
        kwargs["birdUrlTarget"] = birdUrlTarget
        kwargs["cachePath"] = self.__cachePath
        kwargs["ccFileNamePrefix"] = self.__ccFileNamePrefix
        self.__configD = self.__bootstrapConfig(**kwargs)
        return len(self.__configD) >= 3

    def __bootstrapConfig(self, **kwargs):
        """Build on-the-fly default configuration for this wrapper class."""
        # The following few options have no defaults -- and should be specified.
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        cachePath = kwargs.get("cachePath", None)
        ccFileNamePrefix = kwargs.get("ccFileNamePrefix", None)
        logger.info("Bootstrap configuration for prefix %r cc %r bird %r", ccFileNamePrefix, ccUrlTarget, birdUrlTarget)
        # ---
        #  Reasonable values are selected for the remaining options...
        oeFileNamePrefix = "oe-" + ccFileNamePrefix
        try:
            storeConfig = kwargs.get("storeConfig", True)
            molLimit = kwargs.get("molLimit", None)
            useCache = kwargs.get("useCache", False)
            logSizes = kwargs.get("logSizes", False)
            #
            numProc = kwargs.get("numProc", 12)
            maxProc = os.cpu_count()
            numProc = min(numProc, maxProc)
            maxChunkSize = kwargs.get("maxChunkSize", 50)
            #
            logger.debug("+++ >>> Assigning numProc as %d", numProc)
            #
            limitPerceptions = kwargs.get("limitPerceptions", False)
            quietFlag = kwargs.get("quietFlag", True)
            #
            # fpTypeCuttoffD = {"TREE": 0.6, "MACCS": 0.9, "PATH": 0.6, "CIRCULAR": 0.6, "LINGO": 0.9}
            fpTypeCuttoffD = kwargs.get("fpTypeCuttoffD", {"TREE": 0.6, "MACCS": 0.9})
            buildTypeList = kwargs.get("buildTypeList", ["oe-iso-smiles", "oe-smiles", "cactvs-iso-smiles", "cactvs-smiles", "inchi"])
            #
            oesmpKwargs = {
                "ccUrlTarget": ccUrlTarget,
                "birdUrlTarget": birdUrlTarget,
                "cachePath": cachePath,
                "useCache": useCache,
                "ccFileNamePrefix": ccFileNamePrefix,
                "oeFileNamePrefix": oeFileNamePrefix,
                "limitPerceptions": limitPerceptions,
                "minCount": None,
                "maxFpResults": 50,
                "fpTypeCuttoffD": fpTypeCuttoffD,
                "buildTypeList": buildTypeList,
                "screenTypeList": None,
                "quietFlag": quietFlag,
                "numProc": numProc,
                "maxChunkSize": maxChunkSize,
                "molLimit": molLimit,
                "logSizes": logSizes,
                "suppressHydrogens": True,
            }
            ccsiKwargs = {
                "ccUrlTarget": ccUrlTarget,
                "birdUrlTarget": birdUrlTarget,
                "cachePath": cachePath,
                "useCache": useCache,
                "ccFileNamePrefix": ccFileNamePrefix,
                "oeFileNamePrefix": oeFileNamePrefix,
                "limitPerceptions": limitPerceptions,
                "minCount": None,
                "numProc": numProc,
                "quietFlag": quietFlag,
                "maxChunkSize": maxChunkSize,
                "molLimit": None,
                "logSizes": False,
            }
            configD = {"versionNumber": 0.30, "ccsiKwargs": ccsiKwargs, "oesmpKwargs": oesmpKwargs}
            #
            if storeConfig:
                configDirPath = os.path.join(cachePath, "config")
                configFilePath = os.path.join(configDirPath, ccFileNamePrefix + "-config.json")
                logger.info("Saving configuration bootstrap in %r", configFilePath)
                self.__mU.mkdir(configDirPath)
                self.__mU.doExport(configFilePath, configD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return configD

    def readConfig(self, resetCachePath=True):
        """Read a prepared configuration file for the search wrapper class. This will override
        any default configuration settings.

        Args:
             resetCachPath (bool): update cachePath configuration option with the current cachePath setting.

        Returns:
            bool : True for success or False otherwise
        """
        #
        #
        ok = False
        try:
            #
            configFilePath = os.path.join(self.__cachePath, "config", self.__ccFileNamePrefix + "-config.json")
            configD = self.__mU.doImport(configFilePath, fmt="json")
            logger.debug("ConfigD: %r", configD)
            if configD and (len(configD) > 2) and float(configD["versionNumber"]) > 0.2:
                logger.info("Read version %r sections %r from %s", configD["versionNumber"], list(configD.keys()), configFilePath)
                ok = True
                self.__configD = configD
                if resetCachePath:
                    # Allow the configuration to be relocatable.
                    configD["ccsiKwargs"]["cachePath"] = self.__cachePath
                    configD["oesmpKwargs"]["cachePath"] = self.__cachePath
            else:
                logger.error("Reading config file fails from %r", configFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    def buildDependenices(self, ccUrlTarget, birdUrlTarget, **kwargs):
        """Convenience method to build configuration and static dependencies for the chemical search services.

        Args:
            ccUrlTarget (str): path to source concatenated chemical component definition file
            birdUrlTarget (str): path to the source concatenated BIRD definition file

            Other options are propagated to configurations of the wrapped classes in __bootstrapConfig()

        """
        try:
            okT = False
            ok1 = self.setConfig(ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget, **kwargs)
            useCache = kwargs.get("useCache", False)
            ok2 = self.updateChemCompIndex(useCache=useCache)
            ok3 = self.updateSearchIndex(useCache=useCache)
            ok4 = self.updateSearchMoleculeProvider(useCache=useCache)
            okBuild = ok1 and ok2 and ok3 and ok4
            if okBuild:
                fileU = FileUtil()
                dirPathList = [os.path.join(self.__cachePath, subDir) for subDir in ["chem_comp", "oe_mol", "config"]]
                okT = fileU.bundleTarfile(self.__dependTarFilePath, dirPathList, mode="w:gz", recursive=True)
            #
            return okT and okBuild
        except Exception as e:
            logger.exception("Failing build with %r and %r with %s", ccUrlTarget, birdUrlTarget, str(e))
        return False

    def stashDependencies(self, url, dirPath, bundleLabel="A", userName=None, pw=None):
        """Store a copy of the bundled search dependencies remotely -

        Args:
            url (str): URL string for the destination host (e.g. sftp://myserver.net or None for a local file)
            dirPath (str): directory path on the remote resource
            bundleLabel (str, optional): optional label preppended to the stashed dependency bundle artifact (default='A')
            userName (str, optional): optional access information. Defaults to None.
            password (str, optional): optional access information. Defaults to None.

        Returns:
          bool:  True for success or False otherwise

        """
        try:
            ok = False
            fn = self.__makeBundleFileName(self.__dependFileName, bundleLabel=bundleLabel)
            if url and url.startswith("sftp://"):
                sftpU = SftpUtil()
                hostName = url[7:]
                ok = sftpU.connect(hostName, userName, pw=pw, port=22)
                if ok:
                    remotePath = os.path.join("/", dirPath, fn)
                    ok = sftpU.put(self.__dependTarFilePath, remotePath)
            elif not url:
                fileU = FileUtil()
                remotePath = os.path.join(dirPath, fn)
                ok = fileU.put(self.__dependTarFilePath, remotePath)
            else:
                logger.error("Unsupported stash protocol %r", url)
            return ok
        except Exception as e:
            logger.exception("For %r %r failing with %s", url, dirPath, str(e))
        return False

    def __makeBundleFileName(self, rootName, bundleLabel="A"):
        fn = rootName
        try:
            fn = rootName
            fn = "%s-%s" % (bundleLabel.upper(), rootName) if bundleLabel else rootName
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return fn

    def restoreDependencies(self, url, dirPath, bundleLabel="A", userName=None, pw=None):
        """Restore bundled dependencies from remote storage and unbundle these in the
           current local cache directory.

        Args:
            url (str): remote URL
            dirPath (str): remote directory path on the
            bundleLabel (str, optional): optional label preppended to the stashed dependency bundle artifact (default='A')
            userName (str, optional): optional access information. Defaults to None.
            password (str, optional): optional access information. Defaults to None.
        """
        try:
            ok = False
            fileU = FileUtil()
            fn = self.__makeBundleFileName(self.__dependFileName, bundleLabel=bundleLabel)
            if not url:
                remotePath = os.path.join(dirPath, fn)
                ok = fileU.get(remotePath, self.__dependTarFilePath)

            elif url and url.startswith("http://"):
                remotePath = url + os.path.join("/", dirPath, fn)
                ok = fileU.get(remotePath, self.__dependTarFilePath)

            elif url and url.startswith("sftp://"):
                sftpU = SftpUtil()
                ok = sftpU.connect(url[7:], userName, pw=pw, port=22)
                if ok:
                    remotePath = os.path.join(dirPath, fn)
                    ok = sftpU.get(remotePath, self.__dependTarFilePath)
            else:
                logger.error("Unsupported protocol %r", url)
            if ok:
                ok = fileU.unbundleTarfile(self.__dependTarFilePath, dirPath=self.__cachePath)
            return ok
        except Exception as e:
            logger.exception("For %r %r Failing with %s", url, dirPath, str(e))
            ok = False
        return ok

    def updateChemCompIndex(self, useCache=False):
        """Rebuild the basic index of source chemical component and BIRD definitions.
           Update the internal state of this index in the current object instance.

            Resource requirements: 94 sec 1 proc 7GB memory macbook pro

        Args:
            useCache (bool): False to rebuild search index and True to reload

        Returns:
            bool: True for success or false otherwise
        """
        ok = False
        try:
            kwargs = copy.deepcopy(self.__configD["ccsiKwargs"]) if "ccsiKwargs" in self.__configD else None
            if kwargs:
                kwargs["useCache"] = useCache
                ccIdxP = ChemCompIndexProvider(**kwargs)
                ok = ccIdxP.testCache()
                self.__ccIdxP = ccIdxP if ok else None
                logger.info("Chemical component index status %r", ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def getChemCompIndex(self):
        return self.__ccIdxP.getIndex() if self.__ccIdxP else {}

    def getSearchMoleculeProvider(self):
        return self.__oesmP if self.__oesmP else None

    def updateSearchIndex(self, useCache=False):
        """Rebuild the search index from source chemical component and BIRD definitions.
           Update the internal state of this index in the current object instance.

            Resource requirements 771 secs 6 proc macbook pro 7GB memory.

        Args:
            useCache (bool): False to rebuild search index and True to reload

        Returns:
            bool: True for success or false otherwise
        """
        ok = False
        try:
            kwargs = copy.deepcopy(self.__configD["ccsiKwargs"]) if "ccsiKwargs" in self.__configD else None
            if kwargs:
                kwargs["useCache"] = useCache
                siIdxP = ChemCompSearchIndexProvider(**kwargs)
                ok = siIdxP.testCache()
                self.__siIdxP = siIdxP if siIdxP else None
                self.__siIdx = siIdxP.getIndex() if siIdxP and ok else {}
                logger.info("Search index status %r index len %d", ok, len(self.__siIdx) if self.__siIdx else 0)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def updateSearchMoleculeProvider(self, useCache=False):
        """Rebuild the search molecule provider.
           Update the internal state of this object reference in the current object instance.

           Resource requirements: 151 seconds 1 proc  0.5GB memory macbook pro

        Args:
            useCache (bool): False to rebuild molecule store and True to reload

        Returns:
            bool: True for success or false otherwise
        """
        ok = False
        try:
            kwargs = copy.deepcopy(self.__configD["oesmpKwargs"]) if "oesmpKwargs" in self.__configD else None
            if kwargs:
                kwargs["useCache"] = useCache
                oesmP = OeSearchMoleculeProvider(**kwargs)
                ok = oesmP.testCache()
                self.__oesmP = oesmP if oesmP and ok else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reloadSearchDatabase(self):
        """Reload the in-memory search databases from the OE molecule provider.
           Resource requirements: ~90sec load time 0.35 GB memory

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            okmp = self.updateSearchMoleculeProvider(useCache=True)
            if not okmp:
                return ok
            fpTypeCuttoffD = self.__configD["oesmpKwargs"]["fpTypeCuttoffD"] if "fpTypeCuttoffD" in self.__configD["oesmpKwargs"] else {}
            fpTypeList = [k for k, v in fpTypeCuttoffD.items()]
            oesU = OeSearchUtils(self.__oesmP, fpTypeList=fpTypeList)
            ok1 = oesU.testCache()
            self.__oesU = oesU if ok1 else None
            #
            oesubsU = OeSubStructSearchUtils(self.__oesmP)
            ok2 = oesubsU.testCache()
            self.__oesubsU = oesubsU if ok2 else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok1 and ok2

    def searchByDescriptor(self, descriptor, descriptorType, matchOpts="graph-relaxed", searchId=None):
        """Wrapper method for descriptor match and descriptor substructure search methods.

        Args:
            descriptor (str):  molecular descriptor (SMILES, InChI)
            descriptorType (str): descriptor type (SMILES, InChI
            matchOpts (str, optional): graph match criteria (graph-relaxed, graph-relaxed-stereo, graph-strict,
                                       fingerprint-similarity, sub-struct-graph-relaxed, sub-struct-graph-relaxed-stereo,
                                       sub-struct-graph-strict Defaults to "graph-relaxed")
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list, list): status, graph match and finger match lists of type (MatchResults)
                                      -100 descriptor processing error
                                      -200 search execution error
                                         0 search execution success
        """
        if matchOpts.startswith("sub-struct-"):
            return self.subStructSearchByDescriptor(descriptor, descriptorType, matchOpts=matchOpts, searchId=searchId)
        else:
            return self.matchByDescriptor(descriptor, descriptorType, matchOpts=matchOpts, searchId=searchId)

    def matchByDescriptor(self, descriptor, descriptorType, matchOpts="graph-relaxed", searchId=None):
        """Return graph match (w/  finger print pre-filtering) and finger print search results for the
           input desriptor.

        Args:
            descriptor (str):  molecular descriptor (SMILES, InChI)
            descriptorType (str): descriptor type (SMILES, InChI
            matchOpts (str, optional): graph match criteria (graph-relaxed, graph-relaxed-stereo, graph-strict,
                                       fingerprint-similarity, Defaults to "graph-relaxed")
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list, list): status, graph match and finger match lists of type (MatchResults)
                                      -100 descriptor processing error
                                      -200 search execution error
                                         0 search execution success
        """
        ssL = fpL = []
        retStatus = False
        statusCode = -200
        try:
            fpTypeCuttoffD = self.__configD["oesmpKwargs"]["fpTypeCuttoffD"] if "fpTypeCuttoffD" in self.__configD["oesmpKwargs"] else {}
            maxFpResults = self.__configD["oesmpKwargs"]["maxFpResults"] if "maxFpResults" in self.__configD["oesmpKwargs"] else 50
            limitPerceptions = self.__configD["oesmpKwargs"]["limitPerceptions"] if "limitPerceptions" in self.__configD["oesmpKwargs"] else False
            #
            searchId = searchId if searchId else "query"
            messageTag = searchId + ":" + descriptorType
            oeioU = OeIoUtils()
            oeMol = oeioU.descriptorToMol(descriptor, descriptorType, limitPerceptions=limitPerceptions, messageTag=messageTag)
            oeMol = oeioU.suppressHydrogens(oeMol)
            if not oeMol:
                logger.warning("descriptor type %r molecule build fails: %r", descriptorType, descriptor)
                return self.__statusDescriptorError, ssL, fpL
            #
            retStatus, ssL, fpL = self.__oesU.searchSubStructureAndFingerPrint(oeMol, list(fpTypeCuttoffD.items())[:2], maxFpResults, matchOpts=matchOpts)
            statusCode = 0 if retStatus else self.__searchError
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            #
        return statusCode, ssL, fpL

    def subStructSearchByDescriptor(self, descriptor, descriptorType, matchOpts="sub-struct-graph-relaxed", searchId=None):
        """Return graph match (w/  finger print pre-filtering) and finger print search results for the
           input desriptor.

        Args:
            descriptor (str):  molecular descriptor (SMILES, InChI)
            descriptorType (str): descriptor type (SMILES, InChI)
            matchOpts (str, optional): graph match criteria (sub-struct-graph-relaxed, sub-struct-graph-relaxed-stereo,
                                       sub-struct-graph-strict). Defaults to "sub-struct-graph-relaxed".
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list, list): status, substructure search results of type (MatchResults), empty list placeholder
                                      -100 descriptor processing error
                                      -200 search execution error
                                         0 search execution success
        """
        ssL = []
        retStatus = False
        statusCode = -200
        try:
            limitPerceptions = self.__configD["oesmpKwargs"]["limitPerceptions"] if "limitPerceptions" in self.__configD["oesmpKwargs"] else False
            numProc = self.__configD["oesmpKwargs"]["numProc"] if "numProc" in self.__configD["oesmpKwargs"] else 4
            #
            searchId = searchId if searchId else "query"
            messageTag = searchId + ":" + descriptorType
            oeioU = OeIoUtils()
            oeMol = oeioU.descriptorToMol(descriptor, descriptorType, limitPerceptions=limitPerceptions, messageTag=messageTag)
            oeMol = oeioU.suppressHydrogens(oeMol)
            if not oeMol:
                logger.warning("descriptor type %r molecule build fails: %r", descriptorType, descriptor)
                return self.__statusDescriptorError, ssL, []
            #
            ccIdL = self.__oesubsU.prefilterIndex(oeMol, self.__siIdxP, matchOpts=matchOpts)
            retStatus, ssL = self.__oesubsU.searchSubStructure(oeMol, ccIdList=ccIdL, matchOpts=matchOpts, numProc=numProc)
            statusCode = 0 if retStatus else self.__searchError
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            #
        return statusCode, ssL, []

    def matchByFormulaRange(self, elementRangeD, matchSubset=False, searchId=None):
        """Return formula match results for input element range dictionary.

        Args:
            elementRangeD (dict): {'<element_name>: {'min': <int>, 'max': <int>}, ... }
            matchSubset (bool, optional): query for formula subset (default: False)
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list): status, list of chemical component identifiers
        """
        ok = False
        rL = []
        try:
            startTime = time.time()
            searchId = searchId if searchId else "query"
            rL = self.__ccIdxP.matchMolecularFormulaRange(elementRangeD, matchSubset=matchSubset)
            ok = True
            logger.info("%s formula %r matched %d (%.4f seconds)", searchId, elementRangeD, len(rL), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok, rL

    def matchByFormula(self, formula, matchSubset=False, searchId=None):
        """Return formula match results for input molecular formula.

        Args:
            formula (str): molecular formula  (ex. 'C6H6')
            matchSubset (bool, optional): query for formula subset (default: False)
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list): status, list of chemical component identifiers
        """
        ok = False
        rL = []
        try:
            startTime = time.time()
            searchId = searchId if searchId else "query"
            mf = MolecularFormula()
            eD = mf.parseFormula(formula)
            elementRangeD = {k.upper(): {"min": v, "max": v} for k, v in eD.items()}
            rL = self.__ccIdxP.matchMolecularFormulaRange(elementRangeD, matchSubset=matchSubset)
            ok = True
            logger.info("%s formula %r matched %d (%.4f seconds)", searchId, elementRangeD, len(rL), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok, rL

    def status(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Status at %s (up %.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)
