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
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.SingletonClass import SingletonClass

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logger = logging.getLogger(__name__)

MatchResults = namedtuple("MatchResults", "ccId oeMol searchType matchOpts screenType fpType fpScore oeIdx formula", defaults=(None,) * 9)


class ChemCompSearchWrapper(SingletonClass):
    """ Wrapper for chemical component search operations.
    """

    def __init__(self):
        self.__startTime = time.time()
        self.__cachePath = os.environ["CHEM_SEARCH_CACHE_PATH"]
        self.__ccFileNamePrefix = os.environ["CHEM_SEARCH_CC_PREFIX"]
        # ---
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        # ---
        self.__configD = {}
        self.__ccIdxP = None
        self.__siIdx = {}
        self.__oesmP = None
        self.__oesU = None
        # ---
        self.__statusDescriptorError = -100
        self.__searchError = -200
        self.__searchSuccess = 0

    def readConfig(self):
        #
        ok = False
        try:
            self.__cachePath = os.environ.get("CHEM_SEARCH_CACHE_PATH", ".")
            self.__ccFileNamePrefix = os.environ.get("CHEM_SEARCH_CC_PREFIX", "cc-full")
            configFilePath = os.path.join(self.__cachePath, "config", self.__ccFileNamePrefix + "-config.json")
            configD = self.__mU.doImport(configFilePath, fmt="json")
            logger.debug("ConfigD: %r", configD)
            if configD and (len(configD) > 2) and float(configD["versionNumber"]) > 0.1:
                logger.info("Read version %r sections %r from %s", configD["versionNumber"], list(configD.keys()), configFilePath)
                ok = True
                self.__configD = configD
            else:
                logger.error("Reading config file fails from %r", configFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
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
            ok = oesU.testCache()
            self.__oesU = oesU if ok else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def matchByDescriptor(self, descriptor, descriptorType, matchOpts="graph-relaxed", searchId=None):
        """Return graph match (w/  finger print pre-filtering) and finger print search results for the
           input desriptor.

        Args:
            descriptor (str):  molecular descriptor (SMILES, InChI)
            descriptorType (str): descriptor type (SMILES, InChI
            matchOpts (str, optional): graph match criteria (graph-relaxed, graph-relaxed-stereo, graph-strict, fingerprint-similarity). Defaults to "graph-relaxed".
            searchId (str, optional): search identifier for logging. Defaults to None.

        Returns:
            (statusCode, list, list): status, graph match and finger match lists of type (MatchResults)
                                      -200 descriptor processing error
                                      -100 search execution error
                                         0 search execution success
        """
        ssL = fpL = []
        retStatus = False
        try:
            fpTypeCuttoffD = self.__configD["oesmpKwargs"]["fpTypeCuttoffD"] if "fpTypeCuttoffD" in self.__configD["oesmpKwargs"] else {}
            maxFpResults = self.__configD["oesmpKwargs"]["maxFpResults"] if "maxFpResults" in self.__configD["oesmpKwargs"] else 50
            limitPerceptions = self.__configD["oesmpKwargs"]["limitPerceptions"] if "limitPerceptions" in self.__configD["oesmpKwargs"] else False
            #
            searchId = searchId if searchId else "query"
            messageTag = searchId + ":" + descriptorType
            oeioU = OeIoUtils()
            oeMol = oeioU.descriptorToMol(descriptor, descriptorType, limitPerceptions=limitPerceptions, messageTag=messageTag)
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
