##
# File:    ChemCompDepictWrapper.py
# Author:  jdw
# Date:    9-Mar-2020
# Version: 0.001
#
# Updates:
#
##
"""
Wrapper for chemical component depiction operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time

from rcsb.utils.chem.ChemCompSearchWrapper import ChemCompSearchWrapper
from rcsb.utils.chem.OeDepict import OeDepict
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.SingletonClass import SingletonClass

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logger = logging.getLogger(__name__)


class ChemCompDepictWrapper(SingletonClass):
    """ Wrapper for chemical component depiction operations.
    """

    def __init__(self):
        self.__startTime = time.time()
        # ---
        self.__workPath = "."
        self.__mU = MarshalUtil(workPath=self.__workPath)
        self.__configD = None
        self.__cachePath = None
        # ---
        self.__oesU = None
        # ---
        self.__statusDescriptorError = -100
        self.__searchError = -200
        self.__searchSuccess = 0
        self.__imageCount = 0

    def readConfig(self):
        #
        ok = False
        try:
            self.__cachePath = os.environ.get("CHEM_DEPICT_CACHE_PATH", ".")
            configFileName = os.environ.get("CHEM_DEPICT_CONFIG_FILE_NAME", "depict-config.json")
            configFilePath = os.path.join(self.__cachePath, "config", configFileName)
            configD = self.__mU.doImport(configFilePath, fmt="json")
            logger.debug("configD: %r", configD)
            if configD and (len(configD) >= 2) and float(configD["versionNumber"]) > 0.1:
                logger.info("Read version %r sections %r from %s", configD["versionNumber"], list(configD.keys()), configFilePath)
                ok = True
                self.__configD = configD
            else:
                logger.error("Reading config file fails from path %r", configFilePath)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    def setImageCount(self, imageCount):
        self.__imageCount = imageCount

    def getImageCount(self):
        return self.__imageCount

    def __makeImagePath(self):
        imageDirPath = self.__configD["imageDirPath"] if self.__configD and "imageDirPath" in self.__configD else "."
        fileRotateIncrement = self.__configD["fileRotateIncrement"] if self.__configD and "fileRotateIncrement" in self.__configD else 50
        ic = self.__imageCount % fileRotateIncrement
        imagePath = os.path.join(imageDirPath, "image-%s.svg" % ic)
        return imagePath

    def depictMolecule(self, identifier, identifierType, imagePath=None, **kwargs):
        """Create depiction from InChI, SMILES descriptors or PDB identifier.
        """
        try:
            imagePath = imagePath if imagePath else self.__makeImagePath()
            oeio = OeIoUtils()
            if identifierType.lower() in ["smiles"]:
                oeMol = oeio.smilesToMol(identifier)
            elif identifierType.lower() in ["inchi"]:
                oeMol = oeio.inchiToMol(identifier)
            elif identifierType.lower() in ["identifierpdb"]:
                ccsw = ChemCompSearchWrapper()
                oesmP = ccsw.getSearchMoleculeProvider()
                oeMol = oesmP.getMol(identifier)
            #
            ok = self.__depictOne(oeMol, imagePath, **kwargs)
            return imagePath if ok else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def __depictOne(self, oeMol, imagePath, **kwargs):
        """Single

        Args:
            oeMol (object): instance of an OE graph molecule
            imagePath (string): file path for image

        Returns:
            bool: True for success or False otherwise
        """
        try:
            title = kwargs.get("title", None)
            oed = OeDepict()
            oed.setMolTitleList([("Target", oeMol, title)])

            # ---
            bondDisplayWidth = 10.0
            numAtoms = oeMol.NumAtoms()
            if numAtoms > 100 and numAtoms <= 200:
                bondDisplayWidth = 6.0
            elif numAtoms > 200:
                bondDisplayWidth = 4.0
            # ---
            oed.setDisplayOptions(
                imageSizeX=kwargs.get("imageSizeX", 2500),
                imageSizeY=kwargs.get("imageSizeX", 2500),
                labelAtomName=kwargs.get("labelAtomName", False),
                labelAtomCIPStereo=kwargs.get("labelAtomCIPStereo", True),
                labelAtomIndex=kwargs.get("labelAtomIndex", False),
                labelBondIndex=kwargs.get("labelBondIndex", False),
                labelBondCIPStereo=kwargs.get("labelBondCIPStereo", True),
                cellBorders=kwargs.get("cellBorders", True),
                bondDisplayWidth=bondDisplayWidth,
            )
            oed.setGridOptions(rows=1, cols=1, cellBorders=False)
            oed.prepare()
            oed.write(imagePath)
            self.__imageCount += 1
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def status(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Status at %s (up %.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def alignMoleculePair(self, refIdentifier, refIdentifierType, fitIdentifier, fitIdentifierType, imagePath=None, **kwargs):
        """Create aligned depiction for a target molecule InChI, SMILES descriptors or PDB identifier.
        """
        try:
            imagePath = imagePath if imagePath else self.__makeImagePath()
            oeio = OeIoUtils()
            ccsw = ChemCompSearchWrapper()
            oesmP = ccsw.getSearchMoleculeProvider()
            # ---
            if refIdentifierType.lower() in ["smiles"]:
                oeMolRef = oeio.smilesToMol(refIdentifier)
            elif refIdentifierType.lower() in ["inchi"]:
                oeMolRef = oeio.inchiToMol(refIdentifier)
            elif refIdentifierType.lower() in ["identifierpdb"]:
                oeMolRef = oesmP.getMol(refIdentifier)
            #
            if fitIdentifierType.lower() in ["smiles"]:
                oeMolFit = oeio.smilesToMol(fitIdentifier)
            elif fitIdentifierType.lower() in ["inchi"]:
                oeMolFit = oeio.inchiToMol(fitIdentifier)
            elif fitIdentifierType.lower() in ["identifierpdb"]:
                oeMolFit = oesmP.getMol(fitIdentifier)
            # ---
            logger.info("oeMolRef atoms %r", oeMolRef.NumAtoms())
            logger.info("oeMolFit atoms %r", oeMolFit.NumAtoms())

            displayIdRef = "Ref"
            displayIdFit = "Fit"
            ok = self.__depictAlignedPair(oeMolRef, displayIdRef, oeMolFit, displayIdFit, imagePath, **kwargs)
            return imagePath if ok else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def __depictAlignedPair(self, oeMolRef, displayIdRef, oeMolFit, displayIdFit, imagePath, **kwargs):
        """ Depict pairwise MCSS alignment
        """
        try:
            #
            oed = OeDepictMCSAlignPage()
            oed.setSearchType(sType="relaxed")
            #
            oed.setRefMol(oeMolRef, displayIdRef)
            oed.setFitMol(oeMolFit, displayIdFit)
            #
            # imagePath = self.__makeImagePath()
            # ---
            bondDisplayWidth = 10.0
            numAtomsRef = oeMolRef.NumAtoms()
            if numAtomsRef > 100 and numAtomsRef <= 200:
                bondDisplayWidth = 6.0
            elif numAtomsRef > 200:
                bondDisplayWidth = 4.0
            # ---
            oed.setDisplayOptions(
                imageSizeX=kwargs.get("imageSizeX", 2500),
                imageSizeY=kwargs.get("imageSizeX", 2500),
                labelAtomName=kwargs.get("labelAtomName", False),
                labelAtomCIPStereo=kwargs.get("labelAtomCIPStereo", True),
                labelAtomIndex=kwargs.get("labelAtomIndex", False),
                labelBondIndex=kwargs.get("labelBondIndex", False),
                labelBondCIPStereo=kwargs.get("labelBondCIPStereo", True),
                cellBorders=kwargs.get("cellBorders", True),
                bondDisplayWidth=bondDisplayWidth,
                highlightStyleFit=kwargs.get("highlightStyleFit", "ballAndStickInverse"),
            )
            #
            aML = oed.alignPair(imagePath=imagePath)
            logger.info("Aligned atom count %d", len(aML))
            #
            # self.assertGreater(len(aML), 1)
            # if aML:
            #    for (rCC, rAt, tCC, tAt) in aML:
            #        logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False
