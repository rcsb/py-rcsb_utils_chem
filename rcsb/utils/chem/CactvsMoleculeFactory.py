##
# File:    CactvsMoleculeFactory.py
# Author:  jdw
# Date:    25-Jan-2020
# Version: 0.001
#
# Updates:
#
##
"""
Classes to build CACTVS molecule objects data from chemical componenet data.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os

from collections import defaultdict, namedtuple
from pkg_resources import resource_filename, Requirement

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)

ComponentDetails = namedtuple("ComponentDetails", "ccId formula ifCharge")
ComponentAtom = namedtuple("ComponentAtom", "name aType isAromatic isChiral CIP fCharge")
ComponentBond = namedtuple("ComponentBond", "iType isAromatic CIP")
ComponentDescriptors = namedtuple("ComponentDescriptors", "smiles isoSmiles inchi inchiKey")


class CactvsMoleculeFactory(object):
    """ Utility methods for constructing Cactvs molecules from chemical component definition objects.
    """

    def __init__(self, cactvsPythonInterpreterPath, workPath, verbose=False):
        self.__cactvsPythonInterpreterPath = cactvsPythonInterpreterPath
        self.__verbose = verbose
        self.__ccId = None
        self.__atomIdxD = None
        self.__molFilePath = None
        self.__workPath = workPath if workPath else "."
        #

    def annotate(self, jsonPath, aroModel="cactvs"):
        ok = self.__runCactvsPython(self.__molFilePath, jsonPath, aroModel=aroModel)
        return ok

    def __runCactvsPython(self, molFilePath, jsonPath, aroModel="cactvs"):
        fp = resource_filename(Requirement.parse("rcsb.utils.chem"), "rcsb/utils/chem/cactvsAnnotateMol.py")
        logger.info("script path is %r", fp)
        exU = ExecUtils()
        ok = exU.run(self.__cactvsPythonInterpreterPath, [fp, molFilePath, jsonPath, aroModel])
        return ok

    def setFile(self, ccId, filePath, molFormat="mol", atomIdxD=None):
        try:
            _ = molFormat
            self.__ccId = ccId
            self.__atomIdxD = atomIdxD
            self.__molFilePath = filePath

            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", ccId, str(e))
        return False

    def getMoleculeFeatures(self, aroModel="cactvs"):
        """Get the essential features of the constructed obMol for the input component.

         OBConversion object, use SetInAndOutFormat(InCode, OutCode).
         To set a Read Option s, use SetOptions("s", OBConversion::INOPTIONS).

        """
        jsonFilePath = os.path.join(self.__workPath, self.__ccId + "-cactvs.json")
        self.annotate(jsonFilePath, aroModel=aroModel)
        #
        mU = MarshalUtil()
        ctvsMol = mU.doImport(jsonFilePath, fmt="json")
        #
        detailsD = ctvsMol["details"]
        title = detailsD["ccId"]
        molWeight = detailsD["formulaWeight"]
        formula = detailsD["formula"]
        ccId = title
        ifCharge = detailsD["ifCharge"]
        logger.info("%s formula %s charge %d mw %f", title, formula, ifCharge, molWeight)
        #
        desD = ctvsMol["descriptors"]
        inchi = None
        inchiKey = None
        smiles = desD["canSmi"]
        isoSmiles = desD["canIsoSmi"]
        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        #
        #
        aDL = ctvsMol["atoms"]
        typeCounts = defaultdict(int)
        ccAtomD = {}
        ccAtomIdD = {}
        for ii, aD in enumerate(aDL, 1):
            atIdx = aD["idx"]
            aType = aD["type"]
            typeCounts[aType] += 1
            atName = self.__atomIdxD[ii] if ii in self.__atomIdxD else aType + str(typeCounts[aType])
            #
            isAromatic = aD["isAromatic"]
            isChiral = aD["isChiral"]
            iCharge = aD["formalCharge"]
            cipStereo = aD["CIP"]
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
            ccAtomIdD[atIdx] = atName
            logger.debug("%s Atom %s %s %r %r %s", ccId, atName, aType, isAromatic, isChiral, cipStereo)
        #
        ccBondD = {}
        bDL = ctvsMol["bonds"]
        for bD in bDL:
            atI = bD["atI"]
            atJ = bD["atJ"]
            atNameI = ccAtomIdD[atI]
            atNameJ = ccAtomIdD[atJ]
            isAromatic = bD["isAromatic"]
            iType = bD["value_order"]
            cipStereo = bD["CIP"]
            logger.debug("Bond %s %s iType %r cipStereo %r aromatic %r", atNameI, atNameJ, iType, cipStereo, isAromatic)
            #
            ccBondD[(atNameI, atNameJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD
