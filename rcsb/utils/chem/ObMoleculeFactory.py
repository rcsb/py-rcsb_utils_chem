##
# File:    ObMoleculeFactory.py
# Author:  jdw
# Date:    11-Jan-2020
# Version: 0.001
#
# Updates:
#
##
# flake8: noqa
# pylint: skip-file
"""
Classes to build OpenBabel molecule objects from chemical component definition data.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
from collections import defaultdict, namedtuple

from openbabel import openbabel
from openbabel import pybel


logger = logging.getLogger(__name__)

ComponentDetails = namedtuple("ComponentDetails", "ccId formula ifCharge")
ComponentAtom = namedtuple("ComponentAtom", "name aType isAromatic isChiral CIP fCharge")
ComponentBond = namedtuple("ComponentBond", "iType isAromatic CIP")
ComponentDescriptors = namedtuple("ComponentDescriptors", "smiles isoSmiles inchi inchiKey")


class ObMoleculeFactory(object):
    """Utility methods for constructing RdKit molecules from chemical component definition objects."""

    def __init__(self, verbose=False):
        self.__verbose = verbose
        self.__ccId = None
        self.__pybelMol = None
        self.__atomIdxD = None
        #
        self.__obConv = openbabel.OBConversion()
        self.__inputFormatDict = dict([f.split(" -- ") for f in self.__obConv.GetSupportedInputFormat()])
        self.__outputFormatDict = dict([f.split(" -- ") for f in self.__obConv.GetSupportedOutputFormat()])

    def setFile(self, ccId, filePath, molFormat="mol", atomIdxD=None):
        try:
            if molFormat not in self.__inputFormatDict:
                return False
            self.__ccId = ccId
            self.__atomIdxD = atomIdxD
            self.__pybelMol = next(pybel.readfile(molFormat, filePath), None)
            logger.info("pybelMol type %r", self.__pybelMol)
            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", ccId, str(e))
        return False

    def setString(self, ccId, molText, molFormat="mol", atomIdxD=None):
        try:
            if molFormat not in self.__inputFormatDict:
                return False
            self.__ccId = ccId
            self.__atomIdxD = atomIdxD
            self.__pybelMol = pybel.readstring(molFormat, molText)
            logger.info("obMol type %r", self.__pybelMol)
            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", ccId, str(e))
        return False

    def getMoleculeFeatures(self):
        """Get the essential features of the constructed obMol for the input component.

        OBConversion object, use SetInAndOutFormat(InCode, OutCode).
        To set a Read Option s, use SetOptions("s", OBConversion::INOPTIONS).

        """
        title = self.__pybelMol.title
        molWeight = self.__pybelMol.molwt
        formula = self.__pybelMol.formula
        ccId = title
        ifCharge = self.__pybelMol.charge
        logger.info("%s formula %s charge %d mw %f", title, formula, ifCharge, molWeight)
        inchi = self.__pybelMol.write("inchi").strip()
        inchiKey = self.__pybelMol.write("inchikey").strip()
        smiles = self.__pybelMol.write("can", opt={"n": None}).strip()
        isoSmiles = self.__pybelMol.write("can", opt={"i": None, "n": None}).strip()
        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        #
        #
        typeCounts = defaultdict(int)
        ccAtomD = {}
        ccAtomIdD = {}
        for ii, pat in enumerate(self.__pybelMol.atoms, 1):
            at = pat.OBAtom
            atIdx = at.GetIdx()
            # atNo = at.GetAtomicNum()
            aType = at.GetType()
            typeCounts[aType] += 1
            atName = self.__atomIdxD[ii] if ii in self.__atomIdxD else aType + str(typeCounts[aType])
            #
            isAromatic = at.IsAromatic()
            isChiral = at.IsChiral()
            iCharge = at.GetFormalCharge()
            cipStereo = None
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
            ccAtomIdD[atIdx] = atName
            logger.debug("%s Atom %s %s %r %r %s", ccId, atName, aType, isAromatic, isChiral, cipStereo)
        #
        ccBondD = {}
        for bnd in openbabel.OBMolBondIter(self.__pybelMol.OBMol):
            atI = bnd.GetBeginAtomIdx()
            atJ = bnd.GetEndAtomIdx()
            atNameI = ccAtomIdD[atI]
            atNameJ = ccAtomIdD[atJ]
            isAromatic = bnd.IsAromatic()
            iType = bnd.GetBondOrder()
            cipStereo = None
            logger.debug("Bond %s %s iType %r cipStereo %r aromatic %r", atNameI, atNameJ, iType, cipStereo, isAromatic)
            #
            ccBondD[(atNameI, atNameJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD
