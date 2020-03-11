##
# File:    RdMoleculeFactory.py
# Author:  jdw
# Date:    1-Jan-2020
# Version: 0.001
#
# Updates:
#
##
# flake8: noqa
# pylint: skip-file
"""
Classes to build RdKit molecule objects from chemical component definition data.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
from collections import defaultdict, namedtuple

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

logger = logging.getLogger(__name__)

ComponentDetails = namedtuple("ComponentDetails", "ccId formula ifCharge")
ComponentAtom = namedtuple("ComponentAtom", "name aType isAromatic isChiral CIP fCharge")
ComponentBond = namedtuple("ComponentBond", "iType isAromatic CIP")
ComponentDescriptors = namedtuple("ComponentDescriptors", "smiles isoSmiles inchi inchiKey")


class RdMoleculeFactory(object):
    """ Utility methods for constructing RdKit molecules from chemical component definition objects.
    """

    def __init__(self, verbose=False):
        self.__verbose = verbose
        self.__ccId = None
        self.__rdMol = None
        self.__atomIdxD = None

    def setFile(self, ccId, filePath, atomIdxD=None, molFormat="sdf"):
        try:
            self.__ccId = ccId
            self.__atomIdxD = atomIdxD if atomIdxD else {}
            if molFormat in ["sdf", "mol"]:
                self.__rdMol = Chem.MolFromMolFile(filePath)
                return True
            else:
                logger.error("Unsupported format %r", molFormat)
        except Exception as e:
            logger.exception("Failing for %s with %s", ccId, str(e))
        return False

    def updateProperties(self):
        if not self.__rdMol:
            return False
        try:
            self.__rdMol.UpdatePropertyCache()
            for atom in self.__rdMol.GetAtoms():
                atom.UpdatePropertyCache()
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def updateStereoChemAssignments(self):
        if not self.__rdMol:
            return False
        try:
            Chem.AssignAtomChiralTagsFromStructure(self.__rdMol)
            Chem.AssignStereochemistry(self.__rdMol, True, True, True)
            # self.__rdMol.Debug()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def updateAromaticAssignments(self, aromaticModel="default"):
        mD = {
            "simple": Chem.AromaticityModel.AROMATICITY_SIMPLE,
            "default": Chem.AromaticityModel.AROMATICITY_DEFAULT,
            "mdl": Chem.AromaticityModel.AROMATICITY_MDL,
            "rdkit": Chem.AromaticityModel.AROMATICITY_RDKIT,
        }
        if not self.__rdMol or aromaticModel not in mD:
            return False
        try:
            #  self.__rdMol.Debug()
            Chem.SetAromaticity(self.__rdMol, mD[aromaticModel])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return True

    def getMoleculeFeatures(self):
        """Get the essential features of the constructed rdMol for the input component.
        """
        mD = self.__rdMol.GetPropsAsDict()
        logger.debug("mol props %r", mD.items())
        #
        formula = CalcMolFormula(self.__rdMol)
        ccId = self.__ccId
        ifCharge = Chem.rdmolops.GetFormalCharge(self.__rdMol)
        #
        inchiKey = Chem.inchi.MolToInchiKey(self.__rdMol)
        inchi = Chem.inchi.MolToInchi(self.__rdMol)
        smiles = Chem.rdmolfiles.MolToSmiles(self.__rdMol, isomericSmiles=False, canonical=True)
        isoSmiles = Chem.rdmolfiles.MolToSmiles(self.__rdMol, isomericSmiles=True, canonical=True)
        logger.debug("%s formula %s", ccId, formula)
        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        #
        typeCounts = defaultdict(int)
        ccAtomD = {}
        ccAtomIdD = {}
        for ii, at in enumerate(self.__rdMol.GetAtoms(), 1):
            atIdx = at.GetIdx()
            aType = at.GetSymbol()
            typeCounts[aType] += 1
            atName = self.__atomIdxD[ii] if ii in self.__atomIdxD else aType + str(typeCounts[aType])
            # atNo = at.GetAtomicNum()
            isAromatic = at.GetIsAromatic()
            isChiral = at.GetChiralTag() > 0
            iCharge = at.GetFormalCharge()
            # cipStereo = at.GetProp("_CIPCode")
            atD = at.GetPropsAsDict(includePrivate=True, includeComputed=True)
            cipStereo = None
            if "_CIPCode" in atD:
                cipStereo = atD["_CIPCode"]
            if cipStereo and cipStereo not in ["S", "R"]:
                logger.error("%s (%s): Unexpected atom CIP stereo setting %r", ccId, atName, cipStereo)
            #
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
            ccAtomIdD[atIdx] = atName
            # nL = at.GetProp(includePrivate=True, includeComputed=True)
            atD = at.GetPropsAsDict(includePrivate=True, includeComputed=True)
            logger.debug("%s Atom %s %s %r %r %s", ccId, atName, aType, isAromatic, isChiral, cipStereo)
        #
        ccBondD = {}
        for bnd in self.__rdMol.GetBonds():
            atI = bnd.GetBeginAtomIdx()
            atJ = bnd.GetEndAtomIdx()
            atNameI = ccAtomIdD[atI]
            atNameJ = ccAtomIdD[atJ]
            isAromatic = bnd.GetIsAromatic()
            #
            # bType = bnd.GetBondType()
            # iType = 0
            cipStereo = None
            tS = bnd.GetStereo()
            if tS == Chem.rdchem.BondStereo.STEREOE:
                cipStereo = "E"
            elif tS == Chem.rdchem.BondStereo.STEREOZ:
                cipStereo = "Z"

            # bL = bnd.GetPropNames(includePrivate=True, includeComputed=True)
            bD = bnd.GetPropsAsDict(includePrivate=True, includeComputed=True)
            iType = bD["_MolFileBondType"]
            logger.debug("Bond %s %s iType %r cipStereo %r aromatic %r", atNameI, atNameJ, iType, cipStereo, isAromatic)
            #
            if cipStereo and cipStereo not in ["E", "Z"]:
                logger.error("%s (%s %s): Unexpected bond CIP stereo setting %r", ccId, atNameI, atNameJ, cipStereo)
            #
            ccBondD[(atNameI, atNameJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD
