##
# File:    MoleculeAnnotationsCompare.py
# Author:  jdw
# Date:    31-Dec-2019
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to compare molecular feature annotations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
from collections import defaultdict, namedtuple

from rcsb.utils.chem.PdbxChemComp import PdbxChemCompAtomIt, PdbxChemCompBondIt, PdbxChemCompDescriptorIt, PdbxChemCompIt

# from rcsb.utils.chem.PdbxChemCompConstants import PdbxChemCompConstants


logger = logging.getLogger(__name__)

ComponentDetails = namedtuple("ComponentDetails", "ccId formula ifCharge")
ComponentAtom = namedtuple("ComponentAtom", "name aType isAromatic isChiral CIP fCharge")
ComponentBond = namedtuple("ComponentBond", "iType isAromatic CIP")
ComponentDescriptors = namedtuple("ComponentDescriptors", "smiles isoSmiles inchi inchiKey")
ComponentCompare = namedtuple(
    "ComponentCompare",
    "ccId refType refInfo tstType tstInfo difDetails difAromaticAtoms difStereoAtoms difTypeBonds difAromaticBonds difSmiles difIsoSmiles difInchi difInchiKey",
)


class MoleculeAnnotationsCompare(object):
    """Utilities to compare molecular feature annotations."""

    def __init__(self, **kwargs):
        _ = kwargs

    def getChemCompFeatures(self, dataContainer, descriptorProgram="OPENEYE", filterHydrogens=False):
        """Get the essential features of the input chemical component definition."""
        ccIt = iter(PdbxChemCompIt(dataContainer))
        cc = next(ccIt, None)
        formula = cc.getFormulaWithCharge()
        ccId = cc.getId()
        ccName = cc.getName()
        ifCharge = cc.getFormalChargeAsInt()
        isAmbiguous = cc.getAmbiguousFlag() in ["Y", "y"]
        isCurrent = cc.getReleaseStatus() in ["REL"]
        #
        desIt = PdbxChemCompDescriptorIt(dataContainer)
        isoSmiles = smiles = inchi = inchiKey = None
        for des in desIt:
            desType = des.getType().upper()
            desProg = des.getProgram().upper()
            desText = des.getDescriptor().strip()
            if "OPEN" in desProg and desType == "SMILES_CANONICAL" and descriptorProgram == "OPENEYE":
                isoSmiles = desText
            elif "OPEN" in desProg and desType == "SMILES" and descriptorProgram == "OPENEYE":
                smiles = desText
            elif "CAC" in desProg and desType == "SMILES_CANONICAL" and descriptorProgram == "CACTVS":
                isoSmiles = desText
            elif "CAC" in desProg and desType == "SMILES" and descriptorProgram == "CACTVS":
                smiles = desText
            elif desType == "INCHI":
                inchi = desText
            elif desType == "INCHIKEY":
                inchiKey = desText
            logger.debug("type %r prog %r text %r", desType, desProg, desText)
        logger.debug("smiles %r isosmiles %r inchi %r inchikey %r", smiles, isoSmiles, inchi, inchiKey)
        #

        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        #
        atIt = PdbxChemCompAtomIt(dataContainer)
        typeCounts = defaultdict(int)
        ccNameTypeD = {}
        ccAtomD = {}
        sumCharge = 0
        for at in atIt:
            atName = at.getName()
            aType = at.getType().upper()
            ccNameTypeD[atName] = aType
            if filterHydrogens and aType == "H":
                continue
            #
            ccNameTypeD[atName] = aType
            typeCounts[aType] += 1
            isAromatic = at.isAromatic()
            isChiral = at.isChiral()
            cipStereo = at.getCIPStereo()
            #
            iCharge = at.getFormalChargeAsInt()
            sumCharge += iCharge
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
        #
        if sumCharge != ifCharge and ifCharge == 0:
            ifCharge = sumCharge
            formula = cc.getFormulaWithCharge(ifCharge=ifCharge)

        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        ccBondD = {}
        bndIt = PdbxChemCompBondIt(dataContainer)
        for bnd in bndIt:
            atIdI, atIdJ = bnd.getBond()
            if filterHydrogens and (ccNameTypeD[atIdI] == "H" or ccNameTypeD[atIdJ] == "H"):
                continue
            cipStereo = bnd.getStereo()
            isAromatic = bnd.isAromatic()
            iType = bnd.getIntegerType()
            ccBondD[(atIdI, atIdJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"name": ccName, "isCurrent": isCurrent, "isAmbiguous": isAmbiguous, "details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD

    def compare(self, refFD, tstFD, refType="CC", refInfo=None, tstType="CACTVS", tstInfo=None):
        """Compare reference and test molecular features in the input feature dictionaries.

        Args:
            refFD (dict): reference molecule feature dictionary
            tstFD (dict): test molecule feature dictionary
            refType (str, optional): reference molecule type. Defaults to "CC".
            tstType (str, optional): test molecule type. Defaults to "CACTVS".
            refInfo (str, optional): ref molecule additional details. Defaults to None.
            tstInfo (str, optional): test molecule additional details. Defaults to None.


        Returns:
            bool, tuple: flag for no differences, tuple describing comparison details (ComponentCompare)
        """
        compareEachAtom = compareEachBond = False
        difDetails = difAromaticAtoms = difStereoAtoms = difSmiles = difIsoSmiles = difInchi = difInchiKey = False
        difAromaticBonds = difTypeBonds = False
        okR = ok = True
        #
        refDetails = refFD["details"]
        ccId = refDetails.ccId
        #
        isAmbiguous = refFD["isAmbiguous"] if "isAmbiguous" in refFD else False
        isCurrent = refFD["isCurrent"] if "isCurrent" in refFD else False
        refName = refFD["name"] if "name" in refFD else ccId
        #

        #
        tstDetails = tstFD["details"]
        if refDetails != tstDetails:
            logger.error("%s: details differ", ccId)
            logger.error("%s REF: %r", ccId, refDetails)
            logger.error("%s TST: %r", ccId, tstDetails)
            ok = False
        #
        refAtomD = refFD["atoms"]
        refTypeCounts = defaultdict(int)
        refAromaticCount = 0
        refChiralCount = 0
        for atName, atTup in refAtomD.items():
            refTypeCounts[atTup.aType] += 1
            if atTup.isChiral:
                refChiralCount += 1
            if atTup.isAromatic:
                refAromaticCount += 1
        #
        tstAtomD = tstFD["atoms"]
        tstTypeCounts = defaultdict(int)

        tstAromaticCount = 0
        tstChiralCount = 0
        for atName, atTup in tstAtomD.items():
            tstTypeCounts[atTup.aType] += 1
            if atTup.isChiral:
                tstChiralCount += 1
            if atTup.isAromatic:
                tstAromaticCount += 1
        #
        if len(refTypeCounts) != len(tstTypeCounts):
            logger.info("ref %r tst %r", refTypeCounts, tstTypeCounts)
            logger.error("%s: atom types length mismatch", ccId)
            ok = False
        #
        if refChiralCount != tstChiralCount:
            logger.error("%s: chiral atom count missmatch ref: %d tst: %d", ccId, refChiralCount, tstChiralCount)
            ok = False
            difStereoAtoms = True
        #
        if refAromaticCount != tstAromaticCount:
            logger.error("%s: aromatic atom count missmatch ref: %d tst: %d", ccId, refAromaticCount, tstAromaticCount)
            ok = False
            difAromaticAtoms = True
        #
        okR = okR and ok

        for aType in refTypeCounts:
            try:
                ok = refTypeCounts[aType] == tstTypeCounts[aType]
            except Exception:
                ok = False
                okR = okR and ok
            if not ok:
                logger.error("%s:  atom type counts differ for %r", ccId, aType)
        #
        # Atom by atom -
        #
        if compareEachAtom:
            for atName in refAtomD:
                try:
                    ok = refAtomD[atName] == tstAtomD[atName]
                except Exception:
                    ok = False
                    okR = okR and ok
                if not ok:
                    if atName in tstAtomD:
                        logger.error("%s: atom features differ %r: \n -- REF: %r \n -- TST: %r", ccId, atName, refAtomD[atName], tstAtomD[atName])
                    else:
                        logger.error("%s: atom features differ for %r missing atom in tstMol", ccId, atName)
        # ----
        #
        refBondD = refFD["bonds"]
        refBondTypeCounts = defaultdict(int)
        refBondAromaticCount = 0
        for _, bTup in refBondD.items():
            refBondTypeCounts[bTup.iType] += 1
            if bTup.isAromatic:
                refBondAromaticCount += 1
        #
        tstBondD = tstFD["bonds"]
        tstBondTypeCounts = defaultdict(int)
        tstBondAromaticCount = 0
        for _, bTup in tstBondD.items():
            tstBondTypeCounts[bTup.iType] += 1
            if bTup.isAromatic:
                tstBondAromaticCount += 1
        #
        for bType in refBondTypeCounts:
            if bType not in tstBondTypeCounts or refBondTypeCounts[bType] != tstBondTypeCounts[bType]:
                difTypeBonds = True
        #
        if difTypeBonds or len(refBondTypeCounts) != len(tstBondTypeCounts):
            logger.error("%s: bond type mismatch", ccId)
            ok = False
            difTypeBonds = True
        #
        if refBondAromaticCount != tstBondAromaticCount:
            logger.error("%s: aromatic bond count missmatch ref: %d tst: %d", ccId, refBondAromaticCount, tstBondAromaticCount)
            ok = False
            difAromaticBonds = True
        #
        okR = okR and ok
        #
        # Bond by Bond-
        #
        if compareEachBond:
            for ky in refBondD:
                (atNameI, atNameJ) = ky
                try:
                    if (atNameI, atNameJ) in tstBondD:
                        ok = refBondD[(atNameI, atNameJ)] == tstBondD[(atNameI, atNameJ)]
                    elif (atNameJ, atNameI) in tstBondD:
                        ok = refBondD[(atNameI, atNameJ)] == tstBondD[(atNameJ, atNameI)]
                    else:
                        ok = False
                except Exception:
                    ok = False
                #
                okR = okR and ok
                if not ok:
                    if (atNameI, atNameJ) in tstBondD:
                        logger.error(
                            "%s: bond features differ (%r, %r): \n -- REF: %r \n -- TST: %r", ccId, atNameI, atNameJ, refBondD[(atNameI, atNameJ)], tstBondD[(atNameI, atNameJ)]
                        )
                    elif (atNameJ, atNameI) in tstBondD:
                        logger.error(
                            "%s: bond features differ (%r, %r): \n -- REF: %r \n -- TST: %r", ccId, atNameI, atNameJ, refBondD[(atNameI, atNameJ)], tstBondD[(atNameJ, atNameI)]
                        )
                    else:
                        logger.error("%s: bond features differ for (%r, %r) missing bond in tstMol", ccId, atNameI, atNameJ)

        if not okR:
            if not isCurrent or isAmbiguous:
                logger.error("%s: comparison failing - ambiguous flag %r current flag %r name %r", ccId, isAmbiguous, isCurrent, refName)

        #
        # Now do the descriptors
        #
        refDesD = refFD["descriptors"]
        tstDesD = tstFD["descriptors"]
        #
        if refDesD.smiles != tstDesD.smiles:
            logger.error("%s ISOSMILES differ \n -- REF: %r \n -- TST: %r", ccId, refDesD.smiles, tstDesD.smiles)
            difSmiles = True
        if refDesD.isoSmiles != tstDesD.isoSmiles:
            logger.error("%s ISOSMILES differ \n -- REF: %r \n -- TST: %r", ccId, refDesD.isoSmiles, tstDesD.isoSmiles)
            difIsoSmiles = True
        if refDesD.inchi != tstDesD.inchi and tstDesD.inchi:
            logger.error("%s InChIs differ", ccId)
            difInchi = True
        if refDesD.inchiKey != tstDesD.inchiKey and tstDesD.inchiKey:
            logger.error("%s InChI Keys differ", ccId)
            difInchiKey = True
        #
        retCmp = ComponentCompare(
            ccId=ccId,
            refType=refType,
            refInfo=refInfo,
            tstType=tstType,
            tstInfo=tstInfo,
            difDetails=difDetails,
            difAromaticAtoms=difAromaticAtoms,
            difStereoAtoms=difStereoAtoms,
            difTypeBonds=difTypeBonds,
            difAromaticBonds=difAromaticBonds,
            difSmiles=difSmiles,
            difIsoSmiles=difIsoSmiles,
            difInchi=difInchi,
            difInchiKey=difInchiKey,
        )
        return okR, retCmp
