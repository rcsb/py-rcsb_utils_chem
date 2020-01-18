##
# File:    OeMoleculeFactory.py
# Author:  jdw
# Date:    2-Oct-2019
# Version: 0.001
#
# Updates:
# 2-Oct-2019  jdw adapted from OeBuildMol()
##
"""
Classes to build OE molecule objects from chemical component definition data.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"


import logging
from collections import defaultdict, namedtuple

from openeye import oechem

from rcsb.utils.chem.PdbxChemComp import PdbxChemCompAtomIt, PdbxChemCompBondIt, PdbxChemCompDescriptorIt, PdbxChemCompIt
from rcsb.utils.chem.PdbxChemCompConstants import PdbxChemCompConstants


logger = logging.getLogger(__name__)

ComponentDetails = namedtuple("ComponentDetails", "ccId formula ifCharge")
ComponentAtom = namedtuple("ComponentAtom", "name aType isAromatic isChiral CIP fCharge")
ComponentBond = namedtuple("ComponentBond", "iType isAromatic CIP")
ComponentDescriptors = namedtuple("ComponentDescriptors", "smiles isoSmiles inchi inchiKey")


class OeMoleculeFactory(object):
    """ Utility methods for constructing OEGraphMols from chemical component definition objects.
    """

    def __init__(self, verbose=False):
        self.__verbose = verbose
        self.__ccId = None
        self.__oeMol = None
        #
        # dictionary of element counts eD[atno]=count
        self.__eD = {}
        #
        # Source data categories objects from chemical component definitions.
        self.__dataContainer = None
        #
        self.__molXyzL = []
        #

    def __clear(self):
        #
        self.__eD = {}
        self.__molXyzL = []
        self.__oeMol = None

    def setQuiet(self):
        # quiet
        oechem.OEThrow.SetLevel(5)

    def setDebug(self, flag):
        self.__verbose = flag

    def set(self, dataContainer):
        try:
            if dataContainer.exists("chem_comp"):
                ccIt = PdbxChemCompIt(dataContainer)
                for cc in ccIt:
                    ccId = cc.getId()
                self.__ccId = ccId
            else:
                self.__ccId = dataContainer.getName()

            self.__dataContainer = dataContainer
            if self.__dataContainer.exists("chem_comp_atom"):
                return self.__ccId
            else:
                return None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def setOeMol(self, inpOeMol, ccId):
        """  Load this object with an existing oeMOL()
        """
        self.__clear()
        self.__oeMol = oechem.OEGraphMol(inpOeMol)
        self.__ccId = ccId
        # self.getElementCounts()

    def build3D(self, coordType="model", setTitle=True):
        try:
            self.__build3D(coordType=coordType, setTitle=setTitle)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return False

    def __build3D(self, coordType="model", setTitle=True):
        """ Build OE molecule using 3D coordinates and OE stereo perception.
        """
        self.__clear()
        # self.__oeMol=OEGraphMol()
        self.__oeMol = oechem.OEMol()
        #
        if setTitle:
            self.__oeMol.SetTitle(self.__ccId)
        aL = []

        # Atom index dictionary
        aD = {}
        i = 1

        atomIt = PdbxChemCompAtomIt(self.__dataContainer)
        for ccAt in atomIt:

            atName = ccAt.getName()
            aD[atName] = i
            i += 1
            atNo = ccAt.getAtNo()
            if atNo not in self.__eD:
                self.__eD[atNo] = 1
            else:
                self.__eD[atNo] += 1

            atType = ccAt.getType()
            fc = ccAt.getFormalCharge()
            chFlag = ccAt.isChiral()
            arFlag = ccAt.isAromatic()
            isotope = ccAt.getIsotope()
            leavingAtom = ccAt.getLeavingAtomFlag()

            oeAt = self.__oeMol.NewAtom(atNo)
            oeAt.SetName(atName)
            oeAt.SetFormalCharge(fc)
            oeAt.SetStringData("pdbx_leaving_atom_flag", leavingAtom)
            oeAt.SetChiral(chFlag)
            oeAt.SetIsotope(isotope)
            oeAt.SetAromatic(arFlag)
            tChk = oeAt.IsAromatic()
            if tChk != arFlag:
                logger.info("%s inconsistent aromatic setting on %s", self.__ccId, atName)

            cTup = None
            if (coordType == "model") and ccAt.hasModelCoordinates():
                cTup = ccAt.getModelCoordinates()
                self.__oeMol.SetCoords(oeAt, cTup)
            elif (coordType == "ideal") and ccAt.hasIdealCoordinates():
                cTup = ccAt.getIdealCoordinates()
                self.__oeMol.SetCoords(oeAt, cTup)
            else:
                pass

            logger.debug("Atom - %s type %s atno %d isotope %d fc %d (xyz) %r", atName, atType, atNo, isotope, fc, cTup)
            aL.append(oeAt)

        bondIt = PdbxChemCompBondIt(self.__dataContainer)
        for ccBnd in bondIt:
            (at1, at2) = ccBnd.getBond()
            iat1 = aD[at1] - 1
            iat2 = aD[at2] - 1
            iType = ccBnd.getIntegerType()
            logger.debug(" %s %d -- %s %d (%d)", at1, iat1, at2, iat2, iType)
            oeBnd = self.__oeMol.NewBond(aL[iat1], aL[iat2], iType)
            atName1 = oeBnd.GetBgn().GetName()
            atName2 = oeBnd.GetEnd().GetName()
            logger.debug("%s - created bond %s %s", self.__ccId, atName1, atName2)

        #
        # run standard OE perceptions --
        #
        self.updateOePerceptions3D(self.__oeMol)
        # Reset aromatic flags
        self.__transferAromaticFlagsOE()

    def updateOePerceptions3D(self, oeMol, aromaticModel=oechem.OEAroModelOpenEye):
        oeMol.SetDimension(3)
        oechem.OEFindRingAtomsAndBonds(oeMol)
        oechem.OEPerceiveChiral(oeMol)
        #
        # Other aromatic models: oechem.OEAroModelMDL or OEAroModelDaylight or oechem.OEAroModelOpenEye
        #
        if aromaticModel is not None:
            oechem.OEAssignAromaticFlags(oeMol, aromaticModel)
        oechem.OE3DToInternalStereo(oeMol)
        return oeMol

    def updateOePerceptions2D(self, oeMol, aromaticModel=oechem.OEAroModelOpenEye):
        #
        # run standard perceptions --
        oechem.OEFindRingAtomsAndBonds(oeMol)
        oechem.OEPerceiveChiral(oeMol)
        #
        # Other aromatic models: oechem.OEAroModelMDL or OEAroModelDaylight or oechem.OEAroModelOpenEye
        #
        if aromaticModel is not None:
            oechem.OEAssignAromaticFlags(oeMol, aromaticModel)
        return oeMol

    def __transferAromaticFlagsOE(self):
        """Manually assign/transfer aromatic flags from chemical component definition.

           Must follow other perceptions -
        """
        atIt = PdbxChemCompAtomIt(self.__dataContainer)
        ccAtomD = {}
        for at in atIt:
            atName = at.getName()
            isAromatic = at.isAromatic()
            ccAtomD[atName] = isAromatic
        ccBondD = {}
        bndIt = PdbxChemCompBondIt(self.__dataContainer)
        for bnd in bndIt:
            atIdI, atIdJ = bnd.getBond()
            isAromatic = bnd.isAromatic()
            ccBondD[(atIdI, atIdJ)] = isAromatic
            ccBondD[(atIdJ, atIdI)] = isAromatic
        #
        for atom in self.__oeMol.GetAtoms():
            atName = atom.GetName()
            isAromatic = ccAtomD[atName]
            atom.SetAromatic(isAromatic)

        for bond in self.__oeMol.GetBonds():
            atNameI = bond.GetBgn().GetName()
            atNameJ = bond.GetEnd().GetName()
            isAromatic = ccBondD[(atNameI, atNameJ)]
            bond.SetAromatic(isAromatic)

    def updateCIPStereoOE(self, oeMol):
        """ Manually assign OE CIP stereo perceptions for each atom and bond.

            Redundant of oechem.OE3DToInternalStereo(oeMol) -
        """
        for atom in oeMol.GetAtoms():
            cipV = oechem.OEPerceiveCIPStereo(oeMol, atom)
            if cipV in [oechem.OECIPAtomStereo_S, oechem.OECIPAtomStereo_R]:
                logger.info("%s CIP stereo is %r", atom.GetName(), cipV)
                oechem.OESetCIPStereo(oeMol, atom, cipV)
                cipV2 = oechem.OEPerceiveCIPStereo(oeMol, atom)
                logger.info("%s perceive after set CIP stereo (%r) is %r", atom.GetName(), atom.HasStereoSpecified(), cipV2)

        #
        for bond in oeMol.GetBonds():
            if bond.GetOrder() == 2:
                cipV = oechem.OEPerceiveCIPStereo(oeMol, bond)
                if cipV in [oechem.OECIPBondStereo_Z, oechem.OECIPBondStereo_E]:
                    logger.info("%s-%s CIP stereo is %r", bond.GetBgn().GetName(), bond.GetEnd().GetName(), cipV)
                    oechem.OESetCIPStereo(oeMol, bond, cipV)
                    cipV2 = oechem.OEPerceiveCIPStereo(oeMol, bond)
                    logger.info("%s-%s perceive after set CIP stereo (%r) is %r", bond.GetBgn().GetName(), bond.GetEnd().GetName(), bond.HasStereoSpecified(), cipV2)

    def build2D(self, setTitle=True):
        try:
            self.__build2D(setTitle=setTitle)
            return True
        except Exception:
            return False

    def __build2D(self, setTitle=True):
        """  Build molecule using existing assignments of chemical information in the CC definition.
        """
        self.__clear()
        self.__oeMol = oechem.OEGraphMol()
        if setTitle:
            self.__oeMol.SetTitle(self.__ccId)
        aL = []
        i = 1
        # Atom index dictionary
        aD = {}

        atomIt = PdbxChemCompAtomIt(self.__dataContainer)
        for ccAt in atomIt:
            atName = ccAt.getName()
            aD[atName] = i
            i += 1
            atNo = ccAt.getAtNo()
            if atNo not in self.__eD:
                self.__eD[atNo] = 1
            else:
                self.__eD[atNo] += 1
            atType = ccAt.getType()
            fc = ccAt.getFormalCharge()
            chFlag = ccAt.isChiral()
            arFlag = ccAt.isAromatic()
            isotope = ccAt.getIsotope()
            leavingAtom = ccAt.getLeavingAtomFlag()

            oeAt = self.__oeMol.NewAtom(atNo)
            oeAt.SetName(atName)
            oeAt.SetFormalCharge(fc)
            oeAt.SetStringData("pdbx_leaving_atom_flag", leavingAtom)
            oeAt.SetChiral(chFlag)
            oeAt.SetIsotope(isotope)
            oeAt.SetAromatic(arFlag)
            if chFlag:
                st = ccAt.getCIPStereo()
                if st == "S" or st == "R":
                    oeAt.SetStringData("StereoInfo", st)
            logger.debug("Atom - %s type %s atno %d isotope %d fc %d chFlag %r", atName, atType, atNo, isotope, fc, chFlag)
            aL.append(oeAt)

        bondIt = PdbxChemCompBondIt(self.__dataContainer)
        for ccBnd in bondIt:
            (at1, at2) = ccBnd.getBond()
            iat1 = aD[at1] - 1
            iat2 = aD[at2] - 1
            iType = ccBnd.getIntegerType()
            arFlag = ccBnd.isAromatic()
            logger.debug(" %s %d -- %s %d (%d)", at1, iat1, at2, iat2, iType)

            oeBnd = self.__oeMol.NewBond(aL[iat1], aL[iat2], iType)
            oeBnd.SetAromatic(arFlag)
            if arFlag:
                oeBnd.SetIntType(5)
            st = ccBnd.getStereo()
            if st == "E" or st == "Z":
                oeBnd.SetStringData("StereoInfo", st)

        #
        self.__oeMol = self.updateOePerceptions2D(self.__oeMol, aromaticModel=None)
        # run standard perceptions --
        # oechem.OEFindRingAtomsAndBonds(self.__oeMol)
        # oechem.OEPerceiveChiral(self.__oeMol)

        for oeAt in self.__oeMol.GetAtoms():
            st = oeAt.GetStringData("StereoInfo")
            if st == "R":
                oechem.OESetCIPStereo(self.__oeMol, oeAt, oechem.OECIPAtomStereo_R)
            elif st == "S":
                oechem.OESetCIPStereo(self.__oeMol, oeAt, oechem.OECIPAtomStereo_S)

        for oeBnd in self.__oeMol.GetBonds():
            st = oeBnd.GetStringData("StereoInfo")
            if st == "E":
                oechem.OESetCIPStereo(self.__oeMol, oeBnd, oechem.OECIPBondStereo_E)
            elif st == "Z":
                oechem.OESetCIPStereo(self.__oeMol, oeBnd, oechem.OECIPBondStereo_Z)
        if self.__verbose:
            for ii, atm in enumerate(self.__oeMol.GetAtoms()):
                logger.debug("OeBuildMol.build2d - atom  %d %s", ii, atm.GetName())
        # Add new
        self.__transferAromaticFlagsOE()

    def getGraphMolSuppressH(self):
        """ Return the current constructed OE molecule with hydrogens suppressed.
        """
        # OESuppressHydrogens(self.__oeMol, retainPolar=False,retainStereo=True,retainIsotope=True)
        oechem.OESuppressHydrogens(self.__oeMol)
        return self.__oeMol

    def getMol(self):
        """ Return the current constructed OE molecule.
        """
        return oechem.OEMol(self.__oeMol)

    def getCanSMILES(self):
        """ Return the cannonical SMILES string derived from the current OD molecule.
        """
        return oechem.OECreateCanSmiString(self.__oeMol)

    def getIsoSMILES(self):
        """ Return the cannonical stereo SMILES string derived from the current OE molecule.
        """
        return oechem.OECreateIsoSmiString(self.__oeMol)

    def getFormula(self):
        """ Return the Hill order formulat  derived from the current OE molecule.
        """
        return oechem.OEMolecularFormula(self.__oeMol)

    def getInChIKey(self):
        """ Return the InChI key derived from the current OE molecule.
        """
        return oechem.OECreateInChIKey(self.__oeMol)

    def getInChI(self):
        """ Return the InChI string derived from the current OE molecule.
        """
        return oechem.OECreateInChI(self.__oeMol)

    def getTitle(self):
        """ Return the title assigned to the current OE molecule
        """
        return self.__oeMol.GetTitle()

    def getCcId(self):
        """ Return the CC id of this object -
        """
        return self.__ccId

    def getCoords(self):
        """  Return coordinate list if a 3D molecule is built -- otherwise an empty list --

        """
        return self.__molXyzL

    def setSimpleAtomNames(self):
        """
        """
        for atom in self.__oeMol.GetAtoms():
            atom.SetIntType(atom.GetAtomicNum())
            atom.SetType(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()))
        oechem.OETriposAtomNames(self.__oeMol)

    def getElementCounts(self):
        """ Get the dictionary of element counts (eg. eD[iAtNo]=iCount).
        """
        if not self.__eD:
            # calculate from current oeMol
            try:
                self.__eD = {}
                for atom in self.__oeMol.GetAtoms():
                    atNo = atom.GetAtomicNum()
                    if atNo not in self.__eD:
                        self.__eD[atNo] = 1
                    else:
                        self.__eD[atNo] += 1
            except Exception:
                pass

        return self.__eD

    def compare(self):
        """Compare the source and constructed molecule features.
        """
        okR = ok = True
        ccId = self.__ccId
        ccFD = self.__getChemCompFeatures(self.__dataContainer)
        oeFD = self.__getOeMoleculeFeatures()
        isAmbiguous = ccFD["isAmbiguous"]
        isCurrent = ccFD["isCurrent"]
        ccName = ccFD["name"]
        #
        ccDetails = ccFD["details"]
        oeDetails = oeFD["details"]
        if ccDetails != oeDetails:
            logger.error("%s: details differ", ccId)
            logger.error("%s CC: %r", ccId, ccDetails)
            logger.error("%s OE: %r", ccId, oeDetails)
            ok = False
        #
        ccAtomD = ccFD["atoms"]
        ccTypeCounts = defaultdict(int)
        ccAromaticCount = 0
        ccChiralCount = 0
        for atName, atTup in ccAtomD.items():
            ccTypeCounts[atTup.aType] += 1
            if atTup.isChiral:
                ccChiralCount += 1
            if atTup.isAromatic:
                ccAromaticCount += 1
        #
        oeAtomD = oeFD["atoms"]
        oeTypeCounts = defaultdict(int)

        oeAromaticCount = 0
        oeChiralCount = 0
        for atName, atTup in oeAtomD.items():
            oeTypeCounts[atTup.aType] += 1
            if atTup.isChiral:
                oeChiralCount += 1
            if atTup.isAromatic:
                oeAromaticCount += 1
        #
        if len(ccTypeCounts) != len(oeTypeCounts):
            logger.error("%s: atom types length mismatch", ccId)
            ok = False
        #
        if ccChiralCount != oeChiralCount:
            logger.error("%s: chiral atom count missmatch cc: %d oe: %d", ccId, ccChiralCount, oeChiralCount)
            ok = False
        #
        if ccAromaticCount != oeAromaticCount:
            logger.error("%s: aromatic atom count missmatch cc: %d oe: %d", ccId, ccAromaticCount, oeAromaticCount)
            ok = False
        #
        okR = okR and ok

        for aType in ccTypeCounts:
            try:
                ok = ccTypeCounts[aType] == oeTypeCounts[aType]
            except Exception:
                ok = False
                okR = okR and ok
            if not ok:
                logger.error("%s:  atom type counts differ for %r", ccId, aType)
        #
        # Atom by atom -
        #
        for atName in ccAtomD:
            try:
                ok = ccAtomD[atName] == oeAtomD[atName]
            except Exception:
                ok = False
                okR = okR and ok
            if not ok:
                if atName in oeAtomD:
                    logger.error("%s: atom features differ %r: \n -- CC: %r \n -- OE: %r", ccId, atName, ccAtomD[atName], oeAtomD[atName])
                else:
                    logger.error("%s: atom features differ for %r missing atom in oeMol", ccId, atName)
        #
        ccBondD = ccFD["bonds"]
        ccTypeCounts = defaultdict(int)
        ccAromaticCount = 0
        for _, bTup in ccBondD.items():
            ccTypeCounts[bTup.iType] += 1
            if bTup.isAromatic:
                ccAromaticCount += 1
        #
        oeBondD = oeFD["bonds"]
        oeTypeCounts = defaultdict(int)
        oeAromaticCount = 0
        for _, bTup in oeBondD.items():
            oeTypeCounts[bTup.iType] += 1
            if bTup.isAromatic:
                oeAromaticCount += 1
        #
        #
        if len(ccTypeCounts) != len(oeTypeCounts):
            logger.error("%s: bond types length mismatch", ccId)
            ok = False
        #
        if ccAromaticCount != oeAromaticCount:
            logger.error("%s: aromatic bond count missmatch cc: %d oe: %d", ccId, ccAromaticCount, oeAromaticCount)
            ok = False
        #
        okR = okR and ok
        #
        #
        # Bond by Bond-
        #
        for ky in ccBondD:
            (atNameI, atNameJ) = ky
            try:
                if (atNameI, atNameJ) in oeBondD:
                    ok = ccBondD[(atNameI, atNameJ)] == oeBondD[(atNameI, atNameJ)]
                elif (atNameJ, atNameI) in oeBondD:
                    ok = ccBondD[(atNameI, atNameJ)] == oeBondD[(atNameJ, atNameI)]
                else:
                    ok = False
            except Exception:
                ok = False
            #
            okR = okR and ok
            if not ok:
                if (atNameI, atNameJ) in oeBondD:
                    logger.error("%s: bond features differ (%r, %r): \n -- CC: %r \n -- OE: %r", ccId, atNameI, atNameJ, ccBondD[(atNameI, atNameJ)], oeBondD[(atNameI, atNameJ)])
                elif (atNameJ, atNameI) in oeBondD:
                    logger.error("%s: bond features differ (%r, %r): \n -- CC: %r \n -- OE: %r", ccId, atNameI, atNameJ, ccBondD[(atNameI, atNameJ)], oeBondD[(atNameJ, atNameI)])
                else:
                    logger.error("%s: bond features differ for (%r, %r) missing bond in oeMol", ccId, atNameI, atNameJ)

        if not okR:
            if not isCurrent or isAmbiguous:
                logger.error("%s: comparison failing - ambiguous flag %r current flag %r name %r", ccId, isAmbiguous, isCurrent, ccName)

        #
        # Now do the descriptors
        #
        ccDesD = ccFD["descriptors"]
        oeDesD = oeFD["descriptors"]
        #
        if ccDesD.smiles != oeDesD.smiles:
            logger.error("%s SMILES differ", ccId)
        if ccDesD.isoSmiles != oeDesD.isoSmiles:
            logger.error("%s ISOSMILES differ \n -- CC: %r \n -- OE: %r", ccId, ccDesD.isoSmiles, oeDesD.isoSmiles)

        if ccDesD.inchi != oeDesD.inchi:
            logger.error("%s InChIs differ", ccId)
        if ccDesD.inchiKey != oeDesD.inchiKey:
            logger.error("%s InChI Keys differ", ccId)
        #
        return okR

    def __getChemCompFeatures(self, dataContainer):
        """Get the essential features of the input component.
        """
        ccIt = PdbxChemCompIt(dataContainer)
        for cc in ccIt:
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
            if "OPEN" in desProg and desType == "SMILES_CANONICAL":
                isoSmiles = desText
            elif "OPEN" in desProg and desType == "SMILES":
                smiles = desText
            elif desType == "INCHI":
                inchi = desText
            elif desType == "INCHIKEY":
                inchiKey = desText
            logger.debug("type %r prog %r text %r", desType, desProg, desText)
        logger.debug("smiles %r isosmiles %r inchi %r inchikey %r", smiles, isoSmiles, inchi, inchiKey)
        #
        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        #
        atIt = PdbxChemCompAtomIt(dataContainer)
        typeCounts = defaultdict(int)

        ccAtomD = {}
        for at in atIt:
            atName = at.getName()
            aType = at.getType().upper()
            typeCounts[aType] += 1
            isAromatic = at.isAromatic()
            isChiral = at.isChiral()
            cipStereo = at.getCIPStereo()
            #
            iCharge = at.getFormalChargeAsInt()
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
        #
        ccBondD = {}
        bndIt = PdbxChemCompBondIt(dataContainer)
        for bnd in bndIt:
            atIdI, atIdJ = bnd.getBond()
            cipStereo = bnd.getStereo()
            isAromatic = bnd.isAromatic()
            iType = bnd.getIntegerType()
            ccBondD[(atIdI, atIdJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"name": ccName, "isCurrent": isCurrent, "isAmbiguous": isAmbiguous, "details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD

    def __getOeMoleculeFeatures(self):
        """Get the essential features of the constructed OEMol for the input component.
        """
        formula = self.getFormula()
        # ccId = self.__ccId
        ccId = self.getTitle()
        ifCharge = oechem.OENetCharge(self.__oeMol)
        inchi = self.getInChI()
        inchiKey = self.getInChIKey()
        smiles = self.getCanSMILES()
        isoSmiles = self.getIsoSMILES()
        details = ComponentDetails(ccId=ccId, formula=formula, ifCharge=ifCharge)
        descriptors = ComponentDescriptors(smiles=smiles, isoSmiles=isoSmiles, inchi=inchi, inchiKey=inchiKey)
        typeCounts = defaultdict(int)
        ccAtomD = {}
        ccAtomIdD = {}
        for atId, at in enumerate(self.__oeMol.GetAtoms(), 1):
            atName = at.GetName()
            atNo = at.GetAtomicNum()
            aType = PdbxChemCompConstants.periodicTable[atNo - 1]
            typeCounts[aType] += 1
            isAromatic = at.IsAromatic()
            isChiral = at.IsChiral()
            iCharge = at.GetFormalCharge()
            #
            cipStereo = None
            if at.HasStereoSpecified():
                cip = oechem.OEPerceiveCIPStereo(self.__oeMol, at)
                if cip == oechem.OECIPAtomStereo_S:
                    cipStereo = "S"
                if cip == oechem.OECIPAtomStereo_R:
                    cipStereo = "R"
                if cip == oechem.OECIPAtomStereo_NotStereo:
                    cipStereo = "?"
                if isChiral and cip == oechem.OECIPAtomStereo_UnspecStereo:
                    cipStereo = "?"

                #
                if cipStereo and cipStereo not in ["S", "R"]:
                    logger.error("%s (%s): Unexpected atom CIP stereo setting %r", ccId, atName, cipStereo)
            #
            ccAtomD[atName] = ComponentAtom(name=atName, aType=aType, isAromatic=isAromatic, isChiral=isChiral, CIP=cipStereo, fCharge=iCharge)
            ccAtomIdD[atId] = atName
        #
        ccBondD = {}
        for bnd in self.__oeMol.GetBonds():
            atNameI = bnd.GetBgn().GetName()
            atNameJ = bnd.GetEnd().GetName()
            isAromatic = bnd.IsAromatic()
            iType = bnd.GetOrder()
            cipStereo = None
            if iType == 2:
                cip = oechem.OEPerceiveCIPStereo(self.__oeMol, bnd)
                if bnd.HasStereoSpecified():
                    if cip == oechem.OECIPBondStereo_E:
                        cipStereo = "E"
                    if cip == oechem.OECIPBondStereo_Z:
                        cipStereo = "Z"
                    if cip == oechem.OECIPBondStereo_NotStereo:
                        cipStereo = "?"
                    if cip == oechem.OECIPBondStereo_UnspecStereo:
                        cipStereo = "?"

            if cipStereo and cipStereo not in ["E", "Z"]:
                logger.error("%s (%s %s): Unexpected bond CIP stereo setting %r", ccId, atNameI, atNameJ, cipStereo)
            #
            ccBondD[(atNameI, atNameJ)] = ComponentBond(iType=iType, isAromatic=isAromatic, CIP=cipStereo)
        #
        ccD = {"details": details, "descriptors": descriptors, "atoms": ccAtomD, "bonds": ccBondD}
        return ccD
