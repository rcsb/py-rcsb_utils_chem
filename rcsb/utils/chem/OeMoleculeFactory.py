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

from openeye import oechem
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompAtomIt
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompBondIt
from rcsb.utils.chem.PdbxChemComp import PdbxChemCompIt
from rcsb.utils.chem.PdbxChemCompConstants import PdbxChemCompConstants

logger = logging.getLogger(__name__)


class OeMoleculeFactory(object):
    """ Utility methods for constructing OEGraphMols from chemical component definition objects.
    """

    def __init__(self, verbose=False):
        self.__verbose = verbose
        #
        # File system path to the chemical component dictionary definitions in (CVS checkout organization)
        #
        # Internal storage for current OE molecule
        self.__oeMol = None
        #
        # Component identifier
        #
        self.__ccId = None
        #
        # dictionary of element counts eD[atno]=count
        self.__eD = {}
        #
        # Source data categories objects from chemical component definitions.
        self.__dataContainer = None
        #
        self.__molXyzL = []

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

    def __clear(self):
        self.__oeMol = None
        self.__eD = {}

    def serialize(self):
        """ Create a string representing the content of the current OE molecule.   This
            serialization uses the OE internal binary format.
        """
        oms = oechem.oemolostream()
        oms.SetFormat(oechem.OEFormat_OEB)
        oms.openstring()
        oechem.OEWriteMolecule(oms, self.__oeMol)
        logger.debug("SMILES %s", oechem.OECreateCanSmiString(self.__oeMol))
        logger.debug("Atoms = %d", self.__oeMol.NumAtoms())
        return oms.GetString()

    def deserialize(self, oeS):
        """ Reconstruct an OE molecule from the input string serialization (OE binary).

            The deserialized molecule is used to initialize the internal OE molecule
            within this object.

            Returns True for success or False otherwise.
        """
        self.__clear()
        ims = oechem.oemolistream()
        ims.SetFormat(oechem.OEFormat_OEB)
        ims.openstring(oeS)

        nmol = 0
        mList = []
        # for mol in ims.GetOEGraphMols():
        for mol in ims.GetOEMols():
            logger.debug("SMILES %s", oechem.OECreateCanSmiString(mol))
            logger.debug("title  %s", mol.GetTitle())
            logger.debug("atoms  %d", mol.NumAtoms())
            # mList.append(OEGraphMol(mol))
            mList.append(oechem.OEMol(mol))
            nmol += 1
        #
        if nmol >= 1:
            self.__oeMol = mList[0]
            self.__ccId = self.__oeMol.GetTitle()
            #
            logger.debug("mols  %d", nmol)
            logger.debug("id %s", self.__ccId)
            logger.debug("atoms  %d", self.__oeMol.NumAtoms())
            return True
        else:
            return False

    def simpleAtomNames(self):
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

            # if chFlag:
            #    st=ccAt.getCIPStereo()
            #    if st == 'S' or st == 'R':
            #        oeAt.SetStringData("StereoInfo",st)
            # if (self.__verbose):
            #    logger.debug("Atom - %s type %s atno %d isotope %d fc %d chFlag %r\n" % (atName,atType,atNo,isotope,fc,chFlag))
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
            #
            # arFlag = ccBnd.isAromatic()
            logger.debug(" %s %d -- %s %d (%d)", at1, iat1, at2, iat2, iType)
            self.__oeMol.NewBond(aL[iat1], aL[iat2], iType)

            # oeBnd.SetAromatic(arFlag)
            # if arFlag:
            #    oeBnd.SetIntType(5)
            # st=ccBnd.getStereo()
            # if st == 'E' or st =='Z':
            #    oeBnd.SetStringData("StereoInfo",st)
        #
        # run standard perceptions --
        #
        self.updatePerceptions3D()

    def updatePerceptions3D(self):
        self.__oeMol.SetDimension(3)
        oechem.OEFindRingAtomsAndBonds(self.__oeMol)
        oechem.OEPerceiveChiral(self.__oeMol)
        oechem.OEAssignAromaticFlags(self.__oeMol, oechem.OEAroModelOpenEye)
        #
        oechem.OE3DToInternalStereo(self.__oeMol)
        # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
        self.updateCIPStereoOE()

    def updateCIPStereoOE(self):
        """ OE perception of CIP stereo -
        """
        for atom in self.__oeMol.GetAtoms():
            oechem.OEPerceiveCIPStereo(self.__oeMol, atom)

        for bond in self.__oeMol.GetBonds():
            if bond.GetOrder() == 2:
                oechem.OEPerceiveCIPStereo(self.__oeMol, bond)

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
        # run standard perceptions --
        oechem.OEFindRingAtomsAndBonds(self.__oeMol)
        # JDW Turn off chiral perception
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

    def importFile(self, filePath, dType="2D"):
        """  Contruct a OEGraphMol using the content of the input file.  The input
             file must have a file extension recognized by the OE toolkit (e.g. .sdf)
        """
        ifs = oechem.oemolistream()
        if not ifs.open(filePath):
            return False
        # JDW
        self.__oeMol = oechem.OEGraphMol()
        # self.__oeMol = oechem.OEMol()
        oechem.OEReadMolecule(ifs, self.__oeMol)
        #        OETriposAtomNames(self.__oeMol)
        if dType == "2D":
            # run standard perceptions --
            oechem.OEFindRingAtomsAndBonds(self.__oeMol)
            oechem.OEPerceiveChiral(self.__oeMol)

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
        elif dType == "3D":
            # run standard perceptions --
            #
            self.__oeMol.SetDimension(3)
            oechem.OE3DToInternalStereo(self.__oeMol)
            oechem.OEFindRingAtomsAndBonds(self.__oeMol)
            # Other aromatic models: OEAroModelMDL or OEAroModelDaylight
            oechem.OEAssignAromaticFlags(self.__oeMol, oechem.OEAroModelOpenEye)
            self.updateCIPStereoOE()
            oechem.OEAddExplicitHydrogens(self.__oeMol)

        self.__molXyzL = []
        aC = {}
        for ii, atm in enumerate(self.__oeMol.GetAtoms()):
            iAtNum = atm.GetAtomicNum()
            if iAtNum in aC:
                aC[iAtNum] += 1
            else:
                aC[iAtNum] = 1
            atName = PdbxChemCompConstants.periodicTable[iAtNum - 1] + str(aC[iAtNum])
            atm.SetName(atName)
            #
            xyzL = oechem.OEFloatArray(3)
            self.__oeMol.GetCoords(atm, xyzL)
            self.__molXyzL.append((ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), xyzL[0], xyzL[1], xyzL[2]))

        return True

    def importSmiles(self, smiles):
        """  Contruct a OEGraphMol using the input descriptor.
        """
        self.__oeMol = oechem.OEGraphMol()
        oechem.OEParseSmiles(self.__oeMol, smiles)
        oechem.OEFindRingAtomsAndBonds(self.__oeMol)
        oechem.OEPerceiveChiral(self.__oeMol)
        return True

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

    def write(self, filePath, constantMol=False):
        try:
            ofs = oechem.oemolostream()
            ofs.open(filePath)
            logger.info("Writing %s title %s\n", filePath, self.__oeMol.GetTitle())
            if constantMol:
                oechem.OEWriteConstMolecule(ofs, self.__oeMol)
            else:
                oechem.OEWriteMolecule(ofs, self.__oeMol)
            return True
        except Exception as e:
            logger.exception("Failing for %s with %s", filePath, str(e))
        return False
