##
# File: PdbxChemComp.py
# Date: 2-Oct-2019  John Westbrook
#
# Update:
#   2-Oct-2019 jdw adapted from PdbxChemCompPersist()
#   7-Nov-2019 jdw add alternate charge and formula methods
##
"""
A collection of access and iterator classes supporting chemical component dictionary data.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging

from rcsb.utils.chem.PdbxChemCompConstants import PdbxChemCompConstants

logger = logging.getLogger(__name__)


class PdbxCategoryItBase(object):
    """  Base category iterator class.
    """

    def __init__(self, dataCategory, func):
        self.__rL = dataCategory.getRowList() if dataCategory else []
        self.__func = func

    def get(self, index=0):
        try:
            return self.__rL[index]
        except Exception:
            return []

    def __iter__(self):
        return self.forward()

    def __len__(self):
        return len(self.__rL)

    def forward(self):
        # Forward generator
        currentRow = 0
        while currentRow < len(self.__rL):
            row = self.__rL[currentRow]
            currentRow += 1
            yield self.__func(row)

    def reverse(self):
        # The reverse generator
        currentRow = len(self.__rL)
        while currentRow > 0:
            currentRow -= 1
            yield self.__func(self.__rL[currentRow])


class PdbxChemCompIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("chem_comp")
        obj = PdbxChemCompPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompAtomIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("chem_comp_atom")
        obj = PdbxChemCompAtomPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompAtomIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompBondIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("chem_comp_bond")
        obj = PdbxChemCompBondPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompBondIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompDescriptorIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("pdbx_chem_comp_descriptor")
        obj = PdbxChemCompDescriptorPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompDescriptorIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompIdentifierIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("pdbx_chem_comp_identifier")
        obj = PdbxChemCompIdentifierPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompIdentifierIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompAuditIt(PdbxCategoryItBase):
    def __init__(self, dataContainer):
        dataCategory = dataContainer.getObj("pdbx_chem_comp_audit")
        obj = PdbxChemCompAuditPersist([], attributeNameList=dataCategory.getAttributeList() if dataCategory else [])
        super(PdbxChemCompAuditIt, self).__init__(dataCategory, obj.set)


class PdbxChemCompPersist(object):
    """ Accessor methods chemical component attributes.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def getId(self):
        return self.__getAttribute("id")

    def getName(self):
        return self.__getAttribute("name")

    def getType(self):
        return self.__getAttribute("type")

    def getPdbxType(self):
        return self.__getAttribute("pdbx_type")

    def getFormula(self):
        return self.__getAttribute("formula")

    def getSynonyms(self):
        return self.__getAttribute("pdbx_synonyms")

    def getFormalCharge(self):
        return self.__getAttribute("pdbx_formal_charge")

    def getFormalChargeAsInt(self):
        tch = self.getFormalCharge()
        iCharge = int(tch) if tch and tch not in [".", "?"] else 0
        return iCharge

    def getFormulaWithCharge(self, ifCharge=None):
        formula = (self.getFormula()).replace(" ", "")
        fcharge = ifCharge if ifCharge else self.getFormalChargeAsInt()
        if fcharge:
            sign = "+" if fcharge > 0 else "-"
            mag = str(abs(fcharge)) if abs(fcharge) > 1 else ""
            formula = formula + sign + mag
        return formula

    def getModificationDate(self):
        return self.__getAttribute("pdbx_modified_date")

    def getInitialDate(self):
        return self.__getAttribute("pdbx_initial_date")

    def getReleaseStatus(self):
        return self.__getAttribute("pdbx_release_status")

    def getFormulaWeight(self):
        return self.__getAttribute("formula_weight")

    def getSubComponentList(self):
        return self.__getAttribute("pdbx_subcomponent_list")

    def getAmbiguousFlag(self):
        return self.__getAttribute("pdbx_ambiguous_flag")

    def getProcessingSite(self):
        return self.__getAttribute("pdbx_processing_site")

    def getReplacesId(self):
        return self.__getAttribute("pdbx_replaces")

    def getReplacesById(self):
        return self.__getAttribute("pdbx_replaced_by")

    def getNstdParentId(self):
        return self.__getAttribute("mon_nstd_parent_comp_id")

    def getOneLetterCode(self):
        return self.__getAttribute("one_letter_code")

    def getThreeLetterCode(self):
        return self.__getAttribute("three_letter_code")

    def getModelCoordinatesPdbCode(self):
        return self.__getAttribute("pdbx_model_coordinates_db_code")

    def getMissingModelCoordinates(self):
        return self.__getAttribute("pdbx_model_coordinates_missing_flag")

    def getMissingIdealCoordinates(self):
        return self.__getAttribute("pdbx_ideal_coordinates_missing_flag")


class PdbxChemCompAtomPersist(object):
    """ Accessor methods chemical component atom attributes.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def getName(self):
        return self.__getAttribute("atom_id")

    def isChiral(self):
        return self.__getAttribute("pdbx_stereo_config") != "N"

    def getType(self):
        return self.__getAttribute("type_symbol")

    def getLeavingAtomFlag(self):
        return self.__getAttribute("pdbx_leaving_atom_flag")

    def getAtNo(self):
        try:
            tyU = str(self.getType()).upper()
            if (tyU == "D") or (tyU == "T"):
                tyU = "H"
            return PdbxChemCompConstants.periodicTable.index(tyU) + 1
        except Exception:
            # traceback.print_exc(file=self.__lfh)
            return 0

    def getIsotope(self):
        ty = self.getType()
        if ty == "D":
            return 2
        elif ty == "T":
            return 3
        else:
            return 0

    def isAromatic(self):
        return self.__getAttribute("pdbx_aromatic_flag") != "N"

    def getCIPStereo(self):
        tS = self.__getAttribute("pdbx_stereo_config")
        tS = tS if tS not in [".", "?", "N"] else None
        if tS and tS not in ["R", "S"]:
            logger.error("Unexpected CIP stereo type %r in %r", tS, self.__rowData)
        return tS

    def getFormalCharge(self):
        try:
            return int(self.__getAttribute("charge"))
        except Exception:
            return 0

    def getFormalChargeAsInt(self):
        try:
            tch = self.getFormalCharge()
            iCharge = int(tch) if tch and tch not in [".", "?"] else 0
            return iCharge
        except Exception:
            return 0

    def hasModelCoordinates(self):
        xV, yV, zV = self.getModelCoordinates()
        # x=self.__getAttribute('model_Cartn_x')
        # y=self.__getAttribute('model_Cartn_y')
        # z=self.__getAttribute('model_Cartn_z')
        #
        return (xV is not None) and (yV is not None) and (zV is not None)

    def hasIdealCoordinates(self):
        xV, yV, zV = self.getIdealCoordinates()
        # x=self.__getAttribute('pdbx_model_Cartn_x_ideal')
        # y=self.__getAttribute('pdbx_model_Cartn_y_ideal')
        # z=self.__getAttribute('pdbx_model_Cartn_z_ideal')
        #
        return (xV is not None) and (yV is not None) and (zV is not None)

    def getModelCoordinates(self):
        """ Returns (x,y,z)
        """
        try:
            xV = float(self.__getAttribute("model_Cartn_x"))
            yV = float(self.__getAttribute("model_Cartn_y"))
            zV = float(self.__getAttribute("model_Cartn_z"))
            return (xV, yV, zV)
        except Exception:
            return (None, None, None)

    def getIdealCoordinates(self):
        """ Returns (x,y,z)
        """
        try:
            xV = float(self.__getAttribute("pdbx_model_Cartn_x_ideal"))
            yV = float(self.__getAttribute("pdbx_model_Cartn_y_ideal"))
            zV = float(self.__getAttribute("pdbx_model_Cartn_z_ideal"))
            return (xV, yV, zV)
        except Exception:
            return (None, None, None)

    def dump(self):
        logger.info(" --rowData: %r", self.__rowData)


class PdbxChemCompBondPersist(object):
    """ Accessor methods chemical component bond attributes.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def getBond(self):
        """ Returns (atomI,atomJ) atom ids from the atom list.
        """
        return (self.__getAttribute("atom_id_1"), self.__getAttribute("atom_id_2"))

    def getType(self):
        tS = self.__getAttribute("value_order")
        tS = tS if tS not in [".", "?"] else None
        if tS and tS not in ["SING", "DOUB", "TRIP", "QUAD"]:
            logger.error("Unexpected bond type %r in %r", tS, self.__rowData)
        return tS

    def getIntegerType(self):
        bT = self.__getAttribute("value_order")
        if bT == "SING":
            return 1
        elif bT == "DOUB":
            return 2
        elif bT == "TRIP":
            return 3
        elif bT == "QUAD":
            return 4
        else:
            return 0

    def isAromatic(self):
        return self.__getAttribute("pdbx_aromatic_flag") == "Y"

    def getStereo(self):
        tS = self.__getAttribute("pdbx_stereo_config")
        tS = tS if tS not in [".", "?", "N"] else None
        if tS and tS not in ["E", "Z"]:
            logger.error("Unexpected CIP stereo type %r in %r", tS, self.__rowData)
        return tS

    def hasStereo(self):
        return self.__getAttribute("pdbx_stereo_config") != "N"

    def dump(self):
        logger.info(" -- rowData: %r", self.__rowData)


class PdbxChemCompIdentifierPersist(object):
    """ Accessor methods chemical component identifier attributes.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def getIdentifier(self):
        """ Returns the value of the identifier.
        """
        return self.__getAttribute("identifier")

    def getType(self):
        return self.__getAttribute("type")

    def getProgram(self):
        return self.__getAttribute("program")

    def getProgramVersion(self):
        return self.__getAttribute("program_version")

    def dump(self):
        logger.info("-- rowData: %r", self.__rowData)


class PdbxChemCompDescriptorPersist(object):
    """ Accessor methods chemical component descriptor  attributes.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList
        self.__buildTypeD = {
            ("SMILES_CANONICAL", "OPENEYE OETOOLKITS"): "oe-iso-smiles",
            ("SMILES", "OPENEYE OETOOLKITS"): "oe-smiles",
            ("SMILES", "ACDLABS"): "acdlabs-smiles",
            ("SMILES", "CACTVS"): "cactvs-smiles",
            ("SMILES_CANONICAL", "CACTVS"): "cactvs-iso-smiles",
            ("INCHI", "INCHI"): "inchi",
            ("INCHIKEY", "INCHI"): "inchikey",
        }

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def getDescriptor(self):
        """ Returns the value of the descriptor.
        """
        return self.__getAttribute("descriptor")

    def getType(self):
        return self.__getAttribute("type")

    def getProgram(self):
        return self.__getAttribute("program")

    def getProgramVersion(self):
        return self.__getAttribute("program_version")

    def getMolBuildType(self):
        bT = None
        try:
            dType = self.getType()
            dProg = self.getProgram()
            bT = self.__buildTypeD[(dType.upper(), dProg.upper())]
        except Exception as e:
            logger.error("Unknown descriptor build type %r %r %r", dType, dProg, str(e))
        return bT

    def dump(self):
        logger.info(" -- rowData: %r", self.__rowData)


class PdbxChemCompAuditPersist(object):
    """ Accessor methods chemical component audit details.

    """

    def __init__(self, rowData, attributeNameList):
        self.__rowData = rowData
        self.__attributeNameList = attributeNameList

    def __getAttribute(self, name):
        try:
            i = self.__attributeNameList.index(name)
            return self.__rowData[i]
        except Exception:
            return None

    def set(self, rowData=None):
        self.__rowData = rowData
        return self

    def getActionType(self):
        """ Returns the value of the action type.
        """
        return self.__getAttribute("action_type")

    def getDate(self):
        """ Returns the value of audit date.
        """
        return self.__getAttribute("date")

    def getProcessingSite(self):
        """ Returns the value of processing site.
        """
        return self.__getAttribute("processing_site")

    def getAnnotator(self):
        """ Returns the value of audit annotator.
        """
        return self.__getAttribute("annotator")

    def getDetails(self):
        """ Returns the value of audit details.
        """
        return self.__getAttribute("details")
