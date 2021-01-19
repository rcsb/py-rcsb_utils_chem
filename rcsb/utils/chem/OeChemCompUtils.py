##
# File:    OeChemCompUtils.py
# Author:  jdw
# Date:    17-Jan-2021
# Version: 0.001
#
# Updates:
#
#
##
"""
Utilities to build chemical component definitions from OE molecule objects.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import logging
import time

from openeye import oechem
from openeye import oeiupac

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from rcsb.utils.io.MarshalUtil import MarshalUtil


logger = logging.getLogger(__name__)


class OeChemCompUtils(object):
    """Build chemical component definitions from OE molecule objects."""

    def __init__(self):

        self.__oeVersion = oechem.OEToolkitsGetRelease()
        self.__containerList = []

    def write(self, filePath):
        """Write the contents of the current data container of chemical component definitions.

        Args:
            filePath (str): output chemical component definition file path

        Returns:
            (bool): True for success or False otherwise
        """
        try:
            mU = MarshalUtil()
            ok = mU.doExport(filePath, self.__containerList, fmt="mmcif")
        except Exception as e:
            logger.exception("Failing %r with %s", filePath, str(e))
            ok = False
        return ok

    def addOeMol(self, ccId, oeMol, missingModelXyz=True, writeIdealXyz=False):
        """Add the input oeMol to the current PDBx data container as a chemical component definition.

        Args:
            oeMol (obj): instance of OE molecule
            ccId (str): chemical component identifer
            name (str, optional): chemical component name. Defaults to None.
            missingModelXyz (bool, optional): set the missing model coordinate flag. Defaults to True.
            writeIdealXyz (bool, optional): write 3D coordinates in using ideal coordinate data items. Defaults to False.

        Returns:
            (bool): True for success and False otherwise
        """
        try:
            ccIdU = str(ccId).strip().upper()
            curContainer = DataContainer(ccIdU)
            #
            rowD = self.__makeChemCompCategory(ccIdU, oeMol, site="RCSB", missingModelXyz=missingModelXyz)
            aCat = DataCategory("chem_comp", list(rowD.keys()), [rowD])
            curContainer.append(aCat)
            #
            rowDL = self.__makeChemCompAtomCategory(ccIdU, oeMol, writeIdealXyz=writeIdealXyz)
            aCat = DataCategory("chem_comp_atom", list(rowDL[0].keys()), rowDL)
            curContainer.append(aCat)
            #
            rowDL = self.__makeChemCompBondCategory(ccIdU, oeMol)
            aCat = DataCategory("chem_comp_bond", list(rowDL[0].keys()), rowDL)
            curContainer.append(aCat)
            #
            rowDL = self.__makeChemCompDescriptorCategory(ccIdU, oeMol)
            aCat = DataCategory("pdbx_chem_comp_descriptor", list(rowDL[0].keys()), rowDL)
            curContainer.append(aCat)
            #
            rowDL = self.__makeChemCompIdentifierCategory(ccIdU, oeMol)
            aCat = DataCategory("pdbx_chem_comp_identifier", list(rowDL[0].keys()), rowDL)
            curContainer.append(aCat)
            #
            rowD = self.__makeChemCompAuditRow(ccIdU)
            aCat = DataCategory("pdbx_chem_comp_audit", list(rowD.keys()), [rowD])
            curContainer.append(aCat)
            #
            self.__containerList.append(curContainer)
            return True
        except Exception as e:
            logger.exception("Failing %r with %s", ccId, str(e))
            #
        return False

    def __getFormalCharge(self, oeMol):
        fCharge = 0
        for atom in oeMol.GetAtoms():
            fCharge += atom.GetFormalCharge()
        return fCharge

    def __makeChemCompAuditRow(self, ccId, action="CREATE", date=None, processingSite="RCSB", annotator="?", details="?"):
        """
        loop_
        _pdbx_chem_comp_audit.comp_id
        _pdbx_chem_comp_audit.action_type
        _pdbx_chem_comp_audit.date
        _pdbx_chem_comp_audit.processing_site
        _pdbx_chem_comp_audit.annotator
        _pdbx_chem_comp_audit.details
        ARG "Create component"  1999-07-08 RCSB ? ?

        """
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["action_type"] = action
        if date is None:
            date = time.strftime("%Y-%m-%d", time.localtime())
        aRow["date"] = date
        aRow["processing_site"] = processingSite
        aRow["annotator"] = annotator
        aRow["details"] = details
        #
        return aRow

    def __makeChemCompDescriptorCategory(self, ccId, oeMol):
        """
            loop_
            _pdbx_chem_comp_descriptor.comp_id
            _pdbx_chem_comp_descriptor.type
            _pdbx_chem_comp_descriptor.program
            _pdbx_chem_comp_descriptor.program_version
            _pdbx_chem_comp_descriptor.descriptor
            ARG SMILES           ACDLabs              10.04 "O=C(O)C(N)CCCNC(=[NH2+])N"
            ARG SMILES_CANONICAL CACTVS               3.341 "N[C@@H](CCCNC(N)=[NH2+])C(O)=O"
            ARG SMILES           CACTVS               3.341 "N[CH](CCCNC(N)=[NH2+])C(O)=O"
            ARG SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "C(C[C@@H](C(=O)O)N)CNC(=[NH2+])N"
            ARG SMILES           "OpenEye OEToolkits" 1.5.0 "C(CC(C(=O)O)N)CNC(=[NH2+])N"
            ARG InChI            InChI                1.03  "InChI=1S/C6H14N4O2/c7-4(5(11)12)2-1-3-1..... "
            ARG InChIKey         InChI                1.03  ODKSFYDXXFIFQN-BYPYZUCNSA-O
        #
        """
        rowL = []
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "SMILES_CANONICAL"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        aRow["descriptor"] = oechem.OECreateIsoSmiString(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "SMILES"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        aRow["descriptor"] = oechem.OECreateCanSmiString(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "InChI"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        aRow["descriptor"] = oechem.OECreateInChI(oeMol)
        rowL.append(aRow)
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "InChIKey"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        aRow["descriptor"] = oechem.OECreateInChIKey(oeMol)
        rowL.append(aRow)
        #
        return rowL

    def __makeChemCompIdentifierCategory(self, ccId, oeMol):
        """

        loop_
        _pdbx_chem_comp_identifier.comp_id
        _pdbx_chem_comp_identifier.type
        _pdbx_chem_comp_identifier.program
        _pdbx_chem_comp_identifier.program_version
        _pdbx_chem_comp_identifier.identifier
        ATP "SYSTEMATIC NAME" ACDLabs              10.04
        ;adenosine 5'-(tetrahydrogen triphosphate)
        ;
        ATP "SYSTEMATIC NAME" "OpenEye OEToolkits" 1.5.0 "[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-..."
        #
        """
        rowL = []
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "SYSTEMATIC NAME"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        style = oeiupac.OEGetIUPACNamStyle("systematic")
        name = oeiupac.OEToUTF8(oeiupac.OECreateIUPACName(oeMol, style))
        aRow["identifier"] = name
        rowL.append(aRow)
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "COMMON"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        style = oeiupac.OEGetIUPACNamStyle("traditional")
        name = oeiupac.OEToUTF8(oeiupac.OECreateIUPACName(oeMol, style))
        aRow["identifier"] = name
        rowL.append(aRow)
        #
        aRow = {}
        aRow["comp_id"] = ccId
        aRow["type"] = "SYNONYM"
        aRow["program"] = "OpenEye OEToolkits"
        aRow["program_version"] = self.__oeVersion
        style = oeiupac.OEGetIUPACNamStyle("acdname")
        name = oeiupac.OEToUTF8(oeiupac.OECreateIUPACName(oeMol, style))
        aRow["identifier"] = name
        rowL.append(aRow)
        #
        return rowL

    def __makeChemCompCategory(self, ccId, oeMol, site="RCSB", missingModelXyz=False):
        #
        lt = time.strftime("%Y-%m-%d", time.localtime())
        formula = oechem.OEMolecularFormula(oeMol)
        charge = self.__getFormalCharge(oeMol)
        fW = oechem.OECalculateMolecularWeight(oeMol)
        style = oeiupac.OEGetIUPACNamStyle("systematic")
        name = oeiupac.OEToUTF8(oeiupac.OECreateIUPACName(oeMol, style))
        #
        ccRow = {}
        ccRow["id"] = ccId
        if name is not None:
            ccRow["name"] = name
        else:
            ccRow["name"] = "?"
        ccRow["type"] = "NON-POLYMER"
        ccRow["pdbx_type"] = "?"
        if formula is not None:
            ccRow["formula"] = formula
        else:
            ccRow["formula"] = "?"
        ccRow["mon_nstd_parent_comp_id"] = "?"
        # ccRow["pdbx_synonyms"] = "?"
        if charge is not None:
            ccRow["pdbx_formal_charge"] = charge
        else:
            ccRow["pdbx_formal_charge"] = "?"
        ccRow["pdbx_ambiguous_flag"] = "N"
        ccRow["pdbx_initial_date"] = lt
        ccRow["pdbx_modified_date"] = lt
        ccRow["pdbx_release_status"] = "HOLD"
        ccRow["pdbx_replaced_by"] = "?"
        ccRow["pdbx_replaces"] = "?"
        if fW is not None:
            ccRow["formula_weight"] = "%0.3f" % fW
        else:
            ccRow["formula_weight"] = "?"
        ccRow["one_letter_code"] = "?"
        tlc = ccId.split("_")[0]
        ccRow["three_letter_code"] = tlc
        ccRow["pdbx_model_coordinates_details"] = "?"
        ccRow["pdbx_ideal_coordinates_details"] = "?"
        if missingModelXyz:
            ccRow["pdbx_model_coordinates_missing_flag"] = "Y"
        else:
            ccRow["pdbx_model_coordinates_missing_flag"] = "N"
        ccRow["pdbx_model_coordinates_db_code"] = "?"
        ccRow["pdbx_processing_site"] = site
        ccRow["pdbx_subcomponent_list"] = "?"
        return ccRow

    def __makeChemCompAtomCategory(self, ccId, oeMol, writeIdealXyz=True):
        """Populate elements of chemical component definition for atoms and bonds."""
        #
        idCode = ccId
        #
        rowL = []
        for ii, atom in enumerate(oeMol.GetAtoms()):
            #
            atRow = {}
            atRow["comp_id"] = idCode
            atRow["atom_id"] = atom.GetName().strip()
            atRow["alt_atom_id"] = atom.GetName().strip()
            atRow["type_symbol"] = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())
            atRow["charge"] = atom.GetFormalCharge()
            (xC, yC, zC) = oeMol.GetCoords(atom)
            if writeIdealXyz:
                atRow["pdbx_model_Cartn_x_ideal"] = "%0.3f" % xC
                atRow["pdbx_model_Cartn_y_ideal"] = "%0.3f" % yC
                atRow["pdbx_model_Cartn_z_ideal"] = "%0.3f" % zC
            else:
                atRow["model_Cartn_x"] = "%0.3f" % xC
                atRow["model_Cartn_y"] = "%0.3f" % yC
                atRow["model_Cartn_z"] = "%0.3f" % zC
            if atom.GetStringData("pdbx_leaving_atom_flag") in ["Y", "N"]:
                atRow["pdbx_leaving_atom_flag"] = atom.GetStringData("pdbx_leaving_atom_flag")
            else:
                atRow["pdbx_leaving_atom_flag"] = "N"
            if len(atRow["atom_id"]) > 3 or len(atRow["type_symbol"]) == 2:
                atRow["pdbx_align"] = 0
            else:
                atRow["pdbx_align"] = 1

            if atom.IsAromatic():
                atRow["pdbx_aromatic_flag"] = "Y"
            else:
                atRow["pdbx_aromatic_flag"] = "N"
            # oeSt = oechem.OEGetCIPStereo(mol, atom)
            oeSt = None
            # if atom.IsChiral():
            cip = oechem.OEPerceiveCIPStereo(oeMol, atom)
            if atom.HasStereoSpecified():
                if cip == oechem.OECIPAtomStereo_S:
                    oeSt = "S"
                if cip == oechem.OECIPAtomStereo_R:
                    oeSt = "R"
                if cip == oechem.OECIPAtomStereo_NotStereo:
                    oeSt = None
                if cip == oechem.OECIPAtomStereo_UnspecStereo:
                    oeSt = "U"
                # oeSt = oechem.OEGetCIPStereo(oeMol, atom)
                # oeSt = atom.GetCIPStereo()
            if oeSt is not None and (len(oeSt) > 0) and (oeSt == "R" or oeSt == "S"):
                atRow["pdbx_stereo_config"] = oeSt
            else:
                atRow["pdbx_stereo_config"] = "N"
            atRow["pdbx_component_atom_id"] = atom.GetName().strip()
            atRow["pdbx_component_comp_id"] = idCode
            atRow["pdbx_ordinal"] = str(ii + 1)
            rowL.append(atRow)
        #
        return rowL

    def __makeChemCompBondCategory(self, ccId, oeMol):
        idCode = ccId
        rowL = []
        iBondD = {1: "SING", 2: "DOUB", 3: "TRIP", 4: "QUAD", 5: "AROM", 0: "DELO"}

        for ii, bond in enumerate(oeMol.GetBonds()):
            oeOrder = str(bond.GetOrder()).upper()
            b1 = bond.GetBgn().GetName().strip()
            b2 = bond.GetEnd().GetName().strip()
            if (b1 is None) or (len(b1) < 1):
                continue
            if (b2 is None) or (len(b1) < 1):
                continue
            if int(oeOrder) in iBondD:
                bndRow = {}
                bndRow["comp_id"] = idCode
                bndRow["atom_id_1"] = b1
                bndRow["atom_id_2"] = b2
                bndRow["value_order"] = iBondD[int(oeOrder)]

                if bond.IsAromatic():
                    bndRow["pdbx_aromatic_flag"] = "Y"
                else:
                    bndRow["pdbx_aromatic_flag"] = "N"

                oeSt = None
                if bond.GetOrder() == 2:
                    cip = oechem.OEPerceiveCIPStereo(oeMol, bond)
                    if bond.HasStereoSpecified():
                        if cip == oechem.OECIPBondStereo_E:
                            oeSt = "E"
                        if cip == oechem.OECIPBondStereo_Z:
                            oeSt = "Z"
                        if cip == oechem.OECIPBondStereo_NotStereo:
                            oeSt = None
                        if cip == oechem.OECIPBondStereo_UnspecStereo:
                            oeSt = None
                if oeSt is not None and (len(oeSt) > 0) and (oeSt == "E" or oeSt == "Z"):
                    bndRow["pdbx_stereo_config"] = oeSt
                else:
                    bndRow["pdbx_stereo_config"] = "N"

                bndRow["pdbx_ordinal"] = str(ii + 1)
                rowL.append(bndRow)
        return rowL

    def __getElementCounts(self, oeMol):
        """Get the dictionary of element counts (eg. eD[iAtNo]=iCount)."""
        eD = {}
        if len(eD) == 0:
            # calculate from current oeMol
            try:
                eD = {}
                for atom in oeMol.GetAtoms():
                    atNo = atom.GetAtomicNum()
                    if atNo not in eD:
                        eD[atNo] = 1
                    else:
                        eD[atNo] += 1
            except Exception:
                pass
        return eD
