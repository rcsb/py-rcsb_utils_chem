##
# File:    IoUtils.py
# Author:  jdw
# Date:    31-Dec-2019
# Version: 0.001
#
# Updates:
#
##
"""
Utilities to manage general component IO and format conversion operations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging

from rcsb.utils.chem.PdbxChemComp import PdbxChemCompAtomIt, PdbxChemCompBondIt, PdbxChemCompIt
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class IoUtils(object):
    """ Utility methods to manage OE specific IO and format conversion operations.
    """

    def __init__(self, **kwargs):
        self.__dirPath = kwargs.get("dirPath", ".")
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        #

    def getComponentDefinitions(self, ccdFilePath):
        try:
            rdCcObjL = self.__mU.doImport(ccdFilePath, fmt="mmcif")
            logger.info("Read %s with %d definitions", ccdFilePath, len(rdCcObjL))
        except Exception as e:
            logger.exception("Loading %s failing with %s", ccdFilePath, str(e))
        return rdCcObjL

    def makeSdf(self, dataContainer, molBuildType="model-xyz", useAromatic=False):
        sdf = ""
        atomIdxD = {}
        ok = True
        try:
            if not (dataContainer.exists("chem_comp") and dataContainer.exists("chem_comp_atom") and dataContainer.exists("chem_comp_bond")):
                return False, sdf, atomIdxD
            ccIt = PdbxChemCompIt(dataContainer)
            for cc in ccIt:
                ccId = cc.getId()
            #
            now = datetime.datetime.now().strftime("%y%m%d%H%M")
            chargeMap = {1: 3, 2: 2, 3: 1, -1: 5, -2: 6, -3: 7, 0: 0}
            # bondMap = {"arom": 4, "delo": 5, "doub": 2, "pi": 5, "poly": 5, "quad": 8, "sing": 1, "trip": 3}
            oL = []
            atomIt = PdbxChemCompAtomIt(dataContainer)
            numAtom = len(atomIt)
            bondIt = PdbxChemCompBondIt(dataContainer)
            numBond = len(bondIt)
            # Header Line 1 - title -- unformatted text --
            oL.append("%s" % ccId)
            # Header Line 2 - xxNNNNNNNNMMDDYYhhmm
            # oL.append("%2s%8s%10s" % ("  ", "CCTOOLS-", now))
            oL.append("%2s%8s%10s" % ("  ", "CHMUTIL-", now))
            # oL.append("  -ISIS-            %s" % "3D")
            #
            oL.append("")
            chiralFlag = 0
            # Connection table - Counts line ---
            sL = "%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%6s" % (numAtom, numBond, 0, 0, chiralFlag, 0, 0, 0, 0, 0, 999, " V2000")
            oL.append(sL)
            #
            aD = {}
            fcL = []
            isotopeL = []
            iAtom = 1
            atomIt = PdbxChemCompAtomIt(dataContainer)
            for ccAt in atomIt:
                atName = ccAt.getName()
                aD[atName] = iAtom
                atomIdxD[iAtom] = atName
                atNo = ccAt.getAtNo()
                atType = ccAt.getType()
                fc = int(ccAt.getFormalCharge())
                if fc != 0:
                    fcL.append("%4s%4d" % (iAtom, fc))
                # chFlag = ccAt.isChiral()
                # arFlag = ccAt.isAromatic()
                isotope = ccAt.getIsotope()
                if isotope != 0:
                    isotopeL.append("%4s%4d" % (iAtom, isotope))
                cTup = None
                if (molBuildType == "model-xyz") and ccAt.hasModelCoordinates():
                    cTup = ccAt.getModelCoordinates()
                elif (molBuildType == "ideal-xyz") and ccAt.hasIdealCoordinates():
                    cTup = ccAt.getIdealCoordinates()
                else:
                    logger.warning("%s has no %s records", dataContainer.getName(), molBuildType)
                    ok = False
                    cTup = (0.0, 0.0, 0.0)
                iAtom += 1
                #
                aS = "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d" % (cTup[0], cTup[1], cTup[2], atType, 0, chargeMap[fc], 0, 0, 0)
                oL.append(aS)
                logger.debug("Atom - %s type %s atno %d isotope %d fc %d (xyz) %r", atName, atType, atNo, isotope, fc, cTup)
                #

            bondIt = PdbxChemCompBondIt(dataContainer)
            for ccBnd in bondIt:
                (at1, at2) = ccBnd.getBond()

                iType = ccBnd.getIntegerType()
                if useAromatic and ccBnd.isAromatic():
                    iType = 4
                bS = "%3d%3d%3d%3d%3d%3d" % (aD[at1], aD[at2], iType, 0, 0, 0)
                oL.append(bS)
                logger.debug(" %s %d -- %s %d (%d)", at1, aD[at1], at2, aD[at2], iType)
            #
            # Add charge line -
            if fcL:
                oL.append("M  CHG%3s%s" % (len(fcL), "".join(fcL)))
            if isotopeL:
                oL.append("M  ISO%3s%s" % (len(isotopeL), "".join(isotopeL)))
            #
            oL.append("M  END")
            oL.append("$$$$")
            #
            sdf = "\n".join(oL)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False

        return ok, sdf, atomIdxD
