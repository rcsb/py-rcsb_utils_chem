## File: cactvsAnnotateMol.py
# Date: Jan 7, 2020 jdw
#
# $cactvsPath/cspy cactvsAnnotateMol.py <input sdf path> <output json path>
#
# This module is using an external self-contained python library in CACTVS bundle
##
# flake8: noqa
# pylint: skip-file
import json
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


def annotateCactvs(molFilePath, aroModel="cactvs"):
    retD = {}
    try:
        cactvs["aromaticity_model"] = aroModel if aroModel in ["cactvs", "daylight", "tripos"] else "cactvs"
        logger.info("aroModel is %r", aroModel)
        #
        eh = Molfile(molFilePath).read()
        # print(type(eh))
        # print(eh.atoms())
        atomL = []
        iChargeTotal = 0
        for ii, at in enumerate(eh.atoms(), 1):
            aD = {}
            st = at.A_LABEL_STEREO != "undef"
            aD = {
                "idx": ii,
                "label": at.A_LABEL,
                "type": at.A_SYMBOL,
                "CIP": at.A_CIP_STEREO if at.A_CIP_STEREO != "undef" else None,
                "isChiral": st,
                "isAromatic": at.A_ISAROMATIC,
                "formalCharge": at.A_FORMAL_CHARGE,
            }
            iChargeTotal += int(at.A_FORMAL_CHARGE)
            logger.debug("%d %r %r %r %r %r %r", ii, at.A_LABEL, at.A_SYMBOL, at.A_CIP_STEREO, at.A_ISAROMATIC, at.A_FORMAL_CHARGE, at.A_XYZ)
            atomL.append(aD)
        #
        bondL = []
        for ii, bnd in enumerate(eh.bonds(), 1):
            # st = bnd.B_LABEL_STEREO
            aL = bnd.atoms()
            atI = aL[0].A_LABEL
            atJ = aL[1].A_LABEL
            bD = {}
            bD = {
                "idx": ii,
                "label": bnd.B_LABEL,
                "atI": atI,
                "atJ": atJ,
                "value_order": bnd.B_ORDER,
                "CIP": bnd.B_CIP_STEREO if bnd.B_CIP_STEREO != "undef" else None,
                "isAromatic": bnd.B_ISAROMATIC,
            }
            logger.debug("%d %r %r %r %r %r %r ", ii, bnd.B_LABEL, atI, atJ, bnd.B_ORDER, bnd.B_CIP_STEREO, bnd.B_ISAROMATIC)
            bondL.append(bD)

        Prop.Setparam("E_SMILES", {"usearo": True, "usestereo": True, "unique": True})
        canIsoSmi = eh.new("E_SMILES")
        logger.debug("Canonical isomeric smiles %s", canIsoSmi)

        Prop.Setparam("E_SMILES", {"usearo": True, "usestereo": False, "unique": True})
        canSmi = eh.new("E_SMILES")
        logger.debug("Canonical smiles %s", canSmi)
        dD = {"canSmi": canSmi, "canIsoSmi": canIsoSmi}
        #
        # print(dir(eh))
        # print(eh.props())
        #
        detailD = {"ccId": eh.E_MDL_NAME, "formula": eh.E_FORMULA, "ifCharge": iChargeTotal, "formulaWeight": eh.E_WEIGHT}
        #
        retD = {"details": detailD, "descriptors": dD, "atoms": atomL, "bonds": bondL}
    except Exception as e:
        logger.exception("Failing for %s with %s", molFilePath, str(e))
    return retD


def isWritable(filePath):
    if os.path.exists(filePath):
        if os.path.isfile(filePath):
            return os.access(filePath, os.W_OK)
        else:
            return False
    dirPath = os.path.dirname(filePath)
    if not dirPath:
        dirPath = "."
    return os.access(dirPath, os.W_OK)


def main():
    if len(sys.argv) < 2:
        sys.stdout.write("script <inputMolFilePath> <outputAnnotationPath> (%r)\n" % sys.argv)
        return False
    if not os.access(sys.argv[0], os.R_OK):
        sys.stdout.write("script <inputMolFilePath> <outputAnnotationPath> (%r)\n" % sys.argv)
        return False
    #
    if isWritable(sys.argv[1]):
        retD = annotateCactvs(sys.argv[0], aroModel=sys.argv[2])
        if retD:
            with open(sys.argv[1], "w") as ofh:
                json.dump(retD, ofh)
            return True
    return False


if __name__ == "__main__":
    main()
