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
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

#


def annotateCactvs(molFilePath):
    retD = {}
    try:
        cactvs["aromaticity_model"] = "daylight"
        #
        eh = Molfile(molFilePath).read()
        # print(type(eh))
        # print(eh.atoms())
        atomL = []
        for ii, at in enumerate(eh.atoms(), 1):
            aD = {}
            st = at.A_LABEL_STEREO
            aD = {
                "idx": ii,
                "label": at.A_LABEL,
                "type": at.A_SYMBOL,
                "CIP": at.A_CIP_STEREO if at.A_CIP_STEREO != "undef" else None,
                "isAromatic": at.A_ISAROMATIC,
                "formalCharge": at.A_FORMAL_CHARGE,
            }
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
        retD = {"descriptors": dD, "atoms": atomL, "bonds": bondL}
    except Exception as e:
        logger.exception("Failing for %s with %s", molFilePath, str(e))
    return retD


if __name__ == "__main__":
    #
    print(sys.argv)
    retD = annotateCactvs(sys.argv[0])
    with open(sys.argv[1], "w") as ofh:
        json.dump(retD, ofh)
#
