##
# File: OeCommonUtils.py
# Date: 22-Oct-2020
#
##
import logging

from openeye import oechem

logger = logging.getLogger(__name__)


class OeCommonUtils(object):
    @staticmethod
    def getAtomBondExprOpts(matchOpts):
        """Return the OE constants for atom and bond matching criteria.

        Args:
            matchOpts (string): qualitative description of atom and bond matching criteria

        Returns:
            tuple: OE atom and bond comparison expressions
        """
        if matchOpts in ["default", "strict", "graph-strict", "graph-default", "sub-struct-graph-strict"]:
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_Aromaticity
            bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Aromaticity | oechem.OEExprOpts_Chiral
        elif matchOpts in ["relaxed-stereo", "graph-relaxed-stereo", "sub-struct-graph-relaxed-stereo"]:
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_FormalCharge
            bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Chiral
        elif matchOpts in ["relaxed-stereo-sdeq", "graph-relaxed-stereo-sdeq", "sub-struct-graph-relaxed-stereo-sdeq"]:
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_FormalCharge | oechem.OEExprOpts_Aromaticity
            bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_Chiral | oechem.OEExprOpts_EqSingleDouble
        elif matchOpts in ["relaxed", "graph-relaxed", "simple", "sub-struct-graph-relaxed"]:
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge
            bondexpr = oechem.OEExprOpts_BondOrder
        elif matchOpts in ["graph-relaxed-sdeq", "sub-struct-graph-relaxed-sdeq"]:
            atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_FormalCharge | oechem.OEExprOpts_Aromaticity
            bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_EqSingleDouble
        elif matchOpts in ["exact"]:
            atomexpr = oechem.OEExprOpts_ExactAtoms
            bondexpr = oechem.OEExprOpts_ExactBonds
        else:
            atomexpr = oechem.OEExprOpts_DefaultAtoms
            bondexpr = oechem.OEExprOpts_DefaultBonds
            logger.error("Unanticipated match options %r (applying defaults)", matchOpts)
        #
        #
        # oechem.OEExprOpts_EqSingleDouble
        # atomexpr = oechem.OEExprOpts_AtomicNumber|oechem.OEExprOpts_EqAromatic
        # bondexpr = 0
        #
        # atomexpr = oechem.OEExprOpts_AtomicNumber|oechem.OEExprOpts_Aromaticity
        # bondexpr = oechem.OEExprOpts_BondOrder|oechem.OEExprOpts_EqNotAromatic
        #
        return atomexpr, bondexpr
