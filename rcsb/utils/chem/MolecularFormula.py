##
# File: MolecularFormula.py
# Date: 14-Mar-2020
#
# Adapted from recipe a recipe posted by at:
# https://github.com/Zapaan/python-chemical-formula-parser
##

import logging
import re
from collections import Counter

logger = logging.getLogger(__name__)


class MolecularFormula(object):
    def __init__(self):
        self.__regex = r"([A-Z][a-z]*)(\d*)"
        self.__openGroup = "({["
        self.__closeGroup = ")}]"

    def __isBalanced(self, formula):
        """Check if all grouping elements are properly paired by count."""
        cnt = Counter(formula)
        return cnt["["] == cnt["]"] and cnt["{"] == cnt["}"] and cnt["("] == cnt[")"]

    def __toDict(self, tuples):
        """Transform tuples of tuples to a dict of elements."""
        res = dict()
        for atom, num in tuples:
            try:
                res[atom] += int(num or 1)
            except KeyError:
                res[atom] = int(num or 1)
        return res

    def __mergeDict(self, mol1, mol2, weight=1):
        """
        Merge 2 dicts representing molecules. Return a new dict.
        """
        return {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * weight for atom in set(mol1) | set(mol2)}

    def __parse(self, formula):
        """
        Return an element dictionary and length of parsed part.

        Recurse on grouped elements to parse the subpart.
        """
        qL = []
        mol = {}
        ii = 0

        while ii < len(formula):
            # Using a classic loop allow for manipulating the cursor
            token = formula[ii]

            if token in self.__closeGroup:
                # Check for an index for this part
                mtch = re.match(r"\d+", formula[ii + 1 :])
                if mtch:
                    weight = int(mtch.group(0))
                    ii += len(mtch.group(0))
                else:
                    weight = 1

                submol = self.__toDict(re.findall(self.__regex, "".join(qL)))
                return self.__mergeDict(mol, submol, weight), ii

            elif token in self.__openGroup:
                submol, ll = self.__parse(formula[ii + 1 :])
                mol = self.__mergeDict(mol, submol)
                # skip the already read submol
                ii += ll + 1
            else:
                qL.append(token)

            ii += 1
        # merge in all that's left at base level
        return self.__mergeDict(mol, self.__toDict(re.findall(self.__regex, "".join(qL)))), ii

    def parseFormula(self, formula):
        """Parse the formula and return a dictionary of counts of each atom type."""
        rD = self.__parse(formula)[0]
        logger.debug("formula %r result %r", formula, rD)
        return rD
