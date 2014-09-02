import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from text import *
from procutil import *


def run_variational_decomp_fci(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    variational_decomp_fci can be called via :py:func:`~driver.energy`.

    >>> energy('variational_decomp_fci')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Your plugin's PsiMod run sequence goes here
    PsiMod.set_global_option('BASIS', 'sto-3g')
    PsiMod.set_local_option('VARIATIONAL_DECOMP_FCI', 'PRINT', 1)
    energy('scf', **kwargs)
    returnvalue = PsiMod.plugin('variational_decomp_fci.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['variational_decomp_fci'] = run_variational_decomp_fci


def exampleFN():
    # Your Python code goes here
    pass
