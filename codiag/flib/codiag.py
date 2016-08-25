from __future__ import absolute_import, division, print_function

import numpy as np

from . import _flib


def diag2_ii_Abc(a, i, j):
    assert np.isfortran(a)

    c = np.array(0, dtype=float)
    b = np.zeros((3,))
    A = np.zeros((3, 3), order='F')

    _flib.codiag.diag2_ii_abc(a, len(a), i+1, j+1, A, b, c)

    return A, b, c


def _offdiag2_ikki_Abc(a, i, j, k):
    assert np.isfortran(a)

    k = np.asarray(k)
    assert i not in k
    assert j not in k

    c = np.array(0, dtype=float)
    b = np.zeros((3,))
    A = np.zeros((3, 3), order='F')

    _flib.codiag.offdiag2_ikki_abc(a, len(a), i+1, j+1, k+1, A, b, c)

    return A, b, c


def offdiag2_ikki_Abc(a, i, j, k):
    return _offdiag2_ikki_Abc(a, i, j, k)


def offdiag2_jkkj_Abc(a, i, j, k):
    A, b, c = _offdiag2_ikki_Abc(a, i, j, k)

    return A, -b, c
