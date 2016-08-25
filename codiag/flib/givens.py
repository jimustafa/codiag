from __future__ import absolute_import, division, print_function

import numpy as np

from . import _flib


def left_multiply(m, i, j, c, s):
    assert np.isfortran(m)
    _flib.givens.left_multiply(m, i+1, j+1, np.real(c), complex(s))


def right_multiply(m, i, j, c, s):
    assert np.isfortran(m)
    _flib.givens.right_multiply(m, i+1, j+1, np.real(c), complex(s))


def rotate(m, i, j, c, s):
    assert np.isfortran(m)
    _flib.givens.rotate(m, i+1, j+1, np.real(c), complex(s))
