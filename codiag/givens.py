from __future__ import absolute_import, division, print_function

import numpy as np


def random(n, i, j):
    r"""
    Generate a random Givens rotation

    Parameters
    ----------
    n : integer
        dimension of the square Givens rotation :math:`G`
    i : integer
    j : integer

    """
    theta = np.random.random() * 2*np.pi
    phi = np.random.random() * 2*np.pi

    c = np.cos(theta)
    s = np.sin(theta)*np.exp(1j*phi)

    R = np.identity(n, dtype=complex)
    R[i, i] = c
    R[i, j] = s
    R[j, i] = -np.conj(s)
    R[j, j] = c

    return R


def rotation(c, s, n, i, j):
    R = np.identity(n, dtype=complex)
    R[i, i] = c
    R[i, j] = s
    R[j, i] = -1*np.conj(s)
    R[j, j] = c
    return R


def left_multiply(a, i, j, c, s):
    x1 = np.copy(a[..., i, :])
    x2 = np.copy(a[..., j, :])

    a[..., i, :] = c*x1 - s*x2
    a[..., j, :] = np.conj(s)*x1 + c*x2

    return a


def right_multiply(a, i, j, c, s):
    x1 = np.copy(a[..., :, i])
    x2 = np.copy(a[..., :, j])

    a[..., :, i] = c*x1 - np.conj(s)*x2
    a[..., :, j] = s*x1 + c*x2

    return a
