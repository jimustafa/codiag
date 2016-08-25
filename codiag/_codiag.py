from __future__ import absolute_import, division, print_function

import numpy as np


def diag2(a, n):
    r"""
    Compute the sum of the square moduli of the diagonal elements

    Parameters
    ----------
    a : ndarray, shape (..., N, N)
    n : int
        the :math:`n \times n` sub-block
    """
    return np.abs(np.sum(np.abs(np.diagonal(a, axis1=-2, axis2=-1))[..., :n]**2))


def offdiag2(a, n):
    r"""
    Compute the sum of the square moduli of the off-diagonal elements

    Parameters
    ----------
    a : ndarray, shape (..., N, N)
    n : int
        the :math:`n \times n` sub-block
    """
    return np.abs(np.vdot(a[..., :n, :n], a[..., :n, :n]) - diag2(a, n))


def diag2_ii(a, i):
    r"""
    Compute the sum of the square moduli of the diagonal :math:`ii` elements

    Parameters
    ----------
    a : ndarray, shape (..., N, N)
    i : int
        diagonal component
    """
    return np.real(np.sum(np.abs(np.diagonal(a, axis1=-2, axis2=-1))[..., i]**2))


def diag2_iijj(a, i, j):
    r"""
    Compute the sum of the square moduli of the diagonal :math:`ii` and
    :math:`jj` elements

    Parameters
    ----------
    a : ndarray, shape (..., N, N)
    i : int
        diagonal component
    """
    return np.real(np.sum(np.abs(np.diagonal(a, axis1=-2, axis2=-1))[..., (i, j)]**2))


def offdiag2_ij(a, i, j):
    idx = list(np.ogrid[[slice(x) for x in a.shape]][:-2])+[i]+[j]
    return np.real(np.sum(np.abs(a[idx])**2))


def offdiag2_ijji(a, i, j):
    idx1 = list(np.ogrid[[slice(x) for x in a.shape]][:-2])+[i]+[j]
    idx2 = list(np.ogrid[[slice(x) for x in a.shape]][:-2])+[j]+[i]
    return np.real(np.sum(np.abs(a[idx1])**2)) + np.real(np.sum(np.abs(a[idx2])**2))


def v2cs(v):
    v = v / np.linalg.norm(v)
    c = np.sqrt((1 + v[0]) / 2)
    if np.abs(c) < 1e-9:
        # s = np.exp(1j*np.random.random()*2*np.pi)
        s = 1
    else:
        s = (v[1] + 1j*v[2]) / (2*c)

    return c, s


def _diag2_ii_Abc(a, i, j):
    aii = a[:, i, i]
    aij = a[:, i, j]
    aji = a[:, j, i]
    ajj = a[:, j, j]

    aii_ajj = aii+ajj

    c = 1/4 * np.sum(np.real(aii_ajj*aii_ajj.conj()))
    z = np.empty((3, len(a)), dtype=complex)
    z[0] = aii-ajj
    z[1] = -(aij+aji)
    z[2] = 1j*(aij-aji)
    z *= 1/2
    b = np.real(np.sum(np.conj(aii_ajj) * z, axis=1))

    A = np.real(np.dot(z, z.conj().T))

    return A, b, c


def diag2_ii_Abc(a, i, j):
    return _diag2_ii_Abc(a, i, j)


def diag2_iijj_Abc(a, i, j):
    A, _, c = _diag2_ii_Abc(a, i, j)

    return 2*A, None, 2*c


def offdiag2_ijji_Abc(a, i, j):
    aii = a[:, i, i]
    aij = a[:, i, j]
    aji = a[:, j, i]
    ajj = a[:, j, j]

    aii2 = np.real(aii*aii.conj())
    aij2 = np.real(aij*aij.conj())
    aji2 = np.real(aji*aji.conj())
    ajj2 = np.real(ajj*ajj.conj())

    c = 1/4 * np.sum(aij2 + aji2)

    A = np.zeros((3, 3))
    A[0, 0] = 1/2 * np.sum(aij2 + aji2)
    A[1, 1] = 1/2 * np.sum(aii2 + ajj2) - np.sum(np.real(aij*aji.conj()+aii*ajj.conj()))
    A[2, 2] = 1/2 * np.sum(aii2 + ajj2) + np.sum(np.real(aij*aji.conj()-aii*ajj.conj()))

    A[0, 1] = 1/2 * np.sum(np.real((aii-ajj)*(aij.conj()+aji.conj())))
    A[1, 0] = 1/2 * np.sum(np.real((aii-ajj)*(aij.conj()+aji.conj())))

    A[0, 2] = 1/2 * np.sum(np.real((aii-ajj)*np.imag(aij-aji)+(-aij+aji)*np.imag(aii-ajj)))
    A[2, 0] = 1/2 * np.sum(np.real((aii-ajj)*np.imag(aij-aji)+(-aij+aji)*np.imag(aii-ajj)))

    A[1, 2] = np.sum(aji.imag*aij.real-aij.imag*aji.real)
    A[2, 1] = np.sum(aji.imag*aij.real-aij.imag*aji.real)

    return A, None, 2*c


def _offdiag2_ikki_Abc(a, i, j, k):
    k = np.asarray(k)
    assert i not in k
    assert j not in k

    aik = a[:, i, k]
    aki = a[:, k, i]
    ajk = a[:, j, k]
    akj = a[:, k, j]

    aik2 = np.sum(np.real(aik*aik.conj()))
    aki2 = np.sum(np.real(aki*aki.conj()))
    ajk2 = np.sum(np.real(ajk*ajk.conj()))
    akj2 = np.sum(np.real(akj*akj.conj()))

    A = np.zeros((3, 3))

    b = np.zeros((3,))
    b[0] = 1/2 * (aik2 - ajk2 + aki2 - akj2)
    b[1] = -np.sum(np.real(aik*ajk.conj() + aki*akj.conj()))
    b[2] = np.sum(ajk.imag*aik.real - aik.imag*ajk.real - akj.imag*aki.real + aki.imag*akj.real)

    c = 1/2 * (aik2 + ajk2 + aki2 + akj2)

    return A, b, c


def offdiag2_ikki_Abc(a, i, j, k):
    return _offdiag2_ikki_Abc(a, i, j, k)


def offdiag2_jkkj_Abc(a, i, j, k):
    A, b, c = _offdiag2_ikki_Abc(a, i, j, k)

    return A, -b, c
