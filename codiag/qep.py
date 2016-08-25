from __future__ import absolute_import, division, print_function

import numpy as np
import scipy.linalg


def quadeig(A0, A1, A2):
    r"""
    Solve :math:`Q(\lambda)\mathbf{x}=\mathbf{0}`
    where :math:`Q(\lambda) = \lambda^{2}A_{2} + \lambda^{1}A_{1} + A_{0}
    """
    assert A2.shape == A1.shape
    assert A1.shape == A0.shape

    n = A2.shape[0]

    A = np.array(np.bmat([[A1, A0], [-np.eye(n), np.zeros((n, n))]]))
    B = -np.array(np.bmat([[A2, np.zeros((n, n))], [np.zeros((n, n)), np.eye(n)]]))

    eigvals, eigvecs = scipy.linalg.eig(A, b=B)
    eigvecs = eigvecs[:n, :]

    return eigvals, eigvecs


def quadeigvals(A0, A1, A2):
    r"""
    Solve :math:`Q(\lambda)\mathbf{x}=\mathbf{0}`
    where :math:`Q(\lambda) = \lambda^{2}A_{2} + \lambda^{1}A_{1} + A_{0}
    """
    assert A2.shape == A1.shape
    assert A1.shape == A0.shape

    n = A2.shape[0]

    A = np.array(np.bmat([[A1, A0], [-np.eye(n), np.zeros((n, n))]]))
    B = -np.array(np.bmat([[A2, np.zeros((n, n))], [np.zeros((n, n)), np.eye(n)]]))

    eigvals = scipy.linalg.eigvals(A, b=B)

    return eigvals


def quadeigh(A0, A1, A2):
    r"""
    Solve :math:`Q(\lambda)\mathbf{x}=\mathbf{0}`
    where :math:`Q(\lambda) = \lambda^{2}A_{2} + \lambda^{1}A_{1} + A_{0}
    """
    assert A2.shape == A1.shape
    assert A1.shape == A0.shape

    n = A2.shape[0]

    A = np.array(np.bmat([[A1, A0], [A0, np.zeros((n, n))]]))
    B = -np.array(np.bmat([[A2, np.zeros((n, n))], [np.zeros((n, n)), A0]]))

    eigvals, eigvecs = scipy.linalg.eig(A, b=B)
    eigvecs = eigvecs[:n, :]

    return eigvals, eigvecs


def quadeigvalsh(A0, A1, A2):
    r"""
    Solve :math:`Q(\lambda)\mathbf{x}=\mathbf{0}`
    where :math:`Q(\lambda) = \lambda^{2}A_{2} + \lambda^{1}A_{1} + A_{0}
    """
    assert A2.shape == A1.shape
    assert A1.shape == A0.shape

    n = A2.shape[0]

    A = np.array(np.bmat([[A1, A0], [A0, np.zeros((n, n))]]))
    B = -np.array(np.bmat([[A2, np.zeros((n, n))], [np.zeros((n, n)), A0]]))

    eigvals = scipy.linalg.eigvals(A, b=B)

    return eigvals
