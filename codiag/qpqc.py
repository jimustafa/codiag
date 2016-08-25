from __future__ import absolute_import, division, print_function

import numpy as np

from . import qep


def _solve_qp_s2_Q(Q):
    # get the eigenvector corresponding to the largest eigenvalue
    eigvals, eigvecs = np.linalg.eigh(Q)
    eval_imax = np.argmax(eigvals)
    # if eigvals[eval_imax] < 0:
    #     raise Exception

    x = eigvecs[:, eval_imax]
    x /= np.linalg.norm(x)

    return x


def _solve_qp_s2_Qp(Q, p):
    A2 = np.eye(3)
    A1 = -2*Q
    A0 = np.dot(Q, Q) - 1.0/4 * np.outer(p, p)

    eigvals_Q, eigvecs_Q = np.linalg.eigh(Q)
    eigvecs_Q = eigvecs_Q.T

    eigvals = qep.quadeigvals(A0, A1, A2)

    eigvals_real = np.real(eigvals[np.abs(eigvals.imag) < 1e-6])
    eigval_max = np.max(eigvals_real)

    # for (i, eigval) in enumerate(np.sort(eigvals_real)[::-1]):
    if np.any(np.abs(eigvals_Q-eigval_max) < 1e-6):
        print('WARNING')
        iq = np.nonzero(np.abs(eigvals_Q-eigval_max) < 1e-6)[0][0]
        # u = np.dot(np.linalg.pinv(Q-eigval*np.eye(3)), -p/2)
        # if np.sum((np.dot(Q-eigval*np.eye(3), u) - (-p/2))**2) > 1e-6 or np.sum(u**2) > 1:
        u = np.dot(np.linalg.pinv(Q-eigval_max*np.eye(3)), -p/2)
        if np.sum((np.dot(Q-eigval_max*np.eye(3), u) - (-p/2))**2) > 1e-6 or np.sum(u**2) > 1:
            return None
        else:
            if abs(np.sum(u**2) - 1) < 1e-6:
                return u
            else:
                v = eigvecs_Q[iq]
                c = np.sqrt(np.dot(v.conj(), v)/(1-np.dot(u, u)))
                return u + v/c
    else:
        return np.dot(np.linalg.inv(Q - eigval_max*np.eye(3)), -p/2)


def solve_qp_s2(Q, p=None):
    if p is None:
        return _solve_qp_s2_Q(Q)
    else:
        return _solve_qp_s2_Qp(Q, p)


def solve_qp_s2_brute(Q, p=None, ntheta=1025, nphi=1025):
    r"""
    Solve the quadratic program xQx + px on :math:`S^{2}`

    """
    theta = np.linspace(0, np.pi, ntheta)
    phi = np.linspace(0, 2*np.pi, nphi)

    theta, phi = np.meshgrid(phi, theta, indexing='ij')

    theta = theta.ravel()
    phi = phi.ravel()

    x1 = np.cos(theta)
    x2 = np.sin(theta)*np.cos(phi)
    x3 = np.sin(theta)*np.sin(phi)

    x = np.column_stack((x1, x2, x3))

    f = np.einsum('ij,jk,ik->i', x, Q, x)

    if p is not None:
        f += np.einsum('ij,j->i', x, p)

    imax = np.argmax(f)

    return x[imax]
