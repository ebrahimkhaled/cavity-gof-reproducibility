# Shared helpers for the numerical checks of the linear leave-one-out index.
import numpy as np
from numpy.polynomial.hermite_e import hermegauss

GX, GW = hermegauss(32)
GW = GW / GW.sum()

def expit(x):
    return 1.0 / (1.0 + np.exp(-x))

def Lam(u, a, B, tau2):
    """Logistic-normal curve; a, B, tau2 scalars or arrays that broadcast with u."""
    tau = np.sqrt(np.maximum(tau2, 0.0))
    z = (a + B * u + 0 * tau)[..., None] + np.asarray(tau)[..., None] * GX
    return expit(z) @ GW

def Lam_a(u, a, B, tau2):
    tau = np.sqrt(np.maximum(tau2, 0.0))
    z = expit((a + B * u + 0 * tau)[..., None] + np.asarray(tau)[..., None] * GX)
    return (z * (1 - z)) @ GW

def Lam_tau2(u, a, B, tau2):
    tau = np.sqrt(np.maximum(tau2, 0.0))
    z = expit((a + B * u + 0 * tau)[..., None] + np.asarray(tau)[..., None] * GX)
    return (0.5 * z * (1 - z) * (1 - 2 * z)) @ GW

def theta_of(sig1, sig2, sig3, m0, s2):
    """(a, B, tau2) from sigma = (c0'b, beta'St b, b'St b)."""
    B = sig2 / sig3
    return m0 - B * sig1, B, np.maximum(s2 - sig2 ** 2 / sig3, 0.0)

class Design:
    def __init__(self, n, p, lam, eps0=0.0, s=2.0, beta00=0.0, seed=1):
        self.n, self.p, self.lam, self.eps0 = n, p, lam, eps0
        self.b0 = np.full(p, s / np.sqrt(p))
        self.beta00 = beta00
        self.m0, self.s2 = beta00, s * s
        self.D = np.r_[eps0, np.ones(p)] * lam

    def draw(self, rng, m=None):
        m = self.n if m is None else m
        X = rng.standard_normal((m, self.p))
        A = np.c_[np.ones(m), X]
        y = (rng.random(m) < expit(self.beta00 + X @ self.b0)).astype(float)
        return A, y

    def sigma(self, Bm):
        """Rows of Bm are coefficient vectors (intercept first). Sigma = I, mu = 0."""
        return Bm[:, 0], Bm[:, 1:] @ self.b0, np.einsum('ij,ij->i', Bm[:, 1:], Bm[:, 1:])
