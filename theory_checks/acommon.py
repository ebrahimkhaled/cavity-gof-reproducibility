# Independent audit helpers. General design: mu != 0, beta00 != 0, Sigma != I, eps0 > 0, p > n allowed.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
GX, GW = hermegauss(40); GW = GW / GW.sum()

def expit(x):
    return 0.5 * (1 + np.tanh(0.5 * x))

class Des:
    def __init__(s, n, p, lamn, eps0=0.0, sig=2.0, b00=0.0, mu=0.0, ar=0.0, seed=5, generic=True):
        s.n, s.p, s.lam, s.eps0 = n, p, lamn * n, eps0
        r = np.random.default_rng(seed)
        i = np.arange(p)
        Sig = ar ** np.abs(i[:, None] - i[None, :]) if ar > 0 else np.eye(p)
        s.Sig = Sig; s.Lc = np.linalg.cholesky(Sig)
        d = r.standard_normal(p) if generic else np.ones(p)
        s.b0 = d * sig / np.sqrt(d @ Sig @ d)
        m = r.standard_normal(p); s.mu = mu * m / np.linalg.norm(m)
        s.b00 = b00; s.c = np.r_[1.0, s.mu]; s.m0 = b00 + s.mu @ s.b0; s.s2 = sig ** 2
        s.St = np.zeros((p + 1, p + 1)); s.St[1:, 1:] = Sig
        s.Stb = np.r_[0.0, Sig @ s.b0]; s.D = np.r_[eps0, np.ones(p)] * s.lam

    def draw(s, r, m):
        X = s.mu + r.standard_normal((m, s.p)) @ s.Lc.T
        A = np.c_[np.ones(m), X]
        eta = s.b00 + X @ s.b0
        return A, (r.random(m) < expit(eta)).astype(float), eta

    def theta(s, s1, s2, s3):
        B = s2 / s3
        return s.m0 - B * s1, B, np.maximum(s.s2 - s2 ** 2 / s3, 0.0)

def Lam(u, a, B, t2):
    z = a + B * u; tau = np.sqrt(t2)
    out = np.zeros(np.broadcast(z, tau).shape)
    for x, w in zip(GX, GW):
        out += w * expit(z + tau * x)
    return out

def Lders(u, a, B, t2):
    """Lambda, Lambda_a, Lambda_ww/2 by Gauss-Hermite."""
    z = a + B * u; tau = np.sqrt(t2)
    sh = np.broadcast(z, tau).shape
    L = np.zeros(sh); La = np.zeros(sh); L3 = np.zeros(sh)
    for x, w in zip(GX, GW):
        P = expit(z + tau * x); q = P * (1 - P)
        L += w * P; La += w * q; L3 += w * 0.5 * q * (1 - 2 * P)
    return L, La, L3

ZX = np.linspace(-9, 9, 2401); ZW = np.exp(-ZX ** 2 / 2); ZW /= ZW.sum()

def a_of(d, psi, s1, s2, s3):
    """Pair functional E_sigma{F_u (y - xi)} from the definition, numeric derivative of psi; h_sigma is not used."""
    a, B, t2 = d.theta(s1, s2, s3); xi = s1 + np.sqrt(s3) * ZX
    L, La, _ = Lders(xi, a, B, t2); v = L * (1 - L); hh = 1e-4 * np.sqrt(s3)
    dps = (psi(xi + hh) - psi(xi - hh)) / (2 * hh)
    return np.sum(ZW * (dps * v + psi(xi) * B * La * (xi - L)))

def fit_all(d, A, y):
    R = np.linalg.inv(A.T @ A + np.diag(d.D)); b = R @ (A.T @ y)
    return R, b
