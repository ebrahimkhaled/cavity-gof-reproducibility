"""Projected check: R = T - U F^{-1} S (G smoothed cells, MLE score of the 3-parameter curve),
residuals at the exact cavity-law parameters theta_(-i). Compare Cov(R) with the iid projection
Omega0 = E[phi phi' v] - U F^{-1} U' and with the predicted rank-one excess n*tau * a_R a_R'."""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"; os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np

seed, lam, reps = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
outfile = sys.argv[4]
NN, PP = int(sys.argv[5]), int(sys.argv[6])
rng = np.random.default_rng(seed)
n, p, s, b00, G, hb = NN, PP, 2.0, 0.0, 10, 0.1
beta0 = np.full(p, s / np.sqrt(p))
gx, gw = np.polynomial.hermite_e.hermegauss(40); gw = gw / gw.sum()
expit = lambda u: 1.0 / (1.0 + np.exp(-u))
D = np.ones(p + 1); D[0] = 0.0
from statistics import NormalDist
cuts = None

def phis(x):
    up = np.concatenate([np.ones((len(x), 1)), expit((x[:, None] - cuts[None, :]) / hb)], 1)
    lo = np.concatenate([expit((x[:, None] - cuts[None, :]) / hb), np.zeros((len(x), 1))], 1)
    return up - lo   # n x G, smoothed indicators of fixed intervals

def lam_parts(x, a, b, tau):
    U = a[:, None] + b[:, None] * x[:, None] + tau[:, None] * gx[None, :]
    P = expit(U); d1 = P * (1 - P)
    L = P @ gw; da = d1 @ gw; db = x * da; dt = (d1 * gx[None, :]) @ gw
    return L, np.column_stack([da, db, dt])

def fvals(x, y, a, b, tau):
    L, Ld = lam_parts(x, a, b, tau)
    v = L * (1 - L)
    psi = np.column_stack([phis(x), Ld / v[:, None]])      # G cell functions + 3 score functions
    return psi * (y - L)[:, None], psi, L, Ld, v

acc = []
def cavity(rng):
    X = rng.standard_normal((n, p)); Xt = np.column_stack([np.ones(n), X])
    y = (rng.random(n) < expit(b00 + X @ beta0)).astype(float)
    bh = np.zeros(p + 1)
    for it in range(60):
        eta = Xt @ bh; pi = expit(eta); w = pi * (1 - pi)
        M = (Xt * w[:, None]).T @ Xt + lam * np.diag(D)
        step = np.linalg.solve(M, Xt.T @ (y - pi) - lam * D * bh); bh += step
        if np.max(np.abs(step)) < 1e-10: break
    eta = Xt @ bh; pi = expit(eta); w = pi * (1 - pi)
    M = (Xt * w[:, None]).T @ Xt + lam * np.diag(D); Minv = np.linalg.inv(M)
    MX = Xt @ Minv; q = np.einsum("ij,ij->i", MX, Xt); coef = (y - pi) / (1 - w * q)
    Bm = bh[None, :] - MX * coef[:, None]; xi = eta - q * coef
    sl = Bm[:, 1:]; t2 = np.sum(sl ** 2, 1); c = sl @ beta0
    b = c / t2; tau = np.sqrt(np.maximum(s ** 2 - c ** 2 / t2, 1e-8)); a = b00 - b * Bm[:, 0]
    return xi, y, a, b, tau, pi, Minv

xi0 = np.concatenate([cavity(np.random.default_rng(777 + j))[0] for j in range(5)])
cuts = np.quantile(xi0, np.arange(1, G) / G); hb = 0.03 * (np.quantile(xi0, 0.75) - np.quantile(xi0, 0.25))
for r in range(reps):
    xi, y, a, b, tau, pi, Minv = cavity(rng)
    f, psi, L, Ld, v = fvals(xi, y, a, b, tau)
    Z = f.sum(0) / np.sqrt(n)                     # (T_1..T_G, S_1..S_3)
    Uh = psi[:, :G].T @ Ld / n; Fh = (Ld / v[:, None]).T @ Ld / n
    W0 = (psi[:, :G] * v[:, None]).T @ psi[:, :G] / n
    dl = 1e-4
    df = (fvals(xi + dl, y, a, b, tau)[0] - fvals(xi - dl, y, a, b, tau)[0]) / (2 * dl)
    avec = (df * (pi - y)[:, None]).mean(0)
    S = np.diag(D); ntau = n * np.trace(Minv @ S @ Minv @ S)
    acc.append((Z, Uh, Fh, W0, avec, ntau))

np.savez(outfile, Z=np.array([x[0] for x in acc]), U=np.array([x[1] for x in acc]),
         F=np.array([x[2] for x in acc]), W0=np.array([x[3] for x in acc]),
         a=np.array([x[4] for x in acc]), ntau=np.array([x[5] for x in acc]))
