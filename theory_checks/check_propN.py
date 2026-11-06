# Check 3: Proposition on the vector N. c0'N and beta'St N should be O(n^{-1/2}), b'St N - q1 too; E||N||^2 bounded.
# Also the cancellation 0.5*Fbar_uu + Fbar_3 -> 0 and the sizes of Q - q, g - gbar.
import sys, numpy as np
from common import *
from check_expansion import psi_factory, Fam

def one(d, rng):
    n, p = d.n, d.p
    A, y = d.draw(rng, n - 1)
    R = np.linalg.inv(A.T @ A + np.diag(d.D))
    b = R @ (A.T @ y); Rho = A @ R
    ht = np.einsum('ij,ij->i', Rho, A); h = ht / (1 - ht)
    e = (y - A @ b) / (1 - ht)
    Bk = b[None, :] - Rho * e[:, None]; u = y - e
    s1, s2, s3 = d.sigma(Bk)
    t = np.sqrt(np.mean(s3))
    fam = Fam(psi_factory(t), d)
    Fu, Fuu, F1, F2, F3 = fam.all(u, y, s1, s2, s3, 1e-3 * t, 1e-4 * t * t)
    St = np.eye(p + 1); St[0, 0] = 0
    vbar = F1.sum() * np.eye(p + 1)[0] + F2.sum() * np.r_[0, d.b0] + 2 * F3.sum() * (St @ b)
    N = (Rho * (Fu * (1 + h))[:, None]).sum(0) + R @ vbar
    q1 = np.sum(Fu * (1 + h) * e * np.einsum('ij,ij->i', Rho @ St, Rho))
    return np.array([np.r_[1, np.zeros(p)] @ N, np.r_[0, d.b0] @ N, (St @ b) @ N - q1, q1, N @ St @ N, (0.5 * Fuu.mean() + F3.mean())])

for n in [200, 400, 800, 1600, 3200]:
    d = Design(n, n // 4, 0.25 * n)
    rng = np.random.default_rng(n)
    reps = 400 if n <= 800 else (150 if n == 1600 else 40)
    X = np.array([one(d, rng) for _ in range(reps)])
    rms = np.sqrt((X ** 2).mean(0))
    print(n, "sqrt(n)*rms: c'N %.3f beta'N %.3f (b'N-q1) %.3f | q1 rms %.3f  E N'SN %.3f | sqrt(n)*rms(.5Fuu+F3) %.3f" % (np.sqrt(n) * rms[0], np.sqrt(n) * rms[1], np.sqrt(n) * rms[2], rms[3], X[:, 4].mean(), np.sqrt(n) * rms[5]))
