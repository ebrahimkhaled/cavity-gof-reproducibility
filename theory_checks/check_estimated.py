# Check 4: the "second excess". Score covariance with scores at the FIXED theta-bar against scores at the estimated theta(b_hat).
# Usage: python check_estimated.py n p lam reps
import sys, numpy as np
from multiprocessing import Pool
from common import *
from check_pair import D2, fit, pilot

def sc(u, th):
    a, B, t2 = th
    L = Lam(u, a, B, t2); v = L * (1 - L); La = Lam_a(u, a, B, t2)
    return np.stack([La / v, u * La / v, Lam_tau2(u, a, B, t2) / v], -1)

def worker(args):
    d, seed, reps, sig = args
    th0 = theta_of(*sig, d.m0, d.s2)
    rng = np.random.default_rng(seed)
    T0 = np.empty((reps, 3)); T1 = np.empty((reps, 3)); I0 = np.zeros((3, 3)); I1 = np.zeros((3, 3))
    for r in range(reps):
        A, y = d.draw(rng)
        R, b, xi, s1, s2, s3, ht = fit(d, A, y)
        a, B, t2 = theta_of(s1, s2, s3, d.m0, d.s2)
        L = Lam(xi, a, B, t2); v = L * (1 - L)
        thh = theta_of(b[0], b @ d.Stb0, b @ d.St @ b, d.m0, d.s2)
        P0 = sc(xi, th0); P1 = sc(xi, thh)
        T0[r] = P0.T @ (y - L) / np.sqrt(d.n); T1[r] = P1.T @ (y - L) / np.sqrt(d.n)
        I0 += (P0 * v[:, None]).T @ P0 / d.n; I1 += (P1 * v[:, None]).T @ P1 / d.n
    return T0, T1, I0, I1

if __name__ == '__main__':
    n, p, lam, reps = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    d = D2(n, p, lam * n)
    W = 20
    with Pool(W) as pool:
        pm = np.vstack(pool.map(pilot, [(d, 7000 + i, max(4, 400 // W)) for i in range(W)])).mean(0)
        rr = pool.map(worker, [(d, 19000 + i, reps // W, pm[:3]) for i in range(W)])
    T0 = np.vstack([r[0] for r in rr]); T1 = np.vstack([r[1] for r in rr])
    I0 = sum(r[2] for r in rr) / len(T0); I1 = sum(r[3] for r in rr) / len(T0)
    for nm, T, I in (("fixed theta-bar", T0, I0), ("estimated theta(b_hat)", T1, I1)):
        C = T.T @ T / len(T)
        Li = np.linalg.inv(np.linalg.cholesky(I))
        print(n, nm, "ratios", np.round(np.linalg.eigvalsh(Li @ C @ Li.T), 3))
