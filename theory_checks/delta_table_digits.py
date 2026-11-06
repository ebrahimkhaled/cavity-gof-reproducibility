# Table delta under the definition of e_g (cells normalised by w_{g,n}), with more digits. Prints only.
import sys, json, os
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
import check_pair as cp
from scipy.stats import norm
G = cp.G
for line in open(os.path.join(W, "pair_runs.txt"), encoding="utf-8-sig", errors="ignore"):
    line = line.strip()
    if not line.startswith("{"):
        continue
    r = json.loads(line)
    d = cp.D2(r["n"], r["p"], r["lam"] * r["n"], ar=r["ar"], s=r["s"])
    sig = np.array(r["sigma"]); m, t = sig[0], np.sqrt(sig[2])
    cuts = m + t * norm.ppf(np.arange(1, G) / G); hb = 0.03 * 1.349 * t
    xi = m + t * cp.XG
    a0, B0, tau2 = cp.theta_of(*sig, d.m0, d.s2)
    L = cp.Lam(xi, a0, B0, tau2); v = L * (1 - L)
    cells, scores = cp.make_tests(d, sig, cuts, hb)
    wn = np.sqrt((cells(xi) * (v * cp.WG)[:, None]).sum(0))
    avec, Ivec, P, frac, resn, kn, th, Fev = cp.theory(d, sig, r["gamma"], cuts, hb, wn)
    Om = P @ Ivec @ P.T
    Dp = r["gamma"] * (P @ avec)[:, None] * (P @ avec)[None, :]
    print(r["n"], r["p"], r["lam"], r["s"], "trOm %.4f trD %.5g ratio %.5g share %.3g gamma %.4f tau2 %.3f" % (np.trace(Om), np.trace(Dp), np.trace(Dp) / np.trace(Om), frac, r["gamma"], tau2), flush=True)
