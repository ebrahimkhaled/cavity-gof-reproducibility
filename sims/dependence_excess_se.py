"""Projected and unprojected trace excess with a batch-bootstrap standard error (resample whole batch files)."""
import glob, sys
import numpy as np
G = 10

def load(pattern):
    fs = sorted(glob.glob(pattern))
    return [{k: np.load(f)[k] for k in ["Z", "U", "F", "W0"]} for f in fs]

def excess(batches):
    d = {k: np.concatenate([b[k] for b in batches]) for k in ["Z", "U", "F", "W0"]}
    U, F, W0 = d["U"].mean(0), d["F"].mean(0), d["W0"].mean(0)
    sd = np.sqrt(np.diag(W0)); T, S = d["Z"][:, :G], d["Z"][:, G:]
    K = U @ np.linalg.inv(F)
    R = (T - S @ K.T) / sd
    Om0 = (W0 - K @ U.T) / np.outer(sd, sd)
    unproj = np.trace(np.cov((T / sd).T)) - G
    proj = np.trace(np.cov(R.T)) - np.trace(Om0)
    return unproj, proj, len(d["Z"])

rng = np.random.default_rng(1)
for label, pat in [("400", "s400_*.npz"), ("800", "s800_*.npz"), ("1600", "mc_big/t1600_*.npz"),
                   ("3200", "mc_big/t3200_*.npz"), ("6400", "mc_big/t6400_*.npz")]:
    b = load(pat); u, p, R = excess(b)
    boot = np.array([excess([b[j] for j in rng.integers(0, len(b), len(b))])[:2] for _ in range(200)])
    print(f"n={label:>5} R={R:6d}  unprojected {u:.3f} ({boot[:,0].std():.3f})  projected {p:.3f} ({boot[:,1].std():.3f})")
