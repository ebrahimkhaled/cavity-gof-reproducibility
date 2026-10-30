import sys, glob
import numpy as np
files = sorted(glob.glob(sys.argv[1]))
d = {k: np.concatenate([np.load(f)[k] for f in files]) for k in ["Z", "U", "F", "W0", "a", "ntau"]}
G = 10
Z = d["Z"]; R = len(Z)
U = d["U"].mean(0); F = d["F"].mean(0); W0 = d["W0"].mean(0); a = d["a"].mean(0); ntau = d["ntau"].mean()
sd = np.sqrt(np.diag(W0))                  # standardise as e_g does
T, S = Z[:, :G], Z[:, G:]
K = U @ np.linalg.solve(F, np.eye(3))
Rv = (T - S @ K.T) / sd
Om0 = (W0 - K @ U.T) / np.outer(sd, sd)
aR = (a[:G] - K @ a[G:]) / sd
aT = a[:G] / sd
C = np.cov(Rv.T)
CT = np.cov((T / sd).T); W0s = W0 / np.outer(sd, sd)
print(f"files={len(files)} reps={R}  n*tau={ntau:.3f}")
print("UNPROJECTED cells:  trace excess  emp=%.3f  pred=%.3f" % (np.trace(CT - W0s), ntau * aT @ aT))
print("PROJECTED  R:       trace excess  emp=%.3f  pred=%.3f   (trace Omega0=%.3f)" % (np.trace(C - Om0), ntau * aR @ aR, np.trace(Om0)))
u = aR / np.linalg.norm(aR)
print("  directional excess along a_R: emp=%.3f  pred=%.3f" % (u @ (C - Om0) @ u, ntau * aR @ aR))
# sampling SE of directional variance: var(u'R)*sqrt(2/R)
vr = u @ C @ u
print("  se of directional variance ~ %.3f" % (vr * np.sqrt(2 / R)))
ev = np.linalg.eigvalsh(np.linalg.pinv(np.linalg.cholesky(Om0 + 1e-9 * np.eye(G))) @ C @ np.linalg.pinv(np.linalg.cholesky(Om0 + 1e-9 * np.eye(G))).T)
print("  eigenvalues of Om0^-1/2 Cov(R) Om0^-T/2:", np.round(ev, 3))
print("  mean of R:", np.round(Rv.mean(0), 3))
print("  share of |a_T|^2 removed by projection: %.3f" % (1 - (aR @ aR) / (aT @ aT)))
