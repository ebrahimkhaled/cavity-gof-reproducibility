# Scratch checks for the merged section. Nothing here is cited as a file in the text.
import sys, json, os
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
sys.path.insert(0, os.path.join(W, "audit"))

def delta_table():
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
        C = cells(xi)
        w_sq = np.sqrt(((C ** 2) * (v * cp.WG)[:, None]).sum(0))   # as run: E phi^2 v
        w_th = np.sqrt((C * (v * cp.WG)[:, None]).sum(0))          # as in the definition of e_g: E phi v
        out = []
        for wn in (w_sq, w_th):
            avec, Ivec, P, frac, resn, kn, th, Fev = cp.theory(d, sig, r["gamma"], cuts, hb, wn)
            Om = P @ Ivec @ P.T
            Dp = r["gamma"] * (P @ avec)[:, None] * (P @ avec)[None, :]
            out += [np.trace(Om), np.trace(Dp), np.trace(Dp) / np.trace(Om)]
        print(f"n={r['n']} p={r['p']} lam={r['lam']} ar={r['ar']} s={r['s']}: as run trOm {out[0]:.3f} trD {out[1]:.3g} ratio {out[2]:.3g} | stored trD {r['trace_proj_pred']:.3g} | definition trOm {out[3]:.3f} trD {out[4]:.3g} ratio {out[5]:.3g}", flush=True)

def gamma_centred():
    from acommon import Des, fit_all
    for (n, p, lamn, ar, mu, reps) in [(400, 100, 4.16, .7, 0., 40), (800, 200, 4.16, .7, 0., 20), (1600, 400, 4.16, .7, 0., 8), (400, 100, .25, 0., 0., 40), (400, 100, 4.16, .7, 1., 40), (300, 600, 4.16, .7, 0., 30)]:
        d = Des(n, p, lamn, eps0=0., sig=2.0, b00=0., mu=mu, ar=ar); r = np.random.default_rng(1); o = []
        for _ in range(reps):
            A_, y, _ = d.draw(r, n); R, b = fit_all(d, A_, y); H = A_ @ R @ A_.T
            z1 = np.trace(H) / n
            gt = (np.trace(H @ H) - np.sum(np.diag(H) ** 2)) / n / (1 - z1) ** 2
            Hc = H - 1.0 / n
            gc = (np.trace(Hc @ Hc) - np.sum(np.diag(Hc) ** 2)) / n / (1 - z1) ** 2
            RS = R @ d.St; gb = np.trace(RS) + d.c @ R @ d.c; om = np.trace(RS @ RS); gh = n * om / (1 + gb) ** 2
            o.append([gt, gc, gh, np.abs(H.sum(1) - 1).max()])
        m = np.array(o).mean(0)
        print(f"n={n} p={p} lam/n={lamn} ar={ar} mu={mu}: tilde {m[0]:.5f} centred {m[1]:.5f} gamma_hat {m[2]:.5f} gaps {m[0]-m[2]:.5f} {m[1]-m[2]:.5f}  max|H1-1| {m[3]:.1e}", flush=True)

if __name__ == "__main__":
    delta_table()
    gamma_centred()
