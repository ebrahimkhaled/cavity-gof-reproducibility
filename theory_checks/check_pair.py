# Check 2: Monte Carlo pair term K_n of T^0 for FIXED test functions (cells and scores at a fixed theta-bar),
# against gamma_n a a'; then the projected pair term. Usage: python check_pair.py n p lam_over_n reps [ar] [s]
import sys, numpy as np, time, json
from multiprocessing import Pool
from numpy.polynomial.hermite_e import hermegauss
from common import *

G = 10
XG = np.linspace(-8, 8, 1601); WG = np.exp(-XG ** 2 / 2); WG = WG / WG.sum()

class D2(Design):
    def __init__(self, n, p, lam, ar=0.0, s=2.0):
        super().__init__(n, p, lam, eps0=0.0, s=s)
        self.Sig = ar ** np.abs(np.subtract.outer(np.arange(p), np.arange(p))) if ar > 0 else np.eye(p)
        self.Lc = np.linalg.cholesky(self.Sig)
        b = np.tile([.35, -.3, .25, -.2, .15], p // 5 + 1)[:p] if ar > 0 else np.ones(p)
        self.b0 = b * s / np.sqrt(b @ self.Sig @ b)
        self.St = np.zeros((p + 1, p + 1)); self.St[1:, 1:] = self.Sig
        self.Stb0 = np.r_[0, self.Sig @ self.b0]
    def draw(self, rng, m=None):
        m = self.n if m is None else m
        X = rng.standard_normal((m, self.p)) @ self.Lc.T
        A = np.c_[np.ones(m), X]
        y = (rng.random(m) < expit(self.beta00 + X @ self.b0)).astype(float)
        return A, y

def fit(d, A, y):
    R = np.linalg.inv(A.T @ A + np.diag(d.D))
    b = R @ (A.T @ y)
    Rho = A @ R
    ht = np.einsum('ij,ij->i', Rho, A)
    e = (y - A @ b) / (1 - ht)                  # = y_i - xi_i
    xi = y - e
    Bl = b[None, :] - Rho * e[:, None]          # leave-one-out fits
    s1 = Bl[:, 0]; s2 = Bl @ d.Stb0
    s3 = np.einsum('ij,ij->i', Bl @ d.St, Bl)
    return R, b, xi, s1, s2, s3, ht

def pilot(args):
    d, seed, reps = args
    rng = np.random.default_rng(seed); out = []
    for _ in range(reps):
        A, y = d.draw(rng)
        R, b, xi, s1, s2, s3, ht = fit(d, A, y)
        RS = R @ d.St
        gb = np.trace(RS) + R[0, 0]; om = np.trace(RS @ RS)
        H2 = np.sum((A @ R @ A.T) ** 2); Hd = np.sum(ht ** 2); Ht = ht.sum()
        n = d.n
        ghat = (H2 - Hd) / (n * (1 - Ht / n) ** 2)
        out.append([s1.mean(), s2.mean(), s3.mean(), n * om / (1 + gb) ** 2, gb, om * n, ghat])
    return np.array(out)

def make_tests(d, sig, cuts, hb):
    a0, B0, t20 = theta_of(*sig, d.m0, d.s2)
    def cells(u):
        pieces = [np.ones_like(u)]
        c = np.r_[-np.inf, cuts, np.inf]
        out = []
        for g in range(G):
            lo = expit((u - c[g]) / hb) if np.isfinite(c[g]) else 1.0
            hi = expit((u - c[g + 1]) / hb) if np.isfinite(c[g + 1]) else 0.0
            out.append(lo - hi + 0 * u)
        return np.stack(out, -1)
    def scores(u):
        L = Lam(u, a0, B0, t20); v = L * (1 - L)
        La = Lam_a(u, a0, B0, t20)
        return np.stack([La / v, u * La / v, Lam_tau2(u, a0, B0, t20) / v], -1)
    return cells, scores

def worker(args):
    d, seed, reps, sig, cuts, hb, wn = args
    cells, scores = make_tests(d, sig, cuts, hb)
    rng = np.random.default_rng(seed)
    Ts = np.empty((reps, G + 3)); Ib = np.zeros((G + 3, G + 3))
    for r in range(reps):
        A, y = d.draw(rng)
        R, b, xi, s1, s2, s3, ht = fit(d, A, y)
        a, B, t2 = theta_of(s1, s2, s3, d.m0, d.s2)
        L = Lam(xi, a, B, t2)
        Psi = np.c_[cells(xi) / wn, scores(xi)]
        Ts[r] = Psi.T @ (y - L) / np.sqrt(d.n)
        Ib += (Psi * (L * (1 - L))[:, None]).T @ Psi / d.n
    return Ts, Ib

def theory(d, sig, gam, cuts, hb, wn):
    m, c, t2 = sig[0], sig[1], sig[2]
    a0, B0, tau2 = theta_of(m, c, t2, d.m0, d.s2)
    xi = m + np.sqrt(t2) * XG
    L = Lam(xi, a0, B0, tau2); v = L * (1 - L)
    La = Lam_a(xi, a0, B0, tau2); Lp = B0 * La
    h = v * (xi - m) / t2 + Lp * (xi + L - 1)
    cells, scores = make_tests(d, sig, cuts, hb)
    Psi = np.c_[cells(xi) / wn, scores(xi)]
    avec = (Psi * (h * WG)[:, None]).sum(0)
    Ivec = (Psi * (v * WG)[:, None]).T @ Psi
    U = Ivec[:G, G:]; F = Ivec[G:, G:]
    P = np.c_[np.eye(G), -U @ np.linalg.inv(F)]
    # projection of h/v on the span of the scores in L2(v)
    S = np.stack([La / v, xi * La / v, Lam_tau2(xi, a0, B0, tau2) / v], -1)
    k = h / v
    coef = np.linalg.solve((S * (v * WG)[:, None]).T @ S, (S * (v * WG)[:, None]).T @ k)
    res = k - S @ coef
    frac = np.sum(WG * v * res ** 2) / np.sum(WG * v * k ** 2)
    return avec, Ivec, P, frac, np.sum(WG * v * res ** 2), np.sum(WG * v * k ** 2), (a0, B0, tau2), np.linalg.eigvalsh(F)

if __name__ == '__main__':
    n, p, lam, reps = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    ar = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    s = float(sys.argv[6]) if len(sys.argv) > 6 else 2.0
    d = D2(n, p, lam * n, ar=ar, s=s)
    W = 20
    t0 = time.time()
    with Pool(W) as pool:
        pil = np.vstack(pool.map(pilot, [(d, 7000 + i, max(4, 400 // W)) for i in range(W)]))
        pm = pil.mean(0); sig = pm[:3]; gam = pm[3]; ghat = pm[6]
        m, t = sig[0], np.sqrt(sig[2])
        from scipy.stats import norm
        cuts = m + t * norm.ppf(np.arange(1, G) / G)
        hb = 0.03 * 1.349 * t
        # w_g from quadrature
        xi = m + t * XG
        a0, B0, tau2 = theta_of(*sig, d.m0, d.s2)
        Lq = Lam(xi, a0, B0, tau2); vq = Lq * (1 - Lq)
        cells, scores = make_tests(d, sig, cuts, hb)
        wn = np.sqrt(((cells(xi) ** 2) * (vq * WG)[:, None]).sum(0))
        avec, Ivec, P, frac, resn, kn, th, Fev = theory(d, sig, gam, cuts, hb, wn)
        rr = pool.map(worker, [(d, 9000 + i, reps // W, sig, cuts, hb, wn) for i in range(W)])
    Ts = np.vstack([r[0] for r in rr]); Ib = sum(r[1] for r in rr) / len(Ts)
    Cv = (Ts.T @ Ts) / len(Ts)                     # E f = 0 exactly, so no centring
    K = Cv - Ib
    pred = gam * np.outer(avec, avec)
    # bootstrap-free SE of trace pieces through batches
    nb = 20; bs = len(Ts) // nb
    def stats(TT, I):
        C = TT.T @ TT / len(TT); Kx = C - I
        return np.array([np.trace(Kx[:G, :G]), np.trace(P @ Kx @ P.T), avec @ Kx @ avec / (avec @ avec)])
    B = np.array([stats(Ts[i * bs:(i + 1) * bs], Ib) for i in range(nb)])
    se = B.std(0, ddof=1) / np.sqrt(nb)
    Om = P @ Ivec @ P.T
    Dp = gam * (P @ avec)[:, None] * (P @ avec)[None, :]
    print(json.dumps(dict(n=n, p=p, lam=lam, ar=ar, s=s, reps=len(Ts), sigma=list(map(float, sig)), theta=list(map(float, th)), gamma=float(gam), gamma_hat_data=float(ghat),
        trace_cells_MC=float(np.trace(K[:G, :G])), se_tc=float(se[0]), trace_cells_pred=float(np.trace(pred[:G, :G])),
        trace_proj_MC=float(np.trace(P @ K @ P.T)), se_tp=float(se[1]), trace_proj_pred=float(np.trace(Dp)), trace_Omega=float(np.trace(Om)),
        along_a_MC=float(avec @ K @ avec / (avec @ avec)), se_al=float(se[2]), along_a_pred=float(gam * (avec @ avec)),
        frac_unexplained=float(frac), resid_norm=float(resn), k_norm=float(kn), F_eigs=list(map(float, Fev)), time=time.time() - t0)))
    Ssc = Ivec[G:, G:]
    Lc = np.linalg.cholesky(Ib[G:, G:]); Li = np.linalg.inv(Lc)
    print("score cov ratios MC", np.round(np.linalg.eigvalsh(Li @ Cv[G:, G:] @ Li.T), 4))
    print("score cov ratios pred", np.round(np.linalg.eigvalsh(Li @ (Ib[G:, G:] + pred[G:, G:]) @ Li.T), 4))
    print("max abs K-pred", np.abs(K - pred).max(), " max abs pred", np.abs(pred).max())
    np.savez(f"pair_n{n}_p{p}_lam{lam}_ar{ar}_s{s}.npz", K=K, pred=pred, a=avec, P=P, Ib=Ib, Ivec=Ivec)
