# End-to-end check of Theorem null: fit the recalibration root, form e_g, compare with Omega_n + Delta_n.
# Also: does a root exist? (A4); size of the near-null direction; second-excess scaling.
# Usage: python audit_e.py tag n p lamn eps0 sig b00 mu ar reps G
import sys, json, time
from multiprocessing import Pool
from acommon import *

def cells(cuts, hb):
    def phi(u):
        out = []
        ext = [-np.inf] + list(cuts) + [np.inf]
        for lo, hi in zip(ext[:-1], ext[1:]):
            A = 1.0 if lo == -np.inf else expit((u - lo) / hb)
            B = 0.0 if hi == np.inf else expit((u - hi) / hb)
            out.append(A - B)
        return np.stack(out, -1)
    return phi

def scores(u, th):
    a, B, t2 = th
    L, La, L3 = Lders(u, a, B, t2); v = L * (1 - L)
    return L, v, np.stack([La, u * La, L3], -1)      # Lambda, v, Lambda-dot

from numpy.polynomial.hermite_e import hermegauss as _hg
_X, _W = _hg(64); _W = _W / _W.sum()
def _LD(z, a, b, t2):
    w = a + b * z; tau = np.sqrt(t2)
    P = expit(w[:, None] + tau * _X[None, :])
    return P @ _W, (P * (1 - P)) @ _W

def fit_root(xi, y, th0, tmax=30.0):
    """Profile likelihood: for each tau2 solve the (a, b) score equations by Newton, then golden search on tau in [0, sqrt(tmax)].
    Interior root = maximiser strictly inside; boundary = tau2 = 0 (no root of the score equation there); runaway = upper end."""
    m, sd = xi.mean(), xi.std(); z = (xi - m) / sd
    st = [th0[0] + th0[1] * m, th0[1] * sd, th0[2]]
    def prof(tau):
        t2 = tau * tau; sc_ = np.sqrt((1 + .35 * t2) / (1 + .35 * st[2]))
        a, b = st[0] * sc_, st[1] * sc_
        for _ in range(25):
            L, La = _LD(z, a, b, t2); v = np.maximum(L * (1 - L), 1e-12)
            r = (y - L) * La / v; g = np.array([r.sum(), (r * z).sum()])
            q = La * La / v; I = np.array([[q.sum(), (q * z).sum()], [(q * z).sum(), (q * z * z).sum()]])
            dlt = np.linalg.solve(I + 1e-9 * np.eye(2), g)
            nr = np.abs(dlt).max()
            if nr > 1.0:
                dlt = dlt / nr
            a += dlt[0]; b += dlt[1]
            if np.abs(dlt).max() < 1e-7:
                break
        st[0], st[1], st[2] = a, b, t2
        L, _ = _LD(z, a, b, t2); L = np.clip(L, 1e-13, 1 - 1e-13)
        return np.sum(y * np.log(L) + (1 - y) * np.log(1 - L)), a, b
    gr = (5 ** .5 - 1) / 2; lo, hi = 0.0, tmax ** .5
    c, dd = hi - gr * (hi - lo), lo + gr * (hi - lo)
    fc, fd = prof(c)[0], prof(dd)[0]
    for _ in range(22):
        if fc >= fd:
            hi, dd, fd = dd, c, fc; c = hi - gr * (hi - lo); fc = prof(c)[0]
        else:
            lo, c, fc = c, dd, fd; dd = lo + gr * (hi - lo); fd = prof(dd)[0]
    tau = (lo + hi) / 2
    f, a, b = prof(tau)
    f0, a0, b0 = prof(0.0); f1, a1, b1 = prof(tmax ** .5)
    if f0 >= f:
        tau, f, a, b = 0.0, f0, a0, b0
    if f1 >= f:
        tau, f, a, b = tmax ** .5, f1, a1, b1
    t2 = tau * tau
    interior = (t2 > 1e-3) and (t2 < tmax * .99)
    return (a - b * m / sd, b / sd, t2), interior, True

def worker(a):
    cfg, seed, reps = a
    d = Des(*cfg['des'], **cfg['kw']); r = np.random.default_rng(seed)
    phi = cells(cfg['cuts'], cfg['hb']); thb = cfg['thb']; G = len(cfg['cuts']) + 1
    w = np.array(cfg['w']); out = []
    Tp = np.array(cfg['P'])
    for _ in range(reps):
        A_, y, _ = d.draw(r, d.n)
        R, b = fit_all(d, A_, y); Rho = A_ @ R
        ht = np.einsum('ij,ij->i', Rho, A_); ee = (y - A_ @ b) / (1 - ht); xi = y - ee
        # b_(-i) summaries
        s1 = d.c @ b - ee * (Rho @ d.c); s2 = d.Stb @ b - ee * (Rho @ d.Stb)
        Sb = d.St @ b
        rr = np.einsum('ij,ij->i', Rho @ d.St, Rho)
        s3 = b @ Sb - 2 * ee * (Rho @ Sb) + ee ** 2 * rr
        th_i = d.theta(s1, s2, s3)
        Li = Lam(xi, *th_i)
        Ph = phi(xi)
        Lb, vb, Db = scores(xi, thb); Sbar = Db / vb[:, None]
        sn = np.sqrt(d.n)
        T0 = np.r_[Ph.T @ (y - Li) / sn / w, Sbar.T @ (y - Li) / sn]     # T-tilde^0
        pe = Tp @ T0                                                     # projected, deterministic P-bar
        # estimated scores at theta(b_f) (second excess)
        thf = d.theta(d.c @ b, d.Stb @ b, b @ Sb)
        Lf, vf, Df = scores(xi, thf)
        Tf = (Df / vf[:, None]).T @ (y - Li) / sn
        # root
        th, inter, ok = fit_root(xi, y, thb)
        Lh, vh, _ = scores(xi, th)
        eg = (Ph.T @ (y - Lh)) / np.sqrt(Ph.T @ vh)
        out.append(np.r_[eg, pe, T0[G:], Tf, th, float(inter), float(ok)])
    return np.array(out)

if __name__ == '__main__':
    tag, n, p = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    lamn, eps0, sig, b00, mu, ar = map(float, sys.argv[4:10]); reps, G = int(sys.argv[10]), int(sys.argv[11])
    kw = dict(eps0=eps0, sig=sig, b00=b00, mu=mu, ar=ar)
    d = Des(n, p, lamn, **kw); r = np.random.default_rng(3); bs = []; gs = []
    for _ in range(40):
        A_, y_, _ = d.draw(r, n - 1); R, b = fit_all(d, A_, y_); bs.append([d.c @ b, d.Stb @ b, b @ d.St @ b])
        RS = R @ d.St; gb = np.trace(RS) + d.c @ R @ d.c; gs.append(n * np.trace(RS @ RS) / (1 + gb) ** 2)
    s = np.mean(bs, 0); gam = float(np.mean(gs)); thb = tuple(float(x) for x in d.theta(*s)); t = np.sqrt(s[2])
    from statistics import NormalDist
    cuts = [s[0] + t * NormalDist().inv_cdf(g / G) for g in range(1, G)]; hb = 0.05 * t
    phi = cells(cuts, hb)
    xi = s[0] + t * ZX
    L, v, D = scores(xi, thb); S = D / v[:, None]; Ph = phi(xi); W = ZW * v
    F = (S * W[:, None]).T @ S
    w2 = Ph.T @ W; w = np.sqrt(w2)
    U = ((Ph * W[:, None]).T @ S) / w[:, None]
    Pm = np.c_[np.eye(G), -U @ np.linalg.inv(F)]
    Coord = np.c_[Ph / w, S]
    C = (Coord * W[:, None]).T @ Coord
    h = v * (xi - s[0]) / t ** 2 + thb[1] * D[:, 0] * (xi + L - 1)
    av = (Coord * ZW[:, None]).T @ h
    Om = Pm @ C @ Pm.T; pa = Pm @ av; De = gam * np.outer(pa, pa); Vfull = C + gam * np.outer(av, av)
    Fi = np.linalg.inv(F)
    # sandwich: var of root = F^-1 (Var T0(s) restricted) F^-1 + var Z; use information-only value
    info_sd_tau2 = float(np.sqrt(Fi[2, 2] / n))
    cfg = dict(des=(n, p, lamn), kw=kw, cuts=cuts, hb=hb, thb=thb, w=list(w), P=Pm.tolist())
    Wk = 20; t0 = time.time()
    with Pool(Wk) as pool:
        X = np.vstack(pool.map(worker, [(cfg, 77000 + 13 * n + i, reps // Wk) for i in range(Wk)]))
    eg, pe, Ts, Tf, th, inter, ok = X[:, :G], X[:, G:2 * G], X[:, 2 * G:2 * G + 3], X[:, 2 * G + 3:2 * G + 6], X[:, 2 * G + 6:2 * G + 9], X[:, -2], X[:, -1]
    nn = len(X)
    def tr(Y):
        q = (Y ** 2).sum(1)
        return float(q.mean()), float(q.std(ddof=1) / np.sqrt(len(Y)))
    # eigen ratios of score covariance
    Li = np.linalg.inv(np.linalg.cholesky(F))
    r0 = np.linalg.eigvalsh(Li @ (Ts.T @ Ts / nn) @ Li.T); r1 = np.linalg.eigvalsh(Li @ (Tf.T @ Tf / nn) @ Li.T)
    I = inter > 0
    res = dict(tag=tag, n=n, p=p, reps=nn, G=G, theta_bar=[round(x, 3) for x in thb], gamma=round(gam, 4), F_eigs=[float(f'{x:.3g}') for x in np.linalg.eigvalsh(F)],
               info_sd_tau2=round(info_sd_tau2, 3), tr_Om=round(float(np.trace(Om)), 4), tr_De=float(f'{np.trace(De):.3g}'),
               tr_pe=[round(x, 4) for x in tr(pe)], tr_e_all=[round(x, 4) for x in tr(eg)], tr_e_interior=[round(x, 4) for x in tr(eg[I])],
               frac_interior=round(float(I.mean()), 4), frac_tau0=round(float((th[:, 2] <= 1e-3).mean()), 4), frac_top=round(float((th[:, 2] >= 29.7).mean()), 4),
               tau2_q=[round(float(x), 3) for x in np.quantile(th[:, 2], [.05, .25, .5, .75, .95])], sd_tau2_int=round(float(th[I, 2].std()), 3) if I.any() else None,
               ratio_fixed=[round(float(x), 3) for x in r0], ratio_est=[round(float(x), 3) for x in r1], t=round(time.time() - t0, 1))
    print(json.dumps(res), flush=True)
