# Is Delta-hat of Prop est informative at the size of Delta_n? Favourable version: everything at the fixed theta-bar, no root error.
# Usage: python audit_dhat.py n p lamn reps
import sys
from multiprocessing import Pool
from statistics import NormalDist
from acommon import *
from audit_e import cells, scores

def coords(u, thb, phi, w):
    L, v, D = scores(u, thb)
    return np.c_[phi(u) / w, D / v[:, None]], L, v, D

def work(a):
    cfg, seed, reps = a
    d = Des(*cfg['des'], **cfg['kw']); r = np.random.default_rng(seed); thb = cfg['thb']; phi = cells(cfg['cuts'], cfg['hb']); w = np.array(cfg['w']); G = len(w)
    t = cfg['t']; hh = 1e-3 * t; o = []
    for _ in range(reps):
        A_, y, _ = d.draw(r, d.n); R, b = fit_all(d, A_, y); Rho = A_ @ R
        ht = np.einsum('ij,ij->i', Rho, A_); ee = (y - A_ @ b) / (1 - ht); xi = y - ee
        Cc, L, v, D = coords(xi, thb, phi, w)
        Cp = coords(xi + hh, thb, phi, w)[0]; Cm = coords(xi - hh, thb, phi, w)[0]
        dps = (Cp - Cm) / (2 * hh); Lp = thb[1] * D[:, 0]
        Fu = dps * (y - L)[:, None] - Cc * Lp[:, None]
        ah = (Fu * (y - xi)[:, None]).mean(0)
        S = Cc[:, G:]; Fh = (S * v[:, None]).T @ S / d.n
        Uh = (Cc[:, :G] * v[:, None]).T @ S / d.n
        Ph = np.c_[np.eye(G), -Uh @ np.linalg.inv(Fh)]
        H = A_ @ R @ A_.T; z1 = np.trace(H) / d.n; z2 = (np.trace(H @ H) - np.sum(np.diag(H) ** 2)) / d.n
        o.append(z2 / (1 - z1) ** 2 * np.sum((Ph @ ah) ** 2))
    return np.array(o)

if __name__ == '__main__':
    n, p, lamn, reps = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]); G = 10
    kw = dict(eps0=0., sig=2., b00=0., mu=0., ar=0.)
    d = Des(n, p, lamn, **kw); r = np.random.default_rng(3); bs = []; gs = []
    for _ in range(30):
        A_, y_, _ = d.draw(r, n - 1); R, b = fit_all(d, A_, y_); bs.append([d.c @ b, d.Stb @ b, b @ d.St @ b])
        RS = R @ d.St; gb = np.trace(RS) + d.c @ R @ d.c; gs.append(n * np.trace(RS @ RS) / (1 + gb) ** 2)
    s = np.mean(bs, 0); gam = float(np.mean(gs)); thb = tuple(float(x) for x in d.theta(*s)); t = float(np.sqrt(s[2]))
    cuts = [s[0] + t * NormalDist().inv_cdf(g / G) for g in range(1, G)]; hb = .05 * t; phi = cells(cuts, hb)
    xi = s[0] + t * ZX; L, v, D = scores(xi, thb); S = D / v[:, None]; Ph = phi(xi); W = ZW * v; w = np.sqrt(Ph.T @ W)
    F = (S * W[:, None]).T @ S; U = ((Ph * W[:, None]).T @ S) / w[:, None]; Pm = np.c_[np.eye(G), -U @ np.linalg.inv(F)]
    h = v * (xi - s[0]) / t ** 2 + thb[1] * D[:, 0] * (xi + L - 1)
    av = (np.c_[Ph / w, S] * ZW[:, None]).T @ h; trD = gam * np.sum((Pm @ av) ** 2)
    cfg = dict(des=(n, p, lamn), kw=kw, thb=thb, cuts=cuts, hb=hb, w=list(w), t=t)
    with Pool(10) as pool:
        X = np.concatenate(pool.map(work, [(cfg, 500 + i, reps // 10) for i in range(10)]))
    print(f"n={n} tr Delta_n {trD:.3g}  mean tr Delta-hat {X.mean():.3g} (se {X.std()/np.sqrt(len(X)):.2g})  median {np.median(X):.3g}  ratio {X.mean()/trD:.1f}", flush=True)
