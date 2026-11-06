# Brute-force c_1: true A = sum_k d_{k,1} (exact, nonlinear in z_1) by Monte Carlo over z_1. No expansion is used.
# Usage: python audit_c1.py tag n p lamn eps0 sig b00 mu ar D M     |  python audit_c1.py selftest
import sys, json, time
from multiprocessing import Pool
from acommon import *

def mkpsi(lo, hi, hb):
    return lambda u: expit((u - lo) / hb) - expit((u - hi) / hb)

def influence(d, A_, y_, R, b, X1, y1, psis):
    """A(psi) for draws X1 (P x M). Exact closed forms; validated against refits by selftest()."""
    St = d.St; c = d.c; Sb = d.Stb
    Rho = A_ @ R; RS = Rho @ St
    ht = np.einsum('ij,ij->i', Rho, A_); h = ht / (1 - ht); H1 = 1 + h
    e = (y_ - A_ @ b) / (1 - ht); u = y_ - e
    Sbh = St @ b
    cr, br_, hr = Rho @ c, Rho @ Sb, Rho @ Sbh
    rr = np.einsum('ij,ij->i', RS, Rho)
    s1k = c @ b - e * cr; s2k = Sb @ b - e * br_; s3k = b @ Sbh - 2 * e * hr + e * e * rr
    L = Rho @ X1; RX = R @ X1; Q = RS @ RX
    g = np.einsum('ij,ij->j', X1, RX); xi = b @ X1; e1 = y1 - xi
    cR, bR, hR = c @ RX, Sb @ RX, Sbh @ RX
    qq = np.einsum('ij,ij->j', RX, St @ RX)
    Hc = H1[:, None]; ec = e[:, None]
    D = Hc * L * L; ek = (e1 + L * ec) / (1 + g + D)
    un = u[:, None] + Hc * L * ek
    s1n = s1k[:, None] + ek * (cR + Hc * cr[:, None] * L)
    s2n = s2k[:, None] + ek * (bR + Hc * br_[:, None] * L)
    cross = hR + Hc * hr[:, None] * L - ec * Q - ec * Hc * rr[:, None] * L
    s3n = s3k[:, None] + 2 * ek * cross + ek ** 2 * (qq + 2 * Hc * L * Q + Hc ** 2 * rr[:, None] * L * L)
    a0, B0, t0 = d.theta(s1k, s2k, s3k); a1, B1, t1 = d.theta(s1n, s2n, s3n)
    Lold = Lam(u, a0, B0, t0); Lnew = Lam(un, a1, B1, t1)
    yc = y_[:, None]; out = []
    for ps in psis:
        out.append((ps(un) * (yc - Lnew) - (ps(u) * (y_ - Lold))[:, None]).sum(0))
    return out, xi

def selftest():
    d = Des(40, 12, .25, eps0=.3, sig=1.5, b00=-.7, mu=1., ar=.5); r = np.random.default_rng(1)
    A_, y_, _ = d.draw(r, d.n - 1); R, b = fit_all(d, A_, y_)
    X1, y1, _ = d.draw(r, 3); ps = mkpsi(.2, .7, .1)
    out, _ = influence(d, A_, y_, R, b, X1.T, y1, [ps])
    def Psi(k, Aset, yset):
        Rr = np.linalg.inv(Aset.T @ Aset + np.diag(d.D)); bb = Rr @ (Aset.T @ yset)
        s1, s2, s3 = d.c @ bb, d.Stb @ bb, bb @ d.St @ bb
        a, B, t = d.theta(s1, s2, s3); uu = A_[k] @ bb
        return ps(uu) * (y_[k] - Lam(uu, a, B, t))
    for m in range(3):
        tot = 0
        for k in range(d.n - 1):
            keep = np.arange(d.n - 1) != k
            Aw = np.vstack([A_[keep], X1[m]]); yw = np.r_[y_[keep], y1[m]]
            tot += Psi(k, Aw, yw) - Psi(k, A_[keep], y_[keep])
        print("selftest A exact", tot, "closed", out[0][m])

def job(a):
    cfg, seed, D, M = a
    d = Des(*cfg['des'], **cfg['kw']); r = np.random.default_rng(seed)
    p1, p2 = mkpsi(*cfg['p1']), mkpsi(*cfg['p2']); res = []
    for _ in range(D):
        A_, y_, _ = d.draw(r, d.n - 1); R, b = fit_all(d, A_, y_)
        X1, y1, _ = d.draw(r, M); X1 = X1.T
        (A1, A2), xi = influence(d, A_, y_, R, b, X1, y1, [p1, p2])
        s1, s2, s3 = d.c @ b, d.Stb @ b, b @ d.St @ b
        a_, B_, t_ = d.theta(s1, s2, s3); L1 = Lam(xi, a_, B_, t_)
        f1 = p1(xi) * (y1 - L1); f2 = p2(xi) * (y1 - L1)
        RS = R @ d.St; gb = np.trace(RS) + d.c @ R @ d.c; om = np.trace(RS @ RS)
        gam = d.n * om / (1 + gb) ** 2
        aa1, aa2 = a_of(d, p1, s1, s2, s3), a_of(d, p2, s1, s2, s3)
        # E(f|F)=0 exactly, so centring A is a legal control variate
        res.append([np.mean(f1 * (A1 - A1.mean())), np.mean(f1 * (A2 - A2.mean())), np.mean(f2 * (A2 - A2.mean())),
                    gam * aa1 * aa1, gam * aa1 * aa2, gam * aa2 * aa2, f1.mean(), gam])
    return np.array(res)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        selftest(); sys.exit()
    tag, n, p = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    lamn, eps0, sig, b00, mu, ar = map(float, sys.argv[4:10]); D, M = int(sys.argv[10]), int(sys.argv[11])
    kw = dict(eps0=eps0, sig=sig, b00=b00, mu=mu, ar=ar)
    d = Des(n, p, lamn, **kw); r = np.random.default_rng(3); bs = []
    for _ in range(30):
        A_, y_, _ = d.draw(r, n - 1); _, b = fit_all(d, A_, y_); bs.append([d.c @ b, d.Stb @ b, b @ d.St @ b])
    s = np.mean(bs, 0); t = np.sqrt(s[2])
    cfg = dict(des=(n, p, lamn), kw=kw, p1=(s[0] + .2 * t, s[0] + 1.2 * t, .2 * t), p2=(s[0] - 1.0 * t, s[0] + .1 * t, .2 * t))
    W = 20; t0 = time.time()
    with Pool(W) as pool:
        X = np.vstack(pool.map(job, [(cfg, 1000 * n + i, D // W, M) for i in range(W)]))
    m = X.mean(0); nn = len(X)
    dl = X[:, :3] - X[:, 3:6]; dm = dl.mean(0); ds = dl.std(0, ddof=1) / np.sqrt(nn)
    print(json.dumps(dict(tag=tag, n=n, p=p, D=nn, M=M, sigma=[round(float(x), 4) for x in s], theta=[round(float(x), 3) for x in d.theta(*s)],
                          gam=round(float(m[7]), 5), c=[round(float(x), 6) for x in m[:3]], pred=[round(float(x), 6) for x in m[3:6]],
                          diff=[round(float(x), 6) for x in dm], se=[round(float(x), 6) for x in ds],
                          sqrtn_sd_c11=round(float(np.sqrt(n) * X[:, 0].std()), 4), t=round(time.time() - t0, 1))), flush=True)
