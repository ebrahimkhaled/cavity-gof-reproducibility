# Check 1: exact A = sum_k d_{k,1} against the closed expansion of Theorem A, and c_1 against gamma_hat * a * a.
# Run: python check_expansion.py n reps inner
import sys, numpy as np, time
from common import *

def psi_factory(t):
    c0, c1, h = 0.2 * t, 1.2 * t, 0.2 * t      # an off-centre smoothed cell, in index units
    def psi(u):
        return expit((u - c0) / h) - expit((u - c1) / h)
    return psi

class Fam:
    """F(u, y; sigma) and its derivatives by central differences."""
    def __init__(self, psi, d):
        self.psi, self.d = psi, d
    def F(self, u, y, s1, s2, s3):
        a, B, t2 = theta_of(s1, s2, s3, self.d.m0, self.d.s2)
        return self.psi(u) * (y - Lam(u, a, B, t2))
    def all(self, u, y, s1, s2, s3, hu, hs):
        F = self.F
        Fu = (F(u + hu, y, s1, s2, s3) - F(u - hu, y, s1, s2, s3)) / (2 * hu)
        Fuu = (F(u + hu, y, s1, s2, s3) - 2 * F(u, y, s1, s2, s3) + F(u - hu, y, s1, s2, s3)) / hu ** 2
        F1 = (F(u, y, s1 + hs, s2, s3) - F(u, y, s1 - hs, s2, s3)) / (2 * hs)
        F2 = (F(u, y, s1, s2 + hs, s3) - F(u, y, s1, s2 - hs, s3)) / (2 * hs)
        F3 = (F(u, y, s1, s2, s3 + hs) - F(u, y, s1, s2, s3 - hs)) / (2 * hs)
        return Fu, Fuu, F1, F2, F3

def a_of(fam, s1, s2, s3, d, hu=1e-3):
    """a(psi) = E F_u (y - xi) under the law with moments sigma, by double quadrature."""
    a, B, t2 = theta_of(s1, s2, s3, d.m0, d.s2)
    Z = np.linspace(-8, 8, 801); Wz = np.exp(-Z ** 2 / 2); Wz = Wz / Wz.sum()   # dense: the cell is too sharp for Gauss-Hermite
    xi = s1 + np.sqrt(s3) * Z
    L = Lam(xi, a, B, t2)
    ps = fam.psi
    dps = (ps(xi + hu) - ps(xi - hu)) / (2 * hu)
    dL = (Lam(xi + hu, a, B, t2) - Lam(xi - hu, a, B, t2)) / (2 * hu)
    return np.sum(Wz * (dps * L * (1 - L) + ps(xi) * dL * (xi - L)))

def one_outer(d, rng, inner):
    n, p = d.n, d.p
    A, y = d.draw(rng, n - 1)
    M = A.T @ A + np.diag(d.D)
    R = np.linalg.inv(M)
    b = R @ (A.T @ y)
    Rho = A @ R                                   # rows rho_k = R a_k
    ht = np.einsum('ij,ij->i', Rho, A)            # a_k' R a_k = h/(1+h)
    h = ht / (1 - ht)
    e = (y - A @ b) / (1 - ht)                    # y_k - a_k' b^k
    Bk = b[None, :] - Rho * e[:, None]            # rows b^k
    u_k = y - e                                   # a_k' b^k
    s1, s2, s3 = d.sigma(Bk)
    St = np.eye(p + 1); St[0, 0] = 0.0
    tval = np.sqrt(np.mean(s3))
    fam = Fam(psi_factory(tval), d)
    hu, hs = 1e-3 * tval, 1e-4 * tval ** 2
    Fu, Fuu, F1, F2, F3 = fam.all(u_k, y, s1, s2, s3, hu, hs)
    # closed-form coefficients
    RS = R @ St
    gbar = np.trace(RS) + R[0, 0]
    om = np.trace(RS @ RS)
    sb1, sb2, sb3 = d.sigma(b[None, :])
    vbar = F1.sum() * np.eye(p + 1)[0] + F2.sum() * np.r_[0, d.b0] + 2 * F3.sum() * (St @ b)
    N = (Rho * (Fu * (1 + h))[:, None]).sum(0) + R @ vbar
    k1 = Fu * (1 + h) * e
    k2 = Fuu * (1 + h) ** 2
    Rv = []
    A1, y1 = d.draw(rng, inner)
    out_exact = np.empty(inner); out_main = np.empty(inner); f1 = np.empty(inner)
    aB, BB, tB = [float(v[0]) for v in theta_of(sb1, sb2, sb3, d.m0, d.s2)]
    for j in range(inner):
        a1, yy = A1[j], y1[j]
        tau = Rho @ a1
        Ra1 = R @ a1
        g = a1 @ Ra1
        xi1 = a1 @ b
        e1 = yy - xi1
        gk = g + (1 + h) * tau ** 2
        epk = (e1 + tau * e) / (1 + gk)
        Wk = (Ra1[None, :] + Rho * ((1 + h) * tau)[:, None]) * epk[:, None]     # rows w_k
        Bn = Bk + Wk
        t1, t2, t3 = d.sigma(Bn)
        un = u_k + np.einsum('ij,ij->i', A, Wk)
        dk = fam.F(un, y, t1, t2, t3) - fam.F(u_k, y, s1, s2, s3)
        out_exact[j] = dk.sum()
        ep = e1 / (1 + g)
        L = a1 @ N
        Q1 = np.sum(k1 * tau ** 2)
        Q2 = np.sum(k2 * tau ** 2)
        out_main[j] = ep * L + Q1 / (1 + g) + ep ** 2 * (0.5 * Q2 + F3.sum() * om)
        f1[j] = fam.psi(xi1) * (yy - Lam(np.array(xi1), aB, BB, tB).item())
    c1 = np.mean(f1 * out_exact)
    c1m = np.mean(f1 * out_main)
    aval = a_of(fam, sb1[0], sb2[0], sb3[0], d)
    gam = n * om / (1 + gbar) ** 2
    return dict(rms=np.sqrt(np.mean((out_exact - out_main) ** 2)), sdA=np.std(out_exact), c1=c1, c1m=c1m, pred=gam * aval ** 2,
                gam=gam, a=aval, mf=np.mean(f1), se=np.std(f1 * out_exact) / np.sqrt(inner))

if __name__ == '__main__':
    n = int(sys.argv[1]); reps = int(sys.argv[2]); inner = int(sys.argv[3]); lam = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    d = Design(n, n // 4, lam * n, eps0=0.0)
    rng = np.random.default_rng(1000 + n)
    t0 = time.time(); res = [one_outer(d, rng, inner) for _ in range(reps)]
    g = lambda k: np.array([r[k] for r in res])
    print(f"n={n} lam/n={lam} rms(A-main)={np.mean(g('rms')):.4f} sqrt(n)*rms={np.sqrt(n)*np.mean(g('rms')):.3f} sdA={np.mean(g('sdA')):.3f} "
          f"c1 mean={g('c1').mean():.4f} sd over outer={g('c1').std():.4f} (inner se {g('se').mean():.4f}) c1main mean={g('c1m').mean():.4f} "
          f"pred mean={g('pred').mean():.4f} sd={g('pred').std():.4f} gamma={g('gam').mean():.4f} a={g('a').mean():.4f} corr={np.corrcoef(g('c1'),g('pred'))[0,1]:.3f} t={time.time()-t0:.0f}s")
