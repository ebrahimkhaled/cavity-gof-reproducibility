# Check 5: c_1 with the integral over z_1 done by quadrature, piece by piece (no inner Monte Carlo), for a FIXED cell psi.
# c_ana = [EG c'N + a b'St N + x_beta beta'St N]/(1+gbar) + E(f e^2)(q2/2 + Phi3 om)/(1+gbar)^2, against gamma_hat a^2.
import sys, numpy as np
from multiprocessing import Pool
from common import *
from check_expansion import Fam

ZX = np.linspace(-8, 8, 801); ZW = np.exp(-ZX ** 2 / 2); ZW = ZW / ZW.sum()   # dense: the cell is too sharp for 32 Gauss-Hermite nodes
def moments(fam, d, s1, s2, s3):
    a, B, t2 = theta_of(s1, s2, s3, d.m0, d.s2)
    xi = s1 + np.sqrt(s3) * ZX                      # xi nodes
    tau = np.sqrt(max(t2, 0.0))
    eta = (a + B * xi)[:, None] + tau * GX[None, :]  # eta | xi
    P = expit(eta); W = ZW[:, None] * GW[None, :]
    L = Lam(xi, a, B, t2)[:, None]; ps = fam.psi(xi)[:, None]; X = xi[:, None]
    G1 = ps * (1 - L) * (1 - X); G0 = ps * (-L) * (-X)
    EG = np.sum(W * (P * G1 + (1 - P) * G0))
    xb = np.sum(W * P * (1 - P) * (G1 - G0))
    Ef2 = np.sum(W * (P * ps * (1 - L) * (1 - X) ** 2 + (1 - P) * ps * (-L) * X ** 2))
    return EG, xb, Ef2

def a_val(fam, d, s1, s2, s3, hu):
    a, B, t2 = theta_of(s1, s2, s3, d.m0, d.s2)
    xi = s1 + np.sqrt(s3) * ZX
    L = Lam(xi, a, B, t2); ps = fam.psi
    dps = (ps(xi + hu) - ps(xi - hu)) / (2 * hu)
    dL = (Lam(xi + hu, a, B, t2) - Lam(xi - hu, a, B, t2)) / (2 * hu)
    return np.sum(ZW * (dps * L * (1 - L) + ps(xi) * dL * (xi - L)))

def job(args):
    n, lam, seed, reps, tfix = args
    d = Design(n, n // 4, lam * n); p = d.p
    c0, c1_, hb = 0.2 * tfix, 1.2 * tfix, 0.2 * tfix
    fam = Fam(lambda u: expit((u - c0) / hb) - expit((u - c1_) / hb), d)
    rng = np.random.default_rng(seed); out = []
    St = np.eye(p + 1); St[0, 0] = 0; e0 = np.eye(p + 1)[0]; bt = np.r_[0, d.b0]
    for _ in range(reps):
        A, y = d.draw(rng, n - 1)
        R = np.linalg.inv(A.T @ A + np.diag(d.D))
        b = R @ (A.T @ y); Rho = A @ R
        ht = np.einsum('ij,ij->i', Rho, A); h = ht / (1 - ht)
        e = (y - A @ b) / (1 - ht); Bk = b[None, :] - Rho * e[:, None]; u = y - e
        s1, s2, s3 = d.sigma(Bk)
        Fu, Fuu, F1, F2, F3 = fam.all(u, y, s1, s2, s3, 1e-3 * tfix, 1e-4 * tfix ** 2)
        v = F1.sum() * e0 + F2.sum() * bt + 2 * F3.sum() * (St @ b)
        N = (Rho * (Fu * (1 + h))[:, None]).sum(0) + R @ v
        RS = R @ St; gb = np.trace(RS) + R[0, 0]; om = np.trace(RS @ RS)
        rr = np.einsum('ij,ij->i', Rho @ St, Rho)
        q1 = np.sum(Fu * (1 + h) * e * rr); q2 = np.sum(Fuu * (1 + h) ** 2 * rr)
        sb = [float(x[0]) for x in d.sigma(b[None, :])]
        EG, xb, Ef2 = moments(fam, d, *sb)
        aa = a_val(fam, d, *sb, 1e-3 * tfix)
        t1 = EG * (e0 @ N) / (1 + gb); t2 = aa * ((St @ b) @ N) / (1 + gb); t3 = 0.0
        # x_beta part
        a_, B_, tt = theta_of(*sb, d.m0, d.s2)
        xi = sb[0] + np.sqrt(sb[2]) * ZX
        # E[pi'(eta) (G(xi,1)-G(xi,0))]; beta'St N coefficient
        tau = np.sqrt(max(tt, 0.0)); eta = (a_ + B_ * xi)[:, None] + tau * GX[None, :]
        P = expit(eta); L = Lam(xi, a_, B_, tt)[:, None]; ps = fam.psi(xi)[:, None]; X = xi[:, None]
        xbeta = np.sum(ZW[:, None] * GW[None, :] * P * (1 - P) * (ps * (1 - L) * (1 - X) - ps * L * X))
        t3 = xbeta * (bt @ N) / (1 + gb)
        t4 = Ef2 * (0.5 * q2 + F3.sum() * om) / (1 + gb) ** 2
        gam = n * om / (1 + gb) ** 2
        out.append([t1, t2, t3, t4, gam * aa * aa, 0.0, aa * q1 / (1 + gb)])
    return np.array(out)

if __name__ == '__main__':
    lam = 0.25; tfix = 0.2856
    with Pool(20) as pool:
        for n, reps in [(200, 4000), (400, 4000), (800, 2000), (1600, 800), (3200, 200)]:
            X = np.vstack(pool.map(job, [(n, lam, 100 * n + i, reps // 20, tfix) for i in range(20)]))
            m = X.mean(0); se = X.std(0, ddof=1) / np.sqrt(len(X))
            c = X[:, :4].sum(1); dlt = c - X[:, 4]
            print(f"n={n} t(c'N) {m[0]:.5f} t(b'N) {m[1]:.5f} t(beta'N) {m[2]:.5f}({se[2]:.5f}) t(quad) {m[3]:.5f}({se[3]:.5f}) c_ana {c.mean():.5f}({c.std()/np.sqrt(len(c)):.5f}) sd {c.std():.5f} pred {m[4]:.5f} diff {dlt.mean():.5f}({dlt.std()/np.sqrt(len(c)):.5f}) sqrt(n)*diff {np.sqrt(n)*dlt.mean():.3f} xb-a {m[5]:.1e}", flush=True)
