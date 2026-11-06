# gamma-tilde of Prop est against gamma_hat, and the same with the mean direction removed from H (the (c'Rc)^2 term).
from acommon import *
def one(n, p, lamn, ar, eps0, mu, reps, seed=1):
    d = Des(n, p, lamn, eps0=eps0, sig=2.0, b00=0.0, mu=mu, ar=ar); r = np.random.default_rng(seed); o = []
    for _ in range(reps):
        A_, y, _ = d.draw(r, n); R, b = fit_all(d, A_, y); H = A_ @ R @ A_.T
        z1 = np.trace(H) / n; z2 = (np.trace(H @ H) - np.sum(np.diag(H) ** 2)) / n
        gt = z2 / (1 - z1) ** 2
        RS = R @ d.St; gb = np.trace(RS) + d.c @ R @ d.c; om = np.trace(RS @ RS); gh = n * om / (1 + gb) ** 2
        cRc = d.c @ R @ d.c; cRSRc = d.c @ R @ d.St @ R @ d.c
        corr = n * (cRc ** 2 + 2 * cRSRc) / (1 + gb) ** 2
        o.append([gt, gh, corr])
    o = np.array(o); m = o.mean(0)
    print(f"n={n} p={p} lam/n={lamn} ar={ar} eps0={eps0} mu={mu}: gamma_tilde {m[0]:.5f} gamma_n {m[1]:.5f} diff {m[0]-m[1]:.5f} rel {(m[0]-m[1])/m[1]:.3f}  n(c'Rc)^2-type term {m[2]:.5f}  diff-after {m[0]-m[1]-m[2]:.5f}", flush=True)
for a in [(400, 100, 4.16, .7, 0., 0.), (400, 40, 1.04, .7, 0., 0.), (300, 600, 4.16, .7, 0., 0.), (800, 200, 4.16, .7, 0., 0.), (1600, 400, 4.16, .7, 0., 0.),
          (400, 100, .25, 0., 0., 0.), (400, 100, 4.16, .7, 0., 1.)]:
    one(*a, reps=40)
