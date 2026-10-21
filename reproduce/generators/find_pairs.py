"""(Documentation copy in the paper-B MATLAB package; see README.md in this folder.)

Switching-regime pairs for the Perron-seeded closure pilot (paper B, Sec. 4).

Draws sparse pairs with the generator of gen_families.py (both modes irreducible, Euler families defined
at h = 0.1), and keeps those on which an optimized two-segment periodic schedule beats the chattering value
c = min_alpha eta(S_alpha) by at least MARGIN in floating point. This selects candidate switching-regime pairs;
it does not certify the global mixture minimum or a strict inequality lambda_check < c.

Usage: python3 find_pairs.py [per_dim] [max_draws]   -> pairs.json, find_pairs.log
"""
import json, math, os, sys, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
from multiprocessing import Pool
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_families import sample, irreducible, chatter, eta, H_MAX

SEED = 2027
DIMS = (5, 10)
MARGIN = 1e-3


def rate(A1, A2, t1, t2):
    if not (1e-4 < t1 < 40 and 1e-4 < t2 < 40):
        return 1e9
    P = expm(t2 * A2) @ expm(t1 * A1)
    r = max(abs(np.linalg.eigvals(P)))
    return math.log(r) / (t1 + t2) if r > 0 and np.isfinite(r) else 1e9


def best_schedule(A1, A2):
    grid = [(rate(A1, A2, a, b), a, b) for a in (0.2, 0.5, 1, 2, 4) for b in (0.2, 0.5, 1, 2, 4)]
    best = min(grid)
    for _, a, b in sorted(grid)[:3]:
        r = minimize(lambda z: rate(A1, A2, math.exp(z[0]), math.exp(z[1])), [math.log(a), math.log(b)],
                     method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=1500))
        best = min(best, (r.fun, math.exp(r.x[0]), math.exp(r.x[1])))
    return best


def draw(args):
    d, k = args
    rng = np.random.default_rng([SEED, d, k])
    while True:
        A1, A2 = sample(rng, "sparse", d), sample(rng, "sparse", d)
        if not (irreducible(A1) and irreducible(A2)):
            continue
        cA = max(np.max(-np.diag(A1)), np.max(-np.diag(A2)))
        if H_MAX * cA >= 1 or H_MAX * max(eta(A1), eta(A2)) >= 1:
            continue
        break
    c, a = chatter(A1, A2)
    s, t1, t2 = best_schedule(A1, A2)
    return dict(d=d, k=k, A1=A1.tolist(), A2=A2.tolist(), chatter=c, alpha=a, sched=s, t1=t1, t2=t2,
                gain=c - s)


def main():
    per_dim = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    max_draws = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    here = os.path.dirname(os.path.abspath(__file__))
    log = open(os.path.join(here, "find_pairs.log"), "w")
    out = []
    with Pool(10) as pool:
        for d in DIMS:
            T0 = time.time(); found = []; gains = []
            allr = sorted(pool.imap_unordered(draw, [(d, k) for k in range(max_draws)], chunksize=4),
                          key=lambda r: r["k"])
            for r in allr:                        # the first per_dim qualifying draws, in draw order
                gains.append(r["gain"])
                if r["gain"] >= MARGIN and len(found) < per_dim:
                    r["id"] = f"sw{d}_{len(found) + 1}"; found.append(r)
                    msg = (f"{r['id']} (draw {r['k']}): c = {r['chatter']:+.6f} at alpha = {r['alpha']:.3f}; "
                           f"schedule {r['sched']:+.6f} (t1 {r['t1']:.3f}, t2 {r['t2']:.3f}); gain {r['gain']:.2e}")
                    print(msg, flush=True); log.write(msg + "\n"); log.flush()
            g = np.array(gains)
            msg = (f"d = {d}: {len(g)} draws, {np.sum(g >= MARGIN)} with gain >= {MARGIN} "
                   f"({100 * np.mean(g >= MARGIN):.1f}%), {np.sum(g > 1e-6)} with gain > 1e-6; {time.time() - T0:.0f}s")
            print(msg, flush=True); log.write(msg + "\n"); log.flush()
            out += found
    json.dump(out, open(os.path.join(here, "pairs.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
