"""(Documentation copy in the paper-B MATLAB package; see README.md in this folder.)

Random irreducible Metzler pairs for the scaling pilot of the Euler + GM pipeline.

Generator (a 3D generator of the forthcoming companion work, rescaled so the expected off-diagonal row sum does not grow with d):
  dense : off-diagonal entries uniform in [0, 1] * 2/(d-1)          (row sum ~ 1, as for d = 3)
          diagonal uniform in [-3, 1]
  sparse: off-diagonal entry nonzero with probability p_d, then lognormal(0, 1.5) * s_d,
          s_d chosen so the expected row sum is 2.46 (the 3D value 2 * 0.4 * e^{1.125})
          diagonal uniform in [-4, 1]
Every mode is irreducible (resampled otherwise), so cor:gripenberg applies.
Also rejected: h*c_A >= 1 or h*eta(A_i) >= 1 at h = 0.1 (Euler families undefined), and pairs whose
best constant mixture is (near) a single mode, alpha* outside [0.05, 0.95]: there one mode is typically
optimal on its own, the discrete problems have length-one optimal products, and nothing is tested.

Writes families.json with, per pair: id, d, kind, A1, A2, eta(A_i), chattering value
c = min_alpha eta(S_alpha) (an upper bound on lambda_check), alpha*.
"""
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import connected_components

SEED = 2026
DIMS = (5, 10)
PER_KIND = 3
P_SPARSE = {5: 0.5, 10: 0.3}
H_MAX = 0.1


def eta(S):
    return float(np.max(np.linalg.eigvals(S).real))


def irreducible(A):
    off = (A > 0) & ~np.eye(len(A), dtype=bool)
    return connected_components(off, directed=True, connection="strong")[0] == 1


def sample(rng, kind, d):
    if kind == "dense":
        A = rng.uniform(0, 1, (d, d)) * 2 / (d - 1)
        np.fill_diagonal(A, rng.uniform(-3, 1, d))
    else:
        p = P_SPARSE[d]
        s = 2.46 / ((d - 1) * p * np.exp(1.125))
        A = np.where(rng.random((d, d)) < p, np.exp(rng.normal(0, 1.5, (d, d))) * s, 0.0)
        np.fill_diagonal(A, rng.uniform(-4, 1, d))
    return A


def chatter(A1, A2):
    al = np.linspace(0, 1, 2001)
    vals = [eta(a * A1 + (1 - a) * A2) for a in al]
    k = int(np.argmin(vals))
    lo, hi = al[max(k - 1, 0)], al[min(k + 1, len(al) - 1)]
    r = minimize_scalar(lambda a: eta(a * A1 + (1 - a) * A2), bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-12})
    if r.fun < vals[k]:
        return float(r.fun), float(r.x)
    return float(vals[k]), float(al[k])


def main():
    rng = np.random.default_rng(SEED)
    out = []
    rejected = {"dense": 0, "sparse": 0}
    for d in DIMS:
        for kind in ("dense", "sparse"):
            n = 0
            while n < PER_KIND:
                A1, A2 = sample(rng, kind, d), sample(rng, kind, d)
                if not (irreducible(A1) and irreducible(A2)):
                    continue
                cA = max(np.max(-np.diag(A1)), np.max(-np.diag(A2)))
                e1, e2 = eta(A1), eta(A2)
                if H_MAX * cA >= 1 or H_MAX * max(e1, e2) >= 1:
                    continue
                c, a = chatter(A1, A2)
                if not 0.05 <= a <= 0.95:
                    rejected[kind] += 1
                    continue
                n += 1
                out.append({"id": f"d{d}_{kind}{n}", "d": d, "kind": kind, "A1": A1.tolist(),
                            "A2": A2.tolist(), "eta": [e1, e2], "chatter": c, "alpha": a})
                print(f"{out[-1]['id']}: eta = ({e1:+.4f}, {e2:+.4f}), chattering c = {c:+.6f} at alpha = {a:.4f}")
    print("rejected for boundary alpha*:", rejected)
    with open("families.json", "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()
