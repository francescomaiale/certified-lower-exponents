"""Independent Python replay of the certified chattering and schedule bounds: replay() of certify_best_upper.py
(uniformization enclosures, no logarithms), applied to given files.  Extended, in a backward-compatible way, to
schedules with more than two segments (fields 'modes' and 'durations', as written by the MATLAB certifier for the
three-mode example of Ex. 4.10) and to certificates without a chattering part.

Usage: python3 replay_best_upper.py FILE [FILE ...]
"""
import sys, os, json, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common

CASES = common.merged_cases_file()


def unif_upper(A, t, K=200, bits=80):
    """e^{tA} <= e^{-qt} [sum_{k<=K} (qt)^k/k! P^k + tail J], P = I + A/q >= 0 (uniformization), for t / 2^s with
    q t nP / 2^s <= 1, then squared s times (monotone for nonnegative matrices); own dyadic grid."""
    if t <= 0:
        raise ValueError("schedule durations must be positive")
    common.check_family([A])
    G = 2 ** bits
    upg = lambda z: Fr(-((-z.numerator * G) // z.denominator), G)
    d = len(A)
    q = max(Fr(1), max(-A[i][i] for i in range(d)))
    P = [[(Fr(int(i == j)) + A[i][j] / q) for j in range(d)] for i in range(d)]
    if not all(P[i][j] >= 0 for i in range(d) for j in range(d)):
        raise ValueError("uniformization requires a nonnegative matrix")
    nP = max(sum(P[i][j] for j in range(d)) for i in range(d))          # row-sum norm bounds every entry of P^k
    s = 0
    while q * t * nP / 2 ** s > 1:
        s += 1
    E = unif_small(P, q, t / 2 ** s, nP, K, upg)
    for _ in range(s):
        E = [[upg(sum(E[i][l] * E[l][j] for l in range(d))) for j in range(d)] for i in range(d)]
    return E


def unif_small(P, q, t, nP, K, upg):
    d = len(P)
    qt = q * t
    S = [[Fr(0)] * d for _ in range(d)]
    Pk = [[Fr(int(i == j)) for j in range(d)] for i in range(d)]; coef = Fr(1)
    for k in range(K + 1):
        if k:
            Pk = [[upg(sum(Pk[i][l] * P[l][j] for l in range(d))) for j in range(d)] for i in range(d)]
            coef = coef * qt / k
        S = [[S[i][j] + upg(coef * Pk[i][j]) for j in range(d)] for i in range(d)]
    z = qt * nP; term = z ** (K + 1) / math.factorial(K + 1)
    tail = upg(term / (1 - z / (K + 2)))
    assert z / (K + 2) < 1
    emqt = upg(1 / sum(qt ** k / math.factorial(k) for k in range(K + 1)))
    return [[upg(emqt * (S[i][j] + tail)) for j in range(d)] for i in range(d)]


def exp_lower(y, K=80):
    """lower bound for e^y, y rational."""
    if type(K) is not int or K < 1:
        raise ValueError("K must be a positive integer")
    if y >= 0:
        return sum(y ** k / math.factorial(k) for k in range(K + 1))
    z = -y; squarings = 0
    while z >= K + 2:
        z /= 2; squarings += 1
    s = sum(z ** k / math.factorial(k) for k in range(K + 1))
    s += z ** (K + 1) / math.factorial(K + 1) / (1 - z / (K + 2))         # e^z <= s (z < K + 2)
    return (1 / s) ** (2 ** squarings)


def replay_file(f):
    try:
        return _replay_file(f)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError) as exc:
        print(f"replay {os.path.basename(f)}: INVALID ({exc})", flush=True)
        return False


def _replay_file(f):
    R = common.read_json(f); pid = R["id"]
    raw = common.read_json(CASES, parse_float=str, parse_int=str)[pid]
    readings = {"decimal": [[[Fr(v) for v in row] for row in A] for A in raw],
                "binary64": [[[Fr(float(v)) for v in row] for row in A] for A in raw]}
    good = True
    has_c = R.get("U_c") not in (None, [], "")
    has_s = R.get("U_s") not in (None, [], "")
    for name, As in readings.items():
        common.check_family(As)
        rd = R["readings"][name]; d = len(As[0])
        if has_c:
            a = Fr(rd["alpha"]); x = common.vector(rd["x"], d); U = Fr(R["U_c"])
            if len(As) != 2 or not 0 <= a <= 1:
                raise ValueError("mixture requires two modes and 0 <= alpha <= 1")
            # (S_alpha x)_r <= U x_r for every r, cross-multiplied (no division)
            good &= all(x[r] > 0 and sum((a * As[0][r][c] + (1 - a) * As[1][r][c]) * x[c] for c in range(d)) <= U * x[r]
                        for r in range(d))
        if has_s:
            if "modes" in rd or "durations" in rd:
                modes = rd["modes"] if isinstance(rd["modes"], list) else [rd["modes"]]
                ts = [Fr(v) for v in (rd["durations"] if isinstance(rd["durations"], list) else [rd["durations"]])]
            else:
                modes, ts = [0, 1], [Fr(rd["t1"]), Fr(rd["t2"])]
            xs = common.vector(rd["x_sched"], d); Us = Fr(R["U_s"])
            if (not modes or len(modes) != len(ts)
                    or any(type(m) is not int or not 0 <= m < len(As) for m in modes)
                    or any(t <= 0 for t in ts) or any(x <= 0 for x in xs)):
                raise ValueError("invalid schedule or schedule witness")
            Phi = [[Fr(int(i == j)) for j in range(d)] for i in range(d)]
            for m, t in zip(modes, ts):
                E = unif_upper(As[m], t)
                Phi = [[sum(E[i][l] * Phi[l][j] for l in range(d)) for j in range(d)] for i in range(d)]
            rho = max(sum(Phi[i][k] * xs[k] for k in range(d)) / xs[i] for i in range(d))
            good &= all(t > 0 for t in ts) and all(v > 0 for v in xs) and rho <= exp_lower(Us * sum(ts))
    parts = [Fr(R[k]) for k, present in (("U_c", has_c), ("U_s", has_s)) if present]
    good &= bool(parts) and Fr(R["U_best"]) == min(parts)
    if "printed" in R:
        good &= Fr(R["printed"]) >= Fr(R["U_best"])
    print(f"replay {pid:12s}: " + (f"certified c <= {float(Fr(R['U_c'])):+.9f}" if R.get("U_c") else "")
          + (f", schedule <= {float(Fr(R['U_s'])):+.9f}" if R.get("U_s") else "")
          + f": {'VALID' if good else 'INVALID'}", flush=True)
    return good


def main(paths):
    ok = bool(paths) and all([replay_file(f) for f in paths])
    print("replay:", "all VALID" if ok else "SOME INVALID")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
