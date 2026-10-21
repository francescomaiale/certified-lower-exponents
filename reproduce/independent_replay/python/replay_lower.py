"""Independent Python replay of polytope certificates (prop:polytope).

Usage: python3 replay_lower.py FILE [FILE ...]
"""
import sys, os, json, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common

CASES = common.merged_cases_file()


def fmt_dec(x, n):
    s = 10 ** n; m = x * s; assert m.denominator == 1
    m = m.numerator; sign = "-" if m < 0 else ""; m = abs(m)
    return f"{sign}{m // s}.{m % s:0{n}d}"


def replay_one(path):
    try:
        return _replay_one(path)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError, StopIteration) as exc:
        return False, 0, [f"malformed certificate: {exc}"], []


def _replay_one(path):
    """Independent re-check, sharing no code with the certify path: the data are re-read from cases.json as decimal
    strings and as binary64, B_i is rebuilt, (*) is checked by integer cross-multiplication, and t* is recomputed."""
    C = common.read_json(path)
    raw = common.read_json(CASES, parse_float=str, parse_int=str)[C["case"]]
    hq = Fr(C["h"]); t = Fr(C["t"]); orient = C["orientation"]
    if hq <= 0 or orient not in ("N", "T"):
        raise ValueError("positive step and orientation N or T required")
    verts = [[Fr(x) for x in v] for v in C["vertices"]]; d = len(verts[0]); n = len(verts)
    if not d or any(len(v) != d for v in verts):
        raise ValueError("invalid vertices")
    for wt in C["weights"]:
        if (type(wt["mode"]) is not int or not 0 <= wt["mode"] < len(raw)
                or type(wt["vertex"]) is not int or not 0 <= wt["vertex"] < n
                or any(str(int(k)) != k or not 0 <= int(k) < n for k in wt["lambda"])):
            raise ValueError("invalid mode or vertex index")
    msgs = []; ok = True; tstar = None
    pairs = sorted((w["mode"], w["vertex"]) for w in C["weights"])
    if pairs != [(i, j) for i in range(len(raw)) for j in range(n)]:
        ok = False; msgs.append("weights do not cover every (mode, vertex) pair exactly once")
    if not t > 0:
        ok = False; msgs.append("t <= 0")
    if not all(x >= 0 for v in verts for x in v) or not all(sum(v) > 0 for v in verts):
        ok = False; msgs.append("vertex not in K \\ {0}")
    for ver, conv in (("decimal", lambda s: Fr(s)), ("binary64", lambda s: Fr(float(s)))):
        As = [[[conv(x) for x in row] for row in A] for A in raw]
        common.check_family(As)
        if len(As[0]) != d:
            raise ValueError("wrong vertex dimension")
        if not all(hq * abs(A[r][r]) < 1 for A in As for r in range(d)):
            ok = False; msgs.append(f"{ver}: h >= h_0")
        Bs = []
        for A in As:
            B = [[int(r == c) + hq * A[r][c] for c in range(d)] for r in range(d)]
            if orient == "T":
                B = [[B[c][r] for c in range(d)] for r in range(d)]
            if any(x < 0 for row in B for x in row):
                ok = False; msgs.append(f"{ver}: B_i not nonnegative")
            Bs.append(B)
        for wt in C["weights"]:
            i, j = wt["mode"], wt["vertex"]; lam = {int(k): Fr(x) for k, x in wt["lambda"].items()}
            if any(x < 0 for x in lam.values()) or sum(lam.values()) != 1:
                ok = False; msgs.append(f"bad weights at mode {i}, vertex {j}")
            # integer form: L = lcm of all denominators; check L*B v - t L V lam >= 0 componentwise
            Bv = [sum((Bs[i][r][c] * verts[j][c] for c in range(d)), Fr(0)) for r in range(d)]
            Vl = [sum((lam[k] * verts[k][r] for k in lam), Fr(0)) for r in range(d)]
            for r in range(d):
                num_l, den_l = Bv[r].numerator, Bv[r].denominator
                rhs = t * Vl[r]
                if num_l * rhs.denominator < rhs.numerator * den_l:
                    ok = False; msgs.append(f"{ver}: (*) fails at mode {i}, vertex {j}, coordinate {r}")
                if Vl[r] > 0:
                    q = Bv[r] / Vl[r]
                    tstar = q if tstar is None or q < tstar else tstar
    if tstar != Fr(C["t_star"]) or t > tstar:
        ok = False; msgs.append(f"recomputed t* = {float(tstar)} differs from stored {C['t_star_float']} or is < t")
    lam_claim = Fr(C["lambda_lower"])
    if not (t - 1) / hq >= lam_claim:
        ok = False; msgs.append("claimed lambda bound not implied by t")
    # negative controls: the stored weights fail just above t*, and (dual) the polytope is not t-invariant above tau_up
    neg = []
    t_raise = Fr(C["t_star"]) * (1 + Fr(1, 10 ** 12))
    b = C["binding"]; i, j = b["mode"], b["vertex"]
    conv = (lambda s: Fr(s)) if b["data"] == "decimal" else (lambda s: Fr(float(s)))
    A = [[conv(x) for x in row] for row in raw[i]]
    B = [[int(r == c) + hq * A[r][c] for c in range(d)] for r in range(d)]
    if orient == "T":
        B = [[B[c][r] for c in range(d)] for r in range(d)]
    Bv = [sum((B[r][c] * verts[j][c] for c in range(d)), Fr(0)) for r in range(d)]
    wt = next(w for w in C["weights"] if w["mode"] == i and w["vertex"] == j)
    lam = {int(k): Fr(x) for k, x in wt["lambda"].items()}
    Vl = [sum((lam[k] * verts[k][r] for k in lam), Fr(0)) for r in range(d)]
    fails = any(Bv[r] - t_raise * Vl[r] < 0 for r in range(d))
    neg.append(f"weights at t*(1+1e-12): {'fail (expected)' if fails else 'PASS (unexpected)'}")
    ok = ok and fails
    if C["negative_control"]["c"] is not None:
        c = [Fr(x) for x in C["negative_control"]["c"]]
        dual_ok = all(x >= 0 for x in c) and all(sum(c[r] * v[r] for r in range(d)) >= 1 for v in verts)
        tau = sum(c[r] * Bv[r] for r in range(d))
        dual_ok = dual_ok and tau == Fr(C["negative_control"]["tau_up"]) and tau >= Fr(C["t_star"])
        neg.append(f"dual: polytope not t-invariant for t > tau_up = {float(tau):.15f} "
                   f"(tau_up/t* - 1 = {float(tau / Fr(C['t_star']) - 1):.1e}) {'verified' if dual_ok else 'FAILED'}")
        ok = ok and dual_ok
        nxt = 1 + hq * (lam_claim + Fr(1, 10 ** 6))
        neg.append(f"next six-decimal value {fmt_dec(lam_claim + Fr(1, 10 ** 6), 6)}: "
                   f"{'excluded for this polytope' if nxt > tau else 'not excluded (slack between t and t*)'}")
    else:
        neg.append(f"V = {{1}}: t* is the exact min row sum, so any t > t* fails at the binding row")
    return ok, n, msgs, neg


def main(paths):
    allok = bool(paths)
    for p in paths:
        t0 = time.time(); ok, n, msgs, neg = replay_one(p); allok &= ok
        C = common.read_json(p)
        print(f"{os.path.basename(p)}: {'VALID' if ok else 'INVALID'}  lambda_check >= {C['lambda_lower']} "
              f"(t = {float(Fr(C['t'])):.12f}, |V| = {n}, {time.time() - t0:.1f}s)", flush=True)
        for m in msgs + neg:
            print(f"    {m}")
    print(f"replay: {'all certificates VALID' if allok else 'SOME CERTIFICATES INVALID'}")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
