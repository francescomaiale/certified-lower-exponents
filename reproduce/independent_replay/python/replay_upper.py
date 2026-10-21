"""Independent Python replay of implicit-Euler chain certificates (prop:upper).

Usage: python3 replay_upper.py FILE [FILE ...]
"""
import sys, os, json, time, math
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common

_MERGED = common.merged_cases_file()


class _AnySource(dict):
    def __getitem__(self, key):
        return _MERGED


CASES = _AnySource()


def replay_one(path):
    try:
        return _replay_one(path)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError) as exc:
        return False, [f"malformed certificate: {exc}"], "not checked", math.nan


def _replay_one(path):
    """Re-reads the data and checks the certificate using only the matrices I - h A_i (no inverse, no shared code)."""
    C = common.read_json(path); msgs = []; ok = True
    raw = common.read_json(CASES[C["source"]], parse_float=str, parse_int=str)[C["case"]]
    h = Fr(C["h"]); w = common.word(C["word"], len(raw)); k = len(w); P = Fr(C["printed_exact"])
    if h <= 0:
        raise ValueError("h must be positive")
    if k != C["k"] or len(w) == 0:
        ok = False; msgs.append("word length mismatch")
    if P != Fr(C["printed"]):
        ok = False; msgs.append("printed value is not the stored exact value")
    for ver, conv in (("decimal", lambda s: Fr(s)), ("binary64", lambda s: Fr(float(s)))):
        As = [[[conv(x) for x in row] for row in A] for A in raw]; d = len(As[0])
        common.check_family(As)
        Ms = [[[(1 if r == c else 0) - h * A[r][c] for c in range(d)] for r in range(d)] for A in As]
        V = C["versions"][ver]
        if any(M[r][c] > 0 for M in Ms for r in range(d) for c in range(d) if r != c):
            ok = False; msgs.append(f"{ver}: I - hA_i is not a Z-matrix")
        if len(V["z"]) != len(Ms):
            raise ValueError(f"{ver}: need one M-matrix witness per mode")
        for i, (M, zs) in enumerate(zip(Ms, V["z"])):
            z = common.vector(zs, d)
            if not (all(t > 0 for t in z) and all(sum(M[r][c] * z[c] for c in range(d)) > 0 for r in range(d))):
                ok = False; msgs.append(f"{ver}: mode {i}: no strictly positive M-matrix witness")
        xs = [common.vector(x, d) for x in V["chain"]]
        if len(xs) != k + 1 or not all(t > 0 for t in xs[0]):
            raise ValueError(f"{ver}: chain length or x_0 > 0 fails")
        for j in range(1, len(xs)):
            M = Ms[w[j - 1]]
            if any(sum(M[r][c] * xs[j][c] for c in range(d)) < xs[j - 1][r] for r in range(d)):
                ok = False; msgs.append(f"{ver}: chain inequality fails at step {j}"); break
        q = max(xs[-1][r] / xs[0][r] for r in range(d)); r_ = Fr(V["r"])
        if r_ <= 0:
            raise ValueError(f"{ver}: r must be positive")
        if not r_ ** k >= q:
            ok = False; msgs.append(f"{ver}: r^k < q")
        if not (1 - 1 / r_) / h <= P:
            ok = False; msgs.append(f"{ver}: (1 - 1/r)/h exceeds the printed value")
    # Optional diagnostic: a changed letter need not invalidate a chain (e.g. identical modes).
    neg = "not applicable (one mode)"
    if len(raw) >= 2:
        As = [[[Fr(float(x)) for x in row] for row in A] for A in raw]; d = len(As[0])
        M = [[(1 if r == c else 0) - h * As[(w[0] + 1) % len(raw)][r][c] for c in range(d)] for r in range(d)]
        xs = [[Fr(t) for t in x] for x in C["versions"]["binary64"]["chain"][:2]]
        flipped_fails = any(sum(M[r][c] * xs[1][c] for c in range(d)) < xs[0][r] for r in range(d))
        neg = "flipped first letter rejected" if flipped_fails else "flipped first letter NOT rejected"
    # tightness: float value of the word itself
    slack = math.nan
    if ok:
        try:
            Rf = [np.linalg.inv(np.eye(d) - float(h) * np.array(A)) for A in raw_float(raw)]
            Uw = (1 - math.exp(-word_value_independent(Rf, w))) / float(h)
            slack = float(P) - Uw
        except (np.linalg.LinAlgError, ValueError, OverflowError, ZeroDivisionError):
            pass  # Floating-point tightness is not part of the exact proof.
    return ok, msgs, neg, slack


def raw_float(raw):
    return [[[float(x) for x in row] for row in A] for A in raw]


def word_value_independent(R, w):
    lg = 0.0; P = np.eye(R[0].shape[0])
    for a in w:
        P = R[a] @ P; s = np.linalg.norm(P); P /= s; lg += math.log(s)
    return (math.log(max(abs(np.linalg.eigvals(P)))) + lg) / len(w)


def main(paths):
    allok = bool(paths)
    for p in paths:
        t0 = time.time(); ok, msgs, neg, slack = replay_one(p); allok &= ok
        C = common.read_json(p)
        print(f"{os.path.basename(p)}: {'VALID' if ok else 'INVALID'}  lambda_check <= {C['printed']} (k = {C['k']}; "
              f"printed - float value of the word = {slack:.1e}; {neg}; {time.time() - t0:.1f}s)", flush=True)
        for m in msgs:
            print(f"    {m}")
    print(f"replay: {'all certificates VALID' if allok else 'SOME CERTIFICATES INVALID'}")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
