"""Independent Python replay of the invariant-face certificates of Computation 3.8 (certificates/comp38/faces_*.json,
written by the MATLAB function invariant_faces).  New code (the paper's Python 'faces' subcommand printed a log and had
no separate replay); it shares nothing with the MATLAB certifier.

For each orientation, with the binary64 data (transposed for 'T'): the bound for c (x > 0, (S_alpha x)_r <= U x_r), the
list of invariant faces J with their mode sets S(J) recomputed from the zero pattern, and for single-mode faces
x >= 0, x != 0, A_i[J,J] x >= L x (so eta(A_i[J,J]) >= L); the slow flags are L > U.

Usage: python3 replay_faces.py FILE [FILE ...]
"""
import os, sys, json, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common

CASES = common.merged_cases_file()


def as_list(x):
    return x if isinstance(x, list) else [x]


def replay_file(path):
    try:
        return _replay_file(path)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError) as exc:
        print(f"replay {os.path.basename(path)}: INVALID ({exc})", flush=True)
        return False


def _replay_file(path):
    C = common.read_json(path); ok = True; msgs = []
    raw = common.read_json(CASES, parse_float=str, parse_int=str)[C["case"]]
    if not C["orientations"]:
        raise ValueError("no orientations supplied")
    for F in as_list(C["orientations"]):
        As = [[[Fr(float(x)) for x in row] for row in A] for A in raw]; d = len(As[0])
        common.check_family(As)
        if len(As) != 2 or F["orientation"] not in ("A", "N", "T"):
            raise ValueError("expected a pair and orientation A, N or T")
        if F["orientation"] == "T":
            As = [[[A[c][r] for c in range(d)] for r in range(d)] for A in As]
        a = Fr(F["alpha"]); x = [Fr(v) for v in as_list(F["x_c"])]; U = Fr(F["c_upper"])
        if not 0 <= a <= 1 or len(x) != d:
            raise ValueError("invalid mixture weight or witness dimension")
        S = [[a * As[0][r][c] + (1 - a) * As[1][r][c] for c in range(d)] for r in range(d)]
        if not (all(v > 0 for v in x) and all(sum(S[r][c] * x[c] for c in range(d)) <= U * x[r] for r in range(d))):
            ok = False; msgs.append(f"{F['orientation']}: bound for c fails")
        want = set()
        for k in range(1, d):
            for J in itertools.combinations(range(d), k):
                out = [r for r in range(d) if r not in J]
                Si = tuple(i + 1 for i, A in enumerate(As) if all(A[r][j] == 0 for r in out for j in J))
                if Si:
                    want.add((tuple(j + 1 for j in J), Si))
        faces = as_list(F["faces"]) if F["faces"] else []
        have = {(tuple(as_list(f["J"])), tuple(as_list(f["S"]))) for f in faces}
        if want != have:
            raise ValueError(f"{F['orientation']}: invariant faces differ: {sorted(want)} vs {sorted(have)}")
        if len(have) != len(faces):
            raise ValueError("duplicate face")
        allslow = True
        for f in faces:
            J = [j - 1 for j in as_list(f["J"])]; Si = as_list(f["S"]); slow = False
            if len(Si) == 1:
                A = As[Si[0] - 1]; M = [[A[r][c] for c in J] for r in J]
                xj = [Fr(v) for v in as_list(f["x"])]; L = Fr(f["eta_lower"]); n = len(J)
                good = (len(xj) == n and all(v >= 0 for v in xj) and any(v > 0 for v in xj)
                        and all(sum(M[r][c] * xj[c] for c in range(n)) >= L * xj[r] for r in range(n)))
                if not good:
                    ok = False; msgs.append(f"{F['orientation']}: face {as_list(f['J'])}: lower bound fails")
                slow = L > U
            if slow != bool(f["slow"]):
                ok = False; msgs.append(f"{F['orientation']}: face {as_list(f['J'])}: slow flag inconsistent")
            allslow = allslow and slow
        if allslow != bool(F["all_slow"]):
            ok = False; msgs.append(f"{F['orientation']}: all_slow inconsistent")
    print(f"replay faces {C['case']}: {'VALID' if ok else 'INVALID'} {'; '.join(msgs)}", flush=True)
    return ok


def main(paths):
    ok = bool(paths) and all([replay_file(p) for p in paths])
    print("replay:", "all VALID" if ok else "SOME INVALID")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
