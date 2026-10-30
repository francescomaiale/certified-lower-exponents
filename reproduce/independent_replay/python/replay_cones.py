"""Independent Python replay of the invariant cones of Computation 3.9: the 'replay' subcommands of euler_cone.py
(polygons of cases 77 and 677) and euler_cone_boundary.py (the cut cone of case 148), merged and applied to the
certificate files written by the MATLAB certifier (fields case, kind, orientation, h1, vertices, classes).

Exact data (binary64 values), polygon convexity by turn signs, and every image B_i v located in a triangle of the
fan from the first vertex by exact barycentric coordinates.  Generators must be > 0 for kind 'polygon' (prop:cone(ii))
and >= 0 otherwise; the zero patterns of the generators are printed with the stored support classes.

Usage: python3 replay_cones.py FILE [FILE ...]
"""
import os, sys, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common

CASES = common.merged_cases_file()


def replay_file(path):
    try:
        return _replay_file(path)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError) as exc:
        print(f"replay {os.path.basename(path)}: INVALID ({exc})", flush=True)
        return False


def _replay_file(path):
    C = common.read_json(path)
    raw = common.read_json(CASES, parse_float=str, parse_int=str)[C["case"]]
    As = [[[Fr(float(x)) for x in row] for row in A] for A in raw]
    common.check_family(As)
    if len(As[0]) != 3 or C["orientation"] not in ("A", "N", "T"):
        raise ValueError("expected 3-dimensional data and orientation A, N or T")
    if C["orientation"] == "T":
        As = [[[A[c][r] for c in range(3)] for r in range(3)] for A in As]
    h = Fr(C["h1"]); V = [tuple(Fr(x) for x in v) for v in C["vertices"]]
    if h <= 0 or len(V) < 3 or any(len(v) != 3 or any(x < 0 for x in v) or sum(v) <= 0 for v in V):
        raise ValueError("invalid step or cone generators")
    B = [[[Fr(int(i == j)) + h * A[i][j] for j in range(3)] for i in range(3)] for A in As]
    Q = [(v[0] / sum(v), v[1] / sum(v)) for v in V]
    turn = lambda p, q, r: (q[0]-p[0]) * (r[1]-q[1]) - (q[1]-p[1]) * (r[0]-q[0])
    sg = [turn(Q[k], Q[(k+1) % len(Q)], Q[(k+2) % len(Q)]) for k in range(len(Q))]
    convex = all(s >= 0 for s in sg) or all(s <= 0 for s in sg)

    def in_fan(p):
        for k in range(1, len(Q) - 1):
            a, b, c = Q[0], Q[k], Q[k + 1]
            det = (b[0]-a[0]) * (c[1]-a[1]) - (c[0]-a[0]) * (b[1]-a[1])
            if det == 0:
                continue
            l1 = ((p[0]-a[0]) * (c[1]-a[1]) - (c[0]-a[0]) * (p[1]-a[1])) / det
            l2 = ((b[0]-a[0]) * (p[1]-a[1]) - (p[0]-a[0]) * (b[1]-a[1])) / det
            if l1 >= 0 and l2 >= 0 and l1 + l2 <= 1:
                return True
        return False
    imgs = [tuple(sum(Bi[r][c] * v[c] for c in range(3)) for r in range(3)) for Bi in B for v in V]
    if any(any(x < 0 for x in y) or sum(y) <= 0 for y in imgs):
        raise ValueError("an image is not a nonzero nonnegative vector")
    sign_ok = all(x > 0 for v in V for x in v) if C["kind"] == "polygon" else all(x >= 0 for v in V for x in v)
    ok = convex and sign_ok and all(in_fan((y[0] / sum(y), y[1] / sum(y))) for y in imgs)
    supports = {frozenset(i + 1 for i, x in enumerate(v) if x > 0) for v in V}
    classes = set()
    for mask in range(1, 8):
        J = frozenset(i + 1 for i in range(3) if mask & (1 << i))
        if frozenset().union(*(s for s in supports if s <= J)) == J:
            classes.add(",".join(map(str, sorted(J))))
    stored = C["classes"] if isinstance(C["classes"], list) else [C["classes"]]
    ok &= classes == set(stored) and len(stored) == len(classes)
    zero = sorted({tuple(k + 1 for k in range(3) if v[k] == 0) for v in V})
    print(f"replay {C['case']} ({C['kind']}, {C['orientation']}): {len(V)} vertices, {len(imgs)} images, "
          f"{'VALID' if ok else 'INVALID'}; zero patterns of vertices {zero}; stored classes {C['classes']}", flush=True)
    return ok


def main(paths):
    ok = bool(paths) and all([replay_file(p) for p in paths])
    print("replay:", "all VALID" if ok else "SOME INVALID")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
