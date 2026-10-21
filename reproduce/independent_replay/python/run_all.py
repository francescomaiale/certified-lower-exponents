"""Cross-language check: replay every certificate under ../../certificates with the independent Python replays.

Each file is classified by its fields: polytope certificates ('weights') -> replay_lower, implicit Euler chains
('versions') -> replay_upper, chattering/schedule bounds ('readings') -> replay_best_upper, invariant cones ('h1') ->
replay_cones, invariant faces ('orientations') -> replay_faces.  Requires Python 3 with numpy.

Usage: python3 run_all.py [FOLDER ...]       (default: every subfolder of ../../certificates)
       Set OMP_NUM_THREADS=1 etc. when running several replays in parallel.
"""
import sys, os, json, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # the replays run with python3 -I
import common
import replay_lower, replay_upper, replay_best_upper, replay_cones, replay_faces


def kind(C):
    if "weights" in C:
        return "lower"
    if "versions" in C:
        return "upper"
    if "readings" in C:
        return "best_upper"
    if "h1" in C and "vertices" in C:
        return "cone"
    if "orientations" in C:
        return "faces"
    return None


def main(folders):
    if not folders:
        folders = sorted(d for d in glob.glob(os.path.join(common.ROOT, "certificates", "*")) if os.path.isdir(d))
    files = []
    for d in folders:
        found = glob.glob(os.path.join(d, "*.json"))
        if not found:
            print(f"INVALID selection: no certificates in {d}", flush=True)
            return False
        files.extend(found)
    files = sorted(set(files))
    if not files:
        print("INVALID selection: no certificates found", flush=True)
        return False
    groups = {}
    for f in files:
        try:
            k = kind(common.read_json(f))
        except (ValueError, OSError) as exc:
            print(f"INVALID certificate {f}: {exc}", flush=True)
            return False
        if k is None:
            print(f"INVALID certificate {f}: unknown format", flush=True)
            return False
        groups.setdefault(k, []).append(f)
    T0 = time.time(); allok = True; counts = {}
    for k, mod in (("lower", replay_lower), ("upper", replay_upper), ("best_upper", replay_best_upper),
                   ("cone", replay_cones), ("faces", replay_faces)):
        fs = groups.get(k, [])
        if not fs:
            continue
        print(f"==== {k}: {len(fs)} certificates", flush=True)
        ok = mod.main(fs); allok &= ok; counts[k] = (len(fs), ok)
    print(f"==== Python replay of {sum(n for n, _ in counts.values())} certificates: "
          + ", ".join(f"{k} {n} {'VALID' if ok else 'SOME INVALID'}" for k, (n, ok) in counts.items())
          + f"; {time.time() - T0:.0f}s")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:]) else 1)
