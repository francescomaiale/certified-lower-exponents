"""(Documentation copy in the paper-B MATLAB package; see README.md in this folder.)

Larger samples of Obs. 4.7 and 4.8 (journal version): the paper's random streams, continued.

Switching pairs (Obs. 4.8): find_pairs.py draws pair k of dimension d from the stream default_rng([2027, d, k]) and
keeps the qualifying draws in draw order, so the first ten per dimension are the paper's sw5_1..sw5_10 and
sw10_1..sw10_10.  Here: the first 30 qualifying draws for d = 5 and 10, and the first 5 for d = 20.

Chattering pairs (Obs. 4.7): gen_families.py draws from one stream default_rng(2026), three dense and three sparse pairs
for d = 5, then for d = 10.  Here the stream is replayed (the first twelve pairs are checked against
data/scaling_pairs.json) and continued: seven more pairs of each kind for d = 5 and d = 10 (ten of each kind, twenty per
dimension), then three dense and two sparse pairs for d = 20.  For every pair we also record whether an optimized
two-segment schedule beats the numerical estimate of c (evidence for switching, not a proof of mixture optimality).

Sparse pairs for d = 20 need P_SPARSE[20]; we take 0.15, which keeps the expected number of nonzero off-diagonal entries
per row near 3 (0.5 * 4 = 2 for d = 5, 0.3 * 9 = 2.7 for d = 10, 0.15 * 19 = 2.85).

Writes only the new pairs: data/switching_pairs_journal.json and data/scaling_pairs_journal.json (matrices), with
*_info.json (c, alpha*, best schedule found, t1, t2; for the chattering pairs also eta(A_i) and whether a schedule beats c).
Usage: python3 gen_pairs_ext.py [OUTDIR]   (default: ../data; about 20 minutes on 3 cores)
"""
import json, os, sys, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_families as gf          # noqa: E402
import find_pairs as fp            # noqa: E402

gf.P_SPARSE[20] = 0.15
TARGET = {5: 30, 10: 30, 20: 5}
MAX_DRAWS = {5: 1500, 10: 700, 20: 600}
DATA = os.path.join(HERE, "..", "data")


def draw(args):
    gf.P_SPARSE[20] = 0.15
    return fp.draw(args)


def switching():
    out, info = {}, []
    with Pool(3) as pool:
        for d in (5, 10, 20):
            T0 = time.time(); found = 0; k0 = 0
            while found < TARGET[d] and k0 < MAX_DRAWS[d]:
                batch = sorted(pool.map(draw, [(d, k) for k in range(k0, k0 + 60)]), key=lambda r: r["k"])
                for r in batch:
                    if r["gain"] >= fp.MARGIN and found < TARGET[d]:
                        found += 1; pid = f"sw{d}_{found}"
                        out[pid] = [r["A1"], r["A2"]]
                        info.append(dict(id=pid, d=d, k=r["k"], chatter=r["chatter"], alpha=r["alpha"], sched=r["sched"],
                                         t1=r["t1"], t2=r["t2"], gain=r["gain"]))
                k0 += 60
            print(f"switching d = {d}: {found} pairs in the first {k0} draws ({time.time() - T0:.0f} s)", flush=True)
    return out, info


def chattering():
    rng = np.random.default_rng(gf.SEED)
    out, info = {}, []
    plan = [(5, "dense", 3), (5, "sparse", 3), (10, "dense", 3), (10, "sparse", 3),          # the paper's twelve
            (5, "dense", 7), (5, "sparse", 7), (10, "dense", 7), (10, "sparse", 7),          # continued
            (20, "dense", 3), (20, "sparse", 2)]
    count = {}
    for d, kind, m in plan:
        n = 0
        while n < m:
            A1, A2 = gf.sample(rng, kind, d), gf.sample(rng, kind, d)
            if not (gf.irreducible(A1) and gf.irreducible(A2)):
                continue
            cA = max(np.max(-np.diag(A1)), np.max(-np.diag(A2)))
            e1, e2 = gf.eta(A1), gf.eta(A2)
            if gf.H_MAX * cA >= 1 or gf.H_MAX * max(e1, e2) >= 1:
                continue
            c, a = gf.chatter(A1, A2)
            if not 0.05 <= a <= 0.95:
                continue
            n += 1; count[(d, kind)] = count.get((d, kind), 0) + 1
            pid = f"d{d}_{kind}{count[(d, kind)]}"
            s, t1, t2 = fp.best_schedule(A1, A2)
            out[pid] = [A1.tolist(), A2.tolist()]
            info.append(dict(id=pid, d=d, kind=kind, eta=[e1, e2], chatter=c, alpha=a, sched=s, t1=t1, t2=t2,
                             switches=bool(c - s >= fp.MARGIN)))
            print(f"{pid}: c = {c:+.6f} at alpha = {a:.4f}; best schedule {s:+.6f}", flush=True)
    return out, info


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else DATA
    os.makedirs(outdir, exist_ok=True)
    sw, swi = switching()
    old = json.load(open(os.path.join(DATA, "switching_pairs.json")))
    assert all(sw[p] == old[p] for p in old), "the first ten switching pairs per dimension differ from the paper's"
    new = lambda p: p not in old
    json.dump({p: v for p, v in sw.items() if new(p)}, open(os.path.join(outdir, "switching_pairs_journal.json"), "w"))
    json.dump([r for r in swi if new(r["id"])], open(os.path.join(outdir, "switching_pairs_journal_info.json"), "w"), indent=1)
    ch, chi = chattering()
    old = json.load(open(os.path.join(DATA, "scaling_pairs.json")))
    assert all(ch[p] == old[p] for p in old), "the first twelve chattering pairs differ from the paper's"
    new = lambda p: p not in old
    json.dump({p: v for p, v in ch.items() if new(p)}, open(os.path.join(outdir, "scaling_pairs_journal.json"), "w"))
    json.dump([r for r in chi if new(r["id"])], open(os.path.join(outdir, "scaling_pairs_journal_info.json"), "w"), indent=1)
    print("done", flush=True)


if __name__ == "__main__":
    main()
