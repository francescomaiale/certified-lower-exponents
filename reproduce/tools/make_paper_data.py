"""Generate ../../paper_data.m (the data of example_paper.m) from data/ and certificates/ of this folder.

Run after any change of the certificates behind Tables 2 and 3:  python3 tools/make_paper_data.py
"""
import json, os, re
PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # reproduce/
OUT = os.path.join(os.path.dirname(PKG), "paper_data.m")
raw = lambda f: json.load(open(os.path.join(PKG, f)), parse_float=str, parse_int=str)
cert = lambda f: json.load(open(os.path.join(PKG, "certificates", f)))
test, sw = raw("data/test_pairs.json"), raw("data/switching_pairs.json")

def mat(M):
    rows = [" ".join(f'"{x}"' for x in r) for r in M]
    return "[" + "\n     ".join(rows) + "]"

L = []
w = L.append
w("function P = paper_data()")
w("% PAPER_DATA  The pairs of paper B and the recorded settings needed to rebuild its certified bounds.")
w("%   P.pairs.<name>  {A_1, A_2}: string matrices of the decimals printed in the paper (App. A for the test pairs;")
w("%                   the switching pairs of Table 3 are the shortest decimals of their binary64 entries).")
w("%   P.table2        rows of Table 2: pair, h, level (the level at which the paper's search closed), seeds,")
w("%                   from (a larger step whose polytope is transferred, prop:polytope(ii)), upper_word (the")
w("%                   implicit-Euler word of the upper bound), and the printed lower and upper bounds.")
w("%   P.table3        initial closures for Table 3, before the extended refinements marked +:")
w("%                   pair, h, orientation ('N' = B_h, 'T' = transposed family), word [a b] of")
w("%                   the Perron seeds, level, eps, printed lower bound.")
w("%   P.schedules     per switching pair: the two durations of the certified schedule and the printed bound.")
w("%   The levels are floating-point starting points only: every bound is re-certified in exact arithmetic.")
w("")
w("P = struct();")
w("P.pairs = struct();")
names_test = ["planar", "case1", "case1440", "case148"]
for n in names_test + sorted(sw, key=lambda s: (int(s[2:].split('_')[0]), int(s.split('_')[1]))):
    M = test[n] if n in test else sw[n]
    w(f"P.pairs.{n} = {{ ...")
    w(",  ...\n".join(f"    {mat(A)}" for A in M) + "};")
w("")
# ---------------- Table 2
P2 = [("planar","0.1"),("planar","0.05"),("planar","0.025"),("planar","0.0125"),("case1","0.1"),("case1","0.05"),
      ("case1","0.025"),("case1440","0.1"),("case1440","0.05"),("case1440","0.025"),("case148","0.1"),
      ("case148","0.05"),("case148","0.025")]
transfer = {"case148_0.05": "0.1", "case148_0.025": "0.1", "case1440_0.025": "0.05"}
w("% Table 2: pair, h, level, seeds, from, upper_word, lower, upper")
w("T2 = { ...")
for p, h in P2:
    row = f"{p}_{h}"
    U = cert(f"table2/{row}.json")
    if row in transfer:
        src = cert(f"table2/{p}_{transfer[row]}_best.json")
        lev, frm, low = "NaN", transfer[row], src["lambda_lower"]
    else:
        C = cert(f"table2/{row}_best.json")
        lev, frm, low = repr(C["meta"]["t_build"]), "", C["lambda_lower"]
    seeds = "e1e3" if p == "case148" else "eye"
    w(f"    '{p}', '{h}', {lev}, '{seeds}', '{frm}', '{U['word']}', '{low}', '{U['printed']}'; ...")
w("    };")
w("P.table2 = cell2struct(T2, {'pair', 'h', 'level', 'seeds', 'from', 'upper_word', 'lower', 'upper'}, 2);")
w("")
# ---------------- Table 3
w("% Table 3: pair, h, orientation, word, level, eps, lower")
w("T3 = { ...")
for d in (5, 10):
    for k in range(1, 11):
        p = f"sw{d}_{k}"
        for h in ("0.1", "0.05"):
            row = f"{p}_{h}"
            cands = []
            f = os.path.join(PKG, "certificates/switching", f"{row}_best.json")
            if os.path.exists(f): cands.append(json.load(open(f)))
            for f in (os.path.join(PKG, "certificates/table3_transposed", f"{row}_T.json"),):
                if os.path.exists(f): cands.append(json.load(open(f)))
            if not cands:
                w(f"    '{p}', '{h}', '', [0 0], NaN, NaN, '---'; ...")
                continue
            C = max(cands, key=lambda c: float(c["lambda_lower"]))
            m = C["meta"]; wd = m["word"]
            w(f"    '{p}', '{h}', '{C['orientation']}', [{wd[0]} {wd[1]}], {repr(m['t_build'])}, {m['eps']!r}, '{C['lambda_lower']}'; ...")
w("    };")
w("P.table3 = cell2struct(T3, {'pair', 'h', 'orientation', 'word', 'level', 'eps', 'lower'}, 2);")
w("")
w("% Schedules of Table 3 (upper bounds): pair, t1, t2 (A_1 for t1, then A_2 for t2), printed bound")
w("SC = { ...")
for d in (5, 10):
    for k in range(1, 11):
        p = f"sw{d}_{k}"
        S = cert(f"switching/{p}_schedule.json")
        t1, t2 = S["durations"]
        w(f"    '{p}', '{t1}', '{t2}', '{S['printed']}'; ...")
w("    };")
w("P.schedules = cell2struct(SC, {'pair', 't1', 't2', 'upper'}, 2);")
w("end")
open(OUT, "w").write("\n".join(L) + "\n")
print("wrote", OUT, len(L), "lines")
