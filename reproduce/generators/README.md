# Generators of the random pairs (documentation only)

The random pairs of Obs. 4.8 were drawn with numpy's `default_rng`, which MATLAB's random number
generators cannot reproduce. The package therefore ships the drawn pairs as data, and these two Python scripts are
included only to document how they were drawn. They are not needed to reproduce any number of the paper.

| script | seed | output (in this package) |
|---|---|---|
| `gen_families.py` | 2026 | `data/scaling_pairs.json` (matrices), `data/scaling_pairs_info.json` (id, kind, eta(A_i), c, alpha*) |
| `find_pairs.py` | 2027 | `data/switching_pairs.json` (matrices), `data/switching_pairs_info.json` (c, alpha*, best schedule found, t1, t2) |
| `gen_pairs_ext.py` | 2026, 2027 continued | `data/scaling_pairs_journal.json`, `data/switching_pairs_journal.json` and their `*_info.json`: the larger samples of Obs. 4.8 (Tables S5 and S6) |

`find_pairs.py` imports `sample`, `irreducible`, `chatter`, `eta` and `H_MAX` from `gen_families.py` (it expects
both files in the same folder). The matrices in `data/` are the generators' binary64 numbers, written as their
shortest decimal strings, so that each string reads back as the stored number.

Generator of `gen_families.py` (Obs. 4.8): dense pairs with off-diagonal entries uniform in `[0, 2/(d-1)]` and diagonal
uniform in `[-3, 1]`; sparse pairs with each off-diagonal entry nonzero with probability 0.5 (d = 5) or 0.3 (d = 10),
then log-normal with sigma = 1.5, scaled to expected row sum 2.46, and diagonal uniform in `[-4, 1]`. Both modes
irreducible; pairs with the chattering minimizer alpha* outside `[0.05, 0.95]` rejected. `find_pairs.py` (Obs. 4.8)
keeps the sparse pairs on which an optimized two-segment schedule beats c by at least 1e-3.

`gen_pairs_ext.py` (larger samples) replays both streams, checks that they start with the paper's pairs of
`data/scaling_pairs.json` and `data/switching_pairs.json`, and continues them: the first 30 qualifying switching draws
for d = 5 and 10 and the first 5 for d = 20, and seven more dense and seven more sparse chattering pairs for d = 5 and
10, plus three dense and two sparse for d = 20 (off-diagonal probability 0.15 for the sparse ones). It writes only the
new pairs. It imports `gen_families.py` and `find_pairs.py` from this folder.
