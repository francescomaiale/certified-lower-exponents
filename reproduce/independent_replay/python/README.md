# Independent Python replays (cross-language check)

Independent replays derived from the paper's Python scripts, with explicit proof-prerequisite and structure checks, applied to the certificates written by the MATLAB
certifiers. They share no code with MATLAB; they read the pairs from `../../data/*.json` themselves.

| script | replays | origin |
|---|---|---|
| `replay_lower.py` | polytope certificates (Prop. 4.3) | `replay_one` of `certify_lower.py`, with validation |
| `replay_upper.py` | implicit Euler chains (Prop. 4.4) | `replay_one` of `certify_upper.py`, with validation |
| `replay_best_upper.py` | chattering and schedule bounds | `replay` of `certify_best_upper.py` (uniformization, no logarithms), extended to schedules with more than two segments |
| `replay_cones.py` | invariant cones of Computation 3.9 | `replay` of `euler_cone.py` and `euler_cone_boundary.py` |
| `replay_faces.py` | invariant faces of Computation 3.9 | new (the paper's script printed a log only) |
| `run_all.py` | every certificate under `../../certificates` | dispatches by the fields of each file |

Requirements: Python 3 with NumPy (tested with Python 3.12 and the version pinned in the root `requirements.txt`). Usage:

```bash
python3 run_all.py
```

```bash
python3 replay_lower.py ../../certificates/table2/case1_0.1_best.json
```

Run the rejection and scalar-counterexample tests from this directory:

```bash
python3 -m unittest test_replays -v
```

The checkers require positive steps, complete witnesses, valid mixture weights and valid schedules. Unknown formats,
missing folders and empty selections fail. Legacy implicit words are digit strings; integer arrays also support mode
indices above nine. Negative-control and floating-point word-tightness messages in the upper replay are diagnostics;
the exact chain inequalities and the positive M-matrix witnesses establish the bound.

Set `OMP_NUM_THREADS=1` (and the OpenBLAS, MKL and vecLib equivalents) when several replays run in parallel.
