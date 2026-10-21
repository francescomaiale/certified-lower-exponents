# Reproducing the computations of the paper

This folder contains MATLAB drivers, data and certificates for the computations of the paper

> *Euler discretizations, invariant faces and certified brackets for lower exponents of positive switched systems*

(called "the paper" below), except the runs of the algorithm of Guglielmi and Maiale in Appendix B (see the
last section). The original timing wrappers and some final GM seed sets are also absent; the supplied certificates can be checked without them. The code itself is in `../src`. For a family of Metzler matrices `A_1, ..., A_m`, the paper brackets the lower Lyapunov
exponent `lambda_check` of `x' = A_sigma x` by the lower spectral radii of the explicit and implicit Euler families,
`g(h) <= lambda_check <= G(h)`, and certifies every endpoint of the reported continuous-time brackets in exact rational arithmetic. Explicit-word
estimates, mixture minimizations and rate diagnostics remain floating point, as do the Algorithm (L) upper bounds:


* lower bounds by polytopes, `B_i v >= t * sum_u lambda_u u` for every vertex `v` and mode `i` (Prop. 4.3), also at
  every smaller step (Prop. 4.3(ii));
* upper bounds by implicit Euler chains, `(I - h A_{w_j}) x_j >= x_{j-1}` (Prop. 4.4), by an upper bound for the chattering value
  (Collatz-Wielandt at a rational `alpha`) and by periodic schedules (entrywise enclosures of `e^{tA}`);
* the invariant faces and invariant cones of Computation 3.8.

The replayable certificates are written as JSON files and re-checked by a **replay that shares no code with the certifier**
(`replay/`), and by the **independent Python replays** derived from the paper's scripts (`independent_replay/python/`).

## Requirements

* MATLAB R2023a or later (uses `pageeig`); developed and tested with R2026b.
* Optimization Toolbox (`linprog`: `dual-simplex-highs`, with `interior-point` and `interior-point-legacy` as
  fallbacks; R2023a/b uses `dual-simplex` in place of HiGHS). Only R2026b has been exercised in the release review.
* Symbolic Math Toolbox (exact rational arithmetic, `sym`).
* Optional: Python 3 with numpy, for the cross-language replays (`independent_replay/python`, test group (f)).

No other toolbox is needed. The paper used single-threaded runs. Start MATLAB with `-singleCompThread` to reproduce that setting, especially if running several drivers in parallel.

## Quick start

```matlab
cd reproduce
startup_reproduce                 % adds ../src, helpers/, replay/, tables/ and tests/ to the path
run_tests                         % end-to-end and input-validation checks, a few minutes
table2_test_pairs                 % Table 2: certified brackets for the test pairs
r = certified_bracket(load_pair('test', 'case1'), '0.1');   % a certified bracket for any pair and step
```

If MATLAB finds a different Python installation from your shell, set
`run_tests(struct('python_executable', '/path/to/python3'))` to the interpreter in which `requirements.txt` was installed.

From a shell:

```bash
matlab -singleCompThread -batch "startup_reproduce; run_tests"
```

```bash
matlab -singleCompThread -batch "startup_reproduce; run_all('quick')"
```

`run_all('quick')` and `run_all('full')` run the original suite with different workloads. Invoke
`table3_transposed`, `table3_extended` and `larger_samples` separately for the additional experiments; the rate
search also runs in the quick profile and can take about 30 minutes. `summarize_all` refreshes the original bracket
summaries and work-precision figure; extended and larger-sample summaries have their own drivers. `replay_all`
re-checks every certificate. Most certificate-producing drivers skip files that already exist, so runs can be interrupted and
resumed, and several drivers (or row subsets of one driver, option `'rows'`/`'pairs'`) can run in parallel MATLAB
sessions. Summaries go to `results/*.md`, logs to `results/logs/`.

## Layout

| folder | content |
|---|---|
| `../src` | the code: `core`, `lower`, `upper`, `cone`, `glp`, `pipeline` (see the README one level up) |
| `helpers` | data loading (`load_pair`), the recorded search levels (`paper_levels`), `run_lower_row` and `run_upper_row` (one certificate file per row, with its replay), logging |
| `replay` | independent exact replays: `replay_lower`, `replay_upper`, `replay_best_upper`, `replay_faces`, `replay_cone`, and `replay_all` (helpers in `replay/private`, data re-read from `data/`) |
| `tables` | one driver per table, figure and observation (below), `run_all` and `summarize_all` |
| `tests` | `run_tests`, `run_validation_tests` and fixtures (certificates produced by the paper's Python scripts) |
| `data` | the pairs, as JSON (below) |
| `certificates` | the 568 certificates of the paper's bounds, one folder per driver (`appB_*`: Appendix B, see the last section; `table3_extended`, `larger_*`: the extended search and the larger samples) |
| `results` | summaries (`*.md`), Figure 1, and the run outputs of Algorithm (L), Figure 1, the transposed family, the extended search (`table3_extended/`) and the larger samples (`larger_samples/`, with `MANIFEST.txt`) |
| `independent_replay/python` | Independent replays derived from the paper's Python code, with certificate validation, run on the MATLAB certificates |
| `generators` | the Python generators of the random pairs (documentation only) |
| `tools` | `make_paper_data.py`, which writes `../paper_data.m` from `data/` and `certificates/`; `check_rate_theorem.py`, the floating-point checks of Thm. 3.12 (App. E) |

## Data (`data/`)

| file | content |
|---|---|
| `test_pairs.json` | `planar` (Obs. 3.11), `sharp` (Ex. 3.9 with eps = 0.2), `case1`, `case77`, `case148`, `case677`, `case1195`, `case1440` (Appendix A, digit for digit) |
| `scaling_pairs.json`, `scaling_pairs_info.json` | the twelve random pairs of Obs. 4.7 (`d5_dense1`, ..., `d10_sparse3`), and their generator outputs (eta(A_i), c, alpha*) |
| `switching_pairs.json`, `switching_pairs_info.json` | the twenty pairs of Obs. 4.8 (`sw5_1`, ..., `sw10_10`), and c, alpha* and the best two-segment schedule found |
| `switching_pairs_journal.json`, `switching_pairs_journal_info.json` | the larger sample of Obs. 4.8 (Table 10): `sw5_11`, ..., `sw5_30`, `sw10_11`, ..., `sw10_30`, `sw20_1`, ..., `sw20_5` (load_pair shortcut `'switching_journal'`) |
| `scaling_pairs_journal.json`, `scaling_pairs_journal_info.json` | the larger sample of Obs. 4.7 (Table 9): `d5_dense4`, ..., `d10_sparse10`, `d20_dense1`, ..., `d20_sparse2` (shortcut `'scaling_journal'`); the info also records whether a two-segment schedule beats the numerical value of c |
| `bcv13.json` | the three-mode example of Ex. 4.10 |
| `paper_search_results.json` | outputs of the paper's floating-point searches (see "Reproduce and search modes") |

The matrices are stored as JSON numbers, exactly as printed: six-decimal data for the planar pair and case 1, and
the shortest decimal strings of the binary64 entries otherwise. `load_pair` reads them as strings (not doubles), and
`exact_data` forms the two readings used by every certificate: each decimal read as an exact rational, and its exact
binary64 value. MATLAB's random number generators cannot reproduce numpy's draws, so the random pairs are shipped as
data; the generators are in `generators/` for documentation.

## Reproduce and search modes

The paper's lower bounds come from time-budgeted bisections (300 s per closure, 1800 s per row), whose outcome
depends on the machine. The main bracket drivers have two modes.

* `'mode', 'reproduce'` (default): every lower bound is rebuilt by **one** closure at the level where the paper's
  bisection last closed (or, for Algorithm (L), at the recorded shift `nu`), stored in
  `data/paper_search_results.json`. The closure removes the original bisection time budget; vertex counts and final rounding can still depend on the LP solver.
  These recorded levels are only starting points: the certificate is re-derived and checked in exact arithmetic.
* `'mode', 'search'`: the bisection itself, with machine-independent caps (`'pmax'`: vertices processed per closure;
  `'maxprobes'`: bisection steps). Its results may differ from the paper's, since the budgets differ.

Existing certificates are skipped by default. To rebuild a certificate, use a fresh `cdir` or the driver's `force` option. Upper-word searches are rerun when the corresponding certificates are rebuilt.

The extended search and the larger samples (`table3_extended`, `larger_samples`) were run once, with the time budgets
of the paper, and have no reproduce mode: their `'summary'` part rebuilds the printed numbers from the stored
certificates, which the replays check, and their run parts repeat the time-budgeted searches.

## Paper items and drivers

Runtimes are single-threaded, on a 16-core Apple-silicon Mac under load, in reproduce mode.

| paper item | driver | output | runtime | status |
|---|---|---|---|---|
| Table 1, Obs. 3.11 | `table1_sandwich` | `results/table1.md` | 1 min | all 10 rows |
| Computation 3.8 | `comp38_reducible_pairs` | `certificates/comp38`, `results/comp38.md` | 1 min | all statements (see below) |
| Table 2, Rem. 4.5, Obs. 4.6 | `table2_test_pairs` | `certificates/table2`, `results/table2.md` | 15 min (+10 min for case 1 at h = 0.025) | all 13 rows |
| Table 4 (ours), Obs. 4.11 | `table4_best_brackets` | `certificates/table4`, `results/table4.md` | 1 min (+ Table 2 and the scaling driver) | our columns: all 16 pairs |
| Table 4 (Alg. (L)), Obs. 4.11 | `table4_algorithm_l` | `results/glp`, `results/algorithm_l.md` | hours | all printed cells (7 of the 29 runs differ from the earlier Python runs; see below) |
| Obs. 4.7, Table 8 | `table7_scaling` | `certificates/scaling`, `results/table7.md` | hours | all 24 rows |
| Obs. 4.8, Table 3 (B_h) | `table3_switching` | `certificates/switching`, `results/table3.md` | hours | all 20 rows |
| Obs. 4.8, Table 3 (transposed family) | `table3_transposed` | `certificates/table3_transposed`, `results/table3_transposed.md` | 2 hours (7 processes) | all 23 rows |
| Obs. 4.8, Table 3 (extended search, marked +) | `table3_extended` | `certificates/table3_extended`, `results/table3_extended.md` | summary 1 min; runs up to 3.5 hours per row | all 8 bounds marked + and the d = 10 widths |
| Obs. 4.7 and 4.8, Tables 9 and 10 (larger samples) | `larger_samples` | `certificates/larger_switching`, `certificates/larger_scaling`, `results/larger_samples.md`, `results/larger_samples_tables.tex` | summary 1 min; runs about 14 hours (12 to 16 processes) | all 78 rows |
| Obs. 4.9, Table 11 | `table8_rates` | `results/table8.md` | 30 min | all 14 rows |
| Thm. 3.12, Rem. 3.13, App. E (floating-point checks) | `tools/check_rate_theorem.py` (Python, numpy and scipy) | `results/check_rate_theorem.md` | 2 min | Q for all 14 pairs of Table 11; the remainder of Lemma E.1; the identity of Rem. 3.13 |
| Ex. 4.10 | `ex410_bcv` | `certificates/bcv`, `results/ex410.md` | 10 min | all values |
| Obs. 4.12, Fig. 1 | `fig1_work_precision` | `results/wp`, `results/fig1_wp.pdf`, `results/fig1.md` | hours | our 12 points; Alg. (L): 9 of 12 points (see below) |


"Status" counts the paper's printed numbers that the MATLAB run reproduces digit for digit (`results/*.md` compare
every cell with the paper). All certificates written by the drivers pass both the MATLAB replays (`replay_all`) and the
Python replays (`independent_replay/python/run_all.py`).

### What each driver reproduces

* **Table 1** (`table1_sandwich`): g_est, G_up and the exact-schedule rate for both planar pairs and all five steps;
  Obs. 3.11: c = -0.1949763083 at alpha* = 0.3905 and the predicted slope 1.478.
* **Computation 3.8** (`comp38_reducible_pairs`): for cases 77, 148, 677, 1195 and 1440, in both orientations, the
  invariant proper faces with exact exponent bounds and the certified bound for c (`faces_<case>.json`); the cone
  of case 148 (transposed, h1 = 1/5, kappa = 1/20) with its support classes (all nonempty subsets but {3}), touching
  only the slow face {1,3}; the invariant polygons of cases 77 (21 vertices, coordinates >= 0.0796, as in the paper)
  and 677 (see the differences below).
* **Table 2, Rem. 4.5, Obs. 4.6** (`table2_test_pairs`): every lower and upper bound and every width; the transfer
  certificates of Prop. 4.3(ii) (case 148 at h = 0.05, 0.025 from the h = 0.1 polytope; case 1440 at h = 0.025 from
  the h = 0.05 polytope); g(0.1) in [-1.044234, -1.044232] for case 1; the case-148 polytope has 1836 vertices,
  none on the axis of e_2; lambda_min <= (A_2)_22 < -1.368009. The polytopes have 169 to 3131 vertices.
* **Table 4** (`table4_best_brackets`, `table4_algorithm_l`): our lower bounds, the implicit Euler bounds and the
  certified upper bound for c or for a two-segment schedule (c for 14 pairs, schedules for cases 1440 and 148), and the
  widths; the tangent bound on our planar and case-1 polytopes. The Alg. (L) lower bounds are exact (both data
  readings) but, as in the paper, they are not stored as certificates and not replayed independently.
* **Obs. 4.7, Table 8** (`table7_scaling`): the 24 lower and 24 upper certificates for the twelve random pairs
  (d = 5, 10; h = 0.1, 0.05).
* **Obs. 4.8, Table 3** (`table3_switching`, `table3_transposed`): the twenty certified schedules, the
  Perron-seeded lower certificates on B_h (33 rows) and on the transposed family (23 rows: those that closed only at
  eps = 0.1 or not at all on B_h), the floating-point bounds g(h) <= of the best words 1^a 2^b, the widths, and the
  counts of levels that closed. Table 3 of the paper prints, for each row, the better of the two orientations
  (`results/table3_transposed.md`, column "better").
* **Obs. 4.8, Table 3, extended search** (`table3_extended`): the eight d = 10 bounds marked + (six at h = 0.1, two at
  h = 0.05) and the resulting widths, from the 34 runs of the extended search (records in `results/table3_extended/`;
  the certificates of the printed bounds in `certificates/table3_extended`). Its seeds include, on the transposed family
  at h = 0.1, the final vertex sets of the runs of Appendix B, which this package does not contain: pass `'gm_dir'`
  to use them in a new run.
* **Obs. 4.7 and 4.8, Tables 9 and 10** (`larger_samples`): for the 45 further switching pairs, the certified
  schedules and the lower bounds at both steps (both orientations, and the extended search on the d = 10 rows); for
  the 33 further chattering pairs, the certified upper bound for c (a schedule for `d10_sparse8`), the lower bounds
  and the implicit Euler bounds at both steps; the counts and widths quoted in Obs. 4.7 and 4.8. Part `'summary'`
  rebuilds both tables, as LaTeX rows identical to the paper's, from the stored certificates; parts `'switching'`,
  `'chattering'` and `'extended'` rerun the computations (into `results/larger_samples/run/`).
* **Obs. 4.9, Table 11** (`table8_rates`): alpha*, Q and the observed gap and separate slopes for all 14 pairs; the deviations and the
  midpoint accuracy (0.83 h^2); the extrapolation remark for switching pair 4 (-0.4270 above -0.4301).
  `results/rate_search_results.json` retains the original unrounded Python search diagnostics used for the separate
  slopes in the revised table. Its alpha*, Q and gap columns agree with the MATLAB summary to all printed digits.
  Word searches estimate the discrete minima from above; neither the gap nor its midpoint proves optimality.
* **Ex. 4.10** (`ex410_bcv`): a numerical mixture search returning approximately 1/120; the certified schedule bound lambda_check <= -0.299986 (dwell times
  30 in A_3, A_2, A_1); g(h) >= -0.313030 (h = 0.1) and -0.306335 (h = 0.05), approached by the block words to six
  digits; the closures finish at g = -1.01 and not at g = -0.99. The revised paper proves the exact exponent
  lambda_check = -0.3 and lsr(B_h) = ((1-h)(1+0.1*h))^(1/3) for 0 < h < 1.
* **Fig. 1, Obs. 4.12** (`fig1_work_precision`): our certified lower bounds at h = 0.2, 0.1, 0.05, 0.025 for the three
  pairs, identical to the paper's to all printed digits (the word searches and the closures at fixed levels are
  deterministic); the Algorithm (L) points and all times depend on the machine.

## Certificates and replays

| certificate | certifier | statement | MATLAB replay | Python replay |
|---|---|---|---|---|
| polytope (`*_best.json`, `*_transfer.json`, `*_table.json`) | `certify_lower_exact` | `lsr(B_h) >= t`, hence `lambda_check >= (t - 1)/h` (Prop. 4.3) | `replay_lower` | `replay_lower.py` (replay_one of `certify_lower.py`) |
| implicit Euler chain (`<pair>_<h>.json`) | `certify_upper_exact` | `lambda_check <= (1 - 1/r)/h` (Prop. 4.4) | `replay_upper` (uses only `I - h A_i`) | `replay_upper.py` (`certify_upper.py`) |
| chattering / schedule (`*_best_upper.json`, `*_schedule.json`) | `certify_best_upper` | `lambda_check <= c` or `<= log rho(Phi)/T` | `replay_best_upper` (uniformization, no logarithms) | `replay_best_upper.py` (`certify_best_upper.py`) |
| invariant faces (`faces_*.json`) | `invariant_faces` | `lambda_J >= L_J`, `c <= U` (Computation 3.8) | `replay_faces` | `replay_faces.py` |
| invariant cones (`cone_*.json`, `polygon_*.json`) | `certify_cone` (facet normals) | `(I + h_1 A_i) C` in `C` | `replay_cone` (fan triangulation) | `replay_cones.py` (`euler_cone*.py`) |

Each lower certificate stores the exact vertices, for every (mode, vertex) pair the exact weights, the exact level
`t*`, the claimed `t = 1 + h lambda_lower <= t*`, and two negative controls (the binding weights fail at
`t*(1 + 1e-12)`; a dual vector shows that the polytope is not invariant above `tau_up`). The JSON formats are those of
the paper's Python scripts, so each language's replay accepts the other's certificates. To replay one certificate,
or everything:

```matlab
[ok, msgs, neg] = replay_lower('certificates/table2/case1_0.1_best.json')
ok = replay_all                   % every file under certificates/, MATLAB replays
```

```bash
python3 independent_replay/python/run_all.py
```

## Notes on the implementation

* `sym('0.1')` and `str2sym('0.1')` are symbolic floats, not `1/10`. Exact decimals go through `dec2rat`
  (`str2sym(dec2rat('0.1'))` is exactly 1/10) and `sym(x, 'f')` gives the exact binary64 value of a double.
  `str2sym` costs about 0.7 ms per string, so certificates are parsed by `parse_rationals` (and, on the replay side,
  by its own `rat_parse`), which assemble the integers from 15- (9-) digit chunks with vectorized exact arithmetic.
* MATLAB's `linprog` costs 1-5 ms per call. The closure skips LPs whose answer a stored dual vector already implies.
  When both HiGHS solvers fail (tiny entries make the scaled LP ill-conditioned), `lp_antinorm` retries with
  `interior-point-legacy` and another row scaling, and restores a sparse basic solution (crossover) when an interior
  point solution is used.
* MATLAB's LP weights are sometimes less sharp than scipy's (1e-8 relative), which can cost a unit in the sixth
  decimal of a lower bound. `certify_lower_exact` therefore polishes the weights of the near-binding pairs: it solves
  the restricted LP exactly over its bases and keeps the best exactly feasible solution. Only the final exact check
  matters for validity.
* HiGHS occasionally raises an internal error instead of returning a failure flag, or reports "lost feasibility" on an
  LP whose zero vector is feasible (it depends on the last bits of the right-hand side). Every `linprog` call treats an
  error as a failed solve, and in the exact check `lp_antinorm` retries a failed LP with the columns in reverse order
  and with `ConstraintTolerance` 1e-9 before giving up; the closures are unchanged by this.
* The searches use deterministic caps (`pmax`: vertices processed per closure) wherever a result is meant to be
  reproducible; time caps (`probe`, `budget`, `tmax`) are kept for the options that mirror the paper's time budgets.

## Long runs

| run | why it is long | how to run it |
|---|---|---|
| `table7_scaling` | 24 lower rows of 1100 to 2800 vertices | rows in parallel: `table7_scaling('pairs', {...}, 'hs', {...}, 'part', 'lower')`, then `'part', 'summary'` |
| `table3_switching` | 40 rows; d = 10 closures up to 3600 vertices | pairs in parallel, as above |
| `table3_transposed` | 23 rows, up to 6 closures of 300 s each | `'rows'` subsets in parallel |
| `table3_extended` | 34 rows, up to 7 closures of 1800 s each | `'rows'` subsets in parallel |
| `larger_samples` | 234 steps (d = 20 rows up to 1.5 hours), then 60 rows of the extended search | `'pairs'` and `'steps'` subsets in parallel; `'extended'` after `'switching'` |
| `'replay'` of `larger_samples` and `table3_extended` | the MATLAB replay of these 275 certificates (d = 10 and 20 polytopes of up to 5600 vertices) took 12.6 hours on one thread: 6.4 h `larger_switching`, 3.5 h `larger_scaling`, 2.7 h `table3_extended`, up to 38 min for one file | the Python replay takes 17 min; or `replay_lower` on file subsets in parallel (75 min with 14 processes) |
| `table4_algorithm_l` | 30 runs, closures up to 3300 vertices | pairs in parallel |
| `fig1_work_precision` | Algorithm (L) uses its full budget of 1200 s per point | `'methods', {'ours'}` for our points only |


## Differences from the earlier Python computation

The paper's computations were first made with a Python implementation. Where this MATLAB code gives a different value,
the paper prints the MATLAB value, so the drivers reproduce every printed number of the tables. The values that changed:

1. **Algorithm (L) (Table 4, Obs. 4.11).** Its lower bound is the tangent bound of the polytope closed at the shift
   nu, and that polytope depends on the LP solver (MATLAB's HiGHS through `linprog` instead of scipy's). It differed
   in 7 of the 29 runs, changing five cells and two widths of Table 4; the comparison of Obs. 4.11 is unchanged.
2. **Computation 3.8, case 677.** MATLAB's `convhull` treats nearly collinear points differently from Qhull; with the
   inflation factor 5e-5 the driver gets an exactly invariant polygon of 53 vertices with coordinates at least 0.1137.
3. **Computation 3.8, rounding.** The face {1,3} of the transposed case 148 has exponent at least 0.147542 (the
   Python run had rounded 0.147542811 to nearest).
4. **Obs. 3.11(a)** (the slope 1.470) and **Obs. 4.9** (the median deviation 0.15%).

Quantities that depend on the machine and are not part of any claim:

5. **Vertex counts and times.** Polytope sizes and running times depend on the LP solver and the machine. The times
   printed in Sec. 4 are those of the Python implementation; reproduce mode does not bisect.
6. **Fig. 1, Algorithm (L) points.** These runs are limited by a wall-clock budget (1200 s of closures), so their
   results depend on the machine; our 12 points are identical to the paper's.

## Appendix B (not included)

Appendix B reports runs of the Gripenberg-type algorithm of Guglielmi and Maiale (SIAM J. Matrix Anal. Appl. 47
(2026) 77-106). Their code is not part of this repository; it is at
https://github.com/francescomaiale/subradius_computation (the runs of the paper use commit `181d6b5`). The polytope certificates that Appendix B quotes are included, since checking them needs only the code of this repository: `certificates/appB_test` (each lower estimate of Table 5 certified a posteriori, the final vertex sets of the runs, and the closures seeded by them), `appB_fixed` (the fixed-antinorm runs), `appB_switching` (Table 7: vertex sets and seeded closures) and `appB_scaling` (Table 6: seeded closures). `replay_all` and the Python replay check them like every other certificate. The runs themselves were made with the code of Guglielmi and Maiale at the commit above, with a wall-clock limit added.

## Citing and license

See `../CITATION.cff` (to be completed) and `../LICENSE` (the license is still to be chosen).
