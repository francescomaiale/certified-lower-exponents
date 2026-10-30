# Certified brackets for lower exponents of positive switched systems

MATLAB code for the paper *Euler discretizations, invariant faces and certified brackets for lower exponents of
positive switched systems*.

Given Metzler matrices A_1, ..., A_m, the code computes certified bounds

    lower <= lambda <= upper

for the lower Lyapunov exponent lambda of x' = A_sigma(t) x, the best exponential rate that switching can achieve.

- **Lower bound.** A finite polytope invariant under the explicit Euler matrices I + h A_i (Prop. 4.3). It is found
  in floating point, then checked in exact rational arithmetic.
- **Upper bounds.**
  - A word of the implicit Euler matrices (I - h A_i)^(-1), with a chain of vectors (Prop. 4.4).
  - For two modes, also an upper bound from a constant mixture and the best two-segment periodic schedule found,
    certified with exponential enclosures.

Every endpoint of the pipeline's printed continuous-time brackets is checked in exact arithmetic for both readings
of the data: the decimals as written, and the binary64 numbers they round to. The floating-point searches only
propose candidates. The reproduction tables also include explicitly labelled floating-point diagnostics and
comparison columns; these are not certificates of discrete minima or of mixture optimality.

## Requirements

MATLAB with the Optimization Toolbox (`linprog`) and the Symbolic Math Toolbox (exact rationals). Developed and
tested with R2026b; R2023a or later is needed.

## Examples

| script | what it does |
|---|---|
| `example_custom.m` | Certified bounds for your own matrices: edit the data block (matrices, step `h`, time budget) and run it. About a minute for the 3x3 example in the script. |
| `example_paper.m` | Rebuilds Table 2 and the initial Table 3 closures, before the extended refinements marked `+`; compares against their recorded bounds. `quick = true` skips the long rows (roughly 7–20 minutes, depending on the machine). |

Each script adds `src/` to the path itself.

## Contents

| path | content |
|---|---|
| `src/pipeline` | `certified_bracket`: the whole pipeline for one step `h` (searches, then exact certificates). |
| `src/lower` | `polytope_closure` (closure of seed vectors at a level t), `polytope_search` (bisection over t), `perron_seeds`, `lp_antinorm`, and `certify_lower_exact` (the exact check of Prop. 4.3). |
| `src/upper` | `word_search` / `brute_words` (word search), `certify_upper_exact` (Prop. 4.4), `chattering_value`, `best_schedule`, and `certify_best_upper` with its enclosures `exp_upper` / `log_upper`. |
| `src/cone` | Invariant faces and invariant cones (Computation 3.9). |
| `src/glp` | A simple implementation of Algorithm (L) of Guglielmi, Laglia and Protasov, used for the comparison of Obs. 4.9. |
| `src/core` | Exact data readings (`exact_data`), rational helpers, decimal formatting. |
| `paper_data.m` | The pairs of the paper, written as their decimals, and the recorded levels, words and schedules that `example_paper.m` starts from. |
| `reproduce/` | Table and figure drivers, 568 replayable certificates, MATLAB replays and independent Python replays. Algorithm (L) tangent bounds are checked exactly when computed but have no stored replay certificates. See `reproduce/README.md`. |

The reproduction README includes a [numbering guide](reproduce/README.md#paper-numbering-and-stable-filenames)
for historical driver names and archived reports. Current references distinguish the main paper from Tables S1-S6
in its numerical supplement.

## Verification

From this repository's root, Python alone can check all supplied certificates; MATLAB is not needed for replay:

```bash
python3 -m pip install -r requirements.txt
python3 -m unittest discover -s reproduce/independent_replay/python -p 'test_*.py' -v
python3 -I -B reproduce/independent_replay/python/run_all.py
```

The dependency file records the tested NumPy version with Python 3.12. The optional numerical theorem checks and
data generators also use SciPy: install `reproduce/requirements-numerics.txt` for those tools.

To test the MATLAB implementation, including new searches, certificate generation and rejection of altered proofs:

```bash
matlab -singleCompThread -batch "addpath('reproduce'); startup_reproduce; run_tests"
```

The GitHub Actions workflow in `.github/workflows/certificates.yml` runs the Python regression tests and all
certificate replays on pushes and pull requests. MATLAB validation remains a separate local test.

## Related code

Appendix B reports runs of the Gripenberg-type algorithm of Guglielmi and Maiale (SIAM J. Matrix Anal.
Appl. 47 (2026) 77-106). The portable wall-clock wrapper, capped routines and their two upstream helpers are now
bundled in [`reproduce/gm`](reproduce/gm/README.md), based on upstream commit
`181d6b59be56a02dab1a040254b95c905d7f49ae`. That optional comparison component retains the upstream GPL-3.0
licence. The main certification and replay package remains MIT-licensed.

Original Appendix B run records and all 142 final vertex sets are in `reproduce/results/gm_original`, with a
SHA-256 manifest. This includes the seven seeds used automatically by `table3_extended` at h = 0.1.

## Notes

- Times depend on the machine. The lower-bound search of `certified_bracket` has a time budget, so on a slower
  machine it may certify a slightly weaker bound. The recorded levels in `paper_data.m` make `example_paper.m`
  independent of the original bisection budget; the LP solver can still affect vertices and final rounding.
- The explicit Euler matrices must be nonnegative; h < 1 / max |a_jj| is a sufficient restriction. The implicit
  bound additionally requires h * eta(A_i) < 1 for every mode, where eta is the spectral abscissa. The stronger
  condition h * max_i ||A_i||_1 < 1 suffices for both. The pipeline checks both restrictions exactly before searching.
- The repository is about 320 MB, almost all of it the certificates in `reproduce/certificates` (568 files).

## Citing and license

The main software is distributed under the [MIT License](LICENSE). The optional GM comparison code in
`reproduce/gm/` is covered by its own [GPL-3.0 licence](reproduce/gm/LICENSE); see its provenance notice.

See `CITATION.cff` for the citation metadata, which remains to be completed. Before the first public release,
confirm the complete author list and replace the citation placeholders with the repository URL and actual
publication/release metadata. These remaining release details are deferred.
