# Appendix B: GM comparison runs

This optional comparison component contains `gm_fixed_job.m`, `adaptive_capped.m`, `lsr_comp_capped.m`,
`aNorm.m` and `real_antinorm.m`. It is distributed under the GNU General Public License version 3 in
[LICENSE](LICENSE). The main certification/replay code outside this directory retains its MIT licence.

The capped routines adapt `adaptive_subradius_comp.m` and `lsr_comp.m` from
[francescomaiale/subradius_computation](https://github.com/francescomaiale/subradius_computation/tree/181d6b59be56a02dab1a040254b95c905d7f49ae),
commit `181d6b59be56a02dab1a040254b95c905d7f49ae`. The two helpers and licence are copied unchanged from that
commit. The capped routines are the original paper's wrappers: they add a wall-clock check and termination
status, retaining the last completed level's lower estimate on timeout. Their numerical code has not been
changed for packaging. On 22 October 2026 the job driver was made repository-relative, given explicit diagnostic
options, and changed to return its result without changing the caller's working directory or MATLAB path.

No external checkout or CVX installation is required. MATLAB and the Optimization and Symbolic Math Toolboxes
are required, as for the rest of the reproduction package. The upstream LP helper is preserved as used in the
paper and has only been exercised here with MATLAB R2026b. Its numerical estimates are not rational certificates.

From the repository root in MATLAB:

```matlab
addpath('reproduce/gm');
R = gm_fixed_job('test', 'case1', '0.1', 'E', 'adaptive');
R = gm_fixed_job('fixed', 'case1', '0.1', 'E', 'fixed');
```

Groups are `test` (Table 5), `scaling` (Table 6), `switching` (Table 7) and `fixed` (Corollary 4.1).
Use side `E` for explicit or `I` for implicit Euler, and variant `adaptive` or `fixed`.
Defaults preserve the paper's transposed families, coordinate initialization, normalization, delta = h^2/10,
product budget 20000 and caps of 7200, 1800, 1800 and 3600 seconds respectively. Start MATLAB with
`-singleCompThread` for the paper's setting. Caps are checked between batches of products; a long LP call can
exceed them. The product budget is checked between levels and may also be exceeded by the final level.

Outputs go to `reproduce/results/gm_runs/<variant>/<group>/`: a JSON record, a MAT file containing the final
vertex matrix `V`, and any independently checked certificates in `certificates/`. Existing rows are skipped;
`'force', true` regenerates them. `'rdir', folder` changes the base output directory. Explicit runs attempt the
same a posteriori certifications as the original wrapper. Only a successful exact certificate establishes a
rigorous bound; an adaptive estimate, LP tolerance or saved vertex set alone does not.

A short installation check, using a temporary directory, is:

```matlab
R = gm_fixed_job('test', 'case1', '0.1', 'E', 'adaptive', ...
    'rdir', tempname, 'cap', 0, 'products', 2, 'certify', false);
```

`cap`, `products` and `certify` are diagnostic overrides. The result records `paper_protocol=false` whenever
these differ from the paper's defaults. Disabling certification cannot produce a certified result.

`../results/gm_original/` preserves all 142 original run records and final vertex sets byte for byte. Its
`MANIFEST.json` gives SHA-256 hashes and matrix dimensions. It includes the seven `switching/sw10_*_0.1_E_V.mat`
files used by the extended search (indices 1, 3, 4, 6, 8, 9, 10). `table3_extended` loads those seeds by default;
its `gm_dir` option can override their directory or explicitly omit them with an empty string. Fresh timed runs
need not return the same vertices, runtime or final estimate as the original runs.

Run `run_gm_tests` after `addpath('reproduce'); startup_reproduce` for wrapper and saved-seed checks. The Python
regression suite also checks the original asset hashes and the restored Table 9 column.
