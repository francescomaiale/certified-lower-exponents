# Table 3 lower bounds on the transposed family B_h^T (prop:polytope(iii))

The rows of Table 3 (Obs. 4.8) whose lower bound closed only at eps = 0.1 or not at all (---), rerun with our
closure on B_h^T instead of B_h (tables/table3_transposed.m).  Word 1^a 2^b and g_up as in the paper's run; seeds:
Perron vectors of the reversed word for B_h^T; levels t = 1 + h (g_up - eps), eps = 1e-5, 1e-3, 1e-2, 0.1, then
0.15, 0.2, each closure capped at 300 s and 20000 vertices (MATLAB, one thread); the first closed polytope is
certified exactly in both data readings and replayed (MATLAB replay_all and independent_replay/python/run_all.py).
Replays (part 'replay', 2026-10-20 04:46), certificates/table3_transposed (20 files): MATLAB replay_all 20 of 20 VALID
(727 s); independent_replay/python/run_all.py: 20 certificates: lower 20 VALID; 17s.

| pair | h | Table 3, B_h | eps | |V| | time (s) | transposed, B_h^T | eps | |V| | time (s) | gain | better (orientation, eps) | B_h^T closures tried (eps: reason, s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sw10_1 | 0.1 | 0.432612 | 0.1 | 191 | 1887 | 0.382612 | 0.15 | 413 | 1462 + 68 | -0.050000 | 0.432612 (N, 0.1) | 1e-05: time, 418; 0.001: time, 316; 0.01: time, 398; 0.1: time, 311; 0.15: closed, 20 |
| sw10_3 | 0.1 | 0.677933 | 0.1 | 391 | 1516 | 0.677933 | 0.1 | 922 | 1170 + 190 | +0.000000 | 0.677933 (N, 0.1) | 1e-05: time, 491; 0.001: time, 300; 0.01: time, 337; 0.1: closed, 42 |
| sw10_4 | 0.1 | 0.116604 | 0.1 | 529 | 1726 | 0.116604 | 0.1 | 863 | 1072 + 175 | +0.000000 | 0.116604 (N, 0.1) | 1e-05: time, 350; 0.001: time, 300; 0.01: time, 361; 0.1: closed, 60 |
| sw10_6 | 0.1 | 0.535616 | 0.1 | 650 | 1758 | 0.625616 | 0.01 | 2554 | 949 + 1099 | +0.090000 | 0.625616 (T, 0.01) | 1e-05: time, 300; 0.001: time, 300; 0.01: closed, 349 |
| sw10_7 | 0.1 | 1.036617 | 0.1 | 918 | 2054 | 1.136607 | 1e-05 | 2961 | 111 | +0.099990 | 1.136607 (T, 1e-05) | 1e-05: closed, 111 |
| sw10_8 | 0.1 | -0.023001 | 0.1 | 677 | 1682 | -0.023001 | 0.1 | 1086 | 1012 + 239 | +0.000000 | -0.023001 (N, 0.1) | 1e-05: time, 321; 0.001: time, 300; 0.01: time, 325; 0.1: closed, 66 |
| sw10_9 | 0.1 | --- | --- | --- | 2046 | 0.832143 | 0.1 | 2450 | 949 | new | 0.832143 (T, 0.1) | 1e-05: time, 300; 0.001: time, 300; 0.01: time, 300; 0.1: closed, 49 |
| sw10_10 | 0.1 | 0.269425 | 0.1 | 1024 | 1713 | 0.269425 | 0.1 | 1420 | 1241 + 321 | +0.000000 | 0.269425 (N, 0.1) | 1e-05: time, 364; 0.001: time, 361; 0.01: time, 343; 0.1: closed, 173 |
| sw5_1 | 0.05 | 0.366491 | 0.1 | 488 | 1540 | 0.266491 | 0.2 | 1063 | 1684 + 52 | -0.100000 | 0.366491 (N, 0.1) | 1e-05: time, 302; 0.001: time, 301; 0.01: time, 300; 0.1: time, 300; 0.15: time, 307; 0.2: closed, 175 |
| sw5_2 | 0.05 | 0.440719 | 0.01 | 476 | 1205 | 0.440719 | 0.01 | 2228 | 1020 + 271 | +0.000000 | 0.440719 (N, 0.01) | 1e-05: time, 370; 0.001: time, 421; 0.01: closed, 228 |
| sw5_3 | 0.05 | -0.509527 | 0.1 | 448 | 1514 | -0.410527 | 0.001 | 1764 | 428 + 223 | +0.099000 | -0.410527 (T, 0.001) | 1e-05: time, 311; 0.001: closed, 117 |
| sw5_7 | 0.05 | 1.396984 | 0.01 | 1326 | 1424 | 1.405984 | 0.001 | 1545 | 456 + 141 | +0.009000 | 1.405984 (T, 0.001) | 1e-05: time, 372; 0.001: closed, 84 |
| sw5_9 | 0.05 | 0.430659 | 0.01 | 1541 | 1384 | 0.430659 | 0.01 | 5080 | 1021 + 480 | +0.000000 | 0.430659 (N, 0.01) | 1e-05: time, 375; 0.001: time, 423; 0.01: closed, 223 |
| sw10_1 | 0.05 | 0.446053 | 0.1 | 719 | 1525 | 0.346053 | 0.2 | 793 | 2024 + 140 | -0.100000 | 0.446053 (N, 0.1) | 1e-05: time, 410; 0.001: time, 468; 0.01: time, 490; 0.1: time, 300; 0.15: time, 300; 0.2: closed, 55 |
| sw10_2 | 0.05 | 0.553398 | 0.1 | 661 | 1577 | --- | --- | --- | 1978 | --- | 0.553398 (N, 0.1) | 1e-05: time, 408; 0.001: time, 300; 0.01: time, 300; 0.1: time, 369; 0.15: time, 300; 0.2: time, 300 |
| sw10_3 | 0.05 | 0.683448 | 0.1 | 2264 | 1801 | 0.583448 | 0.2 | 752 | 1852 + 130 | -0.100000 | 0.683448 (N, 0.1) | 1e-05: time, 445; 0.001: time, 484; 0.01: time, 300; 0.1: time, 300; 0.15: time, 300; 0.2: closed, 23 |
| sw10_4 | 0.05 | --- | --- | --- | 1815 | 0.081251 | 0.15 | 1848 | 1820 + 564 | new | 0.081251 (T, 0.15) | 1e-05: time, 404; 0.001: time, 303; 0.01: time, 439; 0.1: time, 388; 0.15: closed, 286 |
| sw10_5 | 0.05 | -0.312875 | 0.1 | 1022 | 1567 | -0.412875 | 0.2 | 1152 | 1722 + 259 | -0.100000 | -0.312875 (N, 0.1) | 1e-05: time, 332; 0.001: time, 322; 0.01: time, 369; 0.1: time, 340; 0.15: time, 308; 0.2: closed, 51 |
| sw10_6 | 0.05 | --- | --- | --- | 2334 | 0.536997 | 0.1 | 814 | 1024 + 159 | new | 0.536997 (T, 0.1) | 1e-05: time, 316; 0.001: time, 311; 0.01: time, 300; 0.1: closed, 96 |
| sw10_7 | 0.05 | --- | --- | --- | 1834 | 1.061027 | 0.1 | 479 | 1107 + 77 | new | 1.061027 (T, 0.1) | 1e-05: time, 374; 0.001: time, 416; 0.01: time, 300; 0.1: closed, 16 |
| sw10_8 | 0.05 | --- | --- | --- | 1985 | -0.039419 | 0.15 | 3397 | 1619 + 1853 | new | -0.039419 (T, 0.15) | 1e-05: time, 303; 0.001: time, 304; 0.01: time, 320; 0.1: time, 334; 0.15: closed, 359 |
| sw10_9 | 0.05 | --- | --- | --- | 1855 | 0.738290 | 0.2 | 1748 | 2040 + 557 | new | 0.738290 (T, 0.2) | 1e-05: time, 362; 0.001: time, 444; 0.01: time, 495; 0.1: time, 300; 0.15: time, 300; 0.2: closed, 139 |
| sw10_10 | 0.05 | --- | --- | --- | 1864 | 0.178542 | 0.2 | 1593 | 1797 + 455 | new | 0.178542 (T, 0.2) | 1e-05: time, 317; 0.001: time, 300; 0.01: time, 300; 0.1: time, 319; 0.15: time, 300; 0.2: closed, 261 |

Table 3, B_h: certificates/switching/<pair>_<h>_best.json, eps its level, |V| its vertices; time: the paper's
Python search, all levels tried (data/paper_search_results.json; 600 s per closure).  Transposed, B_h^T: our
certificate, |V| its vertices; time: all closures of the ladder, failed ones included (MATLAB, 300 s cap each; a
cap is checked between batches, so a failed closure can run past 300 s), + exact check and replay.  gain:
transposed minus Table 3 (new where Table 3 has ---; negative where B_h^T closed only at a larger eps).  better:
the larger certified bound, with its orientation (N = B_h, T = B_h^T) and eps.  Reasons: 'time' = 300 s cap,
'size' = 20000 vertices.

## Bracket widths (upper schedule bound minus the best certified lower bound, both h)

| pair | upper | best lower, Table 3 | width, Table 3 | best lower, with B_h^T | width, with B_h^T |
|---|---|---|---|---|---|
| sw5_1 | 0.467683 | 0.465326 | 2.4e-03 | 0.465326 | 2.4e-03 |
| sw5_2 | 0.457887 | 0.442418 | 1.6e-02 | 0.442418 | 1.6e-02 |
| sw5_3 | -0.365833 | -0.457223 | 9.2e-02 | -0.410527 | 4.5e-02 |
| sw5_7 | 1.428992 | 1.396984 | 3.3e-02 | 1.405984 | 2.4e-02 |
| sw5_9 | 0.474726 | 0.430659 | 4.5e-02 | 0.430659 | 4.5e-02 |
| sw10_1 | 0.559043 | 0.446053 | 1.2e-01 | 0.446053 | 1.2e-01 |
| sw10_2 | 0.674041 | 0.629866 | 4.5e-02 | 0.629866 | 4.5e-02 |
| sw10_3 | 0.788865 | 0.683448 | 1.1e-01 | 0.683448 | 1.1e-01 |
| sw10_4 | 0.245147 | 0.116604 | 1.3e-01 | 0.116604 | 1.3e-01 |
| sw10_5 | -0.197966 | -0.228879 | 3.1e-02 | -0.228879 | 3.1e-02 |
| sw10_6 | 0.638392 | 0.535616 | 1.1e-01 | 0.625616 | 1.3e-02 |
| sw10_7 | 1.185699 | 1.036617 | 1.5e-01 | 1.136607 | 5.0e-02 |
| sw10_8 | 0.142660 | -0.023001 | 1.7e-01 | -0.023001 | 1.7e-01 |
| sw10_9 | 0.944065 | --- | --- | 0.832143 | 1.2e-01 |
| sw10_10 | 0.387498 | 0.269425 | 1.2e-01 | 0.269425 | 1.2e-01 |

Per pair, over h = 0.1 and 0.05; rounded up to two digits, as in Table 3 (--- where no lower bound is certified).

## Level that closed (eps = 1e-5 / 1e-3 / 1e-2 / 0.1 / 0.15 / 0.2 / none), rows above

- h = 0.1, d = 10 (8 rows): B_h (Table 3) [0 0 0 7 0 0 1], B_h^T [1 0 1 5 1 0 0]
- h = 0.05, d = 5 (5 rows): B_h (Table 3) [0 0 3 2 0 0 0], B_h^T [0 2 2 0 0 1 0]
- h = 0.05, d = 10 (10 rows): B_h (Table 3) [0 0 0 4 0 0 6], B_h^T [0 0 0 2 2 5 1]
