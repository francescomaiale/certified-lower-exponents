# Table 4 (tab:glp): best certified brackets and Algorithm (L), MATLAB reproduction

Each cell: MATLAB value, then the paper's value in brackets when they differ.

| pair | Alg. (L) lower | ours lower | Alg. (L) upper (float) | implicit Euler | c or schedule | width | all ours match |
|---|---|---|---|---|---|---|---|
| planar | -0.210758 | -0.232706 | -0.194213 | -0.158706 | -0.194976 | 1.6e-02 | yes |
| case1 | -1.028582 | -1.027040 | -1.011626 | -0.999505 | -1.012395 | 1.5e-02 | yes |
| case1440 | 0.334850 | 0.345146 | 0.346401 | 0.346182 | 0.345675^s | 5.3e-04 | yes |
| case148 | --- | -1.247834 | --- | -1.238961 | -1.241336^s | 6.5e-03 | yes |
| d5_dense1 | -0.276163 | -0.222751 | -0.213075 | -0.209005 | -0.213383 | 9.4e-03 | yes |
| d5_dense2 | -0.081527 | -0.059850 | -0.028459 | -0.014798 | -0.029331 | 3.1e-02 | yes |
| d5_dense3 | 0.210187 | 0.222107 | 0.233137 | 0.237887 | 0.232972 | 1.1e-02 | yes |
| d5_sparse1 | -0.204184 | -0.107049 | -0.032087 | -0.001609 | -0.034259 | 7.3e-02 | yes |
| d5_sparse2 | -0.384238 | -0.641285 | -0.102190 | -0.097694 | -0.102660 | 2.9e-01 | yes |
| d5_sparse3 | -0.049779 | -0.071028 | 0.055587 | 0.112334 | 0.053208 | 1.1e-01 | yes |
| d10_dense1 | 0.870186 | 0.459247 | 0.990484 | 0.991569 | 0.990410 | 1.3e-01 | yes |
| d10_dense2 | 0.779913 | 0.525978 | 0.857197 | 0.865203 | 0.856999 | 7.8e-02 | yes |
| d10_dense3 | 0.009496 | -0.229448 | 0.066898 | 0.068764 | 0.066852 | 5.8e-02 | yes |
| d10_sparse1 | -0.064519 | 0.269192^g | 0.318964 | 0.334496 | 0.317685 | 4.9e-02 | yes |
| d10_sparse2 | 0.364290 | 1.356787^g | 1.570804 | 1.583464 | 1.570368 | 2.2e-01 | yes |
| d10_sparse3 | -0.451797 | -0.349697 | -0.179980 | -0.143149 | -0.181089 | 1.7e-01 | yes |

Ours: the best of the steps 0.1 and 0.05 (certificates/table2, certificates/scaling); case 148 at h <= 0.1 by
prop:polytope(ii); ^g: a closure seeded by the final vertex set of a run of [GM] (certificates/appB_scaling, App. B).
c or schedule: certificates/table4/<pair>_best_upper.json (^s: two-segment schedule).  Width =
that column minus the larger lower bound, rounded up to two digits.  Alg. (L): results/glp/ (table4_algorithm_l),
best of tau = 0.1 and 0.05; its upper bounds are floating point, rounded up.
Our certified columns reproduced for 16 of 16 pairs; Alg. (L) columns for 16 of 16.
