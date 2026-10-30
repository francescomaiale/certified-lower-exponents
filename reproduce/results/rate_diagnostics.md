# Table 5 (tab:rates): predicted and observed rate constants, MATLAB reproduction (floating point)

Here g and G denote floating-point word estimates of the discrete minima. The gap and midpoint limits are diagnostics, not certified asymptotics.

| pair | alpha*_1 | predicted Q | (G-g)/h, h = 0.1 | 0.05 | 0.025 | ((g+G)/2 - eta)/h^2 at 0.025 | paper | match |
|---|---|---|---|---|---|---|---|---|
| planar | 0.3905 | 1.4780 | 1.4861 | 1.4800 | 1.4785 | -0.291 | 0.3905, 1.4780, 1.4861, 1.4800, 1.4785 | yes |
| case1 | 0.7890 | 0.5480 | 0.5619 | 0.5507 | 0.5485 | -0.350 | 0.7890, 0.5480, 0.5619, 0.5507, 0.5485 | yes |
| d5_dense1 | 0.8548 | 0.1785 | 0.1785 | 0.1785 | 0.1785 | -0.034 | 0.8548, 0.1785, 0.1785, 0.1785, 0.1785 | yes |
| d5_dense2 | 0.1906 | 0.5888 | 0.5901 | 0.5898 | 0.5889 | -0.076 | 0.1906, 0.5888, 0.5901, 0.5898, 0.5889 | yes |
| d5_dense3 | 0.2520 | 0.1975 | 0.1973 | 0.1974 | 0.1975 | -0.008 | 0.2520, 0.1975, 0.1973, 0.1974, 0.1975 | yes |
| d5_sparse1 | 0.2068 | 1.3307 | 1.3273 | 1.3282 | 1.3304 | -0.204 | 0.2068, 1.3307, 1.3273, 1.3282, 1.3304 | yes |
| d5_sparse2 | 0.8950 | 0.2024 | 0.2026 | 0.2025 | 0.2025 | -0.039 | 0.8950, 0.2024, 0.2026, 0.2025, 0.2025 | yes |
| d5_sparse3 | 0.4736 | 2.3887 | 2.4183 | 2.4010 | 2.3903 | -0.295 | 0.4736, 2.3887, 2.4183, 2.4010, 2.3903 | yes |
| d10_dense1 | 0.2085 | 0.0447 | 0.0445 | 0.0446 | 0.0447 | +0.018 | 0.2085, 0.0447, 0.0445, 0.0446, 0.0447 | yes |
| d10_dense2 | 0.4063 | 0.3178 | 0.3191 | 0.3183 | 0.3179 | +0.100 | 0.4063, 0.3178, 0.3191, 0.3183, 0.3179 | yes |
| d10_dense3 | 0.4197 | 0.0770 | 0.0768 | 0.0770 | 0.0770 | -0.004 | 0.4197, 0.0770, 0.0768, 0.0770, 0.0770 | yes |
| d10_sparse1 | 0.8253 | 0.7134 | 0.6981 | 0.7091 | 0.7123 | -0.369 | 0.8253, 0.7134, 0.6981, 0.7091, 0.7123 | yes |
| d10_sparse2 | 0.3296 | 0.4656 | 0.4732 | 0.4686 | 0.4671 | +0.547 | 0.3296, 0.4656, 0.4732, 0.4686, 0.4671 | yes |
| d10_sparse3 | 0.5997 | 1.5941 | 1.6197 | 1.6006 | 1.5959 | -0.829 | 0.5997, 1.5941, 1.6197, 1.6006, 1.5959 | yes |

Relative deviation |(G-g)/h - Q|/Q: max [2.52 0.64 0.31], median [0.39 0.15 0.03] at h = 0.1, 0.05, 0.025 (paper: max 2.5%, 0.64%, 0.31%;
medians 0.39%, 0.15%, 0.03%).  max |(g+G)/2 - eta| / h^2 at h = 0.025: 0.83 (paper: <= 0.83).
Switching pair 4 (d = 5): 2 L(0.05) - L(0.1) (L: certified polytope lower bound) = -0.4270, above the certified upper bound -0.430082 (paper: -0.4270, -0.4301).
Rows reproduced to the printed digits: 14 of 14.

## Separate slopes from the original unrounded search diagnostics

These floating-point estimates are retained in `rate_search_results.json`; the original columns above agree with that run to all printed digits. Fresh MATLAB runs also report both slopes.

| pair | (c-g_est)/h at h=0.025 | (G_est-c)/h at h=0.025 |
|---|---|---|
| planar | 0.7466 | 0.7320 |
| case1 | 0.2830 | 0.2655 |
| d5_dense1 | 0.0901 | 0.0884 |
| d5_dense2 | 0.2964 | 0.2926 |
| d5_dense3 | 0.0989 | 0.0985 |
| d5_sparse1 | 0.6703 | 0.6601 |
| d5_sparse2 | 0.1022 | 0.1003 |
| d5_sparse3 | 1.2025 | 1.1878 |
| d10_dense1 | 0.0219 | 0.0228 |
| d10_dense2 | 0.1565 | 0.1615 |
| d10_dense3 | 0.0386 | 0.0384 |
| d10_sparse1 | 0.3654 | 0.3469 |
| d10_sparse2 | 0.2199 | 0.2472 |
| d10_sparse3 | 0.8187 | 0.7772 |
