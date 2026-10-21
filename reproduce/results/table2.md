# Table 2 (tab:gripenberg): certified Euler brackets, MATLAB reproduction

| pair | c (float, rounded up) | h | lower (MATLAB) | lower (paper) | upper (MATLAB) | upper (paper) | width | |V| | word length | match |
|---|---|---|---|---|---|---|---|---|---|---|
| planar | -0.194976 | 0.1 | -0.272230 | -0.272230 | -0.123619 | -0.123619 | 1.5e-01 | 230 | 31 | yes |
| planar | -0.194976 | 0.05 | -0.232706 | -0.232706 | -0.158706 | -0.158706 | 7.4e-02 | 249 | 18 | yes |
| planar | -0.194976 | 0.025 | -0.213641 | -0.213641 | -0.176676 | -0.176676 | 3.7e-02 | 190 | 77 | yes |
| planar | -0.194976 | 0.0125 | -0.204261 | -0.204261 | -0.185783 | -0.185783 | 1.9e-02 | 1033 | 141 | yes |
| case1 | -1.012395 | 0.1 | -1.044234 | -1.044234 | -0.988047 | -0.988047 | 5.7e-02 | 225 | 17 | yes |
| case1 | -1.012395 | 0.05 | -1.027040 | -1.027040 | -0.999505 | -0.999505 | 2.8e-02 | 772 | 31 | yes |
| case1 | -1.012395 | 0.025 | -1.019518 | -1.019518 | -1.005757 | -1.005757 | 1.4e-02 | 3131 | 10 | yes |
| case1440 | 0.348177 | 0.1 | 0.344598 | 0.344598 | 0.346673 | 0.346673 | 2.1e-03 | 169 | 24 | yes |
| case1440 | 0.348177 | 0.05 | 0.345146 | 0.345146 | 0.346182 | 0.346182 | 1.1e-03 | 315 | 50 | yes |
| case1440 | 0.348177 | 0.025 | 0.345146^dagger | 0.345146 | 0.345931 | 0.345931 | 7.9e-04 | 315 | 102 | yes |
| case148 | -1.236979 | 0.1 | -1.247834 | -1.247834 | -1.236958 | -1.236958 | 1.1e-02 | 1836 | 27 | yes |
| case148 | -1.236979 | 0.05 | -1.247834^dagger | -1.247834 | -1.238961 | -1.238961 | 8.9e-03 | 1836 | 51 | yes |
| case148 | -1.236979 | 0.025 | -1.247834^dagger | -1.247834 | -1.240094 | -1.240094 | 7.8e-03 | 1836 | 97 | yes |

Width = upper - lower (MATLAB values), rounded up to two digits.  ^dagger: the bound certified at the larger step (case 148: h = 0.1;
case 1440: h = 0.05), which holds here by prop:polytope(ii); the transfer certificate re-checks that polytope
directly at this step (certificates/table2/<row>_transfer.json).

## Remark 4.5 and Obs. 4.6
- case 1, h = 0.1: best explicit word (length 8) gives g(0.1) <= -1.044232 (float, rounded up); certified g(0.1) >= -1.044234.
- case 148, h = 0.1: 1836 vertices, 0 of them on the axis of e_2; min (v_1 + v_3)/1'v = 0.0741.
- case 148: lambda_min <= (A_2)_22 = -1.368010 < -1.368009 (exact binary64 value -1.3680096637401422).
- all rows reproduced to the printed digits: yes (13 of 13)
