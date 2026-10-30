# Computation 3.9 (obs:cone): the reducible test pairs, MATLAB reproduction

## Invariant faces (binary64 data), both orientations

Exponent bounds are exact lower bounds rounded down to 9 decimals; the paper rounds them to 6 decimals.

| case | orientation | certified c <= | invariant proper faces (J: mode, exact eta >=) | verdict | replay |
|---|---|---|---|---|---|
| case77 | A | -0.148186 | {1,3}: A_2, -0.663165228 (not certified slow) | some face not certified slow | VALID |
| case77 | T | -0.148186 | {2}: A_2^T, 0.024892422 (slow) | all invariant faces slow: lambda_min = lambda_check | VALID |
| case148 | A | -1.236979 | {2}: A_2, -1.368009664 (not certified slow); {1,2}: A_1, -1.220814739 (slow) | some face not certified slow | VALID |
| case148 | T | -1.236979 | {3}: A_1^T, -3.669288483 (not certified slow); {1,3}: A_2^T, 0.147542811 (slow) | some face not certified slow | VALID |
| case677 | A | -0.394492 | {2}: A_2, -0.237096535 (slow) | all invariant faces slow: lambda_min = lambda_check | VALID |
| case677 | T | -0.394492 | {1,3}: A_2^T, -1.650614193 (not certified slow) | some face not certified slow | VALID |
| case1195 | A | -0.532919 | {2,3}: A_1, -0.289317810 (slow) | all invariant faces slow: lambda_min = lambda_check | VALID |
| case1195 | T | -0.532919 | {1}: A_1^T, -3.201148067 (not certified slow) | some face not certified slow | VALID |
| case1440 | A | 0.348177 | none | all invariant faces slow: lambda_min = lambda_check | VALID |
| case1440 | T | 0.348177 | none | all invariant faces slow: lambda_min = lambda_check | VALID |

Paper: 1195 (as given) {2,3} under A_1, >= -0.289318, c <= -0.532919; 677 (as given) {2} under A_2, -0.237097,
c <= -0.394492; 77 transposed {2} under A_2^T, 0.024892, c <= -0.148186; 1440 none; 77 as given -0.663165 and
677 transposed -1.650614 (below c); 148: {2} under A_2, about -1.368010, transposed {3} under A_1^T, -3.669288,
and {1,3} under A_2^T, >= 0.147543, c <= -1.236979.

## Case 148, transposed family: a cone that touches only slow faces (cor:eulerface)
- C = {x >= 0 : (19/20)(x1 + x2) >= x3/20}, h1 = 1/5: 4 extreme rays, 0 images outside, invariant: yes; replay VALID 
- support classes of C: 1 | 2 | 1,2 | 1,3 | 2,3 | 1,2,3 (all nonempty subsets except {3})
- invariant faces of the transposed family: 3 | 1,3; touched by C: 1,3
  face {1,3} under A_2^T: exponent >= 0.147542811, slow: yes (c <= -1.236979)

## Cases 77 (as given) and 677 (transposed): invariant polygons of prop:cone(ii), h1 = 0.2
- case77 (A): 21 vertices, min coordinate 0.0796, inflation 0.001, 200 rounds; exact invariance yes; replay VALID 
- case677 (T), inflation 0.0001: 51 vertices, 1 images outside (exact); next factor
- case677 (T): 53 vertices, min coordinate 0.1137, inflation 5e-05, 200 rounds; exact invariance yes; replay VALID 

Paper: 21 and 47 vertices, coordinates at least 0.0796 and 0.1134 (the counts depend on the floating-point hull;
only the exact checks matter).
