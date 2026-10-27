%% This script can be used to certify bounds for the lower exponent of your own positive switched system x' = A_sigma(t) x.
% Edit the data block below and run the script. Every printed bound is certified in exact rational arithmetic,
% for both readings of the data (the decimals as written, and the binary64 numbers they round to).

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));

%% Data
% Metzler matrices (off-diagonal entries >= 0), all of the same size. Strings give exact decimals; plain numbers
% are read as their shortest decimal (0.3 -> 3/10).
A = { ["-1"   "0.5"  "0"
       "0"    "-2"   "1"
       "0.3"  "0"    "-0.5"], ...
      ["-0.5" "0"    "0.2"
       "1"    "-1"   "0"
       "0"    "0.4"  "-2"] };
h = "0.05";          % needs h max|a_jj| < 1 and h eta(A_i) < 1; checked exactly before searching
budget = 60;         % seconds for the search of the lower bound

%% 1. Euler bracket at step h
% lower: a polytope invariant under the explicit Euler matrices I + h A_i (Prop. 4.3)
% upper: a word of the implicit Euler matrices (I - h A_i)^(-1) with a chain of vectors (Prop. 4.4)
res = certified_bracket(A, h, struct('budget', budget));

%% 2. Upper bounds that do not depend on h (two modes)
% Certify an upper bound for c = min_alpha eta(alpha A_1 + (1-alpha) A_2), and the rate of the best
% two-segment schedule found (A_1 for t1, then A_2 for t2). The searches do not certify optimality.
if numel(A) == 2
    D = exact_data(A, h); Af = D.Af;  % supports both numeric matrices and decimal strings
    [rate, t1, t2] = best_schedule(Af);
    R = certify_best_upper(A, struct('schedule', struct('modes', [1 2], 'durations', [t1 t2])));
    fprintf('\nchattering value c = %+.6f (floating point); best schedule t1 = %.4f, t2 = %.4f, rate %+.6f\n', ...
            R.c_float, t1, t2, rate);
    fprintf('certified: lambda <= %s (%s)\n', R.printed, R.best_kind);
    upper = min(res.upper_value, str2double(R.printed));
else
    upper = res.upper_value;
end

%% Summary
fprintf('\n%s <= lambda <= %.6f   (width %.2e)\n', res.lower, upper, upper - res.lower_value);
