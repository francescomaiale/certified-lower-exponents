function [rate, t1, t2] = best_schedule(A)
% BEST_SCHEDULE  Two-segment periodic schedule (A_1 for t1, then A_2 for t2) minimizing log rho(Phi) / (t1 + t2),
% Phi = e^{t2 A_2} e^{t1 A_1}, in floating point (Nelder-Mead in log t1, log t2 from 16 starting points, durations
% in (1e-5, 30)).  A: cell {A_1, A_2} of doubles.  Only the certified
% bound of certify_best_upper matters.
f = @(z) sched_rate(A, exp(z(1)), exp(z(2)));
o = optimset('TolX', 1e-10, 'TolFun', 1e-13, 'MaxIter', 4000, 'MaxFunEvals', 8000, 'Display', 'off');
rate = Inf; t1 = NaN; t2 = NaN;
for a = [0.1 0.3 1 2]
    for b = [0.1 0.3 1 2]
        [z, fv] = fminsearch(f, log([a b]), o);
        if fv < rate, rate = fv; t1 = exp(z(1)); t2 = exp(z(2)); end
    end
end
end


function r = sched_rate(A, t1, t2)
if ~(t1 > 1e-5 && t1 < 30 && t2 > 1e-5 && t2 < 30), r = 1e9; return, end
P = expm(t2 * A{2}) * expm(t1 * A{1});
rho = max(abs(eig(P)));
if rho > 0 && isfinite(rho), r = log(rho) / (t1 + t2); else, r = 1e9; end
end
