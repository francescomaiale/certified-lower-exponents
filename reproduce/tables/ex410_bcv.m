function R = ex410_bcv(varargin)
% EX410_BCV  Ex. 4.10 (ex:bcv): the three-mode benchmark of Blanchini, Colaneri and Valcher (CDC 2013), eps = 1/10
% (port of bcv_benchmark.py).
%   - c = min over the simplex of eta(S_alpha) (floating point; the paper: 1/120);
%   - Computation: the cyclic schedule A_3, A_2, A_1 with dwell time 30 in each mode, certified with exponential
%     enclosures in both data readings (certify_best_upper) and replayed by uniformization: lambda_check <= -0.299986;
%   - the exact identity g(h) = (((1 - h)(1 + h eps))^(1/3) - 1)/h proved in the revised example, and the
%     floating-point rates of the block words B_1^n B_2^n B_3^n, n = 1000;
%   - closures at h = 0.1 seeded by 1, by points of the two-dimensional faces, and by the coordinate vectors, at the
%     levels g = -1.01 (they finish) and g = -0.99 (they do not; every polytope certificate is capped at 1 - h, so no
%     closure above it can finish).  The closures above the cap are stopped after 'tmax' seconds (120, as in the
%     paper) or 'pmax' processed vertices (20000).
%   Certificate: certificates/bcv/bcv13_schedule.json.  Output: results/ex410.md.  Runtime: a few minutes.

p = inputParser; p.addParameter('tmax', 120); p.addParameter('pmax', 20000);
p.addParameter('outdir', ''); p.parse(varargin{:}); o = p.Results;
root = paperB_root();
if ~isempty(o.outdir), root = char(o.outdir); end
say = logger(fullfile(root, 'results', 'ex410.md'));
A = load_pair('bcv', 'bcv13'); Af = cellfun(@str2double, A, 'UniformOutput', false);
say('# Ex. 4.10 (ex:bcv): the three-mode benchmark, MATLAB reproduction');
say('');
% chattering value over the simplex (floating point)
eta = @(al) max(real(eig(al(1) * Af{1} + al(2) * Af{2} + al(3) * Af{3})));
f = @(z) eta(exp(z - max(z)) / sum(exp(z - max(z))));
c = Inf; o2 = optimset('TolX', 1e-10, 'TolFun', 1e-13, 'MaxIter', 4000, 'MaxFunEvals', 8000, 'Display', 'off');
for z0 = [0 0 0; 0 2 4; 4 0 2; 2 4 0; -3 0 3].'
    [z, fv] = fminsearch(f, z0.', o2); c = min(c, fv);
end
for i = 1:3
    for j = i + 1:3
        for t = linspace(0, 1, 2001)
            al = zeros(1, 3); al(i) = t; al(j) = 1 - t; c = min(c, eta(al));
        end
    end
end
say('- numerical mixture search: c approximately %.6f (the paper gives the exact value 1/120 = %.6f).', c, 1 / 120);
say('- the revised paper proves lambda_check = -0.3 and lsr(B_h) = ((1-h)(1+0.1*h))^(1/3) for 0 < h < 1 (exact decimal model).');
% certified schedule
L = paper_levels('bcv');
sch = struct('modes', [3 2 1], 'durations', {cellstr(string(L.durations)).'}, 'rate_float', L.schedule_float);
file = fullfile(root, 'certificates', 'bcv', 'bcv13_schedule.json');
t0 = tic; [C, json] = certify_best_upper(A, struct('id', 'bcv13', 'chatter', false, 'schedule', sch));
write_text(file, json);
[ok, msgs] = replay_best_upper(file);
say('- Computation: schedule A_3, A_2, A_1 with dwell times %s: lambda_check <= %s (exact bound %s, float rate %.7f), %.0fs; replay %s %s', ...
    strjoin(string(L.durations), ', '), C.printed, fmt_dec(parse_rationals(string(C.U_s)), 9, 'ceil'), ...
    L.schedule_float, toc(t0), tf(ok), strjoin(msgs, '; '));
say('  This finite-schedule certificate is consistent with the proved exponent -0.3.');
% Euler bounds
for h = [0.1 0.05]
    lb = (((1 - h) * (1 + h / 10))^(1 / 3) - 1) / h;
    B = cellfun(@(M) eye(3) + h * M, Af, 'UniformOutput', false);
    n = 1000; P = B{1}^n * B{2}^n * B{3}^n;
    r = (max(abs(eig(P)))^(1 / (3 * n)) - 1) / h;
    say('- h = %g: proved formula g(h) = %.6f (rounded); block word B_1^n B_2^n B_3^n, n = 1000: floating-point (rho^(1/3n) - 1)/h = %.6f', ...
        h, lb, r);
end
say('  (paper: -0.313030 at h = 0.1 and -0.306335 at h = 0.05, approached to six digits by these block words)');
% capped closures
h = 0.1; B = cellfun(@(M) eye(3) + h * M, Af, 'UniformOutput', false);
seeds = struct('ones', ones(3, 1), 'faces', [1 1 0; 1 0 1; 0 1 1].', 'coordinates', eye(3));
say('- closures at h = 0.1 (cap at 1 - h, i.e. g = -1):');
for nm = fieldnames(seeds).'
    for g = [-1.01, -0.99]
        [V, st] = polytope_closure(B, 1 + h * g, seeds.(nm{1}), o.tmax, 20000, o.pmax);
        say('  seeds %s, level g = %.2f: %s (|V| = %d, %d processed, %.1fs)', nm{1}, g, st.reason, st.nV, st.processed, st.secs);
    end
end
say('  (paper: the closures finish at g = -1.01, and none finishes at g = -0.99 within 120 s)');
R = C;
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end
