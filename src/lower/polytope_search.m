function S = polytope_search(F, h, t_up, opts)
% POLYTOPE_SEARCH  Bisection on the level t of polytope_closure (the "Search" step of the construction in Sec. 4).
%   S = polytope_search(F, h, t_up, opts)
%     F     cell array of the explicit Euler matrices B_i = I + h A_i (or their transposes), floating point
%     h     step (double), used only for the stopping rule and the printed lambda values (t - 1)/h
%     t_up  upper end of the bisection, e.g. min(1 + h c, product upper bound for lsr(B_h))
%   opts fields (all optional):
%     seeds    d x s starting vertices (default eye(d))
%     t_start  a level to probe first, e.g. the GM lower estimate of lsr(B_h) (default none)
%     probe    time cap of one closure in s (default 10)
%     budget   time budget of the bisection in s (default 60)
%     nmax     vertex cap of one closure (default 6000)
%     pmax     cap on the vertices processed by one closure (default Inf); with pmax set and generous time caps
%              the search is deterministic (the same probes close on every machine)
%     tol      stop when hi - lo <= tol * h (default 2e-7)
%     maxprobes  stop after this many bisection probes (default Inf); with pmax, a budget that does not depend
%              on the machine
%     verbose  print one line per probe (default true)
%   Returns S with fields V (best closed float polytope), t_build (its level), t_fail (lowest failed level, or
%   t_up), t_triv, probes (struct array: lambda, result, nV, processed, secs), secs.
%
%   The lower end is a level that closes: t_start if its
%   closure finishes, otherwise the trivial level t_triv = smallest column sum of the F_i (V = {e_1, ..., e_d}).
%   A failed probe means only that the closure did not finish within its caps.

if nargin < 4, opts = struct(); end
d = size(F{1}, 1);
seeds  = getopt(opts, 'seeds', eye(d));
probe  = getopt(opts, 'probe', 10);
budget = getopt(opts, 'budget', 60);
nmax   = getopt(opts, 'nmax', 6000);
pmax   = getopt(opts, 'pmax', Inf);
tol    = getopt(opts, 'tol', 2e-7);
maxprobes = getopt(opts, 'maxprobes', Inf);
verbose = getopt(opts, 'verbose', true);
t_start = getopt(opts, 't_start', []);

t_triv = Inf;
for i = 1:numel(F), t_triv = min(t_triv, min(sum(F{i}, 1))); end
T0 = tic; probes = struct('lambda', {}, 'result', {}, 'nV', {}, 'processed', {}, 'secs', {});

    function [V, st] = attempt(tb, cap)
        [V, st] = polytope_closure(F, tb, seeds, cap, nmax, pmax);
        why = st.reason;
        probes(end + 1) = struct('lambda', (tb - 1) / h, 'result', why, 'nV', st.nV, 'processed', st.processed, ...
                                 'secs', st.secs);
        if verbose
            fprintf('   probe lambda = %+.9f (t = %.12f): %s, |V| = %d, %d processed, %.1fs\n', (tb - 1) / h, tb, ...
                    why, st.nV, st.processed, st.secs);
        end
    end

V = [];
if ~isempty(t_start) && t_start < t_up
    V = attempt(t_start, min(probe, budget));
end
if ~isempty(V)
    lo = t_start; hi = t_up; Vbest = V;
else
    Vbest = attempt(t_triv, max(probe, 60));
    if isempty(Vbest)
        error('polytope_search:start', ['the closure at the trivial level t = %.12f did not finish; pass ' ...
              'opts.t_start, a level whose closure finishes (needed e.g. for seeds other than eye(d))'], t_triv);
    end
    lo = t_triv; hi = t_up;
    if ~isempty(t_start), hi = min(hi, t_start); end
end
nprobe = 0;
while hi - lo > tol * h && nprobe < maxprobes
    nprobe = nprobe + 1;
    remaining = budget - toc(T0);
    if remaining < 1, break, end
    mid = 0.5 * (lo + hi);
    V = attempt(mid, min(probe, remaining));
    if ~isempty(V)
        lo = mid; Vbest = V;
    else
        hi = mid;
    end
end
S = struct('V', Vbest, 't_build', lo, 't_fail', hi, 't_triv', t_triv, 't_up', t_up, 'probes', probes, ...
           'secs', toc(T0));
end


function v = getopt(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
