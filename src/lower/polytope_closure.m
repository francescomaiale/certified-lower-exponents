function [V, st] = polytope_closure(F, t, seeds, tmax, nmax, pmax, fast)
% POLYTOPE_CLOSURE  Floating-point t-invariant polytope for the family F (closure of the seeds under v -> F_i v / t).
%   [V, st] = polytope_closure(F, t, seeds, tmax, nmax, pmax)
%     F      cell array of nonnegative d x d matrices (B_i = I + h A_i, or their transposes)
%     t      target level
%     seeds  d x s matrix of nonnegative starting vertices (default: eye(d), the coordinate vectors)
%     tmax   time cap in seconds (default Inf), nmax cap on the number of vertices (default 6000),
%     pmax   cap on the number of processed vertices (default Inf; unlike tmax it does not depend on the machine)
%   V is the final d x n vertex matrix, or [] if the closure hit a cap; st has fields processed, nV, secs,
%   reason ('closed', 'time', 'size' or 'work'), lps (linprog calls) and skipped (LPs answered by a cached
%   dual direction).
%   If V is returned, every F_i v / t lies (in floating point) in G = conv V + R^d_+, i.e. V satisfies the
%   polytope certificate (*) of prop:polytope at level t in floating point; certify_lower_exact checks it exactly.
%
%   Closure: a vertex w = F_i v / t (entries below 1e-12 max(w) set to 0, which only
%   shrinks w) is added when lp_antinorm(V, w) < 1 + 1e-10; whenever the vertex count exceeds 1.5 times its value at
%   the last pruning, vertices interior to the hull of the others (a_{V\v}(v) >= 1 + 1e-9) are discarded.
%
%   linprog is slower per call in MATLAB than scipy's HiGHS, so (unless fast = false) LPs whose answer is already
%   implied by a stored dual vector c >= 0 are skipped; this does not change any decision, because
%   a_V(x) <= c'x / min_u c'u for every c >= 0:
%     - add w without an LP if c'w < (1 + 1e-10) min_u c'u for a stored facet direction c (then a_V(w) < 1 + 1e-10);
%     - keep v in a pruning without an LP if its last exposing direction c (the dual of its last pruning LP) still
%       gives c'v < (1 + 1e-9) min_{u ~= v} c'u.

ADD_TOL = 1e-10; ZERO_REL = 1e-12; NF = 500;
d = size(F{1}, 1);
if nargin < 3 || isempty(seeds), seeds = eye(d); end
if nargin < 4 || isempty(tmax), tmax = Inf; end
if nargin < 5 || isempty(nmax), nmax = 6000; end
if nargin < 6 || isempty(pmax), pmax = Inf; end
if nargin < 7 || isempty(fast), fast = true; end
t0 = tic;
S = seeds ./ sum(seeds, 1);
cap = max(1024, 4 * size(S, 2));
Vall = zeros(d, cap); Vall(:, 1:size(S, 2)) = S; nall = size(S, 2);
alive = false(1, cap); alive(1:nall) = true;
E = zeros(d, cap); hasE = false(1, cap);          % exposing direction of each vertex (from its pruning LP)
order = 1:nall;                                   % alive vertices, in insertion order
queue = zeros(1, cap); queue(1:nall) = 1:nall; qh = 1; qt = nall;
Cf = zeros(d, NF); minCf = zeros(NF, 1); nf = 0; pf = 0;   % facet directions (duals of "outside" LPs)
n_prune = numel(order); processed = 0; cnt = [0 0];         % [linprog calls, skipped LPs]
while qh <= qt
    j = queue(qh); qh = qh + 1;
    if ~alive(j), continue, end
    processed = processed + 1;
    M = Vall(:, order);
    for i = 1:numel(F)
        w = F{i} * Vall(:, j) / t;
        w(w < ZERO_REL * max(w)) = 0;
        if any(all(M * (1 + ADD_TOL) <= w, 1))
            continue                                   % w dominates a vertex: inside G
        end
        c = [];
        add = fast && nf > 0 && any(Cf(:, 1:nf).' * w < (1 + ADD_TOL) * minCf(1:nf));
        if add
            cnt(2) = cnt(2) + 1;
        else
            [val, ~, c] = lp_antinorm(M, w); cnt(1) = cnt(1) + 1;
            add = val < 1 + ADD_TOL;
        end
        if add
            if nall == size(Vall, 2)                    % grow the storage
                Vall = [Vall, zeros(d, nall)]; alive = [alive, false(1, nall)]; %#ok<AGROW>
                E = [E, zeros(d, nall)]; hasE = [hasE, false(1, nall)]; %#ok<AGROW>
                queue = [queue, zeros(1, nall)]; %#ok<AGROW>
            end
            nall = nall + 1; Vall(:, nall) = w; alive(nall) = true;
            order(end + 1) = nall; qt = qt + 1; queue(qt) = nall; %#ok<AGROW>
            M = [M, w]; %#ok<AGROW>
            if fast
                minCf(1:nf) = min(minCf(1:nf), Cf(:, 1:nf).' * w);
                if ~isempty(c) && any(c > 0)           % the dual of the LP that added w: a facet direction
                    pf = mod(pf, NF) + 1; nf = max(nf, pf);
                    Cf(:, pf) = c; minCf(pf) = min(c.' * M);
                end
            end
        end
    end
    if numel(order) > 1.5 * n_prune
        [keep, E(:, order), hasE(order), cnt] = prune(M, E(:, order), hasE(order), fast, cnt);
        alive(order(~keep)) = false; order = order(keep); n_prune = numel(order);
        if fast && nf > 0, minCf(1:nf) = min(Cf(:, 1:nf).' * Vall(:, order), [], 2); end
    end
    if qh <= qt && (toc(t0) > tmax || numel(order) > nmax || processed >= pmax)
        if numel(order) > nmax, why = 'size'; elseif processed >= pmax, why = 'work'; else, why = 'time'; end
        V = []; st = struct('processed', processed, 'nV', numel(order), 'secs', toc(t0), 'reason', why, ...
                            'lps', cnt(1), 'skipped', cnt(2));
        return
    end
end
M = Vall(:, order);
[keep, ~, ~, cnt] = prune(M, E(:, order), hasE(order), fast, cnt);
V = M(:, keep);
st = struct('processed', processed, 'nV', size(V, 2), 'secs', toc(t0), 'reason', 'closed', 'lps', cnt(1), ...
            'skipped', cnt(2));
end


function [keep, E, hasE, cnt] = prune(V, E, hasE, fast, cnt)
% Logical mask of the vertices kept: drop v when it lies in conv(others) + R^d_+ with margin 1e-9
% (the antiball is unchanged); large vertices are tried first.
PRUNE_TOL = 1e-9;
n = size(V, 2); keep = true(1, n);
[~, o] = sort(sum(V, 1), 'descend');
for k = o
    others = keep; others(k) = false;
    if ~any(others), continue, end
    if fast && hasE(k)
        ce = E(:, k).';
        if ce * V(:, k) < (1 + PRUNE_TOL) * min(ce * V(:, others))
            cnt(2) = cnt(2) + 1; continue           % still exposed by its last direction: kept
        end
    end
    [val, ~, c] = lp_antinorm(V(:, others), V(:, k)); cnt(1) = cnt(1) + 1;
    if val >= 1 + PRUNE_TOL
        keep(k) = false;
    elseif ~isempty(c)
        E(:, k) = c; hasE(k) = true;
    end
end
end
