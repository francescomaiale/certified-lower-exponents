function [val, mu, c, ub] = lp_antinorm(V, w, strict)
% LP_ANTINORM  Polytope antinorm a_V(w) = max{1'mu : V mu <= w, mu >= 0} by linprog (HiGHS), with its dual.
%   [val, mu, c] = lp_antinorm(V, w) for a d x n vertex matrix V >= 0 and w >= 0.  w lies in G = conv V + R^d_+
%   iff a_V(w) >= 1.  Returns
%     val  a conservative value: mu is scaled down until V mu <= w holds in floating point;
%     mu   the weights (n x 1, >= 0);
%     c    a dual vector (d x 1, >= 0) with V'c >= 1 approximately, so a_V(w) <= c'w / min(V'c);
%     ub   that dual bound c'w / min(V'c) (Inf if unavailable), so ub - val measures the LP's optimality gap.
%   Vertices using a coordinate where w vanishes are excluded, columns are
%   normalized to unit 1-norm and rows scaled by 1/w_r.  The value is compared with the dual bound; if they differ
%   by more than 1e-9 (the dual simplex occasionally stops at a suboptimal basis on near-parallel columns), the LP
%   is re-solved by interior point.  If both fail, the legacy interior-point solver and a second row scaling are
%   tried (tiny entries of w make the scaled matrix ill-conditioned).  Then a column generation on small restricted LPs is tried
%   (n > 60); otherwise val = 0 ("not certified inside") is returned, or an error is raised when strict.

if nargin < 3, strict = false; end
persistent OPTS
if isempty(OPTS)
    simplex = 'dual-simplex-highs'; simplex_tol = 1e-10;
    if verLessThan('matlab', '24.1') %#ok<VERLESSMATLAB>  compatibility with R2023a/b
        simplex = 'dual-simplex'; simplex_tol = 1e-9;  % smallest tolerance accepted by the older solver
    end
    OPTS = {optimoptions('linprog', 'Algorithm', simplex, 'Display', 'off', ...
                         'ConstraintTolerance', simplex_tol, 'OptimalityTolerance', 1e-10), ...
            optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off', ...
                         'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10), ...
            optimoptions('linprog', 'Algorithm', 'interior-point-legacy', 'Display', 'off', ...
                         'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10)};
end
[d, n] = size(V); w = w(:); pos = w > 0;
if any(~pos)
    ok = all(V(~pos, :) <= 0, 1);
else
    ok = true(1, n);
end
cols = find(ok);
val = 0; mu = zeros(n, 1); c = []; ub = Inf;
if isempty(cols), return, end
nrm = sum(V(:, cols), 1);
f = -1 ./ nrm(:); np = nnz(pos); nc = numel(cols);
best = [];
% rows scaled by 1/w_r; if both solvers fail (tiny w_r make the scaled matrix ill-conditioned), the rows are
% scaled by 1/max(w_r, 1e-8 max(w)) instead (same LP, other conditioning)
for sv = 1:2
    if sv == 1, sr = w(pos); else, sr = max(w(pos), 1e-8 * max(w)); end
    Aub = V(pos, cols) ./ nrm ./ sr; bub = w(pos) ./ sr;
    for k = 1:numel(OPTS)
        if k == 3 && ~isempty(best), break, end         % the legacy interior point only if both HiGHS solvers failed
        try                                              % HiGHS occasionally raises an internal error instead of
            [x, ~, flag, ~, lam] = linprog(f, Aub, bub, [], [], zeros(nc, 1), [], OPTS{k});   % returning a flag
        catch
            continue
        end
        if flag ~= 1 || isempty(x), continue, end
        m = zeros(n, 1); m(cols) = max(x, 0) ./ nrm(:);
        Vm = V * m; p = Vm > 0;
        if any(p), sc = min(1, min(w(p) ./ Vm(p))); else, sc = 1; end
        v = sc * sum(m) * (1 - 4e-16 * d);
        cc = zeros(d, 1); cc(pos) = max(lam.ineqlin, 0) ./ sr;
        mm = min(cc' * V(:, cols));
        if mm > 0, ub = (cc' * w) / mm; else, ub = Inf; end
        if isempty(best) || v > best{1}
            best = {v, m, cc, ub};
        end
        if v >= ub * (1 - 1e-9)
            break
        end
    end
    if ~isempty(best), break, end
end
if ~isempty(best)
    [val, mu, c, ub] = best{:};
    if nnz(mu) > d + 1                      % an interior-point solution: restore a sparse (basic) solution
        [val, mu] = crossover(V, w, mu, val, OPTS{1});
    end
    [v1, k1] = single_vertex(V(:, cols), w);    % the best single vertex: a lower bound for the LP value that
    if v1 > val                                  % an inaccurate solver answer can miss (typical in closures,
        val = v1 * (1 - 4e-16 * d); mu = zeros(n, 1); mu(cols(k1)) = val;   % where B v = t u for a vertex u)
    end
    return
end
if nc > 60
    % column generation on restricted LPs of at most 40 + 20 k columns
    Vp = V(pos, cols); wp = w(pos);
    Rt = repmat(wp, 1, nc) ./ Vp; Rt(Vp <= 0) = Inf;
    ratio = min(Rt, [], 1);
    [~, o] = sort(ratio, 'descend');
    work = cols(o(1:40));
    for it = 1:100
        work = sort(unique(work));
        [val, mw, c, ub] = lp_antinorm(V(:, work), w);
        if isempty(c), break, end
        red = c' * V(:, cols); mm = min(red);
        if mm >= 1 - 1e-12 || numel(work) >= nc
            mu = zeros(n, 1); mu(work) = mw;
            if mm > 0, c = c / mm; ub = c' * w; end
            return
        end
        [~, o] = sort(red, 'ascend');
        work = [work, cols(o(1:min(20, nc)))]; %#ok<AGROW>
    end
end
if strict
    [val, mu, c, ub] = strict_retry(V, w, cols, OPTS{1});
    if ~isempty(c), return, end
    error('lp_antinorm:failed', 'LP failed (HiGHS and column generation)');
end
val = 0; mu = zeros(n, 1); c = []; ub = Inf;
end


function [val, mu, c, ub] = strict_retry(V, w, cols, opt)
% Last attempts of the exact check (strict) before giving up.  HiGHS dual simplex occasionally reports "lost
% feasibility" (flag -9) on an LP whose zero vector is feasible, depending on the last bits of w (F*V and F*V(:, j)
% round differently); the same LP then solves with the columns in reverse order or with ConstraintTolerance 1e-9.
% Only the exact check uses these attempts, so the closures, and the vertex sets they produce, are unchanged.
[d, n] = size(V); pos = w > 0; nc = numel(cols);
val = 0; mu = zeros(n, 1); c = []; ub = Inf;
nrm = sum(V(:, cols), 1); sr = w(pos);
Aub = V(pos, cols) ./ nrm ./ sr; bub = ones(nnz(pos), 1); f = -1 ./ nrm(:);
opt2 = optimoptions(opt, 'ConstraintTolerance', 1e-9);
tries = {nc:-1:1, opt; 1:nc, opt2; nc:-1:1, opt2};
for k = 1:size(tries, 1)
    r = tries{k, 1};
    try
        [x, ~, flag, ~, lam] = linprog(f(r), Aub(:, r), bub, [], [], zeros(nc, 1), [], tries{k, 2});
    catch
        continue
    end
    if flag ~= 1 || isempty(x), continue, end
    m = zeros(n, 1); m(cols(r)) = max(x, 0) ./ nrm(r)';
    Vm = V * m; p = Vm > 0;
    if any(p), sc = min(1, min(w(p) ./ Vm(p))); else, sc = 1; end
    val = sc * sum(m) * (1 - 4e-16 * d); mu = m;
    c = zeros(d, 1); c(pos) = max(lam.ineqlin, 0) ./ sr;
    mm = min(c' * V(:, cols));
    if mm > 0, ub = (c' * w) / mm; else, ub = Inf; end
    return
end
end


function [val, mu] = crossover(V, w, mu, val, opt)
% Interior-point solutions put tiny positive weights on every vertex (down to 1e-27), which makes the exact check
% slow and the certificates large.  Re-solve the LP by dual simplex on the columns with significant weight and keep
% that sparse solution if its value is within 1e-12 of val; otherwise drop the weights below 1e-15 max(mu), which
% keeps V mu <= w (V >= 0) and changes the value by less than n * 1e-15.
d = size(V, 1);
[~, o] = sort(mu, 'descend'); S = o(mu(o) > 1e-9 * max(mu));
S = S(1:min(numel(S), 4 * d));
pos = w > 0; nrm = sum(V(:, S), 1);
try
    [x, ~, flag] = linprog(-1 ./ nrm(:), V(pos, S) ./ nrm ./ w(pos), ones(nnz(pos), 1), [], [], zeros(numel(S), 1), [], opt);
catch                                       % an internal HiGHS error: fall back to dropping the small weights
    x = []; flag = 0;
end
if flag == 1 && ~isempty(x)
    m = zeros(size(mu)); m(S) = max(x, 0) ./ nrm(:);
    Vm = V * m; p = Vm > 0;
    if any(p), sc = min(1, min(w(p) ./ Vm(p))); else, sc = 1; end
    v = sc * sum(m) * (1 - 4e-16 * d);
    if v >= val * (1 - 1e-12), val = v; mu = m * sc; return, end
end
m = mu; m(m < 1e-15 * max(m)) = 0;
Vm = V * m; p = Vm > 0;
if any(p), sc = min(1, min(w(p) ./ Vm(p))); else, sc = 1; end
mu = m * sc; val = sc * sum(m) * (1 - 4e-16 * d);
end


function [v, k] = single_vertex(V, w)
% max over columns u of max{s : s u <= w}, the best LP value supported by one vertex
pos = w > 0; R = w(pos) ./ V(pos, :); R(V(pos, :) <= 0) = Inf;
r = min(R, [], 1); r(any(V(~pos, :) > 0, 1)) = 0;
[v, k] = max(r);
if ~isfinite(v), v = 0; end
end
