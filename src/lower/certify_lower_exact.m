function [C, json, info] = certify_lower_exact(V, D, opts)
% CERTIFY_LOWER_EXACT  Exact polytope certificate (prop:polytope) from a closed floating-point polytope.
%   [C, json, info] = certify_lower_exact(V, D, opts)
%     V     d x n vertices, >= 0: doubles (from polytope_closure; their exact binary64 values are the vertices)
%           or sym rationals (e.g. the vertices of an earlier certificate, for prop:polytope(ii))
%     D     exact data from exact_data (both readings of A_i, exact h)
%     opts  name ('pair'), kind ('best'), orientation ('N' or 'T'), meta (struct, stored as is),
%           lambda_claim (decimal string: certify exactly this six-decimal bound instead of floor_6),
%           polish (true): exact basis polishing of the LP weights of the (near-)binding pairs, see below
%   Returns the certificate C (decoded JSON, in the format read by the replays of reproduce/, MATLAB and
%   Python), its JSON text, and info (exact t_star, lambda, timings, polished pairs).
%
%   Certificate: for every vertex v and mode i, weights lambda >= 0 with sum 1 and B_i v >= t V lambda entrywise,
%   B_i = I + h A_i.  Then lsr({B_i}) >= t and lambda_check >= (t - 1)/h.
%   Construction: for each pair (i, v) the LP weights mu of
%   a_V(B_i v) (lp_antinorm) are rationalized (exact binary64 values) and normalized exactly; t_star is the minimum,
%   over pairs, data readings and coordinates, of (B_i v)_r / (V lambda)_r, in exact rational arithmetic.  MATLAB's
%   linprog solutions are sometimes less sharp than scipy's (1e-8 relative), which can cost a unit in the sixth
%   decimal; so the weights of every pair whose level is within 1e-8 of t_star and whose LP has a duality gap
%   above 1e-13 (relative) are polished: the active rows R
%   and the support S of mu define the LP vertex V(R,S) mu_S = (B_i v)_R, which is solved exactly and kept if it
%   is exactly feasible (mu_S >= 0, V mu <= B_i v) and improves the level.  This is repeated until the binding
%   pair is polished, or until the LP dual bounds show that the six-decimal value cannot improve (or, with
%   lambda_claim, until the claim is reached).  Only exactness of the final check matters, not how the weights were
%   found.
%   The claimed bound is lambda = floor_6((t_star - 1)/h) and t = 1 + h lambda <= t_star.  For the binding pair,
%   the LP dual c >= 0, rescaled so that c'u >= 1 for every vertex u, gives tau_up = c'(B_i v): no t > tau_up is
%   supported by this polytope (negative control).

if nargin < 3, opts = struct(); end
name   = getopt(opts, 'name', 'pair');
kind   = getopt(opts, 'kind', 'best');
orient = getopt(opts, 'orientation', 'N');
meta   = getopt(opts, 'meta', struct());
claim  = getopt(opts, 'lambda_claim', '');
do_pol = getopt(opts, 'polish', true);
verbose = getopt(opts, 'verbose', false);
T0 = tic;
d = D.d; m = D.m; n = size(V, 2); I = eye(d);
if ~ismember(string(orient), ["N", "T"])
    error('certify_lower_exact:orientation', 'orientation must be N or T');
end
if size(V, 1) ~= d || n == 0 || ~isreal(V) || any(~isfinite(double(V(:)))) || ...
        any(double(sign(V(:))) < 0) || any(double(sign(sum(V, 1))) <= 0)
    error('certify_lower_exact:vertices', 'vertices must be finite, nonnegative and nonzero, with d rows');
end
check_euler_step(D, 'explicit');
Ff = cell(1, m);
for i = 1:m
    Ff{i} = I + D.hf * D.Af{i};
    if orient == 'T', Ff{i} = Ff{i}.'; end
end
nv = numel(D.ver); Fq = cell(nv, m);
for v = 1:nv
    for i = 1:m
        B = sym(I) + D.hq * D.ver(v).A{i};
        if orient == 'T', B = B.'; end
        Fq{v, i} = B;
    end
end
if isa(V, 'sym'), Vq = V; else, Vq = sym(V, 'f'); end
Vf = double(Vq);
% ---------------------------------------------------------------- LP weights for every (mode, vertex) pair
K = cell(1, m); W = cell(1, m); gap = zeros(m, n); UB = inf(m, n);
for i = 1:m
    FV = Ff{i} * Vf; Ki = cell(1, n); Wi = cell(1, n);
    for j = 1:n
        [val, mu, ~, ub] = lp_antinorm(Vf, FV(:, j), true);
        gap(i, j) = (ub - val) / max(val, realmin); UB(i, j) = ub;
        s = find(mu > 0);
        if isempty(s), error('certify_lower_exact:lp', 'empty LP weights at mode %d, vertex %d', i, j); end
        Ki{j} = s; Wi{j} = mu(s);
    end
    [K{i}, W{i}] = pad(Ki, Wi);
end
t_lp = toc(T0);
if verbose, fprintf('   LP weights for %d pairs: %.1fs\n', m * n, t_lp); end
% exact weights, normalized
Wq = cell(1, m);
for i = 1:m, Wq{i} = sym(W{i}, 'f'); Wq{i} = Wq{i} ./ sum(Wq{i}, 1); end
BV = cell(nv, m);
for v = 1:nv, for i = 1:m, BV{v, i} = Fq{v, i} * Vq; end, end
polished = false(m, n); npol = 0; nfail = 0;
for it = 1:30
    [tstar, bind, Tpair] = levels(Vq, K, Wq, BV, d, n, m, nv);
    if ~do_pol || polished(bind(1), bind(2)) || gap(bind(1), bind(2)) <= 1e-13, break, end
    % stop when polishing cannot change the six decimals: the LP dual bounds cap every pair's level
    lam_now = floor(double((tstar - 1) / D.hq) * 1e6 + 1e-9);
    lam_cap = floor((min(UB(:)) * (1 - 1e-13) - 1) / D.hf * 1e6);
    if isempty(claim) && lam_now >= lam_cap, break, end
    if ~isempty(claim) && double(sign(tstar - (1 + D.hq * parse_rationals(string(claim))))) >= 0, break, end
    near = find(~polished & gap > 1e-13 & Tpair <= double(tstar) * (1 + 1e-8) + 1e-15);
    [~, o] = sort(Tpair(near)); near = near(o(1:min(64, numel(o))));       % the lowest ones first
    if isempty(near), near = sub2ind([m n], bind(1), bind(2)); end
    if verbose
        fprintf('   polish round %d: (t*-1)/h = %.12f, binding pair (%d, %d), %d pairs to polish, %.1fs\n', it, ...
                double((tstar - 1) / D.hq), bind(1), bind(2), numel(near), toc(T0));
    end
    for q = near(:).'
        [i, j] = ind2sub([m n], q);
        polished(i, j) = true;
        lamq = polish(Vq, Vf, Ff{i} * Vf(:, j), BV(:, i), j, K{i}(:, j), W{i}(:, j), Tpair(i, j));
        if isempty(lamq), nfail = nfail + 1; continue, end
        npol = npol + 1;
        [Kc, Wc] = deal(lamq.K(:), lamq.W(:));
        if numel(Kc) > size(K{i}, 1)
            K{i}(end + 1:numel(Kc), :) = 1; W{i}(end + 1:numel(Kc), :) = 0;
            Wq{i}(end + 1:numel(Kc), :) = sym(0);
        end
        K{i}(:, j) = 1; W{i}(:, j) = 0; Wq{i}(:, j) = sym(0);
        K{i}(1:numel(Kc), j) = Kc; W{i}(1:numel(Kc), j) = double(Wc); Wq{i}(1:numel(Kc), j) = Wc;
    end
end
t_ex = toc(T0) - t_lp;
% ---------------------------------------------------------------- dual certificate for the binding pair
i = bind(1); j = bind(2); v = bind(3);
[~, ~, c] = lp_antinorm(Vf, Ff{i} * Vf(:, j), true);
cq = sym(max(c, 0), 'f');
Bv = BV{v, i}(:, j);
for r = 1:d
    if double(sign(Bv(r))) == 0            % coordinates where B_i v vanishes: any weight is free, make it large
        pr = Vf(r, :) > 0;
        if any(pr), cq(r) = exact_max([cq(r); (1 ./ Vq(r, pr)).']); end
    end
end
mm = exact_min((cq.' * Vq).');
if double(sign(mm)) > 0
    cq = cq / mm; tau_up = cq.' * Bv; cstr = cellstr(string(cq(:).'));
else
    tau_up = tstar; cstr = NaN;             % no usable dual (should not happen)
end
% ---------------------------------------------------------------- claimed bound
if isempty(claim)
    lam = floor(((tstar - 1) / D.hq) * sym(10)^6) / sym(10)^6;
else
    lam = parse_rationals(string(claim));
end
t = 1 + D.hq * lam;
if ~isAlways(t > 0), error('certify_lower_exact:level', 'the certified level t must be positive'); end
if double(sign(tstar - t)) < 0
    error('certify_lower_exact:claim', 'claimed bound %s is above the exact level of the polytope (t* - 1)/h = %.9f', ...
          fmt_dec(lam, 6, 'floor'), double((tstar - 1) / D.hq));
end
% ---------------------------------------------------------------- certificate (format of the replays)
vs = string(Vq);
verts = cell(n, 1);
for k = 1:n, verts{k} = cellstr(vs(:, k).'); end
wts = cell(1, m * n);
for i2 = 1:m
    ws = string(Wq{i2});
    for k = 1:n
        sk = W{i2}(:, k) > 0;
        keys = cellstr(string(K{i2}(sk, k) - 1));
        wts{(i2 - 1) * n + k} = struct('mode', i2 - 1, 'vertex', k - 1, ...
            'lambda', containers.Map(keys(:).', cellstr(ws(sk, k)).'));
    end
end
lamstr = fmt_dec(lam, 6);
if orient == 'T', tr = '^T'; else, tr = ''; end
meta.nV = n; meta.lp_secs = round(t_lp, 1); meta.exact_secs = round(toc(T0) - t_lp, 1);
meta.polished_pairs = npol; meta.implementation = 'matlab';
S = struct();
S.row = sprintf('%s_%s', name, D.h_str); S.kind = kind; S.('case') = name; S.h = char(D.h_str);
S.orientation = orient;
S.statement = sprintf('lambda_check >= %s  (lsr(B_h) >= t, B_h = I + h A_i%s)', lamstr, tr);
S.t = char(t); S.lambda_lower = lamstr; S.lambda_lower_exact = char((t - 1) / D.hq);
S.t_star = char(tstar); S.t_star_float = double(tstar);
S.binding = struct('mode', bind(1) - 1, 'vertex', bind(2) - 1, 'data', D.ver(bind(3)).name);
S.negative_control = struct('c', {cstr}, 'tau_up', char(tau_up), 'tau_up_float', double(tau_up));
S.data_versions = {D.ver.name};
S.meta = meta;
S.vertices = verts;
S.weights = wts;
json = jsonencode(S);
C = jsondecode(json);
info = struct('t_star', tstar, 'lambda', lam, 'lambda_star', (tstar - 1) / D.hq, 't', t, 'tau_up', tau_up, ...
              'binding', bind, 'nV', n, 'polished', npol, 'polish_failed', nfail, 'lp_secs', t_lp, ...
              'exact_secs', t_ex, 'secs', toc(T0));
end


function [Kp, Wp] = pad(Ki, Wi)
n = numel(Ki); smax = max(cellfun(@numel, Ki));
Kp = ones(smax, n); Wp = zeros(smax, n);
for j = 1:n
    Kp(1:numel(Ki{j}), j) = Ki{j}; Wp(1:numel(Ki{j}), j) = Wi{j};
end
end


function [tstar, bind, Tpair] = levels(Vq, K, Wq, BV, d, n, m, nv)
% exact t_star (min over pairs, readings, coordinates of (B_i v)_r / (V lambda)_r), the binding pair
% [mode, vertex, reading], and the double levels of all pairs (m x n)
tstar = []; bind = []; Tpair = inf(m, n);
for i = 1:m
    VL = sym(zeros(d, n));
    for p = 1:size(K{i}, 1)
        VL = VL + Vq(:, K{i}(p, :)) .* Wq{i}(p, :);
    end
    Z = double(sign(VL)) == 0;
    for v = 1:nv
        R = BV{v, i} ./ (VL + Z);                 % ratios; Z marks (V lambda)_r = 0 (no constraint)
        bad = any(~Z & double(sign(BV{v, i})) <= 0, 1);
        Rd = double(R); Rd(Z) = Inf; Rd(:, bad) = 0;
        Tpair(i, :) = min(Tpair(i, :), min(Rd, [], 1));
        if any(bad)
            tb = sym(0); jb = find(bad, 1);
        else
            [tb, ib] = exact_min(R(~Z));
            jj = repmat(1:n, d, 1); jj = jj(~Z); jb = jj(ib);
        end
        if isempty(tstar) || double(sign(tb - tstar)) < 0
            tstar = tb; bind = [i, jb, v];
        end
    end
end
end


function out = polish(Vq, Vf, w, BVi, j, Kj, Wj, tcur)
% Exact optimum of the LP of the pair (i, j) restricted to the support S of its float weights: every basic
% solution V(R,S') mu = (B_i v_j)_R with S' a subset of S (all subsets if |S| <= 5, else single vertices and S)
% or one of the five best single vertices of the polytope, and R a set of |S'| rows among the active rows of the
% float solution (all rows for single vertices), is solved exactly; the exactly
% feasible one (mu >= 0, V(:,S') mu <= B_i v_j) of largest level, over all data readings, is returned if it beats
% tcur.  (The float LP solutions differ from the exact optimum by up to 1e-8 relative.)
out = [];
[~, o] = sort(Wj, 'descend'); S = Kj(o(Wj(o) > 0)).';
d = size(Vf, 1);
if numel(S) > d, S = S(1:d); end
if numel(S) <= 5
    subsets = {};
    for k = 1:numel(S), subsets = [subsets, num2cell(combs(S, k), 2).']; end %#ok<AGROW>
else
    subsets = [num2cell(S(:)).', {S}];
end
% single vertices of the whole polytope with the best float levels min_r w_r / v_r
R = w ./ Vf; R(Vf <= 0) = Inf; r = min(R, [], 1); r(any(Vf(w <= 0, :) > 0, 1)) = 0;
[~, o1] = sort(r, 'descend'); subsets = [num2cell(o1(1:min(5, end))), subsets];
mu = zeros(size(Vf, 2), 1); mu(Kj(Wj > 0)) = Wj(Wj > 0);
rows = find(w > 0 & w - Vf * mu <= 1e-9 * w).';      % active rows of the float solution
best = []; tbest = tcur * (1 - 1e-15);
for u = 1:numel(subsets)
    Sp = subsets{u}; k = numel(Sp);
    if k == 1, rr = find(w > 0).'; else, rr = rows; end
    if numel(rr) < k || nchoosek(numel(rr), k) > 400, continue, end
    Rs = combs(rr, k);
    for q = 1:size(Rs, 1)
        R = Rs(q, :);
        if rank(Vf(R, Sp)) < k, continue, end
        for v = 1:numel(BVi)
            b = BVi{v}(:, j);
            x = Vq(R, Sp) \ b(R);
            if any(double(sign(x)) < 0) || any(double(sign(b - Vq(:, Sp) * x)) < 0), continue, end
            if double(sum(x)) < tbest * (1 - 1e-12), continue, end
            lam = x / sum(x); VL = Vq(:, Sp) * lam; P = double(sign(VL)) > 0; tl = [];
            for v2 = 1:numel(BVi)
                b2 = BVi{v2}(:, j);
                if any(double(sign(b2(P))) <= 0), tl = sym(0); break, end
                t2 = exact_min(b2(P) ./ VL(P));
                if isempty(tl) || double(sign(t2 - tl)) < 0, tl = t2; end
            end
            if double(tl) > tbest && (isempty(best) || double(sign(tl - best.t)) > 0)
                best = struct('K', Sp(:), 'W', lam, 't', tl); tbest = double(tl);
            end
        end
    end
end
out = best;
end


function c = combs(v, k)
% all k-subsets of the vector v, one per row (nchoosek, also for a single element)
if numel(v) == 1
    if k == 1, c = v; else, c = zeros(0, k); end
else
    c = nchoosek(v, k);
end
end
