function [alpha, alpha_float, info] = tangent_bound(V, A, versions)
% TANGENT_BOUND  Exact tangent-cone lower bound alpha(Q) of [GLP, Prop. 8] for Q = conv V + K (Obs. 4.11).
%   [alpha, alpha_float] = tangent_bound(V, A) with V a d x n vertex matrix (doubles: their exact binary64
%   values are used) and A the data (cell of decimal-string or numeric matrices).  alpha (sym) is the largest
%   rational alpha certified such that, at every vertex w and for every mode i, A_i w - alpha w lies in
%   cone{v - w : v in V} + K, for both data readings (versions, default {'decimal', 'binary64'}).  Then
%   (I + delta (A_i - alpha I)) Q is contained in Q for small delta > 0, and prop:polytope(i), applied to the shifted
%   family with t = 1, gives lambda_check >= alpha.  Construction:
%     LP per (i, w):  max alpha  s.t.  A_i w - alpha w - sum_v mu_v (v - w) >= 0,  mu >= 0,  |alpha| <= 1e3;
%   then mu is rationalized (exact binary64), scaled down where needed on coordinates with w_r = 0 (where the
%   condition does not involve alpha), and alpha = min over w_r > 0 of (A_i w - comb)_r / w_r, exactly.
%   alpha_float is the floating-point LP minimum.  If every LP solver fails at a vertex, mu = 0 is used there (feasible,
%   so the bound stays valid, only weaker); info.lp_fallbacks counts such vertices.

if nargin < 3, versions = {'decimal', 'binary64'}; end
D = exact_data(A, '1', versions);
d = size(V, 1); n = size(V, 2);
if d ~= D.d || n == 0 || ~isreal(V) || any(~isfinite(double(V(:)))) || ...
        any(double(sign(V(:))) < 0) || any(double(sign(sum(V, 1))) <= 0)
    error('tangent_bound:vertices', 'vertices must be nonzero, finite and nonnegative, with d rows');
end
Vq = sym(V, 'f'); Vf = double(Vq);
simplex = 'dual-simplex-highs';
if verLessThan('matlab', '24.1'), simplex = 'dual-simplex'; end %#ok<VERLESSMATLAB>  compatibility with R2023a/b
opts = {optimoptions('linprog', 'Algorithm', simplex, 'Display', 'off'), ...
        optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off'), ...
        optimoptions('linprog', 'Algorithm', 'interior-point-legacy', 'Display', 'off')};
nfall = 0;
best = []; worst = Inf; T0 = tic;
for i = 1:D.m
    Af = D.Af{i};
    for j = 1:n
        w = Vf(:, j); others = [1:j - 1, j + 1:n];
        Dm = Vf(:, others) - w;
        Aub = [w, Dm]; bub = Af * w;
        f = [-1; zeros(n - 1, 1)];
        flag = 0;
        for k = 1:numel(opts)
            try
                [x, fv, flag] = linprog(f, Aub, bub, [], [], [-1e3; zeros(n - 1, 1)], [1e3; inf(n - 1, 1)], opts{k});
            catch                                  % an internal HiGHS error counts as a failed solve
                flag = 0;
            end
            if flag == 1, break, end
        end
        if flag == 1
            worst = min(worst, -fv);
            sel = find(x(2:end) > 1e-15 * max(x(2:end))); ks = others(sel);   % (drops interior-point dust)
            mu = sym(x(1 + sel), 'f');
        else                                       % every solver failed: mu = 0 is feasible (weaker, still valid)
            nfall = nfall + 1; ks = []; mu = sym(zeros(0, 1));
            p = Vf(:, j) > 0; awf = Af * Vf(:, j); worst = min(worst, min(awf(p) ./ Vf(p, j)));
        end
        wq = Vq(:, j);
        for v = 1:numel(D.ver)
            Aw = D.ver(v).A{i} * wq;
            if isempty(ks), comb = sym(zeros(d, 1)); else, comb = (Vq(:, ks) - repmat(wq, 1, numel(ks))) * mu; end
            s = sym(1);
            for r = 1:d                                % coordinates with w_r = 0: alpha-free rows
                if double(sign(wq(r))) == 0 && double(sign(Aw(r) - comb(r))) < 0
                    neg = Vq(r, ks) * mu;           % (comb(r) > Aw(r) >= 0, so ks is not empty here)
                    if double(sign(Aw(r) / neg - s)) < 0, s = Aw(r) / neg; end
                end
            end
            comb = s * comb;
            p = double(sign(wq)) > 0;
            a = exact_min((Aw(p) - comb(p)) ./ wq(p));
            if isempty(best) || double(sign(a - best)) < 0, best = a; end
        end
    end
end
alpha = best; alpha_float = worst;
info = struct('nV', n, 'secs', toc(T0), 'lp_fallbacks', nfall);
end
