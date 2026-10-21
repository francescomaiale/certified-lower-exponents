function R = glp_algorithm_l(A, tau, opts)
% GLP_ALGORITHM_L  Naive Algorithm (L) of Guglielmi, Laglia and Protasov [GLP, Sec. 4.6] (Obs. 4.11).
%   R = glp_algorithm_l(A, tau, opts)
%     A     cell {A_1, ..., A_m} (decimal strings or numeric); tau dwell time (double)
%     opts  L       exhaustion length (16, 14, 12 for d <= 3, d = 5, d = 10)
%           nu      run a single closure at this shift (reproduction of a recorded search), else search
%           probe   time cap of one closure in s (300); budget  time budget of the search in s (1800)
%           pmax    cap on the vertices processed by one closure (Inf); nmax vertex cap (6000)
%           exact   compute the exact tangent bound (true)
%   Steps: (1) among all products of the e^{tau A_i} of length <= L, the product Pi of length n minimizing
%   rho(Pi)^(1/n) gives the upper bound beta = log rho(Pi) / (n tau) (floating point); (2) the leading eigenvector of
%   Pi is closed under the maps e^{tau A_i} e^{-tau (beta - nu)} (polytope_closure at level 1); the search tries
%   nu = 0, then 0.05, 0.1, 0.2, ... until a closure finishes, then bisects down to 1e-6; (3) the lower bound is the
%   exact tangent bound alpha(Q) of the last closed polytope (tangent_bound), floored to six decimals.
%   R has fields tau, word_len, beta_l, nu, probes, nV, alpha (sym), alpha_str, alpha_float, secs, ...

if nargin < 3, opts = struct(); end
T0 = tic;
D = exact_data(A, '1', {'binary64'});
d = D.d; m = D.m;
E = cellfun(@(M) expm(tau * M), D.Af, 'UniformOutput', false);
Ldef = 12; if d <= 3, Ldef = 16; elseif d == 5, Ldef = 14; end
[bw, bv] = brute_words(E, getopt(opts, 'L', Ldef));
[val, k] = min(bv); word = bw{k};
beta = val / tau;
P = eye(d);
for a = word, P = E{a} * P; P = P / max(abs(P(:))); end
[EV, ev] = eig(P); [~, k] = max(abs(diag(ev)));
v1 = abs(real(EV(:, k))); v1 = v1 / sum(v1);
search_secs = toc(T0);
probe = getopt(opts, 'probe', 300); budget = getopt(opts, 'budget', 1800);
pmax = getopt(opts, 'pmax', Inf); nmax = getopt(opts, 'nmax', 6000);
probes = struct('nu', {}, 'closed', {}, 'nV', {}, 'secs', {});
    function Q = closure(nu, cap)
        F = cellfun(@(M) M * exp(-tau * (beta - nu)), E, 'UniformOutput', false);
        [Q, st] = polytope_closure(F, 1.0, v1, cap, nmax, pmax);
        probes(end + 1) = struct('nu', nu, 'closed', ~isempty(Q), 'nV', st.nV, 'secs', st.secs);
    end
T1 = tic; Qbest = []; nu_ok = [];
if isfield(opts, 'nu') && ~isempty(opts.nu)
    Qbest = closure(opts.nu, Inf);
    if ~isempty(Qbest), nu_ok = opts.nu; end
else
    Q = closure(0, min(probe, budget));
    if ~isempty(Q)
        nu_ok = 0; Qbest = Q;
    else
        nu_bad = 0; nu = 0.05;
        while toc(T1) < budget - 5
            Q = closure(nu, min(probe, budget - toc(T1)));
            if ~isempty(Q), nu_ok = nu; Qbest = Q; break, end
            nu_bad = nu; nu = 2 * nu;
        end
        while ~isempty(nu_ok) && nu_ok - nu_bad > 1e-6 && toc(T1) < budget - 5
            mid = 0.5 * (nu_ok + nu_bad);
            Q = closure(mid, min(probe, budget - toc(T1)));
            if ~isempty(Q), nu_ok = mid; Qbest = Q; else, nu_bad = mid; end
        end
    end
end
closure_secs = toc(T1);
R = struct('tau', tau, 'word', sprintf('%d', word), 'word_len', numel(word), 'beta_l', beta, 'nu', nu_ok, ...
           'probes', probes, 'search_secs', search_secs, 'closure_secs', closure_secs, 'nV', [], ...
           'alpha', [], 'alpha_str', '', 'alpha_float', NaN, 'lp_secs', 0);
if isempty(Qbest) || ~getopt(opts, 'exact', true), R.secs = toc(T0); return, end
T2 = tic;
[a, af] = tangent_bound(Qbest, A);
lam = floor(a * sym(10)^6) / sym(10)^6;
R.nV = size(Qbest, 2); R.alpha = lam; R.alpha_exact = a; R.alpha_str = fmt_dec(lam, 6); R.alpha_float = af;
R.lp_secs = toc(T2); R.secs = toc(T0);
end
