function out = run_lower_row(A, h, spec)
% RUN_LOWER_ROW  One certified lower bound lambda_check >= (t - 1)/h (prop:polytope): closure, exact check,
% certificate file and independent replay.
%   out = run_lower_row(A, h, spec)
%     A     data {A_1, ..., A_m} (decimal strings, e.g. from load_pair); h step (decimal string)
%     spec  fields (all optional except name)
%       name         pair name stored in the certificate (the key of the pair in data/*.json)
%       kind         'best' (default), 'table' or 'transfer'
%       orientation  'N' (default) or 'T'
%       seeds        d x s seeds of the closure (default eye(d))
%       level        reproduce: one closure at this level t (e.g. a level recorded by the paper's search)
%       search       struct of polytope_search options (t_up, t_start, probe, budget, pmax, maxprobes, tol):
%                    bisection instead of a single closure
%       vertices     sym vertex matrix: no closure, certify these vertices (prop:polytope(ii), kind 'transfer')
%       lambda_claim six-decimal string: certify exactly this bound (kind 'table')
%       nmax, pmax, tmax  caps of a single closure (6000, Inf, Inf)
%       file         certificate path (written if set); replay (true)
%       meta         struct stored in the certificate
%   out: name, h, kind, lambda (char), lambda_star (double, exact level), nV, closed (false if the closure did not
%   finish: then no certificate), replay_ok, secs (closure, exact, replay), file, info.

T0 = tic;
D = exact_data(A, h);
check_euler_step(D, 'explicit');
d = D.d; m = D.m;
orient = getopt(spec, 'orientation', 'N');
F = cell(1, m);
for i = 1:m
    F{i} = eye(d) + D.hf * D.Af{i};
    if orient == 'T', F{i} = F{i}.'; end
end
out = struct('name', spec.name, 'h', char(D.h_str), 'kind', getopt(spec, 'kind', 'best'), 'lambda', '', ...
             'lambda_star', NaN, 'nV', NaN, 'closed', false, 'replay_ok', false, 'file', '', ...
             'secs_closure', 0, 'secs_exact', 0, 'secs_replay', 0, 'info', struct());
meta = getopt(spec, 'meta', struct());
if isfield(spec, 'vertices') && ~isempty(spec.vertices)
    V = spec.vertices;
elseif isfield(spec, 'search') && ~isempty(spec.search)
    so = spec.search; so.seeds = getopt(spec, 'seeds', eye(d)); so.verbose = getopt(so, 'verbose', false);
    S = polytope_search(F, D.hf, so.t_up, so);
    V = S.V;
    meta.t_build = S.t_build; meta.t_fail = S.t_fail; meta.probes = numel(S.probes); meta.search_secs = S.secs;
    out.info.search = S;
else
    t = spec.level;
    [V, st] = polytope_closure(F, t, getopt(spec, 'seeds', eye(d)), getopt(spec, 'tmax', Inf), ...
                               getopt(spec, 'nmax', 6000), getopt(spec, 'pmax', Inf));
    meta.t_build = t; meta.closure_secs = round(st.secs, 1); meta.processed = st.processed;
    out.info.closure = st;
    if isempty(V)
        out.secs_closure = toc(T0);
        return
    end
end
out.closed = true; out.secs_closure = toc(T0);
T1 = tic;
copts = struct('name', spec.name, 'kind', out.kind, 'orientation', orient, 'meta', meta, ...
               'lambda_claim', getopt(spec, 'lambda_claim', ''));
[C, json, info] = certify_lower_exact(V, D, copts);
out.secs_exact = toc(T1);
out.lambda = C.lambda_lower; out.lambda_star = double(info.lambda_star); out.nV = info.nV; out.info.exact = info;
file = getopt(spec, 'file', '');
if ~isempty(file), write_text(file, json); out.file = file; end
if getopt(spec, 'replay', true)
    T2 = tic;
    [ok, msgs] = replay_lower(C);
    out.replay_ok = ok; out.info.replay_msgs = msgs; out.secs_replay = toc(T2);
end
end
