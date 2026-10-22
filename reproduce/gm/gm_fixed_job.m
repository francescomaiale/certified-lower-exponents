function R = gm_fixed_job(group, pid, h_str, side, variant, varargin)
% GM_FIXED_JOB  One run of the algorithm of Guglielmi and Maiale (fixed code, github.com/francescomaiale/
% subradius_computation at commit 181d6b5) on an Euler family, for the reruns of paper B, App. B.
%   gm_fixed_job(group, pid, h_str, side, variant)
%     group    'test' (Table 5), 'scaling' (Table 6 / Obs. B.1), 'switching' (Table 7 / Obs. 4.8), 'fixed' (Cor. 4.1)
%     side     'E' (explicit family, lower estimate) or 'I' (implicit family, upper estimate)
%     variant  'adaptive' (adaptive_subradius_comp) or 'fixed' (lsr_comp, the fixed antinorm of Cor. 4.1)
%   Settings as in the paper: transposed families, delta = h^2/10, V = coordinate vectors, at most 2e4 products,
%   normalization by lambda_hat (case 1440: 0.345674608; scaling pairs: c; switching pairs: the certified schedule
%   bound; otherwise none); wall-clock caps 7200 s (test), 3600 s (fixed), 1800 s (scaling, switching).
%   Explicit runs are then certified with our code (paperB_repo): (table) the estimate, rounded to six
%   decimals, by a closure from the coordinate vectors on B_h (test pairs, as Obs. B.1); (gmV) GM's final vertex set
%   by Prop. 4.3 on B_h^T; (gmseed) the closure of GM's vertices at GM's level.  Output in results/<group>/ and
%   certificates/<group>/ in the original runs; packaged outputs are described in README.md.
%   Optional name/value arguments: 'rdir' (base output directory), 'force' (false),
%   'cap' (paper default), 'products' (20000), 'certify' (true).
%   Reduced budgets or certify=false are diagnostics, not paper reproductions.
p = inputParser;
p.addParameter('rdir', ''); p.addParameter('force', false, @(x) islogical(x) && isscalar(x));
p.addParameter('cap', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0));
p.addParameter('products', 20000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && x == fix(x));
p.addParameter('certify', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:}); o = p.Results;
group = validatestring(group, {'test', 'fixed', 'scaling', 'switching'});
side = validatestring(side, {'E', 'I'});
variant = validatestring(variant, {'adaptive', 'fixed'});
h = str2double(h_str);
validateattributes(h, {'double'}, {'scalar', 'real', 'finite', 'positive'});
W = fileparts(mfilename('fullpath')); root = fileparts(W);
oldpath = path; restore = onCleanup(@() path(oldpath)); %#ok<NASGU>
addpath(root); startup_reproduce; addpath(W, '-begin');
base = o.rdir; if isempty(base), base = fullfile(root, 'results', 'gm_runs'); end
rdir = fullfile(base, variant, group); cdir = fullfile(rdir, 'certificates');
if ~isfolder(rdir), mkdir(rdir); end
if ~isfolder(cdir), mkdir(cdir); end
tag = sprintf('%s_%s_%s', pid, h_str, side);
file = fullfile(rdir, [tag '.json']);
if isfile(file) && ~o.force
    R = jsondecode(fileread(file)); fprintf('RESULT %s/%s exists\n', group, tag); return
end

seeds = [];
switch group
    case {'test', 'fixed'}
        As = load_pair('test', pid); lamhat = 0;
        if strcmp(pid, 'case1440'), lamhat = 0.345674608; end
        if strcmp(group, 'test'), tmax = 7200; else, tmax = 3600; end
        d = size(As{1}, 1); seeds = eye(d); if strcmp(pid, 'case148'), seeds = [1 0; 0 0; 0 1]; end
    case 'scaling'
        As = load_pair('scaling', pid);
        info = jsondecode(fileread(fullfile(paperB_root(), 'data', 'scaling_pairs_info.json')));
        lamhat = info(strcmp({info.id}, pid)).chatter; tmax = 1800;
    case 'switching'
        As = load_pair('switching', pid);
        S = jsondecode(fileread(fullfile(paperB_root(), 'certificates', 'switching', [pid '_schedule.json'])));
        lamhat = double(parse_rationals(string(S.U_s))); tmax = 1800;
end
A = cellfun(@str2double, As, 'UniformOutput', false); h = str2double(h_str); d = size(A{1}, 1); I = eye(d);
paper_cap = tmax; if ~isempty(o.cap), tmax = o.cap; end
delta = h^2 / 10; M = o.products;
if side == 'E'
    sc = exp(h * lamhat); F = cellfun(@(X) (I + h * X).', A, 'UniformOutput', false);
else
    sc = 1 / (1 - h * lamhat); F = cellfun(@(X) inv(I - h * X).', A, 'UniformOutput', false);
end
X = [F{:}] / sc;
t0 = tic;
if strcmp(variant, 'fixed')
    [r, mm, ~, status] = lsr_comp_capped(X, [], delta, M, I, [], false, tmax); V = I;
else
    [r, mm, ~, V, ~, status] = adaptive_capped(X, [], delta, M, I, [], false, true, tmax);
end
secs = toc(t0);
lo = r(1) * sc; hi = r(2) * sc;
if side == 'E', bound = (lo - 1) / h; other = (hi - 1) / h; else, bound = (1 - 1 / hi) / h; other = (1 - 1 / lo) / h; end
R = struct('group', group, 'pair', pid, 'h', h_str, 'side', side, 'variant', variant, 'lamhat', lamhat, 'sc', sc, ...
           'lo', lo, 'hi', hi, 'bound', bound, 'other', other, 'depth', mm(2), 'products', mm(3), 'nV', size(V, 2), ...
           'secs', secs, 'status', status, 'delta', delta, 'M', M, 'tmax', tmax, 'code', 'subradius_computation@181d6b5');
fprintf('RESULT %s/%s %s: bound %+.6f (other %+.6f), depth %d, %d products, |V| = %d, %.0f s, %s\n', group, tag, ...
        variant, bound, other, R.depth, R.products, R.nV, secs, status);
R.paper_protocol = tmax == paper_cap && M == 20000 && o.certify;
R.certification_requested = o.certify;
save(fullfile(rdir, [tag '_V.mat']), 'V');

if side == 'E' && o.certify
    if strcmp(group, 'switching')
        L = paper_levels('switching', sprintf('%s_%s', pid, h_str)); R.g_up = L.g_up; R.above_word = bound > L.g_up;
        possible = ~R.above_word;
    else
        possible = true;
    end
    Fu = cellfun(@(Y) (I + h * Y).', A, 'UniformOutput', false);
    % (table) the six-decimal estimate, certified by a closure from the coordinate vectors on B_h (Obs. B.1)
    if strcmp(group, 'test') && possible
        lam = sprintf('%.6f', bound); R.table = struct('claim', lam, 'closed', false, 'nV', NaN);
        margin = 1e-9;
        for attempt = 1:4
            spec = struct('name', pid, 'kind', 'table', 'seeds', seeds, 'level', (1 + h * str2double(lam)) * (1 + margin), ...
                          'lambda_claim', lam, 'file', fullfile(cdir, [tag '_table.json']), 'tmax', 600);
            try
                out = run_lower_row(As, h_str, spec);
                if out.closed, R.table = struct('claim', lam, 'closed', true, 'nV', out.nV, 'replay_ok', out.replay_ok); break, end
            catch err
                fprintf('table certificate, margin %.0e: %s\n', margin, err.message);
            end
            margin = margin * 100;
        end
        fprintf('RESULT %s/%s table certificate of %s: closed %d\n', group, tag, lam, R.table.closed);
    end
    % (gmV) GM's final vertex set, certified on B_h^T (adaptive runs on test and switching pairs)
    if strcmp(variant, 'adaptive') && any(strcmp(group, {'test', 'switching'}))
        out = run_lower_row(As, h_str, struct('name', pid, 'kind', 'gm_vertices', 'orientation', 'T', ...
                            'vertices', sym(V, 'f'), 'file', fullfile(cdir, [tag '_gmV.json']), ...
                            'meta', struct('source', 'GM (fixed code) final vertex set')));
        R.cert_V = struct('lambda', out.lambda, 'nV', out.nV, 'replay_ok', out.replay_ok);
        fprintf('RESULT %s/%s vertex set certifies %s\n', group, tag, out.lambda);
    end
    % (gmseed) the closure of GM's vertices at GM's level, certified on B_h^T
    R.seed = struct('closed', false, 'lambda', '', 'replay_ok', false, 'nV', NaN, 'reason', 'not tried');
    if possible
        stmax = 600; if strcmp(group, 'scaling'), stmax = 300; end
        [Vs, st] = polytope_closure(Fu, lo, V ./ sum(V, 1), stmax, 20000);
        R.seed = struct('closed', ~isempty(Vs), 'lambda', '', 'replay_ok', false, 'nV', st.nV, 'reason', st.reason, ...
                        'secs', st.secs);
        if ~isempty(Vs)
            out = run_lower_row(As, h_str, struct('name', pid, 'kind', 'gm_seeded', 'orientation', 'T', ...
                                'vertices', sym(Vs, 'f'), 'file', fullfile(cdir, [tag '_gmseed.json']), ...
                                'meta', struct('seed', 'GM vertices (fixed code)', 't_build', lo, 'closure_secs', st.secs)));
            R.seed.lambda = out.lambda; R.seed.replay_ok = out.replay_ok; R.seed.nV = out.nV;
        end
        fprintf('RESULT %s/%s seeded closure: %s, certified %s\n', group, tag, st.reason, R.seed.lambda);
    end
end
fid = fopen(file, 'w'); fwrite(fid, jsonencode(R)); fclose(fid);
fprintf('RESULT %s/%s done\n', group, tag);
end
