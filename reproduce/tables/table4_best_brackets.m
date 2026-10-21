function T = table4_best_brackets(varargin)
% TABLE4_BEST_BRACKETS  Table 4 (tab:glp), our columns: best certified lower bound, best implicit Euler bound,
% certified upper bound for c or two-segment schedule, and width; Obs. 4.11 (obs:glp), upper side.
%   T = table4_best_brackets('part', 'certify')   certify c (Collatz-Wielandt at a rational alpha) for all 16 pairs,
%                                                 and the two-segment schedules of cases 1440 and 148 (exponential
%                                                 enclosures), both data readings; replay each certificate
%   T = table4_best_brackets('part', 'summary')   Table 4 from these certificates and those of table2_test_pairs and
%                                                 table7_scaling (lower: best of h = 0.1, 0.05; implicit Euler: best of
%                                                 h = 0.1, 0.05), with the Alg. (L) columns of table4_algorithm_l if
%                                                 present -> results/table4.md
%   Options: 'pairs', 'force', 'replay', 'schedule' ('paper': the rational durations of the paper's optimized
%   schedules (data/paper_search_results.json); 'optimize': re-optimize in floating point (Nelder-Mead) first).
%   Certificates: certificates/table4/<pair>_best_upper.json.

p = inputParser;
p.addParameter('part', 'all'); p.addParameter('pairs', {}); p.addParameter('force', false);
p.addParameter('replay', true); p.addParameter('schedule', 'paper');
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); cdir = fullfile(root, 'certificates', 'table4');
test = {'planar', 'case1', 'case1440', 'case148'};
[~, pilot] = load_pair('scaling');
pairs = o.pairs; if isempty(pairs), pairs = [test, pilot]; end
say = logger(fullfile(root, 'results', 'logs', sprintf('table4_%s.log', o.part)));
if any(strcmp(o.part, {'certify', 'all'}))
    for q = 1:numel(pairs)
        pid = pairs{q}; file = fullfile(cdir, [pid '_best_upper.json']);
        if isfile(file) && ~o.force, say('%s: certificate exists', pid); continue, end
        if ismember(pid, test), A = load_pair('test', pid); else, A = load_pair('scaling', pid); end
        L = paper_levels('best_upper', pid);
        copts = struct('id', pid, 'alpha', L.alpha_float);
        if isfield(L, 't1')
            if strcmp(o.schedule, 'paper')
                copts.schedule = struct('modes', [1 2], 'durations', {{L.t1, L.t2}}, 'rate_float', L.schedule_float);
            else
                Af = cellfun(@str2double, A, 'UniformOutput', false);
                [sf, t1, t2] = best_schedule(Af);
                copts.schedule = struct('modes', [1 2], 'durations', [t1 t2], 'rate_float', sf);
            end
        end
        t0 = tic; [R, json] = certify_best_upper(A, copts); write_text(file, json);
        msg = sprintf('%s: certified c <= %s', pid, fmt_dec(parse_rationals(string(R.U_c)), 9, 'ceil'));
        if ~isempty(R.U_s), msg = [msg sprintf(', schedule <= %s', fmt_dec(parse_rationals(string(R.U_s)), 9, 'ceil'))]; end %#ok<AGROW>
        msg = [msg sprintf('; best %s %s (%.0fs)', R.best_kind, R.printed, toc(t0))]; %#ok<AGROW>
        if o.replay
            t1 = tic; [ok, msgs] = replay_best_upper(file);
            msg = [msg sprintf('; replay %s (%.0fs)', tf(ok), toc(t1))]; %#ok<AGROW>
            if ~ok, msg = [msg ' ' strjoin(msgs, '; ')]; end %#ok<AGROW>
        end
        say('%s', msg);
    end
end
T = [];
if any(strcmp(o.part, {'summary', 'all'})), T = summary(root, test, pilot); end
end


function T = summary(root, test, pilot)
say = logger(fullfile(root, 'results', 'table4.md'));
% paper values (tab:glp): Alg. (L) lower, ours lower, Alg. (L) upper, implicit Euler, c or schedule, width
P = struct( ...
 'planar',      {{'-0.210758', '-0.232706', '-0.194213', '-0.158706', '-0.194976', '1.6e-02'}}, ...
 'case1',       {{'-1.028582', '-1.027040', '-1.011626', '-0.999505', '-1.012395', '1.5e-02'}}, ...
 'case1440',    {{'0.334850', '0.345146', '0.346401', '0.346182', '0.345675', '5.3e-04'}}, ...
 'case148',     {{'---', '-1.247834', '---', '-1.238961', '-1.241336', '6.5e-03'}}, ...
 'd5_dense1',   {{'-0.276163', '-0.222751', '-0.213075', '-0.209005', '-0.213383', '9.4e-03'}}, ...
 'd5_dense2',   {{'-0.081527', '-0.059850', '-0.028459', '-0.014798', '-0.029331', '3.1e-02'}}, ...
 'd5_dense3',   {{'0.210187', '0.222107', '0.233137', '0.237887', '0.232972', '1.1e-02'}}, ...
 'd5_sparse1',  {{'-0.204184', '-0.107049', '-0.032087', '-0.001609', '-0.034259', '7.3e-02'}}, ...
 'd5_sparse2',  {{'-0.384238', '-0.641285', '-0.102190', '-0.097694', '-0.102660', '2.9e-01'}}, ...
 'd5_sparse3',  {{'-0.049779', '-0.071028', '0.055587', '0.112334', '0.053208', '1.1e-01'}}, ...
 'd10_dense1',  {{'0.870186', '0.459247', '0.990484', '0.991569', '0.990410', '1.3e-01'}}, ...
 'd10_dense2',  {{'0.779913', '0.525978', '0.857197', '0.865203', '0.856999', '7.8e-02'}}, ...
 'd10_dense3',  {{'0.009496', '-0.229448', '0.066898', '0.068764', '0.066852', '5.8e-02'}}, ...
 'd10_sparse1', {{'-0.064519', '0.269192', '0.318964', '0.334496', '0.317685', '4.9e-02'}}, ...
 'd10_sparse2', {{'0.364290', '1.356787', '1.570804', '1.583464', '1.570368', '2.2e-01'}}, ...
 'd10_sparse3', {{'-0.451797', '-0.349697', '-0.179980', '-0.143149', '-0.181089', '1.7e-01'}});
say('# Table 4 (tab:glp): best certified brackets and Algorithm (L), MATLAB reproduction');
say('');
say('Each cell: MATLAB value, then the paper''s value in brackets when they differ.');
say('');
say('| pair | Alg. (L) lower | ours lower | Alg. (L) upper (float) | implicit Euler | c or schedule | width | all ours match |');
say('|---|---|---|---|---|---|---|---|');
T = struct('pair', {}, 'match_ours', {}, 'match_glp', {});
for pid = [test, pilot]
    q = pid{1}; pv = P.(q);
    if ismember(q, test), dir = 'table2'; else, dir = 'scaling'; end
    lo = best_of(root, dir, q, '_best', 'lambda_lower', @max);
    if strcmp(q, 'case148'), lo = best_of(root, dir, q, '_transfer', 'lambda_lower', @max, lo); end
    gmark = '';
    if ~ismember(q, test)          % closures seeded by the final vertex sets of the runs of [GM] (App. B)
        g = best_of(root, 'appB_scaling', q, '_E_gmseed', 'lambda_lower', @max);
        if ~strcmp(g, '(missing)') && str2double(g) > str2double(lo), lo = g; gmark = '^g'; end
    end
    up = best_of(root, dir, q, '', 'printed', @min);
    bf = fullfile(root, 'certificates', 'table4', [q '_best_upper.json']);
    cs = '(missing)'; kind = '';
    if isfile(bf), R = jsondecode(fileread(bf)); cs = R.printed; if strcmp(R.best_kind, 'schedule'), kind = '^s'; end, end
    [gl, gu] = glp_best(root, q);
    lmax = max(str2double(lo), str2double(gl));
    wstr = sprintf('%.1e', ceil_sig(str2double(cs) - lmax));
    vals = {gl, lo, gu, up, cs, wstr};
    cells = cell(1, 6);
    for k = 1:6
        if strcmp(vals{k}, pv{k}), cells{k} = vals{k}; else, cells{k} = sprintf('%s [%s]', vals{k}, pv{k}); end
    end
    mo = all(strcmp(vals([2 4 5]), pv([2 4 5])));
    mg = all(strcmp(vals([1 3]), pv([1 3])));
    say('| %s | %s | %s%s | %s | %s | %s%s | %s | %s |', q, cells{1}, cells{2}, gmark, cells{3}, cells{4}, cells{5}, kind, ...
        cells{6}, yesno(mo));
    T(end + 1) = struct('pair', q, 'match_ours', mo, 'match_glp', mg); %#ok<AGROW>
end
say('');
say('Ours: the best of the steps 0.1 and 0.05 (certificates/table2, certificates/scaling); case 148 at h <= 0.1 by');
say('prop:polytope(ii); ^g: a closure seeded by the final vertex set of a run of [GM] (certificates/appB_scaling, App. B).');
say('c or schedule: certificates/table4/<pair>_best_upper.json (^s: two-segment schedule).  Width =');
say('that column minus the larger lower bound, rounded up to two digits.  Alg. (L): results/glp/ (table4_algorithm_l),');
say('best of tau = 0.1 and 0.05; its upper bounds are floating point, rounded up.');
say('Our certified columns reproduced for %d of %d pairs; Alg. (L) columns for %d of %d.', nnz([T.match_ours]), ...
    numel(T), nnz([T.match_glp]), numel(T));
end


function [gl, gu] = glp_best(root, q)
gl = '---'; gu = '---'; lv = -Inf; uv = Inf;
for tau = {'0.1', '0.05'}
    f = fullfile(root, 'results', 'glp', sprintf('glp_%s_%s.json', q, tau{1}));
    if ~isfile(f), continue, end
    R = jsondecode(fileread(f));
    if ~isempty(R.alpha_str) && str2double(R.alpha_str) > lv, lv = str2double(R.alpha_str); gl = R.alpha_str; end
    if R.beta_l < uv, uv = R.beta_l; end
end
if isfinite(uv), gu = fmt_dec(sym(uv, 'f'), 6, 'ceil'); end
end


function v = best_of(root, dir, pid, suffix, field, pick, v0)
v = ''; if nargin > 6, v = v0; end
for h = {'0.1', '0.05', '0.025'}
    f = fullfile(root, 'certificates', dir, sprintf('%s_%s%s.json', pid, h{1}, suffix));
    if ~isfile(f) || (strcmp(h{1}, '0.025') && ~strcmp(suffix, '_transfer')), continue, end
    C = jsondecode(fileread(f)); x = C.(field);
    if strcmp(suffix, '_transfer'), x = C.meta.source_lambda; end
    if isempty(v) || pick(str2double(x), str2double(v)) == str2double(x), v = x; end
end
if isempty(v), v = '(missing)'; end
end


function y = ceil_sig(x)
% round up to two significant digits
e = floor(log10(abs(x))); y = ceil(x / 10^(e - 1) - 1e-9) * 10^(e - 1);
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
