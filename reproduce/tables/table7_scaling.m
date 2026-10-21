function T = table7_scaling(varargin)
% TABLE7_SCALING  Obs. 4.7 (obs:scaling) and Table 8 (tab:scaling): certified brackets for the random pairs in
% dimensions 5 and 10 with apparent chattering behavior (a numerical classification) (data/scaling_pairs.json).
%   T = table7_scaling('mode', 'reproduce')   (default) each lower bound by one closure at the level where the
%                                             paper's 30-minute bisection last closed (data/paper_search_results.json)
%   T = table7_scaling('mode', 'search')      bisection between the trivial level and min(1 + h c, best explicit
%                                             word) with machine-independent caps ('pmax', 'maxprobes')
%   Options: 'cdir' (certificate folder, e.g. a scratch folder for a search-mode run), 'pairs' (cellstr, default all 12), 'hs' ({'0.1', '0.05'}), 'part' ('lower', 'upper', 'summary',
%   'all'), 'force', 'replay'.  Certificates: certificates/scaling/<pair>_<h>_best.json and <pair>_<h>.json.
%   Summary: results/table7.md (Table 8, d = 5) and the d = 10 statements of Obs. 4.7.
%
%   Long: every lower row has 1100 to 2800 vertices; in reproduce mode a d = 5 row takes 5 to 20 minutes and a d = 10
%   row 20 to 60 minutes (closure, exact check, replay), single thread.

p = inputParser;
p.addParameter('mode', 'reproduce'); p.addParameter('pairs', {}); p.addParameter('hs', {'0.1', '0.05'});
p.addParameter('part', 'all'); p.addParameter('force', false); p.addParameter('cdir', ''); p.addParameter('replay', true);
p.addParameter('pmax', 3000); p.addParameter('maxprobes', 12);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); cdir = o.cdir; if isempty(cdir), cdir = fullfile(root, 'certificates', 'scaling'); end
[~, all_pairs] = load_pair('scaling');
pairs = o.pairs; if isempty(pairs), pairs = all_pairs; end
info = jsondecode(fileread(fullfile(root, 'data', 'scaling_pairs_info.json')));
say = logger(fullfile(root, 'results', 'logs', sprintf('table7_%s_%s.log', o.part, strjoin(pairs, '+'))));
if numel(pairs) > 3, say = logger(fullfile(root, 'results', 'logs', sprintf('table7_%s.log', o.part))); end
for q = 1:numel(pairs)
    pid = pairs{q}; A = load_pair('scaling', pid); d = size(A{1}, 1);
    for hh = o.hs
        h = hh{1}; row = sprintf('%s_%s', pid, h);
        if any(strcmp(o.part, {'lower', 'all'}))
            file = fullfile(cdir, [row '_best.json']);
            if ~isfile(file) || o.force
                spec = struct('name', pid, 'kind', 'best', 'file', file, 'replay', o.replay, 'seeds', eye(d));
                if strcmp(o.mode, 'reproduce')
                    spec.level = paper_levels('scaling', row).t_build;
                else
                    hf = str2double(h); Dd = exact_data(A, h);
                    Fb = cellfun(@(M) eye(d) + hf * M, Dd.Af, 'UniformOutput', false);
                    [~, vb] = word_search(Fb, hf, struct('ntop', 2));
                    c = info(strcmp({info.id}, pid)).chatter;
                    spec.search = struct('t_up', min(1 + hf * c, exp(vb)), 'pmax', o.pmax, 'maxprobes', o.maxprobes, ...
                                         'probe', Inf, 'budget', Inf);
                end
                say('%s: lower bound, mode %s ...', row, o.mode);
                out = run_lower_row(A, h, spec);
                if out.closed
                    say('%s: CERTIFIED lambda_check >= %s (exact %+.10f), |V| = %d, closure %.0fs, exact %.0fs, replay %s %.0fs', ...
                        row, out.lambda, out.lambda_star, out.nV, out.secs_closure, out.secs_exact, tf(out.replay_ok), ...
                        out.secs_replay);
                else
                    say('%s: closure did not finish', row);
                end
            end
        end
        if any(strcmp(o.part, {'upper', 'all'}))
            file = fullfile(cdir, [row '.json']);
            if ~isfile(file) || o.force
                out = run_upper_row(A, h, struct('name', pid, 'file', file, 'replay', o.replay, 'source', 'pilot'));
                say('%s: CERTIFIED lambda_check <= %s (k = %d; float %+.9f), %.0fs, replay %s', row, out.upper, out.k, ...
                    out.U_float, out.secs_search + out.secs_exact, tf(out.replay_ok));
            end
        end
    end
end
T = [];
if any(strcmp(o.part, {'summary', 'all'})), T = summary(cdir, root, info); end
end


function T = summary(cdir, root, info)
say = logger(fullfile(root, 'results', 'table7.md'));
% paper values (tab:scaling): lower at h = 0.1, upper at h = 0.1 and 0.05, |V|
P = {'d5_dense1', '-0.222751', '-0.204815', '-0.209005', 1755; 'd5_dense2', '-0.059850', '-0.000609', '-0.014798', 1391;
     'd5_dense3', '0.222107', '0.242747', '0.237887', 1146;   'd5_sparse1', '-0.107049', '0.029982', '-0.001609', 1118;
     'd5_sparse2', '-0.641285', '-0.092940', '-0.097694', 1802; 'd5_sparse3', '-0.071028', '0.170806', '0.112334', 1454};
say('# Table 8 (tab:scaling): random pairs in dimension 5, MATLAB reproduction');
say('');
say('| pair (c) | h | lower (MATLAB) | lower (paper) | g(h) <= (GM run, paper) | upper (MATLAB) | upper (paper) | width | |V| MATLAB / paper | match |');
say('|---|---|---|---|---|---|---|---|---|---|');
T = struct('row', {}, 'match', {});
for k = 1:size(P, 1)
    pid = P{k, 1}; c = info(strcmp({info.id}, pid)).chatter;
    L = read_cert(fullfile(cdir, [pid '_0.1_best.json']), 'lambda_lower');
    nV = read_cert(fullfile(cdir, [pid '_0.1_best.json']), 'nV');
    for h = {'0.1', '0.05'}
        U = read_cert(fullfile(cdir, sprintf('%s_%s.json', pid, h{1})), 'printed');
        Lp = P{k, 2}; Up = P{k, 2 + find(strcmp(h{1}, {'0.1', '0.05'}))};
        dag = ''; if strcmp(h{1}, '0.05'), dag = '^dagger'; end
        gm = paper_levels('scaling', sprintf('%s_%s', pid, h{1})).gm_upper_g;
        match = strcmp(L, Lp) && strcmp(U, Up);
        say('| %s (%s) | %s | %s%s | %s | %s | %s | %s | %.1e | %d / %d | %s |', pid, fmt_dec(sym(c, 'f'), 6, 'ceil'), ...
            h{1}, L, dag, Lp, fmt_dec(sym(gm, 'f'), 6, 'ceil'), U, Up, str2double(U) - str2double(L), nV, P{k, 5}, yesno(match));
        T(end + 1) = struct('row', sprintf('%s_%s', pid, h{1}), 'match', match); %#ok<AGROW>
    end
end
say('');
say('^dagger: the bound certified at h = 0.1, which holds at h = 0.05 by prop:polytope(ii).  The direct h = 0.05');
say('certificates (below) are weaker.  The vertex counts depend on the LP solver and are not part of the claim.');
say('');
say('## All rows (d = 5 and 10), lower and upper bounds');
say('');
say('| row | lower (MATLAB) | lower (paper) | upper (MATLAB) | upper (paper) | |V| | match |');
say('|---|---|---|---|---|---|---|');
for pid = {info.id}
    for h = {'0.1', '0.05'}
        row = sprintf('%s_%s', pid{1}, h{1}); Lv = paper_levels('scaling', row);
        L = read_cert(fullfile(cdir, [row '_best.json']), 'lambda_lower');
        U = read_cert(fullfile(cdir, [row '.json']), 'printed');
        nV = read_cert(fullfile(cdir, [row '_best.json']), 'nV');
        say('| %s | %s | %s | %s | %s | %d | %s |', row, L, Lv.lambda_lower, U, Lv.upper, nV, ...
            yesno(strcmp(L, Lv.lambda_lower) && strcmp(U, Lv.upper)));
    end
end
end


function v = read_cert(f, field)
if ~isfile(f), if strcmp(field, 'nV'), v = NaN; else, v = '(missing)'; end, return, end
C = jsondecode(fileread(f));
if strcmp(field, 'nV'), v = numel(C.vertices); else, v = C.(field); end
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
