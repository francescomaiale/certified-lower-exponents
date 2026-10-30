function T = table3_switching(varargin)
% TABLE3_SWITCHING  Table 3 (tab:switching) and Obs. 4.8 (obs:scaling): random pairs on which switching beats
% chattering (data/switching_pairs.json).
%   Lower bounds (prop:polytope), for h = 0.1 and 0.05: the word w = 1^a 2^b with a, b within 4 of t1/h, t2/h (the
%   best two-segment schedule found) minimizing rho(B_w)^(1/|w|) gives g(h) <= g_up (floating point); the closure is
%   seeded by the Perron vector of B_w and its images along the word, at the levels t = 1 + h (g_up - eps),
%   eps = 1e-5, 1e-3, 1e-2, 1e-1.
%     'mode', 'reproduce'  (default) one closure at the level that closed in the paper's run (none for the rows
%                          marked --- in Table 3);
%     'mode', 'search'     the four levels in turn, each closure capped at 'pmax' processed vertices (default
%                          4000) and 20000 vertices.
%   Upper bounds: the paper's optimized two-segment schedule (rational durations), certified by exponential
%   enclosures in both data readings (certify_best_upper) and replayed by uniformization.
%   Options: 'cdir' (certificate folder, e.g. a scratch folder for a search-mode run), 'pairs', 'hs', 'part' ('lower', 'upper', 'summary', 'all'), 'force', 'replay'.
%   Certificates: certificates/switching/<pair>_<h>_best.json and <pair>_schedule.json; summary results/table3.md.
%
%   Long: the d = 10 rows (500 to 3600 vertices); a few minutes to an hour each.

p = inputParser;
p.addParameter('mode', 'reproduce'); p.addParameter('pairs', {}); p.addParameter('hs', {'0.1', '0.05'});
p.addParameter('part', 'all'); p.addParameter('force', false); p.addParameter('cdir', ''); p.addParameter('replay', true);
p.addParameter('pmax', 4000);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); cdir = o.cdir; if isempty(cdir), cdir = fullfile(root, 'certificates', 'switching'); end
info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_info.json')));
pairs = o.pairs; if isempty(pairs), pairs = {info.id}; end
tag = o.part; if numel(pairs) <= 4, tag = [o.part '_' strjoin(pairs, '+')]; end
say = logger(fullfile(root, 'results', 'logs', sprintf('table3_%s.log', tag)));
EPS = [1e-5, 1e-3, 1e-2, 1e-1];
for q = 1:numel(pairs)
    pid = pairs{q}; I = info(strcmp({info.id}, pid)); A = load_pair('switching', pid);
    if any(strcmp(o.part, {'upper', 'all'}))
        file = fullfile(cdir, [pid '_schedule.json']);
        if ~isfile(file) || o.force
            L = paper_levels('best_upper', pid);
            sch = struct('modes', [1 2], 'durations', {{L.t1, L.t2}}, 'rate_float', L.schedule_float);
            t0 = tic; [R, json] = certify_best_upper(A, struct('id', pid, 'chatter', false, 'schedule', sch));
            write_text(file, json);
            ok = true; if o.replay, ok = replay_best_upper(file); end
            say('%s: schedule (t1 %.4f, t2 %.4f) certified lambda_check <= %s (float %+.7f), %.0fs, replay %s', pid, ...
                R.t1_float, R.t2_float, R.printed, L.schedule_float, toc(t0), tf(ok));
        end
    end
    if ~any(strcmp(o.part, {'lower', 'all'})), continue, end
    for hh = o.hs
        h = hh{1}; row = sprintf('%s_%s', pid, h); file = fullfile(cdir, [row '_best.json']);
        if isfile(file) && ~o.force, continue, end
        D = exact_data(A, h); hf = D.hf; d = D.d;
        B = cellfun(@(M) eye(d) + hf * M, D.Af, 'UniformOutput', false);
        [v, a, b] = best_two_block(B, round(I.t1 / hf), round(I.t2 / hf), 4);
        g_up = (exp(v) - 1) / hf;
        w = [ones(1, a), 2 * ones(1, b)];
        seeds = perron_seeds(B, w);
        Lv = paper_levels('switching', row);
        if strcmp(o.mode, 'reproduce')
            if ~isfield(Lv, 't_build')
                say('%s: word 1^%d 2^%d, g(h) <= %+.6f; no level closed in the paper''s run (---)', row, a, b, g_up);
                continue
            end
            levels = Lv.t_build; epss = Lv.eps; caps = Inf;
        else
            levels = 1 + hf * (g_up - EPS); epss = EPS; caps = o.pmax * ones(1, 4);
        end
        done = false;
        for k = 1:numel(levels)
            spec = struct('name', pid, 'kind', 'best', 'seeds', seeds, 'level', levels(k), 'nmax', 20000, ...
                          'pmax', caps(min(k, end)), 'file', file, 'replay', o.replay, ...
                          'meta', struct('seed', 'perron', 'word', [a b], 'g_up', g_up, 'eps', epss(k)));
            out = run_lower_row(A, h, spec);
            if out.closed
                say('%s: word 1^%d 2^%d, g(h) <= %+.6f; eps %.0e closed: CERTIFIED lambda_check >= g(h) >= %s (exact %+.9f), |V| = %d, closure %.0fs, exact %.0fs, replay %s', ...
                    row, a, b, g_up, epss(k), out.lambda, out.lambda_star, out.nV, out.secs_closure, out.secs_exact, tf(out.replay_ok));
                done = true; break
            end
            say('%s: eps %.0e: closure did not finish (%s, |V| = %d)', row, epss(k), out.info.closure.reason, out.info.closure.nV);
        end
        if ~done, say('%s: word 1^%d 2^%d, g(h) <= %+.6f; no level closed', row, a, b, g_up); end
    end
end
T = [];
if any(strcmp(o.part, {'summary', 'all'})), T = summary(root, cdir, info); end
end


function [v, a, b] = best_two_block(B, a0, b0, r)
% the word 1^a 2^b, a in a0 +- r, b in b0 +- r (a, b >= 1), minimizing log rho(B_2^b B_1^a) / (a + b)
v = Inf;
for aa = max(1, a0 - r):a0 + r
    for bb = max(1, b0 - r):b0 + r
        P = B{2}^bb * B{1}^aa;
        x = log(max(abs(eig(P)))) / (aa + bb);
        if x < v, v = x; a = aa; b = bb; end
    end
end
end


function T = summary(root, cdir, info)
say = logger(fullfile(root, 'results', 'table3.md'));
% paper (tab:switching): c, upper, g<= (0.1), g>= (0.1), g<= (0.05), g>= (0.05), width
P = {'sw5_1', '0.539231', '0.467683', '0.465337', '0.465326', '0.466492', '0.366491', '2.4e-03';
     'sw5_2', '0.459280', '0.457887', '0.442429', '0.442418', '0.450720', '0.440719', '1.6e-02';
     'sw5_3', '-0.333377', '-0.365833', '-0.457212', '-0.457223', '-0.409526', '-0.509527', '9.2e-02';
     'sw5_4', '-0.401497', '-0.430082', '-0.460383', '-0.460394', '-0.443685', '-0.443696', '1.4e-02';
     'sw5_5', '0.587104', '0.548729', '0.531162', '0.531151', '0.540134', '0.539133', '9.6e-03';
     'sw5_6', '0.532754', '0.519922', '0.500647', '0.500636', '0.511408', '0.511397', '8.6e-03';
     'sw5_7', '1.431824', '1.428992', '1.377996', '1.377985', '1.406985', '1.396984', '3.3e-02';
     'sw5_8', '-0.152645', '-0.171842', '-0.210279', '-0.211280', '-0.189358', '-0.189369', '1.8e-02';
     'sw5_9', '0.481345', '0.474726', '0.404665', '0.404654', '0.440660', '0.430659', '4.5e-02';
     'sw5_10', '0.857180', '0.843598', '0.839946', '0.838945', '0.841758', '0.840757', '2.9e-03';
     'sw10_1', '0.567910', '0.559043', '0.532613', '0.432612', '0.546054', '0.446053', '1.2e-01';
     'sw10_2', '0.716928', '0.674041', '0.629877', '0.629866', '0.653399', '0.553398', '4.5e-02';
     'sw10_3', '0.802565', '0.788865', '0.777934', '0.677933', '0.783449', '0.683448', '1.1e-01';
     'sw10_4', '0.253414', '0.245147', '0.216605', '0.116604', '0.231252', '---', '1.3e-01';
     'sw10_5', '-0.168419', '-0.197966', '-0.228868', '-0.228879', '-0.212874', '-0.312875', '3.1e-02';
     'sw10_6', '0.673707', '0.638392', '0.635617', '0.535616', '0.636998', '---', '1.1e-01';
     'sw10_7', '1.200095', '1.185699', '1.136618', '1.036617', '1.161028', '---', '1.5e-01';
     'sw10_8', '0.198883', '0.142660', '0.077000', '-0.023001', '0.110582', '---', '1.7e-01';
     'sw10_9', '0.953083', '0.944065', '0.932144', '---', '0.938291', '---', '---';
     'sw10_10', '0.670254', '0.387498', '0.369426', '0.269425', '0.378543', '---', '1.2e-01'};
say('# Table 3 (tab:switching): switching pairs, MATLAB reproduction');
say('');
say('Each cell: MATLAB value, then the paper''s value in brackets when they differ.');
say('');
say('| d | pair | c | upper | g(0.1) <= | g(0.1) >= | g(0.05) <= | g(0.05) >= | width | match |');
say('|---|---|---|---|---|---|---|---|---|---|');
T = struct('pair', {}, 'match', {});
nres = zeros(4, 4);
for k = 1:size(P, 1)
    pid = P{k, 1}; I = info(strcmp({info.id}, pid)); A = load_pair('switching', pid);
    vals = cell(1, 7);
    vals{1} = fmt_dec(sym(I.chatter, 'f'), 6, 'ceil');
    sf = fullfile(cdir, [pid '_schedule.json']); vals{2} = '(missing)';
    if isfile(sf), R = jsondecode(fileread(sf)); vals{2} = R.printed; end
    lows = [];
    for j = 1:2
        h = {'0.1', '0.05'}; h = h{j}; D = exact_data(A, h);
        B = cellfun(@(M) eye(D.d) + D.hf * M, D.Af, 'UniformOutput', false);
        v = best_two_block(B, round(I.t1 / D.hf), round(I.t2 / D.hf), 4);
        vals{1 + 2 * j} = fmt_dec(sym((exp(v) - 1) / D.hf, 'f'), 6, 'ceil');
        lf = fullfile(cdir, sprintf('%s_%s_best.json', pid, h)); vals{2 + 2 * j} = '---';
        if isfile(lf)
            C = jsondecode(fileread(lf)); vals{2 + 2 * j} = C.lambda_lower; lows(end + 1) = str2double(C.lambda_lower); %#ok<AGROW>
            dd = 1 + (I.d == 10) + 2 * (j == 2);       % rows: d = 5 / 10 at h = 0.1, then at h = 0.05
            nres(dd, :) = nres(dd, :) + (C.meta.eps == [1e-5 1e-3 1e-2 1e-1]);
        end
    end
    if isempty(lows), vals{7} = '---'; else, vals{7} = sprintf('%.1e', ceil_sig(str2double(vals{2}) - max(lows))); end
    cells = vals;
    for j = 1:7, if ~strcmp(vals{j}, P{k, j + 1}), cells{j} = sprintf('%s [%s]', vals{j}, P{k, j + 1}); end, end
    match = all(strcmp(vals, P(k, 2:8)));
    say('| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |', I.d, pid, cells{:}, yesno(match));
    T(end + 1) = struct('pair', pid, 'match', match); %#ok<AGROW>
end
say('');
say('c: chattering value (floating point, rounded up).  upper: certified two-segment schedule, rounded up');
say('(certificates/switching/<pair>_schedule.json).  g(h) <= : floating-point bound of the best word 1^a 2^b, rounded');
say('up.  g(h) >= : certified (certificates/switching/<pair>_<h>_best.json), rounded down; --- no level closed.');
say('Width: upper minus the larger certified lower bound, rounded up to two digits.');
say('Level that closed (certificate within 1e-5 / 1e-3 / 1e-2 / 1e-1 of g(h) <=), numbers of pairs: h = 0.1: d = 5 %s,', ...
    mat2str(nres(1, :)));
say('d = 10 %s; h = 0.05: d = 5 %s, d = 10 %s.  Paper: h = 0.1: d = 5 [8 2 0 0], d = 10 [2 0 0 7];', mat2str(nres(2, :)), ...
    mat2str(nres(3, :)), mat2str(nres(4, :)));
say('h = 0.05: d = 5 [3 2 3 2], d = 10 [0 0 0 4].');
say('Rows reproduced to every printed digit: %d of %d.', nnz([T.match]), numel(T));
end


function y = ceil_sig(x)
e = floor(log10(abs(x))); y = ceil(x / 10^(e - 1) - 1e-9) * 10^(e - 1);
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
