function T = table3_transposed(varargin)
% TABLE3_TRANSPOSED  The weak lower bounds of Table 3 (Obs. 4.8) again, with our closure on the transposed family
% B_h^T (prop:polytope(iii)) instead of B_h (results/table3_transposed.md).
%   Rows: those whose Table 3 bound on B_h closed only at eps = 0.1 or not at all (---): at h = 0.1 the d = 10
%   pairs but sw10_2 and sw10_5; at h = 0.05 all ten d = 10 pairs, and the d = 5 pairs closed at eps >= 1e-2 (sw5_1,
%   sw5_2, sw5_3, sw5_7, sw5_9).  Table 3 of the paper keeps, for each row, the better of the two orientations.
%   The word 1^a 2^b and g_up of the paper's run (data/paper_search_results.json); seeds = Perron
%   vectors of the reversed word for B_h^T; levels t = 1 + h (g_up - eps), eps = 1e-5, 1e-3, 1e-2, 0.1, then 0.15,
%   0.2, each closure capped at 'seed_tmax' seconds and 20000 vertices; the first closed polytope is certified
%   exactly (certify_lower_exact, both data readings) and replayed.
%   Options: 'part' ('run' (default), 'replay' or 'summary'), 'rows' (e.g. {'sw10_1_0.1'}; default all),
%   'seed_tmax' (300), 'force'.
%   'replay': every certificate of certificates/table3_transposed through replay_all and through
%   independent_replay/python/run_all.py (python3 -I -B); the counts go to results/table3_transposed/replays.json
%   and into the summary.
%   Outputs: results/table3_transposed/<pair>_<h>_T.json (the ladder), certificates/table3_transposed/<pair>_<h>_T.json,
%   log results/logs/table3_transposed.log.
%
%   Long: up to 6 x 300 s of closures plus an exact check of up to 30 min per row; the rows took about 2 hours with 7
%   single-threaded MATLAB processes in parallel (one row each).

p = inputParser;
p.addParameter('part', 'run'); p.addParameter('rows', {}); p.addParameter('seed_tmax', 300); p.addParameter('force', false);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); rdir = fullfile(root, 'results', 'table3_transposed');
cdir = fullfile(root, 'certificates', 'table3_transposed');
if ~isfolder(rdir), mkdir(rdir); end
T = [];
switch o.part
    case 'run'
        rows = o.rows; if isempty(rows), rows = all_rows(); end
        for r = rows(:).', run_row(o, rdir, cdir, r{1}); end
    case 'replay'
        T = replays(root, rdir, cdir);
    case 'summary'
        T = summary(root, rdir);
end
end


function rows = all_rows()
rows = [strcat({'sw10_1', 'sw10_3', 'sw10_4', 'sw10_6', 'sw10_7', 'sw10_8', 'sw10_9', 'sw10_10'}, '_0.1'), ...
        strcat(arrayfun(@(k) sprintf('sw10_%d', k), 1:10, 'UniformOutput', false), '_0.05'), ...
        strcat({'sw5_1', 'sw5_2', 'sw5_3', 'sw5_7', 'sw5_9'}, '_0.05')];
end


function [pid, h] = split_row(row)
k = find(row == '_', 1, 'last'); pid = row(1:k - 1); h = row(k + 1:end);
end


function run_row(o, rdir, cdir, row)
say = logger(fullfile(paperB_root(), 'results', 'logs', 'table3_transposed.log'));
[pid, h] = split_row(row);
tag = sprintf('%s_T', row); file = fullfile(rdir, [tag '.json']);
if isfile(file) && ~o.force, return, end
T0 = tic;
As = load_pair('switching', pid); D = exact_data(As, h); d = D.d;
L = paper_levels('switching', row);
BT = cellfun(@(M) (eye(d) + D.hf * M).', D.Af, 'UniformOutput', false);
seeds = perron_seeds(BT, fliplr([ones(1, L.word(1)), 2 * ones(1, L.word(2))]));
R = struct('pair', pid, 'h_str', h, 'd', d, 'orientation', 'T', 'word', L.word(:).', 'g_up', L.g_up, 'tries', [], ...
           'eps', NaN, 'lambda', '', 'lambda_star', NaN, 'nV', NaN, 'replay_ok', false, 'secs_closures', 0, ...
           'secs_exact', 0, 'secs_replay', 0, 'secs_total', 0, 'seed_tmax', o.seed_tmax);
tries = struct('eps', {}, 'closed', {}, 'reason', {}, 'nV', {}, 'processed', {}, 'secs', {});
for e = [1e-5 1e-3 1e-2 0.1 0.15 0.2]
    out = run_lower_row(As, h, struct('name', pid, 'kind', 'best', 'orientation', 'T', 'seeds', seeds, ...
                        'level', 1 + D.hf * (L.g_up - e), 'nmax', 20000, 'tmax', o.seed_tmax, ...
                        'file', fullfile(cdir, [tag '.json']), 'meta', struct('seed', 'perron, reversed word, B_h^T', ...
                        'word', L.word(:).', 'g_up', L.g_up, 'eps', e)));
    st = out.info.closure;
    tries(end + 1) = struct('eps', e, 'closed', out.closed, 'reason', st.reason, 'nV', st.nV, ...
                            'processed', st.processed, 'secs', st.secs); %#ok<AGROW>
    R.secs_closures = R.secs_closures + st.secs;
    say('%s: eps %g: %s, |V| = %d, %.0fs%s', tag, e, st.reason, st.nV, st.secs, ifelse(out.closed, ...
        sprintf('; CERTIFIED lambda_check >= %s (exact %.0fs), replay %s', out.lambda, out.secs_exact, ...
        tf(out.replay_ok)), ''));
    if out.closed
        R.eps = e; R.lambda = out.lambda; R.lambda_star = out.lambda_star; R.nV = out.nV; R.replay_ok = out.replay_ok;
        R.secs_exact = out.secs_exact; R.secs_replay = out.secs_replay;
        break
    end
end
if isnan(R.eps), say('%s: no level closed', tag); end
R.tries = tries; R.secs_total = toc(T0); write_text(file, jsonencode(R));
end


function R = replays(root, rdir, cdir)
t0 = tic; [okm, Tm] = replay_all({'table3_transposed'}); sm = toc(t0);
py = fullfile(root, 'independent_replay', 'python', 'run_all.py');
t0 = tic;
[st, txt] = system(sprintf('OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -I -B "%s" "%s"', py, cdir));
sp = toc(t0);
fprintf('%s', txt);
last = regexp(txt, '==== Python replay of [^\n]*', 'match', 'once');
R = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'n_files', numel(dir(fullfile(cdir, '*.json'))), ...
           'matlab_n', numel(Tm), 'matlab_valid', nnz([Tm.ok]), 'matlab_ok', okm, 'matlab_secs', round(sm), ...
           'python_ok', st == 0, 'python_summary', last, 'python_secs', round(sp));
write_text(fullfile(rdir, 'replays.json'), jsonencode(R));
end


function T = summary(root, rdir)
say = logger(fullfile(root, 'results', 'table3_transposed.md'));
say('# Table 3 lower bounds on the transposed family B_h^T (prop:polytope(iii))');
say('');
say('The rows of Table 3 (Obs. 4.8) whose lower bound closed only at eps = 0.1 or not at all (---), rerun with our');
say('closure on B_h^T instead of B_h (tables/table3_transposed.m).  Word 1^a 2^b and g_up as in the paper''s run; seeds:');
say('Perron vectors of the reversed word for B_h^T; levels t = 1 + h (g_up - eps), eps = 1e-5, 1e-3, 1e-2, 0.1, then');
say('0.15, 0.2, each closure capped at 300 s and 20000 vertices (MATLAB, one thread); the first closed polytope is');
say('certified exactly in both data readings and replayed (MATLAB replay_all and independent_replay/python/run_all.py).');
f = fullfile(rdir, 'replays.json');
if isfile(f)
    P = jsondecode(fileread(f));
    say('Replays (part ''replay'', %s), certificates/table3_transposed (%d files): MATLAB replay_all %d of %d VALID', ...
        P.date, P.n_files, P.matlab_valid, P.matlab_n);
    say('(%d s); independent_replay/python/run_all.py: %s.', P.matlab_secs, ...
        regexprep(P.python_summary, '^==== Python replay of ', ''));
end
say('');
say('| pair | h | Table 3, B_h | eps | |V| | time (s) | transposed, B_h^T | eps | |V| | time (s) | gain | better (orientation, eps) | B_h^T closures tried (eps: reason, s) |');
say('|---|---|---|---|---|---|---|---|---|---|---|---|---|');
info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_info.json')));
rows = all_rows();
[~, ord] = sort(cellfun(@(r) row_key(r, info), rows)); rows = rows(ord);
T = struct('row', {}, 'pair', {}, 'h', {}, 'N', {}, 'epsN', {}, 'Tval', {}, 'epsT', {});
for r = rows
    [pid, h] = split_row(r{1});
    f = fullfile(rdir, [r{1} '_T.json']); src = '';
    if ~isfile(f), continue, end
    C = jsondecode(fileread(f));
    [lo, eN, nN] = table3_value(root, r{1});
    LP = paper_levels('switching', r{1}); trN = LP.tries; if iscell(trN), trN = [trN{:}]; end
    lam = C.lambda; if isempty(lam), lam = '---'; end
    closed = ~isempty(C.eps) && ~isnan(C.eps);         % NaN is stored as null
    if closed, eT = sprintf('%g', C.eps); else, eT = '---'; end
    tr = C.tries; if iscell(tr), tr = [tr{:}]; end
    if isfield(tr, 'reason')
        trs = strjoin(arrayfun(@(t) sprintf('%g: %s, %.0f', t.eps, t.reason, t.secs), tr(:).', 'UniformOutput', false), '; ');
    else
        trs = strjoin(arrayfun(@(t) sprintf('%g: %s, %.0f', t.eps, ifelse(t.closed, 'closed', 'time'), t.secs), ...
                      tr(:).', 'UniformOutput', false), '; ');
    end
    tT = sprintf('%.0f', sum([tr.secs]));
    if isfield(C, 'secs_exact') && closed, tT = sprintf('%.0f + %.0f', sum([tr.secs]), C.secs_exact + C.secs_replay); end
    nv = '---'; if closed, nv = sprintf('%d', C.nV); end
    say('| %s | %s | %s | %s | %s | %.0f | %s%s | %s | %s | %s | %s | %s | %s |', pid, h, lo, eN, nN, sum([trN.secs]), lam, ...
        src, eT, nv, tT, gain(lam, lo), better(lo, eN, lam, eT), trs);
    T(end + 1) = struct('row', r{1}, 'pair', pid, 'h', h, 'N', lo, 'epsN', eN, 'Tval', lam, 'epsT', eT); %#ok<AGROW>
end
say('');
say('Table 3, B_h: certificates/switching/<pair>_<h>_best.json, eps its level, |V| its vertices; time: the paper''s');
say('Python search, all levels tried (data/paper_search_results.json; 600 s per closure).  Transposed, B_h^T: our');
say('certificate, |V| its vertices; time: all closures of the ladder, failed ones included (MATLAB, 300 s cap each; a');
say('cap is checked between batches, so a failed closure can run past 300 s), + exact check and replay.  gain:');
say('transposed minus Table 3 (new where Table 3 has ---; negative where B_h^T closed only at a larger eps).  better:');
say('the larger certified bound, with its orientation (N = B_h, T = B_h^T) and eps.  Reasons: ''time'' = 300 s cap,');
say('''size'' = 20000 vertices.');
% ---------------------------------------------------------------------------------------- widths
say('');
say('## Bracket widths (upper schedule bound minus the best certified lower bound, both h)');
say('');
say('| pair | upper | best lower, Table 3 | width, Table 3 | best lower, with B_h^T | width, with B_h^T |');
say('|---|---|---|---|---|---|');
ids = {info.id}; ids = ids(ismember(ids, {T.pair}));
for pid = ids
    S = jsondecode(fileread(fullfile(root, 'certificates', 'switching', [pid{1} '_schedule.json'])));
    up = str2double(S.printed); lowN = -Inf; lowB = -Inf;
    for h = {'0.1', '0.05'}
        lo = table3_value(root, sprintf('%s_%s', pid{1}, h{1}));
        lowN = max(lowN, str2double(lo)); lowB = max(lowB, str2double(lo));
        k = find(strcmp({T.row}, sprintf('%s_%s', pid{1}, h{1})));
        if ~isempty(k), lowB = max(lowB, str2double(T(k).Tval)); end
    end
    say('| %s | %s | %s | %s | %s | %s |', pid{1}, S.printed, dec6(lowN), wid(up, lowN), dec6(lowB), wid(up, lowB));
end
say('');
say('Per pair, over h = 0.1 and 0.05; rounded up to two digits, as in Table 3 (--- where no lower bound is certified).');
% ---------------------------------------------------------------------------------------- counts
say('');
say('## Level that closed (eps = 1e-5 / 1e-3 / 1e-2 / 0.1 / 0.15 / 0.2 / none), rows above');
say('');
E = [1e-5 1e-3 1e-2 0.1 0.15 0.2];
for h = {'0.1', '0.05'}
    for dd = {'sw5_', 'sw10_'}
        sel = strcmp({T.h}, h{1}) & startsWith({T.pair}, dd{1});
        if ~any(sel), continue, end
        cN = counts({T(sel).epsN}, E); cT = counts({T(sel).epsT}, E);
        say('- h = %s, d = %s (%d rows): B_h (Table 3) %s, B_h^T %s', h{1}, ifelse(strcmp(dd{1}, 'sw5_'), '5', '10'), ...
            nnz(sel), mat2str(cN), mat2str(cT));
    end
end
end


function k = row_key(row, info)
[pid, h] = split_row(row);
d = info(strcmp({info.id}, pid)).d; n = sscanf(pid(find(pid == '_') + 1:end), '%d');
k = 1e4 * strcmp(h, '0.05') + 100 * d + n;
end


function [lo, e, nV] = table3_value(root, row)
lo = '---'; e = '---'; nV = '---';
f = fullfile(root, 'certificates', 'switching', [row '_best.json']);
if isfile(f)
    C = jsondecode(fileread(f)); lo = C.lambda_lower; e = sprintf('%g', C.meta.eps); nV = sprintf('%d', C.meta.nV);
end
end


function s = better(lo, eN, lam, eT)
x = str2double(lo); y = str2double(lam);
if isnan(x) && isnan(y), s = '---';
elseif isnan(y) || x >= y, s = sprintf('%s (N, %s)', lo, eN);
else, s = sprintf('%s (T, %s)', lam, eT);
end
end


function c = counts(es, E)
c = zeros(1, numel(E) + 1);
for k = 1:numel(es)
    j = find(abs(str2double(es{k}) - E) < 1e-12 * max(E));
    if isempty(j), c(end) = c(end) + 1; else, c(j) = c(j) + 1; end
end
end


function s = gain(lam, lo)
if strcmp(lam, '---'), s = '---';
elseif strcmp(lo, '---'), s = 'new';
else, s = sprintf('%+.6f', str2double(lam) - str2double(lo));
end
end


function s = dec6(x)
if isinf(x), s = '---'; else, s = sprintf('%.6f', x); end
end


function s = wid(up, lo)
if isinf(lo), s = '---'; return, end
x = up - lo; e = floor(log10(abs(x))); y = ceil(x / 10^(e - 1) - 1e-9) * 10^(e - 1);
s = sprintf('%.1e', y);
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
