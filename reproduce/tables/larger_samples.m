function T = larger_samples(varargin)
% LARGER_SAMPLES  The larger samples of Obs. 4.8 (Tables S5 and S6): the random streams of the paper continued
% (generators/gen_pairs_ext.py; data/switching_pairs_journal.json and scaling_pairs_journal.json, with *_info.json).
%   'part'
%     'switching'   per pair (default all 45): 'upper', the optimized two-segment schedule (rational durations),
%                   certified in both data readings and replayed; '0.1', '0.05', the lower bound for g(h): the best
%                   word 1^a 2^b with a, b within 4 of t1/h, t2/h gives g_up; Perron seeds; levels
%                   t = 1 + h (g_up - eps), eps = 1e-5, 1e-3, 1e-2, 0.1 with 600, 600, 300, 300 s per closure on B_h;
%                   if no level below 0.1 closed, eps up to 0.2 with 300 s per closure on B_h^T (as Table 3)
%     'chattering'  per pair (default all 33): 'upper', the chattering value certified by Collatz-Wielandt (and a
%                   two-segment schedule if one beats it); '0.1', '0.05', certified_bracket (probe 300 s, budget
%                   1800 s, 20000 vertices: the protocol of Table S4)
%     'extended'    the extended search (helpers/extended_search_row.m, without vertex sets of [GM]) on the d = 10
%                   switching rows of the 'switching' records that stopped at eps >= 1e-2 or found none, both
%                   orientations; at h = 0.05 the best certified h = 0.1 polytope of the orientation (switching or
%                   extended) is a seed and a baseline.  Run it after 'switching'; the h = 0.1 rows go first.
%     'replay'      certificates/larger_switching and larger_scaling: MATLAB replay_all (about 10 hours on one
%                   thread) and the Python replay (15 minutes)
%     'summary'     (default) Tables S5 and S6 and the counts of Obs. 4.8 from the stored certificates and run
%                   records: results/larger_samples.md and results/larger_samples_tables.tex
%   Options: 'pairs' (cell of names), 'steps' ({'upper', '0.1', '0.05'}), 'rdir' (records and certificates of a run;
%   default results/larger_samples/run, so that a run never overwrites the stored ones), 'cap' (extended search,
%   seconds per closure, 1800), 'force'.
%   Every closed polytope is certified exactly in both data readings and replayed (run_lower_row).
%   Long: the paper's runs took about 14 hours with 12 to 16 single-threaded MATLAB processes, one pair and step each
%   (d = 20 rows up to 1.5 hours, the extended search up to 3.5 hours per row); split the work with 'pairs'/'steps'.
%   The stored records are in results/larger_samples/{switching,scaling,extended}; MANIFEST.txt lists the files.

p = inputParser;
p.addParameter('part', 'summary'); p.addParameter('pairs', {}); p.addParameter('steps', {'upper', '0.1', '0.05'});
p.addParameter('rdir', ''); p.addParameter('cap', 1800); p.addParameter('force', false);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); T = [];
rdir = o.rdir; if isempty(rdir), rdir = fullfile(root, 'results', 'larger_samples', 'run'); end
switch o.part
    case 'switching'
        info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_journal_info.json')));
        pairs = o.pairs; if isempty(pairs), pairs = {info.id}; end
        for q = pairs(:).', for s = o.steps(:).', switching_row(o, root, rdir, q{1}, s{1}, info); end, end
    case 'chattering'
        info = jsondecode(fileread(fullfile(root, 'data', 'scaling_pairs_journal_info.json')));
        pairs = o.pairs; if isempty(pairs), pairs = {info.id}; end
        for q = pairs(:).', for s = o.steps(:).', chattering_row(o, root, rdir, q{1}, s{1}, info); end, end
    case 'extended'
        extended(o, root, rdir);
    case 'replay'
        T = [replays(root, 'larger_switching'), replays(root, 'larger_scaling')];
    case 'summary'
        T = summary(root);
end
end


% ============================================================================================== runs
function switching_row(o, root, rdir, pid, step, info)
say = logger(fullfile(root, 'results', 'logs', 'larger_samples.log'));
d0 = fullfile(rdir, 'switching'); cdir = fullfile(rdir, 'certificates');
if ~isfolder(d0), mkdir(d0); end
if ~isfolder(cdir), mkdir(cdir); end
file = fullfile(d0, sprintf('%s_%s.json', pid, step));
if isfile(file) && ~o.force, return, end
A = load_pair('switching_journal', pid); I = info(strcmp({info.id}, pid)); T0 = tic;
if strcmp(step, 'upper')
    cf = fullfile(cdir, sprintf('%s_schedule.json', pid));
    [C, json] = certify_best_upper(A, struct('id', pid, 'chatter', false, ...
                'schedule', struct('modes', [1 2], 'durations', [I.t1 I.t2])));
    write_text(cf, json); ok = replay_best_upper(cf);
    R = struct('pair', pid, 'c_float', I.chatter, 'schedule_float', I.sched, 'upper', C.printed, 'replay_ok', ok, ...
               'secs', toc(T0));
    say('%s schedule: lambda_check <= %s (float %+.6f, c %+.6f), replay %d', pid, C.printed, I.sched, I.chatter, ok);
else
    D = exact_data(A, step); d = D.d; h = D.hf;
    B = cellfun(@(M) eye(d) + h * M, D.Af, 'UniformOutput', false);
    [v, a, b] = best_two_block(B, round(I.t1 / h), round(I.t2 / h), 4);
    g_up = (exp(v) - 1) / h; w = [ones(1, a), 2 * ones(1, b)];
    R = struct('pair', pid, 'h', step, 'd', d, 'word', [a b], 'g_up', g_up, 'N', [], 'T', [], 'secs', 0);
    R.N = ladder(A, step, B, perron_seeds(B, w), g_up, h, [1e-5 1e-3 1e-2 0.1], [600 600 300 300], 'N', pid, cdir);
    if isnan(R.N.eps) || R.N.eps >= 0.1
        BT = cellfun(@(M) M.', B, 'UniformOutput', false);
        R.T = ladder(A, step, BT, perron_seeds(BT, fliplr(w)), g_up, h, [1e-5 1e-3 1e-2 0.1 0.15 0.2], ...
                     300 * ones(1, 6), 'T', pid, cdir);
    end
    R.secs = toc(T0);
    say('%s h = %s: word 1^%d 2^%d, g(h) <= %+.6f; N: %s (eps %g)', pid, step, a, b, g_up, R.N.lambda, R.N.eps);
end
write_text(file, jsonencode(R));
end


function L = ladder(A, h_str, F, seeds, g_up, h, epss, caps, orient, pid, cdir)
L = struct('eps', NaN, 'lambda', '', 'nV', NaN, 'replay_ok', false, 'tries', []);
tries = struct('eps', {}, 'closed', {}, 'nV', {}, 'secs', {}, 'reason', {});
for k = 1:numel(epss)
    t = 1 + h * (g_up - epss(k));
    [V, st] = polytope_closure(F, t, seeds, caps(k), 20000);
    tries(end + 1) = struct('eps', epss(k), 'closed', ~isempty(V), 'nV', st.nV, 'secs', st.secs, 'reason', st.reason); %#ok<AGROW>
    if ~isempty(V)
        out = run_lower_row(A, h_str, struct('name', pid, 'kind', 'best', 'orientation', orient, 'vertices', sym(V, 'f'), ...
                            'file', fullfile(cdir, sprintf('%s_%s_%s.json', pid, h_str, orient)), ...
                            'meta', struct('seed', 'perron', 'g_up', g_up, 'eps', epss(k), 't_build', t, 'closure_secs', st.secs)));
        L.eps = epss(k); L.lambda = out.lambda; L.nV = out.nV; L.replay_ok = out.replay_ok;
        break
    end
end
L.tries = tries;
end


function [v, a, b] = best_two_block(B, a0, b0, r)
v = Inf; a = NaN; b = NaN;
for aa = max(1, a0 - r):a0 + r
    for bb = max(1, b0 - r):b0 + r
        x = log(max(abs(eig(B{2}^bb * B{1}^aa)))) / (aa + bb);
        if x < v, v = x; a = aa; b = bb; end
    end
end
end


function chattering_row(o, root, rdir, pid, step, info)
say = logger(fullfile(root, 'results', 'logs', 'larger_samples.log'));
d0 = fullfile(rdir, 'scaling'); cdir = fullfile(rdir, 'certificates');
if ~isfolder(d0), mkdir(d0); end
if ~isfolder(cdir), mkdir(cdir); end
file = fullfile(d0, sprintf('%s_%s.json', pid, step));
if isfile(file) && ~o.force, return, end
A = load_pair('scaling_journal', pid); I = info(strcmp({info.id}, pid)); T0 = tic;
if strcmp(step, 'upper')
    opts = struct('id', pid, 'chatter', true);
    if I.switches, opts.schedule = struct('modes', [1 2], 'durations', [I.t1 I.t2]); end
    cf = fullfile(cdir, sprintf('%s_best_upper.json', pid));
    [C, json] = certify_best_upper(A, opts); write_text(cf, json); ok = replay_best_upper(cf);
    R = struct('pair', pid, 'c_float', I.chatter, 'upper', C.printed, 'kind', C.best_kind, 'replay_ok', ok, 'secs', toc(T0));
    say('%s best upper: lambda_check <= %s (%s), replay %d', pid, C.printed, C.best_kind, ok);
else
    r = certified_bracket(A, step, struct('name', pid, 'probe', 300, 'budget', 1800, 'nmax', 20000, 'outdir', cdir, ...
                                          'source', 'matlab', 'verbose', false));
    okL = replay_lower(r.lower_cert, A); okU = replay_upper(r.upper_cert, A);
    R = struct('pair', pid, 'h', step, 'lower', r.lower, 'upper', r.upper, 'nV', r.search.nV, 't_build', r.search.t_build, ...
               'probes', numel(r.search.probes), 'word', r.search.word, 'explicit_word', r.search.explicit_word, ...
               'lsr_product_bound', r.search.lsr_product_bound, 'replay_lower', okL, 'replay_upper', okU, 'secs', toc(T0));
    say('%s h = %s: %s <= lambda_check <= %s (|V| = %d), replays %d/%d', pid, step, r.lower, r.upper, r.search.nV, okL, okU);
end
write_text(file, jsonencode(R));
end


function extended(o, root, rdir)
say = logger(fullfile(root, 'results', 'logs', 'larger_samples.log'));
d0 = fullfile(rdir, 'extended'); cdir = fullfile(rdir, 'certificates'); sw = fullfile(rdir, 'switching');
if ~isfolder(d0), mkdir(d0); end
info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_journal_info.json')));
pairs = o.pairs; if isempty(pairs), pairs = {info([info.d] == 10).id}; end
for h = {'0.1', '0.05'}
    for q = pairs(:).'
        f = fullfile(sw, sprintf('%s_%s.json', q{1}, h{1}));
        if ~isfile(f), continue, end
        R1 = jsondecode(fileread(f)); lam = -Inf;
        for oo = 'NT'
            if isstruct(R1.(oo)) && ~isempty(R1.(oo).lambda), lam = max(lam, str2double(R1.(oo).lambda)); end
        end
        if R1.g_up - lam < 0.01 - 1e-9, continue, end              % resolved to eps < 1e-2: not in the search
        I = info(strcmp({info.id}, q{1})); A = load_pair('switching_journal', q{1});
        for oo = 'NT'
            tag = sprintf('%s_%s_%s', q{1}, h{1}, oo); file = fullfile(d0, [tag '.json']);
            if isfile(file) && ~o.force, continue, end
            src01 = '';
            if strcmp(h{1}, '0.05'), src01 = best_cert_01(cdir, q{1}, oo); end
            R = extended_search_row(A, q{1}, h{1}, oo, R1.g_up, R1.g_up - lam, [I.t1 I.t2], [], src01, ...
                                    struct('cap', o.cap, 'cdir', cdir, 'tag', tag));
            R.item1_lower = lam;
            write_text(file, jsonencode(R));
            say('%s extended: best %s (eps %g), %.0f s', tag, R.best_lambda, R.best_eps, R.secs);
        end
    end
end
end


function f = best_cert_01(cdir, pid, orient)
c = [{fullfile(cdir, sprintf('%s_0.1_%s.json', pid, orient))}, ...
     arrayfun(@(x) fullfile(x.folder, x.name), dir(fullfile(cdir, sprintf('%s_0.1_%s_eps*.json', pid, orient))).', ...
              'UniformOutput', false)];
f = ''; best = -Inf;
for k = 1:numel(c)
    if ~isfile(c{k}), continue, end
    C = jsondecode(fileread(c{k}));
    if str2double(C.lambda_lower) > best, best = str2double(C.lambda_lower); f = c{k}; end
end
end


% ============================================================================================== replay
function R = replays(root, folder)
t0 = tic; [okm, Tm] = replay_all({folder}); sm = toc(t0);
py = fullfile(root, 'independent_replay', 'python', 'run_all.py');
t0 = tic;
[st, txt] = system(sprintf('OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -I -B "%s" "%s"', py, ...
                   fullfile(root, 'certificates', folder)));
sp = toc(t0);
R = struct('folder', folder, 'n_files', numel(dir(fullfile(root, 'certificates', folder, '*.json'))), ...
           'matlab_n', numel(Tm), 'matlab_valid', nnz([Tm.ok]), 'matlab_ok', okm, 'matlab_secs', round(sm), ...
           'python_ok', st == 0, 'python_summary', regexp(txt, '==== Python replay of [^\n]*', 'match', 'once'), ...
           'python_secs', round(sp));
fprintf('%s: MATLAB %d of %d VALID (%d s); Python: %s (%d s)\n', folder, R.matlab_valid, R.matlab_n, R.matlab_secs, ...
        R.python_summary, R.python_secs);
end


% ============================================================================================== summary
function T = summary(root)
cs = fullfile(root, 'certificates', 'larger_switching'); cc = fullfile(root, 'certificates', 'larger_scaling');
rs = fullfile(root, 'results', 'larger_samples');
say = logger(fullfile(root, 'results', 'larger_samples.md'));
tex = {};
say('# Larger samples of Obs. 4.8 (Tables S5 and S6)');
say('');
say('From the stored certificates (certificates/larger_switching, larger_scaling) and run records');
say('(results/larger_samples).  Switching pairs: direct certificates at each step, as in Table 3; T = transposed family,');
say('+ = extended search.  Chattering pairs: at h = 0.05 the larger of the direct certificate and the h = 0.1 bound');
say('(dagger, prop:polytope(ii)), as in Table S4.  Widths: upper (schedule, or c) minus the larger lower bound of the two');
say('steps, rounded up to two digits.  LaTeX rows: results/larger_samples_tables.tex.');
% ------------------------------------------------------------------------- switching
info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_journal_info.json')));
num = arrayfun(@(I) str2double(I.id(find(I.id == '_') + 1:end)), info);
[~, ord] = sortrows([[info.d].', num(:)]); info = info(ord);
say(''); say('## Switching pairs (Table S6)'); say('');
say('| d | pair | c | upper | g(0.1) <= | g(0.1) >= | g(0.05) <= | g(0.05) >= | width |');
say('|---|---|---|---|---|---|---|---|---|');
S = struct('d', {}, 'eps01', {}, 'eps005', {}, 'width', {});
tex{end + 1} = '% Table S6 (tab:switchingnew)';
for q = 1:numel(info)
    I = info(q); pid = I.id; k = str2double(pid(find(pid == '_') + 1:end));
    U = jsondecode(fileread(fullfile(cs, [pid '_schedule.json']))); up = str2double(U.printed);
    cells = {}; lows = -Inf; ep = [NaN NaN]; hh = {'0.1', '0.05'};
    for j = 1:2
        R = jsondecode(fileread(fullfile(rs, 'switching', sprintf('%s_%s.json', pid, hh{j}))));
        [lam, mark] = best_switching(cs, pid, hh{j});
        cells{end + 1} = up6(R.g_up); %#ok<AGROW>
        if isinf(lam), cells{end + 1} = '---'; else, cells{end + 1} = sprintf('%.6f%s', lam, mk(mark)); ep(j) = R.g_up - lam; end %#ok<AGROW>
        lows = max(lows, lam);
        tg{2 * j - 1} = texnum(ceil6(R.g_up), ''); %#ok<AGROW>
        if isinf(lam), tg{2 * j} = '---'; else, tg{2 * j} = texnum(lam, mark); end %#ok<AGROW>
    end
    w = up - lows;
    say('| %d | %d | %s | %s | %s | %s | %s | %s | %s |', I.d, k, up6(I.chatter), U.printed, cells{:}, wtxt(w));
    tex{end + 1} = sprintf('%d & %d & %s & %s & %s & %s & %s & %s & %s\\\\', I.d, k, texnum(ceil6(I.chatter), ''), ...
                           texnum(up, ''), tg{:}, wtex(w)); %#ok<AGROW>
    S(end + 1) = struct('d', I.d, 'eps01', ep(1), 'eps005', ep(2), 'width', w); %#ok<AGROW>
end
say('');
for dd = [5 10 20]
    X = S([S.d] == dd);
    say('- d = %d (%d pairs): h = 0.1 %s; h = 0.05 %s; widths %s to %s (median %.2g)', dd, numel(X), buckets([X.eps01]), ...
        buckets([X.eps005]), wtxt(min([X.width])), wtxt(max([X.width])), median([X.width]));
end
% ------------------------------------------------------------------------- chattering
info = jsondecode(fileread(fullfile(root, 'data', 'scaling_pairs_journal_info.json')));
key = arrayfun(@(I) 1e4 * I.d + 100 * strcmp(I.kind, 'sparse') + str2double(regexprep(I.id, '.*[a-z]', '')), info);
[~, ord] = sort(key); info = info(ord);
W = jsondecode(fileread(fullfile(root, 'results', 'larger_samples', 'scaling_word_estimates.json')));
say(''); say('## Chattering pairs (Table S5)'); say('');
say('Explicit-word estimates are floating point, from scaling_word_estimates.json (original log excerpts); rounded upward to six decimals.');
say('');
say('| d | pair | c | lower 0.1 | g(h) <= 0.1 | upper 0.1 | lower 0.05 | upper 0.05 | width | vertices (0.1) |');
say('|---|---|---|---|---|---|---|---|---|---|');
tex{end + 1} = '% Table S5 (tab:scalingnew)';
C = struct('d', {}, 'width', {});
for q = 1:numel(info)
    I = info(q); pid = I.id; k = regexprep(pid, '.*[a-z]', '');
    B = jsondecode(fileread(fullfile(cc, [pid '_best_upper.json']))); c = str2double(B.printed);
    cm = ''; if ~strcmp(B.best_kind, 'chattering'), cm = 's'; end
    L1 = jsondecode(fileread(fullfile(cc, [pid '_0.1_best.json']))); lo1 = str2double(L1.lambda_lower);
    f5 = fullfile(cc, [pid '_0.05_best.json']); lo5 = -Inf;
    if isfile(f5), L5 = jsondecode(fileread(f5)); lo5 = str2double(L5.lambda_lower); end
    dag = lo1 > lo5; lo5b = max(lo1, lo5);
    U1 = jsondecode(fileread(fullfile(cc, [pid '_0.1.json']))); U5 = jsondecode(fileread(fullfile(cc, [pid '_0.05.json'])));
    W1 = W.records(strcmp({W.records.pair}, pid) & strcmp({W.records.h}, '0.1'));
    assert(isscalar(W1), 'larger_samples:wordEstimate', 'Missing or duplicated word estimate for %s', pid);
    gup = ceil6(str2double(W1.g_up_logged));
    w = c - max(lo1, lo5);
    say('| %d | %s %s | %s%s | %.6f | %.6f | %s | %.6f%s | %s | %s | %d |', I.d, I.kind, k, B.printed, mk(cm), lo1, gup, U1.printed, ...
        lo5b, mk(ifelse(dag, 'dagger', '')), U5.printed, wtxt(w), nverts(L1));
    tex{end + 1} = sprintf('%d & %s %s & %s & %s & %s & %s & %s & %s & %s\\\\', I.d, I.kind, k, texnum(c, cm), texnum(lo1, ''), texnum(gup, ''), ...
                           texnum(str2double(U1.printed), ''), texnum(lo5b, ifelse(dag, '\dagger', '')), ...
                           texnum(str2double(U5.printed), ''), wtex(w)); %#ok<AGROW>
    C(end + 1) = struct('d', I.d, 'width', w); %#ok<AGROW>
end
say('');
for dd = [5 10 20]
    X = C([C.d] == dd);
    say('- d = %d (%d pairs): widths %s to %s (median %.2g)', dd, numel(X), wtxt(min([X.width])), wtxt(max([X.width])), ...
        median([X.width]));
end
write_text(fullfile(root, 'results', 'larger_samples_tables.tex'), strjoin(tex, newline));
T = struct('switching', S, 'chattering', C);
end


function [lam, mark] = best_switching(cs, pid, h)
lam = -Inf; mark = '';
for f = dir(fullfile(cs, sprintf('%s_%s_*.json', pid, h))).'
    X = jsondecode(fileread(fullfile(f.folder, f.name)));
    if str2double(X.lambda_lower) > lam + 1e-12
        lam = str2double(X.lambda_lower);
        mark = [ifelse(X.orientation == 'T', '\top', '') ifelse(contains(f.name, '_eps'), '+', '')];
    end
end
end


function s = buckets(e)
E = [1.5e-5 1.5e-3 0.0105 0.0205 0.0505 0.1005 0.1505 0.2005];          % levels, plus the six-decimal rounding
lab = {'1e-5', '1e-3', '1e-2', '0.02', '0.05', '0.1', '0.15', '0.2'};
c = zeros(1, numel(E) + 1);
for x = e
    j = find(x <= E, 1);
    if isnan(x) || isempty(j), c(end) = c(end) + 1; else, c(j) = c(j) + 1; end
end
parts = {};
for j = 1:numel(E), if c(j), parts{end + 1} = sprintf('%s: %d', lab{j}, c(j)); end, end %#ok<AGROW>
if c(end), parts{end + 1} = sprintf('none: %d', c(end)); end
s = strjoin(parts, ', ');
end


function n = nverts(C)
if isfield(C, 'meta') && isfield(C.meta, 'nV'), n = C.meta.nV;
elseif iscell(C.vertices), n = numel(C.vertices);
else, n = size(C.vertices, 1); end
end


function y = ceil6(x)
y = ceil(x * 1e6 - 1e-9) / 1e6;
end


function s = up6(x)
s = sprintf('%.6f', ceil6(x));
end


function s = texnum(x, mark)
if x < 0, s = sprintf('$%.6f', x); else, s = sprintf('$\\phantom{-}%.6f', x); end
if ~isempty(mark), s = [s '^{' mark '}']; end
s = [s '$'];
end


function [m, e] = sig2up(x)
e = floor(log10(x)); m = ceil(x / 10^(e - 1) - 1e-9);
if m >= 100, m = floor(m / 10); e = e + 1; end
end


function s = wtex(x)
[m, e] = sig2up(x);
if e == 0, s = sprintf('$%.1f$', m / 10); else, s = sprintf('$%.1f\\times10^{%d}$', m / 10, e); end
end


function s = wtxt(x)
[m, e] = sig2up(x); s = sprintf('%.1fe%d', m / 10, e);
end


function s = mk(m)
if isempty(m), s = ''; else, s = ['^' strrep(m, '\top', 'T')]; end
end


function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
