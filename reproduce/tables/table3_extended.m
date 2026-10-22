function T = table3_extended(varargin)
% TABLE3_EXTENDED  The extended search of Obs. 4.8 on the d = 10 rows of Table 3 (bounds marked + there).
%   Rows: the d = 10 pairs at h = 0.1 whose Table 3 bound (better orientation, before this search) is at eps >= 1e-2
%   (sw10_1, 3, 4, 6, 8, 9, 10), and all ten at h = 0.05; both orientations each (34 runs).
%   Search (helpers/extended_search_row.m): seeds from the five best words 1^a 2^b, on B_h^T at h = 0.1 also the
%   final vertex set of the run of [GM] (Appendix B), at h = 0.05 also the certified h = 0.1 polytope of Table 3 if it
%   has this orientation; levels eps = 0.15, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-5 below the current one, 1800 s and
%   20000 vertices per closure; every closed polytope certified exactly in both data readings and replayed.
%   The original GM vertex sets are bundled in results/gm_original/switching and used by default.
%   'gm_dir' overrides that directory (<pair>_0.1_E_V.mat, variable V, d x n); an explicit '' skips GM seeds.
%   'part'  'run'      the rows ('rows', e.g. {'sw10_1_0.1_N'}; default all 34); records and certificates go to
%                      results/table3_extended/run/ (option 'rdir'), so that a run never overwrites the stored ones
%           'replay'   certificates/table3_extended (the printed bounds): MATLAB replay_all (about 2.7 hours on one
%                      thread) and the Python replay (2 minutes)
%           'summary'  (default) the d = 10 rows of Table 3 with the extended search, from the stored records
%                      (results/table3_extended/*.json) and certificates: results/table3_extended.md
%   Options: 'cap' (1800), 'force'.
%   Long: the paper's runs took up to 3.5 hours per row (34 single-threaded MATLAB processes, 12 at a time).

p = inputParser;
p.addParameter('part', 'summary'); p.addParameter('rows', {}); p.addParameter('gm_dir', 'bundled'); p.addParameter('cap', 1800);
p.addParameter('rdir', ''); p.addParameter('force', false);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); T = [];
if strcmp(o.gm_dir, 'bundled'), o.gm_dir = fullfile(root, 'results', 'gm_original', 'switching'); end
switch o.part
    case 'run'
        rdir = o.rdir; if isempty(rdir), rdir = fullfile(root, 'results', 'table3_extended', 'run'); end
        cdir = fullfile(rdir, 'certificates');
        if ~isfolder(cdir), mkdir(cdir); end
        rows = o.rows; if isempty(rows), rows = all_rows(); end
        for r = rows(:).', run_row(o, root, rdir, cdir, r{1}); end
    case 'replay'
        T = replays(root, 'table3_extended');
    case 'summary'
        T = summary(root);
end
end


function rows = all_rows()
rows = {};
for k = [1 3 4 6 8 9 10], for o = 'NT', rows{end + 1} = sprintf('sw10_%d_0.1_%s', k, o); end, end %#ok<AGROW>
for k = 1:10, for o = 'NT', rows{end + 1} = sprintf('sw10_%d_0.05_%s', k, o); end, end %#ok<AGROW>
end


function run_row(o, root, rdir, cdir, tag)
say = logger(fullfile(root, 'results', 'logs', 'table3_extended.log'));
parts = strsplit(tag, '_'); pid = [parts{1} '_' parts{2}]; h = parts{3}; orient = parts{4}; row = [pid '_' h];
file = fullfile(rdir, [tag '.json']);
if isfile(file) && ~o.force, return, end
A = load_pair('switching', pid);
info = jsondecode(fileread(fullfile(root, 'data', 'switching_pairs_info.json'))); I = info(strcmp({info.id}, pid));
L = paper_levels('switching', row);
P = table3_row(pid, h); eps_cur = L.g_up - str2double(P.lower);
extra = [];
if orient == 'T' && strcmp(h, '0.1') && ~isempty(o.gm_dir)
    g = fullfile(o.gm_dir, sprintf('%s_E_V.mat', row));
    assert(isfile(g), 'table3_extended:missingSeed', 'Missing GM seed: %s', g);
    S = load(g, 'V'); extra = S.V;
    assert(isreal(extra) && size(extra, 1) == size(A{1}, 1) && ~isempty(extra) && ...
           all(isfinite(extra(:))) && all(extra(:) >= 0) && all(sum(extra, 1) > 0), ...
           'table3_extended:invalidSeed', 'Invalid GM vertex set: %s', g);
end
src01 = '';
if strcmp(h, '0.05')
    P01 = table3_row(pid, '0.1');
    if P01.orientation == orient
        if orient == 'N', src01 = fullfile(root, 'certificates', 'switching', sprintf('%s_0.1_best.json', pid));
        else, src01 = fullfile(root, 'certificates', 'table3_transposed', sprintf('%s_0.1_T.json', pid)); end
        if ~isfile(src01), src01 = ''; end
    end
end
R = extended_search_row(A, pid, h, orient, L.g_up, eps_cur, [I.t1 I.t2], extra, src01, ...
                        struct('cap', o.cap, 'cdir', cdir, 'tag', tag, 'extra_label', 'final vertex set of [GM]'));
R.paper_lower = P.lower; R.paper_orientation = P.orientation;
write_text(file, jsonencode(R));
say('%s: start eps %.3g; best %s (eps %g), %.0f s; Table 3 before: %s (%s)', tag, eps_cur, R.best_lambda, R.best_eps, ...
    R.secs, P.lower, P.orientation);
end


function P = table3_row(pid, h)
if isempty(which('paper_data')), addpath(fileparts(paperB_root())); end   % paper_data.m: top level of the package
Q = paper_data(); T = Q.table3;
P = T(strcmp({T.pair}, pid) & strcmp({T.h}, h));
end


function T = summary(root)
say = logger(fullfile(root, 'results', 'table3_extended.md'));
say('# Table 3, d = 10: the extended search of Obs. 4.8');
say('');
say('Per row: the Table 3 bound before the extended search (better orientation; N = B_h, T = B_h^T), the best bound of');
say('the extended search (records results/table3_extended/*.json; certificates of the printed values in');
say('certificates/table3_extended), the printed value and the width (schedule bound minus the larger lower bound of the');
say('two steps, rounded up to two digits).  At h = 0.05 the re-certified h = 0.1 polytope (prop:polytope(ii)) served as');
say('a seed and a baseline only: Table 3 reports direct certificates at each step.');
say('');
say('| pair | h | before | eps | extended, N | extended, T | printed | eps | width |');
say('|---|---|---|---|---|---|---|---|---|');
rdir = fullfile(root, 'results', 'table3_extended');
T = struct('pair', {}, 'h', {}, 'lower', {}, 'mark', {});
for k = 1:10
    pid = sprintf('sw10_%d', k);
    S = jsondecode(fileread(fullfile(root, 'certificates', 'switching', [pid '_schedule.json'])));
    lows = -Inf; lines = {};
    for h = {'0.1', '0.05'}
        P = table3_row(pid, h{1}); L = paper_levels('switching', [pid '_' h{1}]);
        best = str2double(P.lower); mark = ifelse(P.orientation == 'T', 'T', ''); ext = {'---', '---'};
        for q = 1:2
            o = 'NT'; f = fullfile(rdir, sprintf('%s_%s_%s.json', pid, h{1}, o(q)));
            if ~isfile(f), continue, end
            R = jsondecode(fileread(f));
            if isempty(R.best_lambda), ext{q} = '(none)'; continue, end
            ext{q} = sprintf('%s (eps %g)', R.best_lambda, R.best_eps);
            if str2double(R.best_lambda) > best + 1e-12
                best = str2double(R.best_lambda); mark = [ifelse(o(q) == 'T', 'T', '') '+'];
                c = dir(fullfile(root, 'certificates', 'table3_extended', sprintf('%s_%s_%s_eps*.json', pid, h{1}, o(q))));
                ok = any(arrayfun(@(f) strcmp(jsondecode(fileread(fullfile(f.folder, f.name))).lambda_lower, R.best_lambda), c));
                if ~ok, mark = [mark ' (certificate missing)']; end %#ok<AGROW>
            end
        end
        lows = max(lows, best);
        lines{end + 1} = {pid, h{1}, P.lower, L.g_up - str2double(P.lower), ext{1}, ext{2}, best, mark, L.g_up - best}; %#ok<AGROW>
        T(end + 1) = struct('pair', pid, 'h', h{1}, 'lower', sprintf('%.6f', best), 'mark', mark); %#ok<AGROW>
    end
    w = wid(str2double(S.printed) - lows);
    for q = 1:numel(lines)
        x = lines{q};
        say('| %s | %s | %s | %.2g | %s | %s | %.6f%s | %.2g | %s |', x{1}, x{2}, x{3}, x{4}, x{5}, x{6}, x{7}, ...
            ifelse(isempty(x{8}), '', ['^' x{8}]), x{9}, ifelse(q == 1, w, ''));
    end
end
end


function R = replays(root, folder)
t0 = tic; [okm, Tm] = replay_all({folder}); sm = toc(t0);
py = fullfile(root, 'independent_replay', 'python', 'run_all.py');
t0 = tic;
[st, txt] = system(sprintf('OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -I -B "%s" "%s"', py, ...
                   fullfile(root, 'certificates', folder)));
sp = toc(t0);
fprintf('%s', txt);
R = struct('folder', folder, 'n_files', numel(dir(fullfile(root, 'certificates', folder, '*.json'))), ...
           'matlab_n', numel(Tm), 'matlab_valid', nnz([Tm.ok]), 'matlab_ok', okm, 'matlab_secs', round(sm), ...
           'python_ok', st == 0, 'python_summary', regexp(txt, '==== Python replay of [^\n]*', 'match', 'once'), ...
           'python_secs', round(sp));
fprintf('%s: MATLAB %d of %d VALID; Python: %s\n', folder, R.matlab_valid, R.matlab_n, R.python_summary);
end


function s = wid(x)
e = floor(log10(abs(x))); y = ceil(x / 10^(e - 1) - 1e-9) * 10^(e - 1);
s = sprintf('%.1e', y);
end


function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
