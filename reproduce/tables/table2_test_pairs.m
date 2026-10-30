function T = table2_test_pairs(varargin)
% TABLE2_TEST_PAIRS  Table 2 (tab:gripenberg), Remark 4.5 (rem:certified) and Obs. 4.7: certified Euler brackets
% for the planar pair and cases 1, 1440 and 148.
%   T = table2_test_pairs('mode', 'reproduce')     (default) every lower bound is rebuilt by one closure at the level
%                                                   where the paper's bisection last closed
%                                                   (data/paper_search_results.json), then checked exactly
%   T = table2_test_pairs('mode', 'search')         bisection with machine-independent caps instead
%                                                   ('pmax', 'maxprobes'; see polytope_search)
%   Further options: 'cdir' (certificate folder, e.g. a scratch folder for a search-mode run), 'rows' (cellstr subset of the row names below, default all), 'part' ('lower', 'upper',
%   'transfer', 'summary' or 'all'), 'force' (recompute existing certificates), 'replay' (true).
%   Certificates: certificates/table2/<row>_{best,transfer}.json (lower) and <row>.json (upper).
%   The summary (results/table2.md) compares every printed digit with the paper.
%
%   Long rows (reproduce mode, single thread, 16-core Mac): case1_0.025 (3130 vertices in the paper; about an hour),
%   case148_0.1 and planar_0.0125 (minutes); the others take seconds to a minute.

p = inputParser;
p.addParameter('mode', 'reproduce'); p.addParameter('rows', {}); p.addParameter('part', 'all');
p.addParameter('force', false); p.addParameter('cdir', ''); p.addParameter('replay', true);
p.addParameter('pmax', 4000); p.addParameter('maxprobes', 30);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); cdir = o.cdir; if isempty(cdir), cdir = fullfile(root, 'certificates', 'table2'); end
% paper values: lower, upper (tab:gripenberg); c rounded up
P = {'planar_0.1',    '-0.272230', '-0.123619'; 'planar_0.05',   '-0.232706', '-0.158706';
     'planar_0.025',  '-0.213641', '-0.176676'; 'planar_0.0125', '-0.204261', '-0.185783';
     'case1_0.1',     '-1.044234', '-0.988047'; 'case1_0.05',    '-1.027040', '-0.999505';
     'case1_0.025',   '-1.019518', '-1.005757';
     'case1440_0.1',  '0.344598',  '0.346673';  'case1440_0.05', '0.345146',  '0.346182';
     'case1440_0.025', '0.345146', '0.345931';
     'case148_0.1',   '-1.247834', '-1.236958'; 'case148_0.05',  '-1.247834', '-1.238961';
     'case148_0.025', '-1.247834', '-1.240094'};
direct = {'planar_0.1', 'planar_0.05', 'planar_0.025', 'planar_0.0125', 'case1_0.1', 'case1_0.05', ...
          'case1_0.025', 'case1440_0.1', 'case1440_0.05', 'case148_0.1'};
transfer = {'case148_0.05', 'case148_0.1'; 'case148_0.025', 'case148_0.1'; 'case1440_0.025', 'case1440_0.05'};
rows = o.rows; if isempty(rows), rows = P(:, 1).'; end
part = o.part;
say = logger(fullfile(root, 'results', 'logs', sprintf('table2_%s_%s.log', part, strjoin(rows, '+'))));
if numel(rows) > 4, say = logger(fullfile(root, 'results', 'logs', sprintf('table2_%s.log', part))); end

% ------------------------------------------------------------------------------------------- lower, direct
if any(strcmp(part, {'lower', 'all'}))
    for r = intersect(rows, direct, 'stable')
        row = r{1}; [pair, h] = split_row(row);
        file = fullfile(cdir, [row '_best.json']);
        if isfile(file) && ~o.force, say('%s: lower certificate exists (%s)', row, file); continue, end
        A = load_pair('test', pair); d = size(A{1}, 1);
        spec = struct('name', pair, 'kind', 'best', 'file', file, 'replay', o.replay, 'seeds', eye(d));
        if strcmp(pair, 'case148'), spec.seeds = [1 0; 0 0; 0 1]; end   % e_1, e_3: no vertex on the axis of e_2
        L = paper_levels('table2', row);
        if strcmp(o.mode, 'reproduce')
            spec.level = L.t_build;
        else
            hf = str2double(h); Dd = exact_data(A, h);
            Fb = cellfun(@(M) eye(d) + hf * M, Dd.Af, 'UniformOutput', false);
            [~, vb] = word_search(Fb, hf, struct('ntop', 2));
            t_up = min([1 + hf * chattering_value(Dd.Af), exp(vb)]);
            spec.search = struct('t_up', t_up, 'pmax', o.pmax, 'maxprobes', o.maxprobes, 'probe', Inf, 'budget', Inf);
            if strcmp(pair, 'case148'), spec.search.t_start = 1 + hf * -3.669288; end
        end
        say('%s: lower bound, mode %s ...', row, o.mode);
        out = run_lower_row(A, h, spec);
        if ~out.closed, say('%s: closure did not finish', row); continue, end
        say('%s: CERTIFIED lambda_check >= %s (exact (t*-1)/h = %+.10f), |V| = %d, closure %.0fs, exact %.0fs, replay %s %.0fs', ...
            row, out.lambda, out.lambda_star, out.nV, out.secs_closure, out.secs_exact, tf(out.replay_ok), out.secs_replay);
    end
end
% ------------------------------------------------------------------------------------------- lower, transfers
if any(strcmp(part, {'transfer', 'all'}))
    for k = 1:size(transfer, 1)
        row = transfer{k, 1}; src = transfer{k, 2};
        if ~ismember(row, rows), continue, end
        file = fullfile(cdir, [row '_transfer.json']);
        if isfile(file) && ~o.force, say('%s: transfer certificate exists', row); continue, end
        sfile = fullfile(cdir, [src '_best.json']);
        if ~isfile(sfile), say('%s: source certificate %s missing; run the lower part first', row, src); continue, end
        [pair, h] = split_row(row); A = load_pair('test', pair);
        S = jsondecode(fileread(sfile));
        spec = struct('name', pair, 'kind', 'transfer', 'vertices', cert_vertices(S), 'file', file, ...
                      'replay', o.replay, 'meta', struct('source', [src '_best'], 'source_lambda', S.lambda_lower));
        out = run_lower_row(A, h, spec);
        say('%s: TRANSFER of the polytope of %s (lambda >= %s): lambda_check >= %s at h = %s (exact %+.10f), |V| = %d, %.0fs, replay %s', ...
            row, src, S.lambda_lower, out.lambda, h, out.lambda_star, out.nV, out.secs_exact, tf(out.replay_ok));
    end
end
% ------------------------------------------------------------------------------------------- upper
if any(strcmp(part, {'upper', 'all'}))
    for r = rows
        row = r{1}; [pair, h] = split_row(row);
        file = fullfile(cdir, [row '.json']);
        if isfile(file) && ~o.force, say('%s: upper certificate exists', row); continue, end
        A = load_pair('test', pair);
        out = run_upper_row(A, h, struct('name', pair, 'file', file, 'replay', o.replay, 'source', 'table2'));
        say('%s: CERTIFIED lambda_check <= %s (word length %d; float %+.9f), %.0fs + %.0fs, replay %s', row, out.upper, ...
            out.k, out.U_float, out.secs_search, out.secs_exact, tf(out.replay_ok));
    end
end
% ------------------------------------------------------------------------------------------- summary
T = [];
if any(strcmp(part, {'summary', 'all'})), T = summary(P, cdir, root); end
end


function T = summary(P, cdir, root)
say = logger(fullfile(root, 'results', 'table2.md'));
say('# Table 2 (tab:gripenberg): certified Euler brackets, MATLAB reproduction');
say('');
say('| pair | c (float, rounded up) | h | lower (MATLAB) | lower (paper) | upper (MATLAB) | upper (paper) | width | |V| | word length | match |');
say('|---|---|---|---|---|---|---|---|---|---|---|');
T = struct('row', {}, 'lower', {}, 'upper', {}, 'match', {});
cs = struct();
for k = 1:size(P, 1)
    row = P{k, 1}; [pair, h] = split_row(row);
    if ~isfield(cs, pair)
        A = load_pair('test', pair); cs.(pair) = chattering_value(cellfun(@str2double, A, 'UniformOutput', false));
    end
    lf = fullfile(cdir, [row '_best.json']); tfile = fullfile(cdir, [row '_transfer.json']);
    lo = ''; nV = NaN; dag = '';
    if isfile(lf)
        C = jsondecode(fileread(lf)); lo = C.lambda_lower; nV = numel(C.vertices);
    end
    if isfile(tfile)
        C = jsondecode(fileread(tfile));
        lo2 = C.meta.source_lambda;       % prop:polytope(ii): the bound of the larger step holds here
        if str2double(C.lambda_lower) < str2double(lo2)
            warning('table2:transfer', '%s: transfer certificate below its source bound', row);
        end
        if isempty(lo) || str2double(lo2) > str2double(lo), lo = lo2; nV = numel(C.vertices); dag = '^dagger'; end
    end
    up = ''; k2 = NaN;
    uf = fullfile(cdir, [row '.json']);
    if isfile(uf), U = jsondecode(fileread(uf)); up = U.printed; k2 = U.k; end
    match = strcmp(lo, P{k, 2}) && strcmp(up, P{k, 3});
    w = str2double(up) - str2double(lo);
    say('| %s | %s | %s | %s%s | %s | %s | %s | %.1e | %d | %d | %s |', pair, fmt_dec(sym(cs.(pair), 'f'), 6, 'ceil'), h, ...
        lo, dag, P{k, 2}, up, P{k, 3}, ceil_sig(w), nV, k2, yesno(match));
    T(end + 1) = struct('row', row, 'lower', lo, 'upper', up, 'match', match); %#ok<AGROW>
end
say('');
say('Width = upper - lower (MATLAB values), rounded up to two digits.  ^dagger: the bound certified at the larger step (case 148: h = 0.1;');
say('case 1440: h = 0.05), which holds here by prop:polytope(ii); the transfer certificate re-checks that polytope');
say('directly at this step (certificates/table2/<row>_transfer.json).');
% Remark 4.5 details
say('');
say('## Remark 4.5 and Obs. 4.7');
A = load_pair('test', 'case1'); D = exact_data(A, '0.1');
B = cellfun(@(M) eye(3) + 0.1 * M, D.Af, 'UniformOutput', false);
[wb, vb] = word_search(B, 0.1, struct('ntop', 2));
say('- case 1, h = 0.1: best explicit word (length %d) gives g(0.1) <= %s (float, rounded up); certified g(0.1) >= %s.', ...
    numel(wb), fmt_dec(sym((exp(vb) - 1) / 0.1, 'f'), 6, 'ceil'), T(5).lower);
f148 = fullfile(cdir, 'case148_0.1_best.json');
if isfile(f148)
    V = double(cert_vertices(f148));
    onaxis = nnz(V(1, :) == 0 & V(3, :) == 0);
    say('- case 148, h = 0.1: %d vertices, %d of them on the axis of e_2; min (v_1 + v_3)/1''v = %.4f.', size(V, 2), ...
        onaxis, min((V(1, :) + V(3, :)) ./ sum(V, 1)));
end
A = load_pair('test', 'case148'); a22 = sym(str2double(A{2}(2, 2)), 'f');
say('- case 148: lambda_min <= (A_2)_22 = %s < -1.368009 (exact binary64 value %.16f).', fmt_dec(a22, 6, 'floor'), double(a22));
say('- all rows reproduced to the printed digits: %s (%d of %d)', yesno(all([T.match])), nnz([T.match]), numel(T));
end


function [pair, h] = split_row(row)
k = find(row == '_', 1, 'last'); pair = row(1:k - 1); h = row(k + 1:end);
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end


function y = ceil_sig(x)
% round up to two significant digits
e = floor(log10(abs(x))); y = ceil(x / 10^(e - 1) - 1e-9) * 10^(e - 1);
end
