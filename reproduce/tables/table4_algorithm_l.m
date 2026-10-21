function T = table4_algorithm_l(varargin)
% TABLE4_ALGORITHM_L  Obs. 4.11 (obs:glp) and the Alg. (L) columns of Table 4: the naive Algorithm (L) of
% Guglielmi, Laglia and Protasov on the 15 pairs (all but case 148), tau = 0.1 and 0.05 (glp_algorithm_l), and the
% tangent bound on our polytopes (planar pair at every step, case 1 at h = 0.1).
%   'mode', 'reproduce'  (default) a single closure at the shift nu found by the paper's search
%                        (data/paper_search_results.json), then the exact tangent bound
%   'mode', 'search'     the search of glp_algorithm_l with time caps 'probe' (300 s) and 'budget' (1800 s), or
%                        machine-independent caps 'pmax'
%   Options: 'pairs', 'taus' ({'0.1', '0.05'}), 'part' ('glp', 'tangent', 'all'), 'force'.
%   Output: results/glp/glp_<pair>_<tau>.json and tangent_<pair>_<h>.json (read by table4_best_brackets), and
%   results/algorithm_l.md.  The lower bounds alpha are exact (both data readings); the upper bounds beta are
%   floating point, as in the paper.  Not replayed independently (as in the paper).
%   Long: closures of 400 to 3300 vertices; d = 10 rows take 10 to 60 minutes each in reproduce mode.

p = inputParser;
p.addParameter('mode', 'reproduce'); p.addParameter('pairs', {}); p.addParameter('taus', {'0.1', '0.05'});
p.addParameter('part', 'all'); p.addParameter('force', false);
p.addParameter('probe', 300); p.addParameter('budget', 1800); p.addParameter('pmax', Inf);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); rdir = fullfile(root, 'results', 'glp');
[~, pil] = load_pair('scaling');
pairs = o.pairs; if isempty(pairs), pairs = [{'planar', 'case1', 'case1440'}, pil]; end
tag = o.part; if numel(pairs) <= 4, tag = [o.part '_' strjoin(pairs, '+')]; end
say = logger(fullfile(root, 'results', 'logs', sprintf('algorithm_l_%s.log', tag)));
T = struct('pair', {}, 'tau', {}, 'alpha', {}, 'beta', {});
if any(strcmp(o.part, {'glp', 'all'}))
    for q = 1:numel(pairs)
        pid = pairs{q};
        if ismember(pid, {'planar', 'case1', 'case1440'}), A = load_pair('test', pid); else, A = load_pair('scaling', pid); end
        for tt = o.taus
            tau = tt{1}; file = fullfile(rdir, sprintf('glp_%s_%s.json', pid, tau));
            if isfile(file) && ~o.force, continue, end
            L = paper_levels('glp', sprintf('%s_%s', pid, tau));
            gopts = struct('probe', o.probe, 'budget', o.budget, 'pmax', o.pmax);
            if strcmp(o.mode, 'reproduce')
                if isempty(L.nu)                     % no closure finished within the paper's budget: nothing to rebuild
                    R = struct('tau', str2double(tau), 'word', '', 'word_len', L.word_len, 'beta_l', L.beta_l, 'nu', [], ...
                               'probes', [], 'nV', [], 'alpha', [], 'alpha_str', '', 'alpha_float', NaN, 'secs', 0, ...
                               'id', pid, 'mode', o.mode, 'paper', struct('alpha', L.alpha, 'nu', L.nu, 'nV', L.nV), ...
                               'note', 'no closure finished within the paper''s budget; not rerun in reproduce mode');
                    write_text(file, jsonencode(R));
                    say('%s tau = %s: no closure in the paper''s run (---); not rerun', pid, tau);
                    continue
                end
                gopts.nu = L.nu;
            end
            R = glp_algorithm_l(A, str2double(tau), gopts);
            R.id = pid; R.mode = o.mode; R.paper = struct('alpha', L.alpha, 'nu', L.nu, 'nV', L.nV);
            if isa(R.alpha, 'sym'), R.alpha = char(R.alpha); R.alpha_exact = char(R.alpha_exact); end
            write_text(file, jsonencode(R));
            say('%s tau = %s: candidate length %d, beta = %+.6f; nu = %.6g; |V| = %d; alpha >= %s (paper %s, |V| %d); %.0fs', ...
                pid, tau, R.word_len, R.beta_l, R.nu, R.nV, R.alpha_str, L.alpha, L.nV, R.secs);
        end
    end
end
if any(strcmp(o.part, {'tangent', 'all'}))
    for r = {'planar_0.1', 'planar_0.05', 'planar_0.025', 'planar_0.0125', 'case1_0.1'}
        cf = fullfile(root, 'certificates', 'table2', [r{1} '_best.json']);
        file = fullfile(rdir, ['tangent_' r{1} '.json']);
        if ~isfile(cf) || (isfile(file) && ~o.force), continue, end
        k = find(r{1} == '_', 1, 'last'); pid = r{1}(1:k - 1);
        C = jsondecode(fileread(cf)); V = double(cert_vertices(C));
        t0 = tic; [a, af] = tangent_bound(V, load_pair('test', pid));
        lam = floor(a * sym(10)^6) / sym(10)^6;
        res = struct('row', r{1}, 'cert_lower', C.lambda_lower, 'alpha', fmt_dec(lam, 6), 'alpha_float', af, ...
                     'nV', size(V, 2), 'secs', toc(t0));
        write_text(file, jsonencode(res));
        say('tangent bound on our polytope %s: (t-1)/h = %s, tangent alpha >= %s (float %+.6f), |V| = %d', r{1}, ...
            C.lambda_lower, res.alpha, af, size(V, 2));
    end
end
% summary
say2 = logger(fullfile(root, 'results', 'algorithm_l.md'));
say2('# Algorithm (L) (obs:glp), MATLAB reproduction');
say2('');
say2('| pair | tau | candidate length | beta (float) | nu | |V| | alpha (exact, floored) | paper alpha | match |');
say2('|---|---|---|---|---|---|---|---|---|');
fs = dir(fullfile(rdir, 'glp_*.json'));
for f = fs.'
    R = jsondecode(fileread(fullfile(f.folder, f.name)));
    pa = ''; if isfield(R, 'paper') && ~isempty(R.paper.alpha), pa = R.paper.alpha; end
    say2('| %s | %g | %d | %+.6f | %.6g | %d | %s | %s | %s |', R.id, R.tau, R.word_len, R.beta_l, R.nu, R.nV, ...
         ifempty(R.alpha_str, '---'), ifempty(pa, '---'), yesno(strcmp(R.alpha_str, pa)));
    T(end + 1) = struct('pair', R.id, 'tau', R.tau, 'alpha', R.alpha_str, 'beta', R.beta_l); %#ok<AGROW>
end
fs = dir(fullfile(rdir, 'tangent_*.json'));
say2('');
for f = fs.'
    R = jsondecode(fileread(fullfile(f.folder, f.name)));
    say2('- tangent bound on our polytope %s: %s, against its level (t-1)/h = %s (|V| = %d)', R.row, R.alpha, R.cert_lower, R.nV);
end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end


function s = ifempty(s, d)
if isempty(s), s = d; end
end
