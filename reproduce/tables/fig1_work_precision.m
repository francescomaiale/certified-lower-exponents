function T = fig1_work_precision(varargin)
% FIG1_WORK_PRECISION  Obs. 4.9 (obs:glp) and Figure 2 (fig:wp): certified error U - (certified lower bound)
% against computing time, for our Euler closure and the naive Algorithm (L), on three pairs (port of wp.py and
% wp_plot.py).
%   Pairs and fixed certified upper bounds U: case 1440 (two-segment schedule, 0.345674607), pair 4 in dimension 5 of
%   Table 3 (sw5_4, two-segment schedule, -0.4300830) and case 1 (chattering value, -1.012395841).
%   Ours, for h in 'hs' (0.2, 0.1, 0.05, 0.025): word search on B_h (word_search), closure seeded by the Perron orbit of
%   the best word at the levels 1 + h (g_up - eps), eps = 1e-5, 1e-3, 1e-2, 1e-1, capped at 'caps' seconds (900, 600,
%   300, 300), exact check (certify_lower_exact); the lower bound is the exact level (t* - 1)/h.
%   Algorithm (L), tau = h: glp_algorithm_l with 'glp_budget' seconds of closures (1200) and 300 s per closure.
%   Times include the searches, every closure (failed ones too) and the exact checks; they depend on the machine.
%   Options: 'pairs', 'hs', 'methods' ({'ours', 'glp'}), 'part' ('run', 'plot', 'all'), 'force'.
%   Output: results/wp/<method>_<pair>_<h>.json, results/fig1_wp.pdf/.png, results/fig1.md.
%   Long: the Algorithm (L) runs use their full budget (12 runs x 20 minutes); ours at h = 0.025 take minutes.

p = inputParser;
p.addParameter('pairs', {'case1440', 'sw5_4', 'case1'}); p.addParameter('hs', {'0.2', '0.1', '0.05', '0.025'});
p.addParameter('methods', {'ours', 'glp'}); p.addParameter('part', 'all'); p.addParameter('force', false);
p.addParameter('caps', [900 600 300 300]); p.addParameter('glp_budget', 1200);
p.parse(varargin{:}); o = p.Results;
root = paperB_root(); rdir = fullfile(root, 'results', 'wp');
U = struct('case1440', 0.345674607, 'sw5_4', -0.4300830, 'case1', -1.012395841);
EPS = [1e-5 1e-3 1e-2 1e-1];
tag = sprintf('%s_%s_%s', o.part, strjoin(o.pairs, '+'), strjoin(o.methods, '+'));
say = logger(fullfile(root, 'results', 'logs', ['fig1_' tag '.log']));
if any(strcmp(o.part, {'run', 'all'}))
    for q = 1:numel(o.pairs)
        pid = o.pairs{q};
        if startsWith(pid, 'sw'), A = load_pair('switching', pid); else, A = load_pair('test', pid); end
        for hh = o.hs
            h = hh{1};
            for mm = o.methods
                file = fullfile(rdir, sprintf('%s_%s_%s.json', mm{1}, pid, h));
                if isfile(file) && ~o.force, continue, end
                if strcmp(mm{1}, 'ours'), R = ours(A, pid, h, EPS, o.caps); else, R = glp(A, pid, h, o.glp_budget); end
                write_text(file, jsonencode(R));
                if isempty(R.lower)
                    say('%-4s %-8s h=%-6s: no certificate, %.0fs', mm{1}, pid, h, R.secs);
                else
                    say('%-4s %-8s h=%-6s: lower %.9f, error %.3e, %.0fs, word length %d', mm{1}, pid, h, R.lower, ...
                        U.(pid) - R.lower, R.secs, R.word_len);
                end
            end
        end
    end
end
T = [];
if any(strcmp(o.part, {'plot', 'all'})), T = plot_wp(root, rdir, U); end
end


function R = ours(A, pid, h, EPS, caps)
T0 = tic;
D = exact_data(A, h); hf = D.hf; d = D.d;
B = cellfun(@(M) eye(d) + hf * M, D.Af, 'UniformOutput', false);
[w, v] = word_search(B, hf);
g_up = (exp(v) - 1) / hf;
seeds = perron_seeds(B, w);
R = struct('method', 'ours', 'id', pid, 'h', h, 'word_len', numel(w), 'g_up', g_up, 'tries', [], 'eps', [], ...
           'lower', [], 'lower_exact', '', 'nV', [], 'secs', 0);
tries = struct('eps', {}, 'closed', {}, 'nV', {}, 'secs', {});
for k = 1:numel(EPS)
    [V, st] = polytope_closure(B, 1 + hf * (g_up - EPS(k)), seeds, caps(k), 20000);
    tries(end + 1) = struct('eps', EPS(k), 'closed', ~isempty(V), 'nV', st.nV, 'secs', st.secs); %#ok<AGROW>
    if ~isempty(V)
        [~, ~, info] = certify_lower_exact(V, D, struct('name', pid));
        R.eps = EPS(k); R.lower = double(info.lambda_star); R.lower_exact = char(info.lambda_star); R.nV = info.nV;
        break
    end
end
R.tries = tries; R.secs = toc(T0);
end


function R = glp(A, pid, h, budget)
G = glp_algorithm_l(A, str2double(h), struct('probe', 300, 'budget', budget));
lo = []; if ~isempty(G.alpha), lo = double(G.alpha); end
R = struct('method', 'glp', 'id', pid, 'h', h, 'lower', lo, 'nu', G.nu, 'word_len', G.word_len, 'beta_l', G.beta_l, ...
           'nV', G.nV, 'secs', G.search_secs + G.closure_secs + G.lp_secs);
end


function T = plot_wp(root, rdir, U)
say = logger(fullfile(root, 'results', 'fig1.md'));
say('# Figure 2 (fig:wp): work and precision, MATLAB reproduction (times depend on the machine)');
say('');
say('| pair | method | h or tau | certified lower | error U - lower | time (s) |');
say('|---|---|---|---|---|---|');
pairs = fieldnames(U).'; H = {'0.2', '0.1', '0.05', '0.025'};
titles = struct('case1440', 'case 1440 (d = 3, switching)', 'sw5_4', 'pair 4, d = 5 (switching)', ...
                'case1', 'case 1 (d = 3, chattering)');
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 16 6], 'Color', 'w');
try, fig.Theme = 'light'; catch, end            % R2025a+ may default to the dark theme
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
style = struct('ours', {{'-o', [0.165 0.471 0.839], 'Euler closure (ours)'}}, ...
               'glp', {{'--s', [0.922 0.408 0.204], 'Algorithm (L), naive'}});
T = [];
for q = 1:numel(pairs)
    ax = nexttile(tl); hold(ax, 'on');
    for mm = {'ours', 'glp'}
        pts = zeros(0, 2); lab = {};
        for k = 1:numel(H)
            f = fullfile(rdir, sprintf('%s_%s_%s.json', mm{1}, pairs{q}, H{k}));
            if ~isfile(f), continue, end
            R = jsondecode(fileread(f));
            if isempty(R.lower)
                say('| %s | %s | %s | --- | --- | %.0f |', pairs{q}, mm{1}, H{k}, R.secs); continue
            end
            err = U.(pairs{q}) - R.lower;
            say('| %s | %s | %s | %.9f | %.3e | %.0f |', pairs{q}, mm{1}, H{k}, R.lower, err, R.secs);
            pts(end + 1, :) = [max(R.secs, 0.1), err]; lab{end + 1} = H{k}; %#ok<AGROW>
        end
        if isempty(pts), continue, end
        s = style.(mm{1});
        loglog(ax, pts(:, 1), pts(:, 2), s{1}, 'Color', s{2}, 'MarkerFaceColor', s{2}, 'LineWidth', 1.2, ...
               'MarkerSize', 4, 'DisplayName', s{3});
        text(ax, pts(:, 1) * 1.15, pts(:, 2), lab, 'FontSize', 6, 'Color', [0.32 0.32 0.31]);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.07 4000], 'XTick', [0.1 1 10 100 1000], 'FontSize', 7);
    grid(ax, 'on'); ax.XMinorGrid = 'off'; ax.YMinorGrid = 'off'; ax.Box = 'off';
    title(ax, titles.(pairs{q}), 'FontSize', 8, 'FontWeight', 'normal'); xlabel(ax, 'time (s)');
    if q == 1, ylabel(ax, 'certified error U - lower bound'); end
    if q == 1, legend(ax, 'Location', 'southwest', 'FontSize', 6, 'Box', 'off'); end
end
exportgraphics(fig, fullfile(root, 'results', 'fig1_wp.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(root, 'results', 'fig1_wp.png'), 'Resolution', 200);
close(fig);
say('');
say('Figure: results/fig1_wp.pdf (.png).  Times below 0.1 s are drawn at 0.1 s.');
end
