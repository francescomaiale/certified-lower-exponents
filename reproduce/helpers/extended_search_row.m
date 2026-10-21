function R = extended_search_row(A, pid, h, orient, g_up, eps_cur, t12, extra, src01, o)
% EXTENDED_SEARCH_ROW  The extended search of Obs. 4.8 for one row (pair, step h, orientation) of a switching pair.
%   R = extended_search_row(A, pid, h, orient, g_up, eps_cur, t12, extra, src01, o)
%     A        data {A_1, A_2} (decimal strings, load_pair); pid its name; h the step (decimal string)
%     orient   'N' (B_h = I + h A_i) or 'T' (the transposed family, prop:polytope(iii))
%     g_up     floating-point upper bound for g(h) (best two-block word); eps_cur the current eps of the row
%              (g_up minus the best certified bound so far, Inf if none)
%     t12      [t1 t2], the durations of the optimized two-segment schedule
%     extra    d x s matrix of further seeds (e.g. the final vertex set of a run of [GM]), or []; o.extra_label
%     src01    for h = 0.05: a certificate file of a certified h = 0.1 polytope of this orientation, or ''.  Its
%              vertices are seeds, and it is re-certified at h = 0.05 (prop:polytope(ii)) as a baseline.
%     o        cap (seconds per closure, 1800), nmax (20000), cdir (certificate folder), tag (file prefix)
%   Seeds: the Perron vector of B_w and its images along w for the five best words w = 1^a 2^b with a, b within 4 of
%   t1/h, t2/h (reversed words on B_h^T), plus extra and the vertices of src01, all normalized to sum 1.  Levels
%   t = 1 + h (g_up - eps), eps = 0.15, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-5, those below eps_cur (by more than 5e-6:
%   certified bounds are rounded down to six decimals), from the easiest; the ladder stops at the first level that
%   does not close.  Every closed polytope is certified exactly in both data readings and replayed (run_lower_row):
%   certificates <cdir>/<tag>_eps<eps>.json and <tag>_transfer.json.
%   R: words, n_seeds, seeds (labels), transfer (lambda, source, replay_ok), tries (eps, closed, nV, secs, reason,
%   lambda), best_lambda, best_eps, secs.

T0 = tic;
D = exact_data(A, h); d = D.d; hf = D.hf;
B = cellfun(@(M) eye(d) + hf * M, D.Af, 'UniformOutput', false);
a0 = round(t12(1) / hf); b0 = round(t12(2) / hf); cand = [];
for aa = max(1, a0 - 4):a0 + 4
    for bb = max(1, b0 - 4):b0 + 4
        cand(end + 1, :) = [log(max(abs(eig(B{2}^bb * B{1}^aa)))) / (aa + bb), aa, bb]; %#ok<AGROW>
    end
end
cand = sortrows(cand, 1); words = cand(1:5, 2:3);
if orient == 'T', F = cellfun(@(M) M.', B, 'UniformOutput', false); else, F = B; end
seeds = zeros(d, 0);
for k = 1:size(words, 1)
    w = [ones(1, words(k, 1)), 2 * ones(1, words(k, 2))];
    if orient == 'T', w = fliplr(w); end
    seeds = [seeds, perron_seeds(F, w)]; %#ok<AGROW>
end
src = {sprintf('Perron orbits of 5 words (best 1^%d 2^%d)', words(1, 1), words(1, 2))};
if ~isempty(extra)
    seeds = [seeds, extra ./ sum(extra, 1)]; src{end + 1} = sprintf('%s (%d)', getopt(o, 'extra_label', 'extra seeds'), size(extra, 2));
end
R = struct('pair', pid, 'h', h, 'orientation', orient, 'g_up', g_up, 'eps_start', eps_cur, 'words', words, ...
           'cap', getopt(o, 'cap', 1800), 'transfer', [], 'tries', [], 'best_lambda', '', 'best_eps', NaN);
cdir = o.cdir; tag = o.tag;
if ~isempty(src01)
    C = jsondecode(fileread(src01)); Vq = cert_vertices(C); Vd = double(Vq);
    seeds = [seeds, Vd ./ sum(Vd, 1)]; src{end + 1} = sprintf('certified h = 0.1 polytope (%d)', size(Vd, 2));
    out = run_lower_row(A, h, struct('name', pid, 'kind', 'transfer', 'orientation', orient, 'vertices', Vq, ...
                        'file', fullfile(cdir, sprintf('%s_transfer.json', tag)), ...
                        'meta', struct('source', file_name(src01), 'source_lambda', C.lambda_lower)));
    R.transfer = struct('lambda', out.lambda, 'source', file_name(src01), 'source_lambda', C.lambda_lower, ...
                        'replay_ok', out.replay_ok);
    eps_cur = min(eps_cur, g_up - str2double(out.lambda));
end
R.seeds = strjoin(src, '; '); R.n_seeds = size(seeds, 2);
tries = struct('eps', {}, 'closed', {}, 'nV', {}, 'secs', {}, 'reason', {}, 'lambda', {});
for e = [0.15 0.1 0.05 0.02 0.01 1e-3 1e-5]
    if e >= eps_cur - 5e-6, continue, end
    t = 1 + hf * (g_up - e);
    [V, st] = polytope_closure(F, t, seeds, R.cap, getopt(o, 'nmax', 20000));
    lam = '';
    if ~isempty(V)
        out = run_lower_row(A, h, struct('name', pid, 'kind', 'best', 'orientation', orient, 'vertices', sym(V, 'f'), ...
                            'file', fullfile(cdir, sprintf('%s_eps%g.json', tag, e)), ...
                            'meta', struct('seeds', R.seeds, 'eps', e, 't_build', t, 'closure_secs', st.secs)));
        lam = out.lambda; R.best_lambda = lam; R.best_eps = e;
    end
    tries(end + 1) = struct('eps', e, 'closed', ~isempty(V), 'nV', st.nV, 'secs', st.secs, 'reason', st.reason, ...
                            'lambda', lam); %#ok<AGROW>
    if isempty(V), break, end
end
R.tries = tries; R.secs = toc(T0);
end


function f = file_name(p)
[~, n, x] = fileparts(p); f = [n x];
end
