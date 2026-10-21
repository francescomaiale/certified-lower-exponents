function T = table8_rates(varargin)
% TABLE8_RATES  Obs. 4.9 (obs:rates) and Table 11 (tab:rates): predicted and observed rate constants of the Euler
% sandwich for fourteen pairs with apparent optimal chattering (port of rate_constants.py); floating point.
%   Conjectured asymptotic diagnostic: at the chattering point alpha*, S = S_alpha*, eta = eta(S),
%     u'S = eta u', S v = eta v, u'v = 1, Q = u'(sum_i alpha*_i A_i^2) v - eta^2;
%     (G - g)/h -> Q and ((g + G)/2 - eta)/h -> 0 are not consequences of the one-sided rate theorem alone.
%   Observed: g(h) and G(h) estimated by the best words of B_h and R_h (word_search: all short words, two-segment
%   words, balanced words p/q with q <= min(4/h, 320), local search), for h = 0.1, 0.05, 0.025.
%   Also the extrapolation remark for pair 4 in dimension 5 of Table 3 (needs certificates/switching).
%   Options: 'pairs' (default: planar, case1 and the twelve scaling pairs).  Output: results/table8.md.
%   Runtime: 10 to 30 minutes (the d = 10 searches at h = 0.025 dominate).

p = inputParser; p.addParameter('pairs', {}); p.parse(varargin{:});
root = paperB_root();
[~, pil] = load_pair('scaling');
pairs = p.Results.pairs; if isempty(pairs), pairs = [{'planar', 'case1'}, pil]; end
P = struct('planar', [0.3905 1.4780 1.4861 1.4800 1.4785], 'case1', [0.7890 0.5480 0.5619 0.5507 0.5485], ...
    'd5_dense1', [0.8548 0.1785 0.1785 0.1785 0.1785], 'd5_dense2', [0.1906 0.5888 0.5901 0.5898 0.5889], ...
    'd5_dense3', [0.2520 0.1975 0.1973 0.1974 0.1975], 'd5_sparse1', [0.2068 1.3307 1.3273 1.3282 1.3304], ...
    'd5_sparse2', [0.8950 0.2024 0.2026 0.2025 0.2025], 'd5_sparse3', [0.4736 2.3887 2.4183 2.4010 2.3903], ...
    'd10_dense1', [0.2085 0.0447 0.0445 0.0446 0.0447], 'd10_dense2', [0.4063 0.3178 0.3191 0.3183 0.3179], ...
    'd10_dense3', [0.4197 0.0770 0.0768 0.0770 0.0770], 'd10_sparse1', [0.8253 0.7134 0.6981 0.7091 0.7123], ...
    'd10_sparse2', [0.3296 0.4656 0.4732 0.4686 0.4671], 'd10_sparse3', [0.5997 1.5941 1.6197 1.6006 1.5959]);
H = [0.1 0.05 0.025];
say = logger(fullfile(root, 'results', 'table8.md'));
say('# Table 11 (tab:rates): predicted and observed rate constants, MATLAB reproduction (floating point)');
say('');
say('| pair | alpha*_1 | predicted Q | (G_est-g_est)/h, h = 0.1 | 0.05 | 0.025 | (c-g_est)/h at 0.025 | (G_est-c)/h at 0.025 | ((g_est+G_est)/2 - c)/h^2 at 0.025 | paper (original columns) | match |');
say('|---|---|---|---|---|---|---|---|---|---|---|');
T = struct('pair', {}, 'vals', {}, 'dev', {}, 'mid', {}, 'separate_slopes', {}, 'match', {});
for q = 1:numel(pairs)
    pid = pairs{q};
    if ismember(pid, {'planar', 'case1'}), A = load_pair('test', pid); else, A = load_pair('scaling', pid); end
    A = cellfun(@str2double, A, 'UniformOutput', false); d = size(A{1}, 1); I = eye(d);
    eta = @(a) max(real(eig(a * A{1} + (1 - a) * A{2})));
    a = fminbnd(eta, 0, 1, optimset('TolX', 1e-13));
    S = a * A{1} + (1 - a) * A{2};
    [Vr, er] = eig(S); [e, i] = max(real(diag(er))); v = abs(real(Vr(:, i)));
    [Vl, el] = eig(S.'); [~, i] = max(real(diag(el))); u = abs(real(Vl(:, i))); u = u / (u.' * v);
    Q = u.' * (a * A{1}^2 + (1 - a) * A{2}^2) * v - e^2;
    sl = zeros(1, 3); mid = NaN;
    for k = 1:3
        h = H(k);
        B = cellfun(@(M) I + h * M, A, 'UniformOutput', false);
        R = cellfun(@(M) inv(I - h * M), A, 'UniformOutput', false);
        [~, vb] = word_search(B, h); [~, vr] = word_search(R, h);
        g = (exp(vb) - 1) / h; G = (1 - exp(-vr)) / h;
        sl(k) = (G - g) / h;
        if k == 3, mid = ((g + G) / 2 - e) / h^2; end
    end
    vals = [a, Q, sl];
    s = arrayfun(@(x) sprintf('%.4f', x), vals, 'UniformOutput', false);
    pv = arrayfun(@(x) sprintf('%.4f', x), P.(pid), 'UniformOutput', false);
    match = isequal(s, pv);
    separate = [(e - g) / H(3), (G - e) / H(3)];
    say('| %s | %s | %s | %s | %s | %s | %.4f | %.4f | %+.3f | %s | %s |', ...
        pid, s{:}, separate(1), separate(2), mid, strjoin(pv, ', '), yesno(match));
    T(end + 1) = struct('pair', pid, 'vals', vals, 'dev', abs(sl - Q) / Q, 'mid', mid, ...
        'separate_slopes', separate, 'match', match); %#ok<AGROW>
end
D = reshape([T.dev], 3, []).';
say('');
say('Relative deviation |(G-g)/h - Q|/Q: max %s, median %s at h = 0.1, 0.05, 0.025 (paper: max 2.5%%, 0.64%%, 0.31%%;', ...
    mat2str(round(100 * max(D, [], 1), 2)), mat2str(round(100 * median(D, 1), 2)));
say('medians 0.39%%, 0.15%%, 0.03%%).  max |(g+G)/2 - eta| / h^2 at h = 0.025: %.2f (paper: <= 0.83).', max(abs([T.mid])));
f4 = @(h) fullfile(root, 'certificates', 'switching', sprintf('sw5_4_%s_best.json', h));
fs = fullfile(root, 'certificates', 'switching', 'sw5_4_schedule.json');
if isfile(f4('0.1')) && isfile(f4('0.05')) && isfile(fs)
    g1 = str2double(jsondecode(fileread(f4('0.1'))).lambda_lower);
    g2 = str2double(jsondecode(fileread(f4('0.05'))).lambda_lower);
    U = jsondecode(fileread(fs)).printed;
    say('Switching pair 4 (d = 5): 2 L(0.05) - L(0.1) = %.4f, where L is the certified polytope lower bound; this exceeds the certified upper bound %s (paper: -0.4270, -0.4301).', ...
        2 * g2 - g1, U);
end
say('Rows reproduced to the printed digits: %d of %d.', nnz([T.match]), numel(T));
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
