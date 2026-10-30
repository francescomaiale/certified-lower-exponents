function T = table1_sandwich()
% TABLE1_SANDWICH  Table 1 (tab:sandwich) and Obs. 3.13 (obs:sandwich): the Euler sandwich in floating point for the
% symmetric pair of Ex. 3.10 (eps = 0.2) and the asymmetric planar pair (port of sandwich_explore.py).
%   For h = 0.2, 0.1, 0.05, 0.025, 0.0125 and all words of length <= 16:
%     g_est = (min_w rho(B_w)^(1/|w|) - 1)/h     (an estimate of g(h) from products, not a bound),
%     G_up  = (1 - 1/min_w rho(R_w)^(1/|w|))/h   (the bound of thm:implicit(i) for the best word),
%     exact schedule: the rate log rho(prod e^{h A_i}) / (h |w|) of that best implicit word as a schedule.
%   For the symmetric pair, g(h) and G(h) are known exactly (Ex. 3.10) and printed for comparison.  Also the slope
%   (G - g)/h against the prediction h u'(sum alpha_i (A_i - S)^2) v of Obs. 3.13.  Output: results/table1.md.
%   Runtime: about a minute.

root = paperB_root();
say = logger(fullfile(root, 'results', 'table1.md'));
hs = [0.2 0.1 0.05 0.025 0.0125];
paper = {[0.098979 0.301021 0.201335; 0.149874 0.250126 0.200333; 0.174984 0.225016 0.200083; ...
          0.187498 0.212502 0.200021; 0.193750 0.206250 0.200005], ...
         [-0.358479 -0.056019 -0.188068; -0.271832 -0.123611 -0.192983; -0.232157 -0.158550 -0.194214; ...
          -0.213186 -0.176426 -0.194521; -0.203857 -0.185482 -0.194598]};
pairs = {'sharp', 'planar'};
say('# Table 1 (tab:sandwich): Euler sandwich in floating point, words of length <= 16, MATLAB reproduction');
say('');
say('| pair | h | g_est | G_up | exact schedule | paper (g_est, G_up, schedule) | match |');
say('|---|---|---|---|---|---|---|');
T = struct('pair', {}, 'h', {}, 'vals', {}, 'match', {});
for q = 1:2
    A = cellfun(@str2double, load_pair('test', pairs{q}), 'UniformOutput', false);
    [c, a] = chattering_value(A);
    for k = 1:numel(hs)
        h = hs(k); I = eye(2);
        B = cellfun(@(M) I + h * M, A, 'UniformOutput', false);
        R = cellfun(@(M) inv(I - h * M), A, 'UniformOutput', false);
        [bw, bv] = brute_words(B, 16); [~, kb] = min(bv);
        [rw, rv] = brute_words(R, 16); [~, kr] = min(rv);
        g_est = (exp(bv(kb)) - 1) / h; G_up = (1 - exp(-rv(kr))) / h;
        w = rw{kr}; P = I;
        for x = w, P = expm(h * A{x}) * P; end
        sch = log(max(abs(eig(P)))) / (h * numel(w));
        vals = [g_est, G_up, sch];
        s = arrayfun(@(x) sprintf('%.6f', x), vals, 'UniformOutput', false);
        ps = arrayfun(@(x) sprintf('%.6f', x), paper{q}(k, :), 'UniformOutput', false);
        match = isequal(s, ps);
        say('| %s | %g | %s | %s | %s | %s | %s |', pairs{q}, h, s{:}, strjoin(ps, ', '), yesno(match));
        T(end + 1) = struct('pair', pairs{q}, 'h', h, 'vals', vals, 'match', match); %#ok<AGROW>
    end
    if q == 1
        say('| sharp, exact | | g(h) = 0.2 - (1 - sqrt(1 - h^2))/h, G(h) = 0.2 + (1 - sqrt(1 - h^2))/h: %s | | | | |', ...
            strjoin(arrayfun(@(h) sprintf('%.6f / %.6f', 0.2 - (1 - sqrt(1 - h^2)) / h, 0.2 + (1 - sqrt(1 - h^2)) / h), ...
            hs, 'UniformOutput', false), '; '));
    else
        S = a * A{1} + (1 - a) * A{2};
        [Vr, er] = eig(S); [~, i] = max(real(diag(er))); v = abs(real(Vr(:, i)));
        [Vl, el] = eig(S.'); [~, i] = max(real(diag(el))); u = abs(real(Vl(:, i))); u = u / (u.' * v);
        pred = u.' * (a * (A{1} - S)^2 + (1 - a) * (A{2} - S)^2) * v;
        obs = (T(end).vals(2) - T(end).vals(1)) / hs(end);
        say('');
        say('Obs. 3.13: planar pair, c = min_alpha eta(S_alpha) = %.10f at alpha* = %.4f (paper: -0.1949763083, 0.3905);', ...
            c, a);
        say('predicted slope of G - g: %.3f (paper 1.478); observed (G_up - g_est)/h at h = %g: %.3f (paper 1.472).', ...
            pred, hs(end), obs);
    end
end
say('');
say('Rows reproduced to the printed digits: %d of %d.', nnz([T.match]), numel(T));
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
