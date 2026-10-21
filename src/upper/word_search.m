function [w, v, info] = word_search(F, h, opts)
% WORD_SEARCH  Floating-point search for a word w minimizing log rho(F_w) / |w| over a family F.
%   [w, v, info] = word_search(F, h, opts)
%     F     cell array of nonnegative d x d matrices, e.g. R_i = (I - h A_i)^{-1} (or B_i = I + h A_i)
%     h     step, used only for the default search sizes
%     opts  L (all words of length <= L; default 16 for d <= 3, 14 for d = 5, and 12 otherwise), n (words 1^a 2^b with a, b <= n; default
%           min(ceil(4/h), 400)), Q (balanced words p/q with q <= Q; default min(ceil(4/h), 320)), ntop (local
%           searches from the ntop best candidates; default 6), sweeps (40), kmax (1200)
%   Returns the word w (row vector of mode indices 1..m; the first letter acts first, F_w = F_{w_k} ... F_{w_1}),
%   its value v = log rho(F_w) / |w|, and info (candidates, secs).  For F = R_h, (1 - exp(-v))/h is the upper
%   bound of thm:implicit(i); for F = B_h, exp(v) is an upper bound for lsr(B_h).  Only the result matters:
%   certify_upper_exact certifies the word exactly.
%
%   Search: all words up to length L (best per length), two-segment words 1^a 2^b,
%   balanced (Christoffel) words of every frequency p/q, then single-letter flip / insert / delete local search
%   from the best candidates.  Two-segment and balanced words are formed for the first two modes.

if nargin < 3, opts = struct(); end
T0 = tic;
d = size(F{1}, 1);
Ldef = 12; if d <= 3, Ldef = 16; elseif d == 5, Ldef = 14; end
% Preserve the pair defaults; cap the final enumeration batch for larger families.
if numel(F) > 2, Ldef = max(1, floor(log(2^Ldef) / log(numel(F)))); end
L     = getopt(opts, 'L', Ldef);
n     = getopt(opts, 'n', min(ceil(4 / h), 400));
Q     = getopt(opts, 'Q', min(ceil(4 / h), 320));
ntop  = getopt(opts, 'ntop', 6);
sweeps = getopt(opts, 'sweeps', 40);
kmax  = getopt(opts, 'kmax', 1200);
cw = {}; cv = [];
[bw, bv] = brute_words(F, L);
cw = [cw, bw]; cv = [cv, bv];
if numel(F) >= 2
    [tw, tv] = two_segment(F, n); cw{end + 1} = tw; cv(end + 1) = tv;
    [qw, qv] = balanced(F, Q); cw{end + 1} = qw; cv(end + 1) = qv;
end
% unique words, sorted by value, best ntop
keys = cellfun(@(x) sprintf('%d,', x), cw, 'UniformOutput', false);
[~, first] = unique(keys, 'stable');
cw = cw(first); cv = cv(first);
[cv, o] = sort(cv); cw = cw(o);
cw = cw(1:min(ntop, numel(cw))); cv = cv(1:numel(cw));
iw = cell(1, numel(cw)); iv = zeros(1, numel(cw));
for k = 1:numel(cw)
    [iv(k), iw{k}] = local_search(F, cw{k}, sweeps, kmax);
end
allw = [iw, cw]; allv = [iv, cv];
[v, k] = min(allv); w = allw{k};
info = struct('candidates', {cw}, 'candidate_values', cv, 'secs', toc(T0));
end


function l = batch_logrho(P)
% log spectral radius of every page of P (d x d x N), as a row
l = reshape(log(max(abs(pageeig(P)), [], 1)), 1, []);
end


function [P, lg] = powers(Fi, n)
d = size(Fi, 1); P = zeros(d, d, n + 1); lg = zeros(1, n + 1); P(:, :, 1) = eye(d);
for a = 1:n
    M = Fi * P(:, :, a); s = max(abs(M(:))); P(:, :, a + 1) = M / s; lg(a + 1) = lg(a) + log(s);
end
end


function [bw, bv] = two_segment(F, n)
% best word 1^a 2^b, a, b <= n
[P1, l1] = powers(F{1}, n); [P2, l2] = powers(F{2}, n);
bv = Inf; bw = [];
b = 1:n;
for a = 1:n
    v = (batch_logrho(pagemtimes(P2(:, :, 2:end), P1(:, :, a + 1))) + l1(a + 1) + l2(2:end)) ./ (a + b);
    [vj, j] = min(v);
    if vj < bv
        bv = vj; bw = [ones(1, a), 2 * ones(1, j)];
    end
end
end


function [bw, bv] = balanced(F, Q)
% best balanced (Christoffel) word of frequency p/q, 2 <= q <= Q, gcd(p, q) = 1
d = size(F{1}, 1); bv = Inf; bw = [];
for q = 2:Q
    ps = find(gcd(1:q - 1, q) == 1);
    j = 0:q - 1;
    letters = floor((j + 1) .* ps(:) / q) - floor(j .* ps(:) / q);   % numel(ps) x q, 1 -> mode 1
    np = numel(ps);
    P = repmat(eye(d), 1, 1, np); lg = zeros(1, np);
    for jj = 1:q
        M = repmat(F{2}, 1, 1, np);
        sel = letters(:, jj) == 1;
        M(:, :, sel) = repmat(F{1}, 1, 1, nnz(sel));
        P = pagemtimes(M, P);
        if mod(jj - 1, 16) == 15 || jj == q
            s = reshape(max(abs(P), [], [1 2]), 1, []); P = P ./ reshape(s, 1, 1, []); lg = lg + log(s);
        end
    end
    v = (batch_logrho(P) + lg) / q;
    [vi, i] = min(v);
    if vi < bv
        bv = vi; bw = 2 - letters(i, :);
    end
end
end


function [val, w] = local_search(F, w, sweeps, kmax)
% single-letter flips, insertions and deletions, best improvement first
m = numel(F); d = size(F{1}, 1);
val = word_value(F, w);
for sw = 1:sweeps
    k = numel(w);
    pre = zeros(d, d, k + 1); lpre = zeros(1, k + 1); pre(:, :, 1) = eye(d);
    for j = 1:k
        M = F{w(j)} * pre(:, :, j); s = max(abs(M(:))); pre(:, :, j + 1) = M / s; lpre(j + 1) = lpre(j) + log(s);
    end
    suf = zeros(d, d, k + 1); lsuf = zeros(1, k + 1); suf(:, :, k + 1) = eye(d);
    for j = k:-1:1
        M = suf(:, :, j + 1) * F{w(j)}; s = max(abs(M(:))); suf(:, :, j) = M / s; lsuf(j) = lsuf(j + 1) + log(s);
    end
    cv = []; cm = zeros(0, 3);                      % moves: [type position letter], 1 flip, 2 insert, 3 delete
    for a = 1:m - 1                                 % flips: letter w(j) -> the a-th other letter
        Rf = zeros(d, d, k); nl = zeros(1, k);
        for j = 1:k
            others = setdiff(1:m, w(j)); nl(j) = others(a); Rf(:, :, j) = F{nl(j)};
        end
        v = (batch_logrho(pagemtimes(suf(:, :, 2:k + 1), pagemtimes(Rf, pre(:, :, 1:k)))) + lsuf(2:k + 1) + lpre(1:k)) / k;
        cv = [cv, v]; cm = [cm; ones(k, 1), (1:k).', nl.']; %#ok<AGROW>
    end
    if k + 1 <= kmax
        for a = 1:m                                 % insertions: new letter a before position j (j = 1..k+1)
            v = (batch_logrho(pagemtimes(suf, pagemtimes(F{a}, pre))) + lsuf + lpre) / (k + 1);
            cv = [cv, v]; cm = [cm; 2 * ones(k + 1, 1), (1:k + 1).', a * ones(k + 1, 1)]; %#ok<AGROW>
        end
    end
    if k > 1                                        % deletions of letter j
        v = (batch_logrho(pagemtimes(suf(:, :, 2:k + 1), pre(:, :, 1:k))) + lsuf(2:k + 1) + lpre(1:k)) / (k - 1);
        cv = [cv, v]; cm = [cm; 3 * ones(k, 1), (1:k).', zeros(k, 1)]; %#ok<AGROW>
    end
    if isempty(cv), break, end
    [bv, b] = min(cv);
    if bv >= val - 1e-15, break, end
    mv = cm(b, :);
    switch mv(1)
        case 1, w(mv(2)) = mv(3);
        case 2, w = [w(1:mv(2) - 1), mv(3), w(mv(2):end)];
        case 3, w(mv(2)) = [];
    end
    val = word_value(F, w);
end
end


function v = getopt(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
