function [bw, bv] = brute_words(F, L)
% BRUTE_WORDS  Best word of every length 1..L by exhaustive enumeration (floating point).
%   [bw, bv] = brute_words(F, L): bw{k} is a word of length k (row of modes 1..m, the first letter acts first)
%   minimizing log rho(F_w) / k, and bv(k) that value.
m = numel(F); d = size(F{1}, 1);
P = cat(3, F{:}); W = (1:m).'; lg = zeros(1, m);
bw = cell(1, L); bv = zeros(1, L);
for k = 1:L
    v = (batch_logrho(P) + lg) / k;
    [bv(k), j] = min(v); bw{k} = W(j, :);
    if k == L, break, end
    N = size(P, 3);
    P2 = zeros(d, d, m * N); W2 = zeros(m * N, k + 1);
    for a = 1:m
        P2(:, :, (a - 1) * N + (1:N)) = pagemtimes(F{a}, P);
        W2((a - 1) * N + (1:N), :) = [W, a * ones(N, 1)];
    end
    lg = repmat(lg, 1, m);
    s = reshape(max(abs(P2), [], [1 2]), 1, []);
    P = P2 ./ reshape(s, 1, 1, []); lg = lg + log(s); W = W2;
end
end


function l = batch_logrho(P)
l = reshape(log(max(abs(pageeig(P)), [], 1)), 1, []);
end
