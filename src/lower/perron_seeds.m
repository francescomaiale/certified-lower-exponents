function S = perron_seeds(F, w)
% PERRON_SEEDS  Seeds of the closure from a word (Obs. 4.8): the Perron vector of F_w and its images along w.
%   S = perron_seeds(F, w) for a family F (cell of nonnegative matrices) and a word w (row of modes, the first
%   letter acts first) returns d x |w| columns: v (the Perron vector of F_w = F_{w_k} ... F_{w_1}, normalized to sum
%   1), then F_{w_1} v, F_{w_2} F_{w_1} v, ..., up to the image under the first |w| - 1 letters, each normalized.
d = size(F{1}, 1); P = eye(d);
for a = w, P = F{a} * P; P = P / max(abs(P(:))); end
[EV, ev] = eig(P); [~, k] = max(abs(diag(ev)));
v = abs(real(EV(:, k))); v = v / sum(v);
S = zeros(d, numel(w)); S(:, 1) = v;
for j = 1:numel(w) - 1
    v = F{w(j)} * v; v = max(v, 0); v = v / sum(v); S(:, j + 1) = v;
end
end
