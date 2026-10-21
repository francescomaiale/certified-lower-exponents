function v = word_value(F, w)
% WORD_VALUE  log rho(F_w) / |w| for a word w (row of modes; the first letter acts first), floating point.
P = eye(size(F{1}, 1)); lg = 0;
for a = w
    P = F{a} * P; s = max(abs(P(:))); P = P / s; lg = lg + log(s);
end
v = (log(max(abs(eig(P)))) + lg) / numel(w);
end
