function [c, alpha] = chattering_value(A)
% CHATTERING_VALUE  c = min over alpha in [0, 1] of eta(alpha A_1 + (1 - alpha) A_2), an upper bound for lambda_check.
%   [c, alpha] = chattering_value({A1, A2}); eta is the spectral abscissa (the Perron root of a Metzler matrix).
%   Floating point (grid of 401 points, then fminbnd); used only for the upper end 1 + h c of the bisection.
%   For m ~= 2 modes it returns min_i eta(A_i).

eta = @(M) max(real(eig(M)));
if numel(A) ~= 2
    c = Inf; alpha = NaN;
    for i = 1:numel(A), c = min(c, eta(A{i})); end
    return
end
f = @(a) eta(a * A{1} + (1 - a) * A{2});
g = linspace(0, 1, 401); v = arrayfun(f, g);
[c, j] = min(v); alpha = g(j);
lo = g(max(j - 1, 1)); hi = g(min(j + 1, numel(g)));
if hi > lo
    [a2, c2] = fminbnd(f, lo, hi, optimset('TolX', 1e-13));
    if c2 < c, c = c2; alpha = a2; end
end
end
