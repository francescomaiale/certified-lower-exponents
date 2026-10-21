function L = log_upper(r, n)
% LOG_UPPER  Rational upper bound for log(r), r > 0 a sym rational.
%   r = 2^b y with y in [1, 2); log y = 2 artanh(z), z = (y - 1)/(y + 1) in [0, 1/3), bounded above by the first n
%   terms of the series plus the remainder 2 z^(2n+1) / ((2n+1)(1 - z^2)); log 2 = 2 artanh(1/3) is bounded below
%   or above by the same series as the sign of b requires.
if nargin < 2, n = 60; end
if ~isscalar(r) || ~isreal(r) || ~isfinite(double(r)) || ~isAlways(r > 0)
    error('log_upper:input', 'r must be finite and positive');
end
if ~isscalar(n) || ~isreal(n) || ~isfinite(n) || n < 1 || n ~= round(n)
    error('log_upper:terms', 'n must be a positive integer');
end
b = 0;
while double(sign(r - 2)) >= 0, r = r / 2; b = b + 1; end
while double(sign(r - 1)) < 0, r = r * 2; b = b - 1; end
[~, ly] = artanh_bounds((r - 1) / (r + 1), n);
[l2lo, l2hi] = artanh_bounds(sym(1) / 3, n);
if b >= 0, L = ly + b * l2hi; else, L = ly + b * l2lo; end
end


function [lo, hi] = artanh_bounds(z, n)
% 2 artanh(z) for 0 <= z < 1: partial sum (lower bound) and partial sum plus remainder (upper bound)
s = sym(0); zk = z; z2 = z * z;
for i = 0:n - 1
    s = s + zk / (2 * i + 1);
    zk = zk * z2;
end
lo = 2 * s;
hi = lo + 2 * zk / ((2 * n + 1) * (1 - z2));
end
