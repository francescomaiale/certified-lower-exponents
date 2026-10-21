function E = exp_upper(A, t, K)
% EXP_UPPER  Entrywise rational upper enclosure of e^{tA} for a Metzler matrix A (exact sym) and rational t > 0.
%   E = exp_upper(A, t, K) with A, t sym (K = 40 Taylor terms by default); the Construction
%   of Sec. 4:
%     N = A + cI >= 0 with c = max(0, max_i -a_ii); X = up(N t / 2^j) with j the least integer such that
%     ||N t||_1 / 2^j <= 1/2 and r = ||X||_1; S = sum_{k <= K} X^k/k! (every term rounded upward), plus the tail bound
%     r^(K+1)/(K+1)! / (1 - r/(K+2)) added to every entry; j squarings (rounded upward); finally multiplied by
%     1/sum_{k <= Kd} (ct)^k/k! >= e^{-ct}.  up() rounds upward to the grid 2^-90.  Every step only increases the
%   entries of a nonnegative matrix, so E >= e^{tA} entrywise.

if nargin < 3, K = 40; end
if ~isscalar(t) || ~isreal(t) || ~isfinite(double(t)) || ~isAlways(t > 0)
    error('exp_upper:time', 't must be finite and positive');
end
if ~isscalar(K) || ~isreal(K) || ~isfinite(K) || K < 1 || K ~= round(K)
    error('exp_upper:terms', 'K must be a positive integer');
end
if isempty(A) || size(A, 1) ~= size(A, 2) || ~isreal(A) || any(~isfinite(double(A(:))))
    error('exp_upper:matrix', 'A must be a finite real square matrix');
end
G = sym(2)^90; up = @(x) ceil(x * G) / G;
d = size(A, 1); I = sym(eye(d));
c = xmax_sym([sym(0); -diag(A)]);
N = A + c * I;
if any(double(sign(N(:))) < 0), error('exp_upper:metzler', 'A is not Metzler'); end
nrm = xmax_sym(sum(N, 1).') * t;                 % max column sum, times t
j = 0;
while double(sign(nrm / sym(2)^j - sym(1) / 2)) > 0, j = j + 1; end
X = up(N * t / sym(2)^j);
r = xmax_sym(sum(X, 1).');                      % ||X||_1 bounds every entry of X^k by r^k
if ~isAlways(r < K + 2), error('exp_upper:tail', 'Taylor tail bound not applicable'); end
S = I; term = I;
for k = 1:K
    term = up(term * X / k);
    S = S + term;
end
tail = up(r^(K + 1) / factorial(sym(K + 1)) / (1 - r / (K + 2)));
E = up(S + tail);
for s = 1:j
    E = up(E * E);
end
x = c * t;
Kd = max(K, ceil(3 * double(x)) + 30);
ps = sym(0); term = sym(1);
for k = 0:Kd
    if k > 0, term = term * x / k; end
    ps = ps + term;
end
damp = up(1 / ps);
E = up(damp * E);
end


function m = xmax_sym(v)
% exact maximum of a vector of symbolic rationals
v = v(:); m = v(1);
for k = 2:numel(v)
    if double(sign(v(k) - m)) > 0, m = v(k); end
end
end
