function r = limit_denominator(x, N)
% LIMIT_DENOMINATOR  Closest rational p/q to the exact value of the double x with q <= N, as a sym.
%   r = limit_denominator(x, N), as Python's Fraction(x).limit_denominator(N) (continued fractions on the exact
%   binary64 value of x, then the better of the last convergent and the best semiconvergent).
X = sym(x, 'f'); N = sym(N);
p0 = sym(0); q0 = sym(1); p1 = sym(1); q1 = sym(0);
[n, d] = numden(X);
if double(sign(d - N)) <= 0, r = X; return, end
while true
    a = floor(n / d);
    q2 = q0 + a * q1;
    if double(sign(q2 - N)) > 0, break, end
    [p0, q0, p1, q1] = deal(p1, q1, p0 + a * p1, q2);
    [n, d] = deal(d, n - a * d);
    if double(d) == 0, break, end
end
k = floor((N - q0) / q1);
b1 = (p0 + k * p1) / (q0 + k * q1); b2 = p1 / q1;
if double(sign(abs(b2 - X) - abs(b1 - X))) <= 0, r = b2; else, r = b1; end
end
