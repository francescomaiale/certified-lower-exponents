function [m, i] = exact_min(x)
% EXACT_MIN  Exact minimum of a vector of symbolic rationals, and its (first) index.
%   Candidates are the entries whose double values lie within 1e-9 of the smallest double value (double() is
%   accurate to about 1e-16, so the exact minimum is among them); they are reduced by a vectorized knockout
%   tournament of exact comparisons.  Ties go to the smaller index.
x = x(:); xd = double(x);
mn = min(xd);
cand = find(xd <= mn + 1e-9 * max(1, abs(mn)));
if isempty(cand), cand = (1:numel(x)).'; end
v = x(cand); id = cand;
while numel(v) > 1
    n2 = 2 * floor(numel(v) / 2);
    a = v(1:2:n2); b = v(2:2:n2);
    keepa = double(sign(a - b)) <= 0;
    ia = id(1:2:n2); ib = id(2:2:n2);
    w = b; w(keepa) = a(keepa); iw = ib; iw(keepa) = ia(keepa);
    if n2 < numel(v), w = [w; v(end)]; iw = [iw; id(end)]; end %#ok<AGROW>
    v = w; id = iw;
end
m = v; i = id;
end
