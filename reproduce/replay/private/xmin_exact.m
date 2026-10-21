function [m, i] = xmin_exact(x)
% XMIN_EXACT  (replay side) exact minimum of a vector of symbolic rationals and an index attaining it.
%   Entries whose double value is within 1e-9 (relative) of the smallest are compared exactly, in rounds of
%   vectorized pairwise comparisons.
x = x(:); xd = double(x); lo = min(xd);
idx = find(xd - lo <= 1e-9 * max(abs(lo), 1));
while numel(idx) > 1
    h = floor(numel(idx) / 2);
    p = idx(1:h); q = idx(h + 1:2 * h);
    s = double(sign(x(q) - x(p)));            % q replaces p where x(q) < x(p)
    p(s < 0) = q(s < 0);
    idx = [p; idx(2 * h + 1:end)];
end
i = idx; m = x(i);
end
