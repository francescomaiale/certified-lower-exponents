function [m, i] = exact_max(x)
% EXACT_MAX  Exact maximum of a vector of symbolic rationals, and its index (see exact_min).
[m, i] = exact_min(-x); m = -m;
end
