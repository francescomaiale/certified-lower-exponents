function [C, info] = certify_cone(V, A, orient, h1, opts)
% CERTIFY_CONE  Exact check that a polyhedral cone in R^3_+ is invariant under the Euler maps I + h1 A_i
% (prop:cone, cor:eulerface; Computation 3.9).
%   [C, info] = certify_cone(V, A, orient, h1, opts)
%     V       3 x n generators of the cone (sym, or doubles taken as exact binary64 values), any order
%     A       cell {A_1, ..., A_m} (decimal strings or numeric); the binary64 values are used, as in the paper
%     orient  'A' (family as given) or 'T' (transposed family)
%     h1      step, sym or decimal string (e.g. '1/5' or '0.2')
%     opts    name ('pair'), kind ('polygon' or 'cut cone'), require_positive (false)
%   The generators are ordered counterclockwise around their centroid (floating point; only the order is
%   used).  For consecutive generators the facet normal n_k = v_k x v_{k+1}, oriented so that n_k . v >= 0 at the
%   centroid, defines cone(V) = {x : n_k . x >= 0 for all k}.  Checked exactly: every n_k . v_j >= 0 (convexity,
%   so every v_j is a generator of this cone), every n_k . (B_i v_j) >= 0 (B_i cone(V) in cone(V)), V >= 0, and,
%   if require_positive, V > 0 (then C \ {0} lies in int K, prop:cone(ii)).
%   The support classes Sigma(C) = {supp x : x in C \ {0}} are read off exactly from the zero patterns of the
%   generators (a point of C has x_k = 0 iff it is a nonnegative combination of generators with x_k = 0).
%   replay_cone re-checks the certificate by a different method (fan triangulation, barycentric coordinates).

if nargin < 5, opts = struct(); end
if ~isa(h1, 'sym'), h1 = parse_rationals(string(h1)); end
if ~isa(V, 'sym'), V = sym(V, 'f'); end
if ~isscalar(h1) || ~isAlways(h1 > 0) || ~isfinite(double(h1)) || ~ismember(string(orient), ["A", "T"])
    error('certify_cone:input', 'positive finite step and orientation A or T required');
end
if size(V, 1) ~= 3 || size(V, 2) < 3 || any(~isfinite(double(V(:)))) || ...
        any(double(sign(V(:))) < 0) || any(double(sign(sum(V, 1))) <= 0) || rank(V) < 3
    error('certify_cone:generators', 'need nonzero nonnegative generators spanning R^3');
end
kind = getopt(opts, 'kind', 'polygon');
m = numel(A); d = 3;
Bq = cell(1, m);
for i = 1:m
    if isnumeric(A{i}), Af = A{i}; else, Af = str2double(A{i}); end
    if ~isequal(size(Af), [3 3]) || ~isreal(Af) || any(~isfinite(Af(:))) || any(Af(~eye(3)) < 0)
        error('certify_cone:matrix', 'need finite 3-dimensional Metzler data');
    end
    Ai = sym(Af, 'f'); if orient == 'T', Ai = Ai.'; end
    Bq{i} = sym(eye(d)) + h1 * Ai;
end
% counterclockwise order of the projected generators around their centroid
P = double(V ./ sum(V, 1)); cen = mean(P, 2);
[~, o] = sort(atan2(P(2, :) - cen(2), P(1, :) - cen(1)));
V = V(:, o); n = size(V, 2);
cq = sum(V ./ sum(V, 1), 2) / n;                       % exact centroid of the normalized generators
N = sym(zeros(n, 3));
for k = 1:n
    a = V(:, k); b = V(:, mod(k, n) + 1);
    nk = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)];
    if double(sign(nk * cq)) < 0, nk = -nk; end
    N(k, :) = nk;
end
convex = all(double(sign(N * V)) >= 0, 'all') && all(any(double(sign(N)) ~= 0, 2));
outside = 0;
for i = 1:m
    outside = outside + nnz(any(double(sign(N * (Bq{i} * V))) < 0, 1));
end
nonneg = all(double(sign(V)) >= 0, 'all');
pos = all(double(sign(V)) > 0, 'all');
ok = m > 0 && convex && outside == 0 && nonneg && ...
     (~(getopt(opts, 'require_positive', false) || strcmp(kind, 'polygon')) || pos);
% support classes: unions of the supports of the generators within each zero pattern (supports as bit masks)
supp = double(sign(V)) ~= 0; masks = [];
for Z = {[], 1, 2, 3, [1 2], [1 3], [2 3]}
    W = find(all(~supp(Z{1}, :), 1));
    U = [];
    for w = W
        s = sum(2 .^ (find(supp(:, w)) - 1));
        U = unique([U, s, bitor(U, s)]);
    end
    masks = unique([masks, U]);
end
classes = arrayfun(@(b) find(bitget(b, 1:3)), masks, 'UniformOutput', false);
[~, o] = sort(cellfun(@(c) numel(c) * 1000 + polyval(c, 10), classes)); classes = classes(o);
classes = cellfun(@(c) strjoin(string(c), ','), classes, 'UniformOutput', false);   % e.g. '1,3'
Pn = double(V ./ sum(V, 1));
C = struct('case', getopt(opts, 'name', 'pair'), 'kind', kind, 'orientation', orient, ...
           'h1', char(h1), 'ok', ok, 'n_vertices', n, 'classes', {classes}, ...
           'min_coordinates', min(Pn, [], 2).', 'vertices', {arrayfun(@(k) cellstr(string(V(:, k))).', 1:n, ...
           'UniformOutput', false)});
info = struct('convex', convex, 'images_outside', outside, 'nonneg', nonneg, 'positive', pos);
end
