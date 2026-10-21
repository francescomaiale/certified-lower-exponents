function [P, info] = orbit_polygon(A, orient, h1, eps_infl, variant)
% ORBIT_POLYGON  Floating-point polygon in the simplex whose cone is invariant under the Euler maps I + h1 A_i
% (prop:cone(iii)), for 3 x 3 families.
%   [P, info] = orbit_polygon(A, orient, h1, eps_infl)
%     A         cell {A_1, A_2} (decimal strings or numeric), orient 'A' or 'T', h1 step (double)
%     eps_infl  relative inflation per round (1e-3 for case 77, 1e-4 for case 677, as in the paper)
%   reach: the hull of the projective orbit of 1, iterating P <- hull(P u B_1 P u B_2 P) until it is stable
%   (same vertex count, sorted coordinates within 1e-12 + 1e-5 |x|, numpy's allclose; at most 4000 rounds).  Then, until every image lies inside with
%   margin 1e-12: inflate the vertices about their centroid by 1 + eps_infl, renormalize, add the images and take
%   the hull again (at most 200 rounds).  P is 3 x n (columns normalized to sum 1), counterclockwise in (x1, x2).
%   variant 'paper' (default) is this order (inflate, add images, hull, test), as in the paper; variant
%   'robust' takes the hull of the polygon and its images first and then inflates and tests, so that a vertex at a
%   fixed direction of some B_i (whose image lies on the boundary) is pushed out too.
%   Only certify_cone's exact check matters.

m = numel(A); d = 3; Bs = cell(1, m);
for i = 1:m
    if isnumeric(A{i}), Ai = A{i}; else, Ai = str2double(A{i}); end
    if orient == 'T', Ai = Ai.'; end
    Bs{i} = eye(d) + h1 * Ai;
end
P = ones(d, 1) / d; conv = false;
for k = 1:4000
    Q = images(P, Bs); Pn = hull_vertices(Q);
    if size(Pn, 2) == size(P, 2) && all(abs(sort(Pn, 2) - sort(P, 2)) <= 1e-12 + 1e-5 * abs(sort(P, 2)), 'all')
        conv = true; break
    end
    P = Pn;
end
info = struct('reach_rounds', k, 'reach_converged', conv, 'reach_vertices', size(P, 2), ...
              'reach_min_coordinate', min(P, [], 'all'));
if nargin < 5, variant = 'paper'; end
for it = 1:200
    if strcmp(variant, 'paper')
        c = mean(P, 2); P = c + (1 + eps_infl) * (P - c); P = P ./ sum(P, 1);
        P = hull_vertices(images(P, Bs));
    else
        P = hull_vertices(images(P, Bs));
        c = mean(P, 2); P = c + (1 + eps_infl) * (P - c); P = P ./ sum(P, 1);
    end
    if margin(P, Bs) > 1e-12, break, end
end
info.variant = variant; info.inflation_rounds = it; info.margin = margin(P, Bs); info.n = size(P, 2); info.min_coordinate = min(P, [], 'all');
end


function Q = images(P, Bs)
Q = P;
for i = 1:numel(Bs), Q = [Q, Bs{i} * P]; end %#ok<AGROW>
Q = Q ./ sum(Q, 1);
end


function H = hull_vertices(Q)
% vertices of the convex hull of the points (columns), counterclockwise in the coordinates (x1, x2)
k = convhull(Q(1, :).', Q(2, :).');
H = Q(:, k(1:end - 1));
end


function w = margin(P, Bs)
c = mean(P, 2); n = size(P, 2); N = zeros(n, 3);
for k = 1:n
    v = cross(P(:, k), P(:, mod(k, n) + 1)).';
    if v * c < 0, v = -v; end
    N(k, :) = v;
end
w = Inf;
for i = 1:numel(Bs)
    Y = Bs{i} * P; Y = Y ./ sum(Y, 1);
    w = min(w, min(N * Y, [], 'all'));
end
end
