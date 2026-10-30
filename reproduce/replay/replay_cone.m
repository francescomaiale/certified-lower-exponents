function [ok, msgs, info] = replay_cone(C, A)
% REPLAY_CONE  Independent exact re-check of an invariant-cone certificate (certify_cone; Computation 3.9).
%   [ok, msgs, info] = replay_cone(C, A)
%     C  certificate (.json file, JSON text or decoded struct); A data (optional; default from data/*.json)
%   Shares no code with certify_cone (port of the 'replay' subcommands of euler_cone.py and euler_cone_boundary.py):
%   the binary64 data are re-read, B_i = I + h1 A_i (transposed for orientation 'T'), the generators are projected
%   to the simplex coordinates (x1, x2), and it checks exactly that the stored polygon is convex (all turns of one
%   sign), that every projected image B_i v lies in a triangle of the fan from the first vertex (barycentric
%   coordinates >= 0 with sum <= 1), that the generators are >= 0 (> 0 for kind 'polygon'), and that the stored
%   support classes are those of the cone: J is a class iff the generators with support inside J cover J.

T0 = tic;
if ischar(C) || isstring(C)
    if isfile(C), C = jsondecode(fileread(C)); else, C = jsondecode(char(C)); end
end
if nargin < 2 || isempty(A), A = replay_data(C.('case')); end
msgs = {}; m = numel(A); %#ok<*AGROW>  diagnostics and support classes are short variable-length lists
h1 = rat_parse(string(C.h1));
info = struct();
if xsgn(h1) <= 0 || ~ismember(string(C.orientation), ["A", "N", "T"])
    ok = false; msgs{end + 1} = 'positive step and valid orientation required'; return
end
vs = C.vertices; if ~iscell(vs), vs = num2cell(vs, 2); end
n = numel(vs); V = sym(zeros(3, n));
for k = 1:n, V(:, k) = rat_parse(string(vs{k}(:)).'); end
if n < 3 || any(xsgn(V(:)) < 0) || any(xsgn(sum(V, 1)) <= 0)
    ok = false; msgs{end + 1} = 'invalid cone generators'; return
end
B = cell(1, m);
for i = 1:m
    if isnumeric(A{i}), Ai = sym(A{i}, 'f'); else, Ai = sym(str2double(string(A{i})), 'f'); end
    if ~isequal(size(Ai), [3 3]) || any(xsgn(Ai(~eye(3))) < 0)
        ok = false; msgs{end + 1} = 'need 3-dimensional Metzler data'; return
    end
    if strcmp(C.orientation, 'T'), Ai = Ai.'; end
    B{i} = sym(eye(3)) + h1 * Ai;
end
Q = [V(1, :) ./ sum(V, 1); V(2, :) ./ sum(V, 1)];
turn = sym(zeros(1, n));
for k = 1:n
    p = Q(:, k); q = Q(:, mod(k, n) + 1); r = Q(:, mod(k + 1, n) + 1);
    turn(k) = (q(1) - p(1)) * (r(2) - q(2)) - (q(2) - p(2)) * (r(1) - q(1));
end
ts = xsgn(turn); convex = all(ts >= 0) || all(ts <= 0);
if ~convex, msgs{end + 1} = 'polygon not convex'; end
nonneg = all(xsgn(V(:)) >= 0); pos = all(xsgn(V(:)) > 0);
if ~nonneg, msgs{end + 1} = 'a generator is not >= 0'; end
if strcmp(C.kind, 'polygon') && ~pos, msgs{end + 1} = 'a generator of a polygon cone is not > 0'; end
nout = 0;
for i = 1:m
    Y = B{i} * V;
    if any(xsgn(Y(:)) < 0) || any(xsgn(sum(Y, 1)) <= 0)
        ok = false; msgs{end + 1} = 'an image is not nonzero and nonnegative'; return
    end
    Y = [Y(1, :) ./ sum(Y, 1); Y(2, :) ./ sum(Y, 1)];
    for j = 1:n
        if ~in_fan(Q, Y(:, j)), nout = nout + 1; end
    end
end
if nout > 0, msgs{end + 1} = sprintf('%d images outside the polygon', nout); end
% support classes
Z = xsgn(V) ~= 0; cls = strings(0);
for J = 1:7
    Jb = logical(bitget(J, 1:3)).';
    inside = all(~Z(~Jb, :), 1);
    if any(inside) && isequal(any(Z(:, inside), 2), Jb)
        cls(end + 1) = strjoin(string(find(Jb).'), ',');
    end
end
stored = sort(string(C.classes(:)).');
if ~isequal(sort(cls), stored), msgs{end + 1} = 'stored support classes differ from the recomputed ones'; end
ok = isempty(msgs);
info = struct('n', n, 'images', m * n, 'classes', {cellstr(sort(cls))}, 'secs', toc(T0));
end


function in = in_fan(Q, p)
in = false; n = size(Q, 2); a = Q(:, 1);
for k = 2:n - 1
    b = Q(:, k); c = Q(:, k + 1);
    dt = (b(1) - a(1)) * (c(2) - a(2)) - (c(1) - a(1)) * (b(2) - a(2));
    if xsgn(dt) == 0, continue, end
    l1 = ((p(1) - a(1)) * (c(2) - a(2)) - (c(1) - a(1)) * (p(2) - a(2))) / dt;
    l2 = ((b(1) - a(1)) * (p(2) - a(2)) - (p(1) - a(1)) * (b(2) - a(2))) / dt;
    if xsgn(l1) >= 0 && xsgn(l2) >= 0 && xsgn(1 - l1 - l2) >= 0, in = true; return, end
end
end
