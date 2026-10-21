function [ok, msgs, neg, info] = replay_lower(C, A)
% REPLAY_LOWER  Independent exact re-check of a polytope certificate (prop:polytope), with negative controls.
%   [ok, msgs, neg, info] = replay_lower(C, A)
%     C  certificate: a .json file (MATLAB or Python, the format of write_cert in certify_lower.py), JSON text,
%        or the decoded struct returned by certify_lower_exact
%     A  the data: cell array {A_1, ..., A_m} of numeric matrices or of decimal-string matrices (read_cases_json)
%   ok is true iff every check passes; msgs lists failed checks, neg the negative controls, info the exact t*.
%
%   Shares no code with the certify path (only MATLAB builtins): the data are re-read in both readings (exact
%   decimals and exact binary64 values), h is the exact decimal C.h, B_i = I + h A_i is rebuilt (transposed for
%   orientation 'T'), and in exact rational arithmetic it checks
%     - h |a_rr| < 1 and B_i >= 0; every vertex is >= 0 and nonzero; t > 0;
%     - the weights cover every (mode, vertex) pair exactly once, are >= 0 and sum to 1;
%     - B_i v - t V lambda >= 0 for every pair and both readings;
%     - t* recomputed equals the stored t_star, and t <= t*; (t - 1)/h >= the claimed lambda_lower;
%   and the negative controls (port of replay_one in certify_lower.py):
%     - the stored weights of the binding pair fail at t*(1 + 1e-12);
%     - the dual c: c >= 0, c'u >= 1 for every vertex u, tau_up = c'(B_i v) equals the stored value and is >= t*,
%       so the polytope is not t-invariant for t > tau_up.

T0 = tic;
if ischar(C) || isstring(C)
    if isfile(C), C = jsondecode(fileread(C)); else, C = jsondecode(char(C)); end
end
if nargin < 2 || isempty(A), A = replay_data(C.('case')); end   % data re-read from data/*.json
msgs = {}; neg = {}; ok = true;
% ---- parse every rational of the certificate at once (rat_parse, exact)
verts = C.vertices;
if ~iscell(verts{1}) && ~isstring(verts{1}), verts = {verts}; end
n = numel(verts); d = numel(verts{1});
vstr = strings(d, n);
for k = 1:n, vstr(:, k) = string(verts{k}(:)); end
wts = C.weights;
if iscell(wts), wts = [wts{:}]; end
nw = numel(wts); wm = zeros(1, nw); wv = zeros(1, nw); kk = cell(1, nw); ws = cell(1, nw);
for e = 1:nw
    wm(e) = wts(e).mode + 1; wv(e) = wts(e).vertex + 1;
    f = fieldnames(wts(e).lambda);
    kk{e} = cellfun(@(s) str2double(regexprep(s, '^x', '')), f).' + 1;
    ws{e} = string(struct2cell(wts(e).lambda)).';
end
cnull = isempty(C.negative_control.c) || (isnumeric(C.negative_control.c) && all(isnan(C.negative_control.c)));
if cnull, cs = strings(1, 0); else, cs = string(C.negative_control.c(:)).'; end
allw = [ws{:}];
Q = rat_parse([vstr(:).', allw, string(C.t), string(C.t_star), string(C.negative_control.tau_up), ...
               string(C.lambda_lower), string(C.h), cs]);
o = 0;
V = reshape(Q(o + 1:o + d * n), d, n); o = o + d * n;
Wall = Q(o + 1:o + numel(allw)); o = o + numel(allw);
t = Q(o + 1); t_star_stored = Q(o + 2); tau_stored = Q(o + 3); lam_claim = Q(o + 4); hq = Q(o + 5); o = o + 5;
cq = Q(o + 1:end).';
orient = C.orientation;
info = struct();
if sgn(hq) <= 0 || ~ismember(string(orient), ["N", "T"])
    ok = false; msgs{end + 1} = 'positive step and orientation N or T required'; return
end
% ---- data, both readings
raw = A; m = numel(raw);
dstr = cell(1, m);
for i = 1:m
    if isnumeric(raw{i}), dstr{i} = short_repr(raw{i}); else, dstr{i} = string(raw{i}); end
end
readings = {'decimal', 'binary64'};
% ---- structural checks
cnt = zeros(m, n);
if any(~isfinite(wm) | ~isfinite(wv) | wm ~= round(wm) | wv ~= round(wv) | wm < 1 | wm > m | wv < 1 | wv > n)
    ok = false; msgs{end + 1} = 'weights refer to a mode or vertex that does not exist';
    return
else
    for e = 1:nw, cnt(wm(e), wv(e)) = cnt(wm(e), wv(e)) + 1; end
end
if any(cnt(:) ~= 1)
    ok = false; msgs{end + 1} = 'weights do not cover every (mode, vertex) pair exactly once';
end
if ~(sgn(t) > 0), ok = false; msgs{end + 1} = 't <= 0'; end
if any(sgn(V(:)) < 0) || any(sgn(sum(V, 1)) <= 0)
    ok = false; msgs{end + 1} = 'vertex not in K \ {0}';
end
% weights per mode as padded (smax x n) index / value arrays
Kp = cell(1, m); Wp = cell(1, m); pos = 0;
for i = 1:m
    sm = 1;
    for e = find(wm == i), sm = max(sm, numel(kk{e})); end
    Kp{i} = ones(sm, n); Wp{i} = sym(zeros(sm, n));
end
for e = 1:nw
    s = numel(kk{e}); i = wm(e); j = wv(e);
    Kp{i}(1:s, j) = kk{e}(:); Wp{i}(1:s, j) = Wall(pos + 1:pos + s).'; pos = pos + s;
end
if any(sgn(Wall) < 0) || any(cellfun(@(x) any(~isfinite(x(:)) | x(:) ~= round(x(:)) | x(:) < 1 | x(:) > n), Kp))
    ok = false; msgs{end + 1} = 'negative weight or weight on a nonexistent vertex';
    return
end
VL = cell(1, m);
for i = 1:m
    if any(sgn(sum(Wp{i}, 1) - 1) ~= 0)
        ok = false; msgs{end + 1} = sprintf('weights of mode %d do not sum to 1', i - 1);
    end
    VL{i} = sym(zeros(d, n));
    for p = 1:size(Kp{i}, 1)
        VL{i} = VL{i} + V(:, min(max(Kp{i}(p, :), 1), n)) .* Wp{i}(p, :);
    end
end
% ---- the certificate inequality, both readings
tstar = []; Bs = struct();
for r = 1:numel(readings)
    rd = readings{r};
    for i = 1:m
        if strcmp(rd, 'decimal')
            Ai = reshape(rat_parse(dstr{i}(:).'), d, d);
        else
            Ai = sym(str2double(dstr{i}), 'f');
        end
        if any(sgn(hq * abs(diag(Ai)) - 1) >= 0)
            ok = false; msgs{end + 1} = sprintf('%s: h >= h_0', rd); %#ok<*AGROW>
        end
        B = sym(eye(d)) + hq * Ai;
        if orient == 'T', B = B.'; end
        if any(sgn(B(:)) < 0)
            ok = false; msgs{end + 1} = sprintf('%s: B_%d not nonnegative', rd, i - 1);
        end
        Bs.(rd){i} = B;
        BV = B * V;
        Dm = BV - t * VL{i};
        [rr, jj] = find(sgn(Dm) < 0, 1);
        if ~isempty(rr)
            ok = false; msgs{end + 1} = sprintf('%s: (*) fails at mode %d, vertex %d, coordinate %d', rd, i - 1, jj - 1, rr - 1);
        end
        P = sgn(VL{i}) > 0;
        if any(P(:))
            q = xmin(BV(P) ./ VL{i}(P));
            if isempty(tstar) || sgn(q - tstar) < 0, tstar = q; end
        end
    end
end
if isempty(tstar) || sgn(tstar - t_star_stored) ~= 0 || sgn(t - tstar) > 0
    ok = false;
    msgs{end + 1} = sprintf('recomputed t* = %.15g differs from stored %.15g or is < t', double(tstar), double(t_star_stored));
end
if sgn((t - 1) / hq - lam_claim) < 0
    ok = false; msgs{end + 1} = 'claimed lambda bound not implied by t';
end
% ---- negative controls
b = C.binding; i = b.mode + 1; j = b.vertex + 1;
if i < 1 || i > m || j < 1 || j > n || ~isfield(Bs, b.data)
    ok = false; msgs{end + 1} = 'binding pair does not exist';
else
    B = Bs.(b.data){i};
    Bv = B * V(:, j);
    t_raise = t_star_stored * (1 + sym(1) / sym(10)^12);
    fails = any(sgn(Bv - t_raise * VL{i}(:, j)) < 0);
    if fails, s = 'fail (expected)'; else, s = 'PASS (unexpected)'; end
    neg{end + 1} = sprintf('weights at t*(1+1e-12): %s', s);
    ok = ok && fails;
    if ~cnull
        dual_ok = numel(cq) == d && all(sgn(cq) >= 0) && all(sgn(cq.' * V - 1) >= 0);
        if dual_ok
            tau = cq.' * Bv;
            dual_ok = sgn(tau - tau_stored) == 0 && sgn(tau - t_star_stored) >= 0;
        else
            tau = tau_stored;
        end
        if dual_ok, s = 'verified'; else, s = 'FAILED'; end
        neg{end + 1} = sprintf('dual: polytope not t-invariant for t > tau_up = %.15f (tau_up/t* - 1 = %.1e) %s', ...
                               double(tau), double(tau / t_star_stored - 1), s);
        ok = ok && dual_ok;
        nxt = 1 + hq * (lam_claim + sym(1) / sym(10)^6);
        if sgn(nxt - tau) > 0, s = 'excluded for this polytope'; else, s = 'not excluded (slack between t and t*)'; end
        neg{end + 1} = sprintf('next six-decimal value: %s', s);
    else
        neg{end + 1} = 'no dual stored (trivial certificate): t* is the exact binding ratio, any t > t* fails there';
    end
end
info = struct('n', n, 'd', d, 't', t, 't_star', tstar, 'lambda_lower', C.lambda_lower, 'secs', toc(T0));
end


function [m, i] = xmin(x)
% exact minimum (vectorized comparisons, replay/private/xmin_exact)
[m, i] = xmin_exact(x);
end


function s = sgn(x)
% exact sign of symbolic rationals, as doubles
s = double(sign(x));
end


function s = short_repr(x)
s = strings(size(x));
for k = 1:numel(x)
    for p = 1:17
        t = sprintf('%.*g', p, x(k));
        if str2double(t) == x(k), break, end
    end
    s(k) = string(t);
end
end
