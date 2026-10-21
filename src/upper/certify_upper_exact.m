function [C, json, info] = certify_upper_exact(w, D, opts)
% CERTIFY_UPPER_EXACT  Inverse-free chain certificate (prop:upper) for the implicit Euler bound of a word.
%   [C, json, info] = certify_upper_exact(w, D, opts)
%     w     word (row vector of modes 1..m, first letter acts first), e.g. from word_search on R_h
%     D     exact data from exact_data
%     opts  name ('pair'), source ('matlab'), U_float_search (float value of the word, stored), meta (struct)
%   Returns the certificate C (decoded JSON struct, in the format read by the replays of reproduce/), its JSON text
%   and info (exact U per reading, printed value, timings).
%
%   Certificate, checked in exact rational arithmetic for each reading of the data:
%     - for every mode i a vector z_i > 0 with (I - h A_i) z_i > 0: then I - h A_i is a nonsingular M-matrix,
%       so R_i = (I - h A_i)^{-1} >= 0 and h eta(A_i) < 1;
%     - a chain x_0 > 0, x_1, ..., x_k with (I - h A_{w_j}) x_j >= x_{j-1}: then x_k >= R_w x_0 and
%       rho(R_w) <= q = max_r (x_k)_r / (x_0)_r (Collatz-Wielandt);
%     - a rational r > 0 with r^k >= q: then lambda_check <= (1 - 1/r)/h <= P, the printed six-decimal value
%       (rounded up; the maximum over the readings).
%   Construction: x_0 is the rationalized (power-iterated, floored) Perron
%   vector of R_w; x_j is (1 + 2^-40) R_{w_j} x_{j-1}, with R_{w_j} the exact inverse, rounded up to a multiple
%   of 2^-100 min(x_{j-1}); z_i is the rationalized (I - h A_i) \ 1; r is a binary64 value of q^(1/k), increased
%   until r^k >= q.

if nargin < 3, opts = struct(); end
name = getopt(opts, 'name', 'pair');
source = getopt(opts, 'source', 'matlab');
Ufs = getopt(opts, 'U_float_search', NaN);
meta = getopt(opts, 'meta', struct());
EPS = sym(2)^-40; BITS = 100;
T0 = tic;
d = D.d; m = D.m; k = numel(w); I = eye(d);
w = w(:).';
if isempty(w) || ~isnumeric(w) || ~isreal(w) || any(~isfinite(w) | w < 1 | w > m | w ~= round(w))
    error('certify_upper_exact:word', 'word must contain integer mode indices 1..m');
end
check_euler_step(D, 'implicit');
% x_0: Perron vector of R_w (floating point, power-iterated), floored and rationalized
Rf = cell(1, m);
for i = 1:m, Rf{i} = inv(I - D.hf * D.Af{i}); end
P = I;
for a = w, P = Rf{a} * P; P = P / max(abs(P(:))); end
[EV, ev] = eig(P); [~, ib] = max(abs(diag(ev)));
x = abs(real(EV(:, ib))); x = x / max(x); x = max(x, 1e-9);
for it = 1:50
    x = P * x; x = x / max(x); x = max(x, 1e-12);
end
x0 = sym(x, 'f');
vers = struct(); Ps = sym([]); Us = sym([]);
for v = 1:numel(D.ver)
    Mq = cell(1, m); Rq = cell(1, m); z = cell(1, m);
    for i = 1:m
        Mq{i} = sym(I) - D.hq * D.ver(v).A{i};
        z{i} = pos_vector(double(Mq{i}), Mq{i});
        Rq{i} = inv(Mq{i});
    end
    X = sym(zeros(d, k + 1)); X(:, 1) = x0; xp = x0;
    for j = 1:k
        a = w(j);
        y = (1 + EPS) * (Rq{a} * xp);
        s = BITS - floor(log2(min(double(xp))));
        xn = ceil(y * sym(2)^s) / sym(2)^s;
        if any(double(sign(Mq{a} * xn - xp)) < 0)
            error('certify_upper_exact:chain', 'chain inequality failed at step %d', j);
        end
        X(:, j + 1) = xn; xp = xn;
    end
    q = xmax(X(:, k + 1) ./ x0);
    r = sym(double(q)^(1 / k) * (1 + 4e-16), 'f');
    while double(sign(r^k - q)) < 0
        r = r * (1 + sym(2)^-50);
    end
    U = (1 - 1 / r) / D.hq;
    P6 = ceil(U * sym(10)^6) / sym(10)^6;
    Ps(v) = P6; Us(v) = U;
    zs = cell(m, 1);
    for i = 1:m, zs{i} = cellstr(string(z{i}).'); end
    Xs = string(X); ch = cell(k + 1, 1);
    for j = 1:k + 1, ch{j} = cellstr(Xs(:, j).'); end
    vers.(D.ver(v).name) = struct('z', {zs}, 'chain', {ch}, 'q', char(q), 'r', char(r), 'U', char(U), ...
                                  'U_float', double(U));
end
Pmax = xmax(Ps);
meta.implementation = 'matlab';
S = struct();
S.row = sprintf('%s_%s', name, D.h_str); S.source = source; S.('case') = name; S.h = char(D.h_str);
if all(w <= 9), S.word = sprintf('%d', w); else, S.word = num2cell(w); end
S.k = k; S.versions = vers;
S.printed = fmt_dec(Pmax, 6); S.printed_exact = char(Pmax);
S.U_float_search = Ufs; S.secs = round(toc(T0), 1); S.meta = meta;
json = jsonencode(S);
C = jsondecode(json);
info = struct('U', Us, 'printed', S.printed, 'printed_exact', Pmax, 'k', k, 'secs', toc(T0));
end


function z = pos_vector(Mf, Mq)
% rational z > 0 with M z > 0 for the Z-matrix M = I - h A (exact check)
d = size(Mf, 1);
cands = {Mf \ ones(d, 1)};
[EV, ev] = eig(inv(Mf)); [~, ib] = max(real(diag(ev)));
pv = abs(real(EV(:, ib))); pv = pv / max(pv);
for fl = [1e-6, 1e-3, 1e-1, 1], cands{end + 1} = max(pv, fl); end %#ok<AGROW>
for c = 1:numel(cands)
    zf = cands{c};
    if ~all(isfinite(zf)) || any(zf <= 0), continue, end
    z = sym(zf / max(zf), 'f');
    if all(double(sign(z)) > 0) && all(double(sign(Mq * z)) > 0)
        return
    end
end
error('certify_upper_exact:mmatrix', 'no positive vector found: I - hA is not a nonsingular M-matrix?');
end


function [m, i] = xmin(x)
% exact minimum of a vector of symbolic rationals: candidates from the double values, then exact comparisons
x = x(:); xd = double(x);
mn = min(xd); cand = find(xd <= mn + 1e-9 * max(1, abs(mn)));
if isempty(cand), cand = (1:numel(x)).'; end
i = cand(1); m = x(i);
for k = cand(2:end).'
    if double(sign(x(k) - m)) < 0, m = x(k); i = k; end
end
end


function [m, i] = xmax(x)
% exact maximum of a vector of symbolic rationals
[m, i] = xmin(-x); m = -m;
end


function v = getopt(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
