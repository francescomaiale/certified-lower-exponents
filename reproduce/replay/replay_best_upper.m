function [ok, msgs, info] = replay_best_upper(R, A)
% REPLAY_BEST_UPPER  Independent exact re-check of a chattering / schedule certificate (certify_best_upper).
%   [ok, msgs, info] = replay_best_upper(R, A)
%     R  certificate: .json file, JSON text or decoded struct (MATLAB or Python certify_best_upper.py format)
%     A  data {A_1, ..., A_m} (optional; default: re-read from data/*.json by the name R.id)
%   Shares no code with the certify path and takes no logarithm (port of replay() in certify_best_upper.py):
%     - chattering: for both data readings, x > 0 and (S_alpha x)_r <= U_c x_r for every r, cross-multiplied, with
%       S_alpha = alpha A_1 + (1 - alpha) A_2; then c <= U_c by Collatz-Wielandt;
%     - schedule: each e^{tA} is enclosed from above by uniformization, e^{tA} = e^{-qt} sum_k (qt)^k/k! P^k with
%       P = I + A/q >= 0, applied to t/2^s with the tail bounded by the row-sum norm of P, then squared s times
%       (on its own dyadic grid 2^-80); Phi is the product of the enclosures, rho = max_i (Phi x)_i / x_i, and
%       rho <= exp_lower(U_s T), a rational lower bound for e^{U_s T}, is checked; then
%       log rho(Phi) / T <= U_s.
%     - U_best = min(U_c, U_s) and the printed value is >= U_best.

T0 = tic;
if ischar(R) || isstring(R)
    if isfile(R), R = jsondecode(fileread(R)); else, R = jsondecode(char(R)); end
end
if nargin < 2 || isempty(A), A = replay_data(R.id); end
msgs = {}; ok = true; m = numel(A); d = size(A{1}, 1);
names = {'decimal', 'binary64'};
has_c = isfield(R, 'U_c') && ~isempty(R.U_c);
has_s = isfield(R, 'U_s') && ~isempty(R.U_s);
if has_c, Uc = rat_parse(string(R.U_c)); end
if has_s, Us = rat_parse(string(R.U_s)); end
rho_max = [];
for v = 1:2
    if ~isfield(R.readings, names{v})
        ok = false; msgs{end + 1} = sprintf('%s: no certificate for this reading', names{v}); continue %#ok<*AGROW>
    end
    rd = R.readings.(names{v});
    As = cell(1, m);
    for i = 1:m
        if isnumeric(A{i}), s = shortest_decimal(A{i});
        else, s = string(A{i}); end
        if v == 1, As{i} = reshape(rat_parse(s(:).'), d, d);
        else, As{i} = sym(str2double(s), 'f'); end
        if any(xsgn(As{i}(~eye(d))) < 0)
            ok = false; msgs{end + 1} = 'data must be Metzler';
            info = struct(); return
        end
    end
    if has_c
        a = rat_parse(string(rd.alpha)); x = rat_parse(string(rd.x(:)).').';
        if m ~= 2 || xsgn(a) < 0 || xsgn(a - 1) > 0 || numel(x) ~= d
            ok = false; msgs{end + 1} = sprintf('%s: invalid mixture or witness', names{v}); continue
        end
        S = a * As{1} + (1 - a) * As{2};
        if any(xsgn(x) <= 0) || any(xsgn(S * x - Uc * x) > 0)
            ok = false; msgs{end + 1} = sprintf('%s: Collatz-Wielandt bound for c fails', names{v});
        end
    end
    if has_s
        if isfield(rd, 'modes') || isfield(rd, 'durations')
            modes = rd.modes(:).' + 1; ts = rat_parse(string(rd.durations(:)).');
        else
            modes = [1 2]; ts = rat_parse([string(rd.t1), string(rd.t2)]);
        end
        xs = rat_parse(string(rd.x_sched(:)).').';
        if isempty(modes) || numel(modes) ~= numel(ts) || any(xsgn(ts) <= 0) || ...
                any(~isfinite(modes) | modes ~= round(modes) | modes < 1 | modes > m) || ...
                numel(xs) ~= d || any(xsgn(xs) <= 0)
            ok = false; msgs{end + 1} = sprintf('%s: invalid schedule', names{v}); continue
        end
        Phi = sym(eye(d));
        for k = 1:numel(modes)
            Phi = unif_upper(As{modes(k)}, ts(k)) * Phi;
        end
        rho = xmax_exact((Phi * xs) ./ xs);
        T = sum(ts);
        if any(xsgn(xs) <= 0) || xsgn(rho - exp_lower(Us * T)) > 0
            ok = false; msgs{end + 1} = sprintf('%s: schedule bound fails (rho > exp(U_s T))', names{v});
        end
        rho_max = [rho_max, double(rho)];
    end
end
best = [];
if has_c, best = Uc; end
if has_s && (isempty(best) || xsgn(Us - best) < 0), best = Us; end
if isempty(best)
    ok = false; msgs{end + 1} = 'no bound in the certificate';
else
    if xsgn(rat_parse(string(R.U_best)) - best) ~= 0
        ok = false; msgs{end + 1} = 'U_best is not min(U_c, U_s)';
    end
    if isfield(R, 'printed') && xsgn(rat_parse(string(R.printed)) - best) < 0
        ok = false; msgs{end + 1} = 'printed value below U_best';
    end
end
info = struct('U_best', double(best), 'rho', rho_max, 'secs', toc(T0));
end


function E = unif_upper(A, t)
% e^{tA} <= e^{-qt} [sum_{k<=K} (qt)^k/k! P^k + tail], P = I + A/q, for t/2^s with q t ||P||_inf / 2^s <= 1,
% then squared s times (nonnegative matrices: monotone); grid 2^-80, K = 200
K = 200; G = sym(2)^80; upg = @(z) ceil(z * G) / G;
d = size(A, 1); I = sym(eye(d));
q = sym(1);
for i = 1:d, if xsgn(-A(i, i) - q) > 0, q = -A(i, i); end, end
P = I + A / q;
if any(xsgn(P(:)) < 0), error('replay_best_upper:metzler', 'A is not Metzler'); end
rs = sum(P, 2); nP = rs(1);
for i = 2:d, if xsgn(rs(i) - nP) > 0, nP = rs(i); end, end
s = 0;
while xsgn(q * t * nP / sym(2)^s - 1) > 0, s = s + 1; end
tt = t / sym(2)^s; qt = q * tt;
S = sym(zeros(d)); Pk = I; coef = sym(1); esum = sym(0);
for k = 0:K
    if k > 0
        Pk = upg(Pk * P); coef = coef * qt / k;
    end
    S = S + upg(coef * Pk); esum = esum + coef;
end
z = qt * nP;
if xsgn(z / (K + 2) - 1) >= 0, error('replay_best_upper:tail', 'tail bound not applicable'); end
tail = upg(z^(K + 1) / factorial(sym(K + 1)) / (1 - z / (K + 2)));
E = upg(upg(1 / esum) * (S + tail));
for k = 1:s
    E = upg(E * E);
end
end


function L = exp_lower(y)
% rational lower bound for e^y (80 terms; for y < 0 the reciprocal of an upper bound for e^{-y})
K = 80; s = sym(0); term = sym(1); z = abs(y); squarings = 0;
if xsgn(y) < 0
    while xsgn(z - (K + 2)) >= 0, z = z / 2; squarings = squarings + 1; end
end
for k = 0:K
    if k > 0, term = term * z / k; end
    s = s + term;
end
if xsgn(y) >= 0
    L = s;
else
    s = s + z^(K + 1) / factorial(sym(K + 1)) / (1 - z / (K + 2));
    L = (1 / s)^(2^squarings);
end
end


function s = shortest_decimal(x)
s = strings(size(x));
for k = 1:numel(x)
    for p = 1:17
        t = sprintf('%.*g', p, x(k));
        if str2double(t) == x(k), break, end
    end
    s(k) = string(t);
end
end


function m = xmax_exact(x)
m = -xmin_exact(-x);
end
