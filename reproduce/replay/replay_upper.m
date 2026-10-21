function [ok, msgs, neg, info] = replay_upper(C, A)
% REPLAY_UPPER  Independent exact re-check of a chain certificate (prop:upper), using only the matrices I - h A_i.
%   [ok, msgs, neg, info] = replay_upper(C, A)
%     C  certificate: a .json file (MATLAB or Python, the format of certify() in certify_upper.py), JSON text,
%        or the decoded struct returned by certify_upper_exact
%     A  the data: cell array {A_1, ..., A_m} of numeric matrices or of decimal-string matrices (read_cases_json)
%   ok is true iff every check passes; msgs lists failed checks, neg the negative control, info.slack the
%   printed value minus the floating-point value (1 - rho(R_w)^(-1/k))/h of the word.
%
%   Shares no code with the certify path and inverts no matrix: for both readings of the data (exact decimals and
%   exact binary64 values) and the exact decimal h, with M_i = I - h A_i, it checks in exact rational arithmetic
%     - M_i is a Z-matrix, z_i > 0 and M_i z_i > 0 (so R_i >= 0 and h eta(A_i) < 1);
%     - the chain has k + 1 vectors, x_0 > 0, and M_{w_j} x_j >= x_{j-1} for j = 1..k;
%     - r > 0, r^k >= q = max_r (x_k)_r / (x_0)_r, and (1 - 1/r)/h <= the printed value P;
%     - P is the stored exact value of the printed six decimals.
%   Optional diagnostic: test whether changing the first letter invalidates the chain.
%   This need not fail for identical modes and is not part of the proof.

T0 = tic;
if ischar(C) || isstring(C)
    if isfile(C), C = jsondecode(fileread(C)); else, C = jsondecode(char(C)); end
end
if nargin < 2 || isempty(A), A = replay_data(C.('case')); end   % data re-read from data/*.json
msgs = {}; neg = {}; ok = true;
raw = A; m = numel(raw); d = size(raw{1}, 1);
dstr = cell(1, m);
for i = 1:m
    if isnumeric(raw{i}), dstr{i} = short_repr(raw{i}); else, dstr{i} = string(raw{i}); end
end
if isnumeric(C.word), w = C.word(:).'; else, w = double(char(C.word)) - double('0'); end
k = numel(w);
hq = rat_parse(string(C.h)); P = rat_parse(string(C.printed_exact));
info = struct('k', k, 'printed', C.printed, 'slack', NaN, 'U_word', NaN, 'secs', toc(T0));
if sgn(hq) <= 0 || k ~= C.k || k == 0 || any(~isfinite(w) | w ~= round(w) | w < 1 | w > m)
    ok = false; msgs{end + 1} = 'word length mismatch or invalid letter';
    if sgn(hq) <= 0, msgs{end + 1} = 'h must be positive'; end
    return
end
if sgn(P - rat_parse(string(C.printed))) ~= 0
    ok = false; msgs{end + 1} = 'printed value is not the stored exact value';
end
readings = {'decimal', 'binary64'}; Ms = struct();
for rr = 1:numel(readings)
    rd = readings{rr};
    if ~isfield(C.versions, rd)
        ok = false; msgs{end + 1} = sprintf('%s: no certificate for this reading', rd); continue
    end
    Vr = C.versions.(rd);
    if numel(Vr.z) ~= m
        ok = false; msgs{end + 1} = sprintf('%s: need one witness per mode', rd); return
    end
    for i = 1:m
        if strcmp(rd, 'decimal')
            Ai = reshape(rat_parse(dstr{i}(:).'), d, d);
        else
            Ai = sym(str2double(dstr{i}), 'f');
        end
        Ms.(rd){i} = sym(eye(d)) - hq * Ai;
    end
    for i = 1:m
        M = Ms.(rd){i}; off = M(~eye(d));
        if any(sgn(off) > 0)
            ok = false; msgs{end + 1} = sprintf('%s: I - hA_%d is not a Z-matrix', rd, i - 1); %#ok<*AGROW>
        end
        if numel(Vr.z) < i
            ok = false; msgs{end + 1} = sprintf('%s: mode %d: no witness z', rd, i - 1); continue
        end
        z = rat_parse(vector_strings(Vr.z{i}).').';
        if numel(z) ~= d || any(sgn(z) <= 0) || any(sgn(M * z) <= 0)
            ok = false; msgs{end + 1} = sprintf('%s: mode %d: no strictly positive M-matrix witness', rd, i - 1);
        end
    end
    ch = Vr.chain; nc = numel(ch);
    if nc ~= k + 1 || any(cellfun(@(x) numel(vector_strings(x)), ch) ~= d)
        ok = false; msgs{end + 1} = sprintf('%s: wrong chain length or dimension', rd); return
    end
    cs = strings(d, nc);
    for j = 1:nc, cs(:, j) = vector_strings(ch{j}); end
    X = reshape(rat_parse(cs(:).'), d, nc);
    if nc ~= k + 1 || any(sgn(X(:, 1)) <= 0)
        ok = false; msgs{end + 1} = sprintf('%s: chain length or x_0 > 0 fails', rd); continue
    end
    for i = 1:m                                   % M_{w_j} x_j - x_{j-1} >= 0, all steps j with w_j = i
        J = find(w == i);
        if isempty(J), continue, end
        Dm = Ms.(rd){i} * X(:, J + 1) - X(:, J);
        [~, jj] = find(sgn(Dm) < 0, 1);
        if ~isempty(jj)
            ok = false; msgs{end + 1} = sprintf('%s: chain inequality fails at step %d', rd, J(jj));
        end
    end
    q = xmax(X(:, end) ./ X(:, 1)); r = rat_parse(string(Vr.r));
    if ~(sgn(r) > 0) || sgn(r^k - q) < 0
        ok = false; msgs{end + 1} = sprintf('%s: r <= 0 or r^k < q', rd);
    end
    if sgn(r) > 0 && sgn((1 - 1 / r) / hq - P) > 0
        ok = false; msgs{end + 1} = sprintf('%s: (1 - 1/r)/h exceeds the printed value', rd);
    end
end
if ~ok, return, end
% negative control: the stored binary64 chain fails its first step once the first letter is flipped
if k >= 1 && m >= 2 && isfield(Ms, 'binary64') && isfield(C.versions, 'binary64') && numel(C.versions.binary64.chain) >= 2
    alt = 1 + mod(w(1), m);                       % 1 <-> 2 for pairs
    ch = C.versions.binary64.chain;
    x0 = rat_parse(vector_strings(ch{1}).').'; x1 = rat_parse(vector_strings(ch{2}).').';
    fails = any(sgn(Ms.binary64{alt} * x1 - x0) < 0);
    if fails, neg{end + 1} = 'flipped first letter rejected'; else, neg{end + 1} = 'flipped first letter NOT rejected'; end
end
% tightness: float value of the word itself
hf = double(hq); Pf = eye(d); lg = 0;
for a = w
    Ra = inv(eye(d) - hf * str2double(dstr{a}));
    Pf = Ra * Pf; s = norm(Pf); Pf = Pf / s; lg = lg + log(s);
end
Uw = (1 - exp(-(log(max(abs(eig(Pf)))) + lg) / k)) / hf;
info = struct('k', k, 'printed', C.printed, 'slack', double(P) - Uw, 'U_word', Uw, 'secs', toc(T0));
end


function [m, i] = xmin(x)
% exact minimum (vectorized comparisons, replay/private/xmin_exact)
[m, i] = xmin_exact(x);
end


function [m, i] = xmax(x)
% exact maximum of a vector of symbolic rationals
[m, i] = xmin(-x); m = -m;
end


function s = sgn(x)
s = double(sign(x));
end


function s = vector_strings(x)
% MATLAB can decode a scalar witness as char, including all digits of a rational.
if ischar(x), s = string(x); else, s = string(x(:)); end
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
