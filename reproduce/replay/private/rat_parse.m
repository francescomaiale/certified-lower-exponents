function Q = rat_parse(s)
% RAT_PARSE  (replay side) Exact sym rationals from strings: integers, 'p/q', decimals ('-0.25', '1e-05').
%   Independent of the certify path.  Digit strings are split into blocks of 9 digits from the right and
%   combined with exact symbolic arithmetic; this avoids str2sym, which costs about 0.7 ms per string.
s = string(s); sz = size(s); s = strtrim(s(:)); n = numel(s);
if n == 0, Q = sym(zeros(sz)); return, end
sg = ones(n, 1); nd = strings(n, 1); dd = strings(n, 1); ex = zeros(n, 1);
for k = 1:n
    t = char(s(k));
    if ~isempty(t) && (t(1) == '-' || t(1) == '+')
        if t(1) == '-', sg(k) = -1; end
        t = t(2:end);
    end
    sl = find(t == '/', 1);
    if ~isempty(sl)
        nd(k) = t(1:sl - 1); dd(k) = t(sl + 1:end);
    else
        ie = find(t == 'e' | t == 'E', 1);              % exponent split by hand (no optional regexp groups)
        if isempty(ie), mant = t; xp = 0; else, mant = t(1:ie - 1); xp = str2double(t(ie + 1:end)); end
        if ~isfinite(xp) || xp ~= round(xp), error('replay:parse', 'not a number: "%s"', s(k)); end
        dot = find(mant == '.');
        if numel(dot) > 1, error('replay:parse', 'not a number: "%s"', s(k)); end
        if isempty(dot), ip = mant; fp = ''; else, ip = mant(1:dot - 1); fp = mant(dot + 1:end); end
        if isempty(ip) && isempty(fp), error('replay:parse', 'not a number: "%s"', s(k)); end
        nd(k) = [ip fp]; dd(k) = "1"; ex(k) = -numel(fp) + xp;
    end
    if isempty(char(nd(k))) || isempty(char(dd(k))) || any(~isstrprop(char(nd(k)), 'digit')) || ...
            any(~isstrprop(char(dd(k)), 'digit'))
        error('replay:parse', 'not a number: "%s"', s(k));
    end
end
N = digits_to_sym(nd); D = digits_to_sym(dd);
if any(isAlways(D == 0)), error('replay:parse', 'zero denominator'); end
Q = sym(sg) .* N ./ D;
pe = ex > 0; ne = ex < 0;
if any(pe), Q(pe) = Q(pe) .* sym(10) .^ sym(ex(pe)); end
if any(ne), Q(ne) = Q(ne) ./ sym(10) .^ sym(-ex(ne)); end
Q = reshape(Q, sz);
end

function X = digits_to_sym(t)
L = strlength(t); nb = max(1, ceil(max(L) / 9));
X = sym(zeros(numel(t), 1));
blocks = zeros(numel(t), nb);                  % blocks(:, 1) least significant
for k = 1:numel(t)
    c = char(t(k));
    for b = 1:ceil(numel(c) / 9)
        hi = numel(c) - 9 * (b - 1); lo = max(1, hi - 8);
        blocks(k, b) = sum((double(c(lo:hi)) - 48) .* 10 .^ (hi - lo:-1:0));
    end
end
for b = nb:-1:1
    X = X * sym(1e9) + sym(blocks(:, b));
end
end
