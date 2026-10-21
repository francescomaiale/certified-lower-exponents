function Q = parse_rationals(s)
% PARSE_RATIONALS  Exact symbolic rationals from strings 'p', 'p/q' or decimals such as '-0.25' or '1e-05'.
%   Q = parse_rationals(s) for a string array (or cellstr, or char) s returns a sym array of the same size.
%   str2sym costs about 0.7 ms per string, which dominates for certificates with 10^5 entries.  Here every
%   integer is cut into chunks of 15 digits, each chunk is an exact double, and the integers are assembled with
%   vectorized symbolic operations (exact).  Decimals are first mapped to 'p/q' by dec2rat.

s = string(s); sz = size(s);
if isempty(s), Q = sym(zeros(sz)); return, end
s = strtrim(s(:));
dec = ~contains(s, "/");
if any(dec), s(dec) = dec2rat(s(dec)); end
neg = startsWith(s, "-"); s = regexprep(s, '^[+-]', '');
hasq = contains(s, "/");
num = s; den = repmat("1", size(s));
if any(hasq)
    pq = split(s(hasq), "/", 2);
    num(hasq) = pq(:, 1); den(hasq) = pq(:, 2);
end
if any(~cellfun(@(c) all(c >= '0' & c <= '9') && ~isempty(c), cellstr([num; den])))
    error('parse_rationals:parse', 'not a rational number string');
end
N = int_from_digits(num); Dn = int_from_digits(den);
if any(isAlways(Dn == 0)), error('parse_rationals:parse', 'zero denominator'); end
N(neg) = -N(neg);
Q = reshape(N ./ Dn, sz);
end


function X = int_from_digits(t)
% exact nonnegative integers from digit strings, 15 digits per chunk
L = strlength(t); nc = max(1, ceil(max(L) / 15));
pad = arrayfun(@(k) repmat('0', 1, nc * 15 - L(k)), (1:numel(t)).', 'UniformOutput', false);
C = char(strcat(string(pad), t));                       % numel(t) x (15 nc), right aligned
X = sym(zeros(numel(t), 1)); base = sym(10)^15;
for c = 1:nc
    blk = C(:, (c - 1) * 15 + (1:15));
    v = (double(blk) - 48) * (10 .^ (14:-1:0)).';       % exact: < 10^15 < 2^53
    X = X * base + sym(v);
end
end
