function s = fmt_dec(x, n, mode)
% FMT_DEC  Exact n-decimal string of a symbolic rational, rounded down or up.
%   s = fmt_dec(x, n, 'floor') is the largest multiple of 10^-n that is <= x, as a char such as '-0.272230';
%   fmt_dec(x, n, 'ceil') the smallest multiple that is >= x.  fmt_dec(x, n) requires x * 10^n to be an
%   integer and prints it exactly.

if nargin < 3, mode = 'exact'; end
y = x * sym(10)^n;
switch mode
    case 'floor', y = floor(y);
    case 'ceil',  y = ceil(y);
    otherwise
        if ~isAlways(y == round(y)), error('fmt_dec:notdecimal', 'value is not a multiple of 10^-%d', n); end
end
m = char(y); sgn = '';
if m(1) == '-', sgn = '-'; m = m(2:end); end
m = [repmat('0', 1, max(0, n + 1 - numel(m))) m];
s = [sgn m(1:end - n) '.' m(end - n + 1:end)];
end
