function r = dec2rat(s)
% DEC2RAT  Exact rational form 'p/q' of decimal strings, for str2sym.
%   r = dec2rat(s) maps a string array (or char, or cellstr) of decimals such as '-2.163042579047693', '0.3',
%   '1e-05' or '7' to a string array of the same size holding 'p' or 'p/q' with the exact value of each decimal.
%   str2sym('0.1') is a symbolic float, not 1/10, so decimals go through this function: str2sym(dec2rat('0.1'))
%   is exactly 1/10.  Rational strings 'p/q' are passed through unchanged.

s = string(s); r = s;
for k = 1:numel(s)
    t = strtrim(s(k));
    if contains(t, "/")
        r(k) = t; continue
    end
    % (MATLAB's regexp drops a capture nested in an optional non-capturing group, so the exponent group is a
    % plain capture including the letter e)
    tok = regexp(char(t), '^([+-]?)(\d*)\.?(\d*)([eE][+-]?\d+)?$', 'tokens', 'once');
    if ~isempty(tok), tok(end + 1:4) = {''}; if ~isempty(tok{4}), tok{4} = tok{4}(2:end); end, end
    if isempty(tok) || (isempty(tok{2}) && isempty(tok{3}))
        error('dec2rat:parse', 'not a decimal number: "%s"', t);
    end
    sgn = tok{1}; if strcmp(sgn, '+'), sgn = ''; end
    digs = [tok{2} tok{3}];
    e = -numel(tok{3});
    if ~isempty(tok{4}), e = e + str2double(tok{4}); end
    digs = regexprep(digs, '^0+', '');
    if isempty(digs)
        r(k) = "0"; continue
    end
    if e >= 0
        r(k) = string([sgn digs repmat('0', 1, e)]);
    else
        r(k) = string([sgn digs '/1' repmat('0', 1, -e)]);
    end
end
end
