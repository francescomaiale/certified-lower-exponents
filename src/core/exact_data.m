function D = exact_data(A, h, versions)
% EXACT_DATA  Floating-point and exact rational readings of the data (A_1, ..., A_m; h).
%   D = exact_data(A, h) with A a cell array of d x d matrices and h the step.
%   Each A{i} is either numeric or a string/cellstr matrix of decimals (e.g. read from cases.json with
%   read_cases_json).  h is numeric or a decimal string.
%
%   Two exact readings of the data are formed:
%     'decimal'   every entry is the exact decimal of its string; for numeric input, the string is the shortest
%                 decimal that rounds to the double (Python's repr), e.g. 0.3 -> 3/10;
%     'binary64'  every entry is the exact binary64 value of that decimal, e.g. 0.3 -> 5404319552844595/2^54.
%   The step h is the exact decimal in both readings (h = 0.1 -> 1/10).  versions selects a subset, e.g.
%   {'binary64'}.  Every certificate is checked against every reading in D.ver.
%
%   Fields: d, m, Af (cell, double), Astr (cell, string decimals), h_str, hq (sym), hf (double),
%           ver (struct array with fields name, A (cell of sym)).

if nargin < 3 || isempty(versions), versions = {'decimal', 'binary64'}; end
if ~iscell(A) || isempty(A), error('exact_data:input', 'A must be a nonempty cell array of matrices'); end
m = numel(A); d = size(A{1}, 1);
if d == 0, error('exact_data:input', 'matrices must be nonempty'); end
D.d = d; D.m = m;
D.Astr = cell(1, m); D.Af = cell(1, m);
for i = 1:m
    if isnumeric(A{i})
        if ~isreal(A{i}) || any(~isfinite(A{i}(:))), error('exact_data:input', 'matrices must be finite and real'); end
        D.Astr{i} = short_repr(A{i});
    else
        D.Astr{i} = string(A{i});
    end
    D.Af{i} = str2double(D.Astr{i});
    if ~isreal(D.Af{i}) || any(~isfinite(D.Af{i}(:))) || ~isequal(size(D.Af{i}), [d d])
        error('exact_data:input', 'A{%d} is not a %d x %d matrix of numbers', i, d, d);
    end
end
if isnumeric(h)
    if ~isscalar(h) || ~isreal(h) || ~isfinite(h), error('exact_data:input', 'h must be a finite real scalar'); end
    D.h_str = short_repr(h);
else
    D.h_str = string(h);
end
D.hf = str2double(D.h_str);
if ~isscalar(D.hf) || ~isreal(D.hf) || ~isfinite(D.hf) || D.hf <= 0
    error('exact_data:input', 'h must be a finite positive scalar');
end
D.hq = str2sym(dec2rat(D.h_str));
D.ver = struct('name', {}, 'A', {});
for v = 1:numel(versions)
    Aq = cell(1, m);
    for i = 1:m
        switch versions{v}
            case 'decimal'
                Aq{i} = reshape(str2sym(cellstr(dec2rat(D.Astr{i}))), d, d);
            case 'binary64'
                Aq{i} = sym(D.Af{i}, 'f');
            otherwise
                error('exact_data:version', 'unknown data version %s', versions{v});
        end
        if any(double(sign(Aq{i}(~eye(d)))) < 0)
            error('exact_data:metzler', 'A{%d} is not Metzler in the %s reading', i, versions{v});
        end
    end
    D.ver(v).name = versions{v}; D.ver(v).A = Aq;
end
end


function s = short_repr(x)
% Shortest decimal string that rounds to x (as Python's repr).
s = strings(size(x));
for k = 1:numel(x)
    for p = 1:17
        t = sprintf('%.*g', p, x(k));
        if str2double(t) == x(k), break, end
    end
    s(k) = string(t);
end
end
