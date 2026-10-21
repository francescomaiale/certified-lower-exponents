function s = short_repr(x)
% SHORT_REPR  Shortest decimal strings that round to the doubles x (as Python's repr), e.g. 0.1 -> "0.1".
s = strings(size(x));
for k = 1:numel(x)
    for p = 1:17
        t = sprintf('%.*g', p, x(k));
        if str2double(t) == x(k), break, end
    end
    s(k) = string(t);
end
end
