function v = getopt(s, f, v)
% GETOPT  v = getopt(s, f, default): the field f of the options struct s if present and nonempty, else default.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
