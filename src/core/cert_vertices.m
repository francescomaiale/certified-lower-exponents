function V = cert_vertices(C)
% CERT_VERTICES  The vertices of a polytope certificate (file, JSON text or decoded struct) as a d x n sym matrix.
if ischar(C) || isstring(C)
    if isfile(C), C = jsondecode(fileread(C)); else, C = jsondecode(char(C)); end
end
vs = C.vertices; n = numel(vs);
if ~iscell(vs), vs = num2cell(vs, 2); n = numel(vs); end
d = numel(vs{1}); s = strings(d, n);
for k = 1:n, s(:, k) = string(vs{k}(:)); end
V = parse_rationals(s);
end
