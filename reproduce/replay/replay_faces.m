function [ok, msgs, info] = replay_faces(C, A)
% REPLAY_FACES  Independent exact re-check of an invariant-face certificate (invariant_faces; Computation 3.8).
%   [ok, msgs, info] = replay_faces(C, A)
%     C  certificate (.json file, JSON text or decoded struct) with fields case and orientations (one entry per
%        orientation, as written by comp38_reducible_pairs); A data (optional; default from data/*.json)
%   Shares no code with invariant_faces.  For each orientation, with the binary64 data (transposed for 'T'):
%     - the bound for c: x > 0 and (S_alpha x)_r <= U x_r for every r (Collatz-Wielandt, cross-multiplied);
%     - the list of faces J with S(J) nonempty, and S(J) itself, are recomputed from the zero pattern of the data
%       and compared with the stored list;
%     - for each face with a single mode i: x >= 0, x ~= 0 and A_i[J,J] x >= L x entrywise, so eta >= L;
%     - the flag slow equals (L > U), and all_slow is the conjunction of the flags.

T0 = tic;
if ischar(C) || isstring(C)
    if isfile(C), C = jsondecode(fileread(C)); else, C = jsondecode(char(C)); end
end
if nargin < 2 || isempty(A), A = replay_data(C.('case')); end
msgs = {}; m = numel(A); d = size(A{1}, 1);
OR = C.orientations; if iscell(OR), OR = [OR{:}]; end
info = struct();
if isempty(OR) || m ~= 2
    ok = false; msgs{end + 1} = 'need a matrix pair and at least one orientation'; return
end
nf = 0;
for o = 1:numel(OR)
    F = OR(o); tag = sprintf('orientation %s', F.orientation);
    if ~ismember(string(F.orientation), ["A", "N", "T"])
        ok = false; msgs{end + 1} = 'invalid orientation'; return
    end
    As = cell(1, m);
    for i = 1:m
        if isnumeric(A{i}), As{i} = sym(A{i}, 'f'); else, As{i} = sym(str2double(string(A{i})), 'f'); end
        if any(xsgn(As{i}(~eye(d))) < 0), ok = false; msgs{end + 1} = 'data must be Metzler'; return, end
        if strcmp(F.orientation, 'T'), As{i} = As{i}.'; end
    end
    a = rat_parse(string(F.alpha)); x = rat_parse(string(F.x_c(:)).').'; U = rat_parse(string(F.c_upper));
    if xsgn(a) < 0 || xsgn(a - 1) > 0 || numel(x) ~= d
        ok = false; msgs{end + 1} = 'invalid mixture or witness'; return
    end
    S = a * As{1} + (1 - a) * As{2};
    if any(xsgn(x) <= 0) || any(xsgn(S * x - U * x) > 0)
        msgs{end + 1} = sprintf('%s: Collatz-Wielandt bound for c fails', tag); %#ok<*AGROW>
    end
    % faces with S(J) nonempty, recomputed
    want = {};
    for mask = 1:2^d - 2
        J = find(bitget(mask, 1:d)); out = setdiff(1:d, J);
        Si = find(cellfun(@(M) all(xsgn(M(out, J)) == 0, 'all'), As));
        if ~isempty(Si), want{end + 1} = sprintf('%s|%s', mat2str(J), mat2str(Si)); end
    end
    faces = F.faces; if iscell(faces), faces = [faces{:}]; end
    have = arrayfun(@(f) sprintf('%s|%s', mat2str(f.J(:).'), mat2str(f.S(:).')), faces, 'UniformOutput', false);
    if ~isequal(sort(want(:)), sort(have(:)))
        msgs{end + 1} = sprintf('%s: stored invariant faces differ from the recomputed ones', tag);
        ok = false; return
    end
    allslow = true;
    for f = faces(:).'
        nf = nf + 1; J = f.J(:).'; slow = false;
        if numel(f.S) == 1
            M = As{f.S}(J, J); xj = rat_parse(string(f.x(:)).').'; L = rat_parse(string(f.eta_lower));
            if numel(xj) ~= numel(J) || any(xsgn(xj) < 0) || all(xsgn(xj) == 0) || any(xsgn(M * xj - L * xj) < 0)
                msgs{end + 1} = sprintf('%s: face %s: Collatz-Wielandt lower bound fails', tag, mat2str(J));
            end
            slow = xsgn(L - U) > 0;
        end
        if slow ~= logical(f.slow)
            msgs{end + 1} = sprintf('%s: face %s: slow flag inconsistent', tag, mat2str(J));
        end
        allslow = allslow && slow;
    end
    if allslow ~= logical(F.all_slow)
        msgs{end + 1} = sprintf('%s: all_slow inconsistent', tag);
    end
end
ok = isempty(msgs);
info = struct('faces', nf, 'secs', toc(T0));
end
