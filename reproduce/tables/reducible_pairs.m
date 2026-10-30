function T = reducible_pairs(varargin)
% REDUCIBLE_PAIRS  Computation 3.9 (obs:cone): the reducible test pairs, cases 77, 148, 677, 1195, 1440.
%   - Invariant faces (both orientations): every proper face F_J invariant under some mode, with an exact lower
%     bound for its exponent eta(A_i[J,J]) (Collatz-Wielandt with a nonnegative eigenvector), against a certified
%     upper bound for the chattering value c >= lambda_check (Collatz-Wielandt at a rational alpha); binary64 data.
%     Certificates: certificates/reducible_pairs/faces_<case>.json, replayed by replay_faces.
%   - Case 148, transposed family, h1 = 1/5: the cone C = {x >= 0 : (1 - kappa)(x1 + x2) >= kappa x3}, kappa = 1/20,
%     with extreme rays e1, e2, (kappa, 0, 1 - kappa), (0, kappa, 1 - kappa), is invariant under both Euler maps, and
%     its support classes are all nonempty subsets of {1,2,3} except {3}.  certificates/reducible_pairs/cone_case148_T.json.
%   - Cases 77 (as given) and 677 (transposed), h1 = 0.2: invariant polygons of prop:cone(ii), from the hull of the
%     projective Euler orbit of 1, inflated (eps 1e-3 and 1e-4 as in the paper; if MATLAB's hull then misses exact
%     invariance by an image on the boundary, the next factor of a short list).  certificates/reducible_pairs/polygon_case77_A.json and
%     polygon_case677_T.json.
%   Every cone certificate is checked exactly (facet normals) and replayed by a different method (replay_cone: fan
%   triangulation and barycentric coordinates).  Summary: results/reducible_pairs.md.  Runtime: about a minute.

p = inputParser; p.addParameter('replay', true); p.parse(varargin{:}); o = p.Results;
root = paperB_root(); cdir = fullfile(root, 'certificates', 'reducible_pairs');
say = logger(fullfile(root, 'results', 'reducible_pairs.md'));
say('# Computation 3.9 (obs:cone): the reducible test pairs, MATLAB reproduction');
say('');
say('## Invariant faces (binary64 data), both orientations');
say('');
say('Exponent bounds are exact lower bounds rounded down to 9 decimals; the paper rounds them to 6 decimals.');
say('');
say('| case | orientation | certified c <= | invariant proper faces (J: mode, exact eta >=) | verdict | replay |');
say('|---|---|---|---|---|---|');
T = struct();
for cs = {'case77', 'case148', 'case677', 'case1195', 'case1440'}
    pid = cs{1}; A = load_pair('test', pid);
    F = struct('case', pid, 'orientations', {{invariant_faces(A, 'A'), invariant_faces(A, 'T')}});
    file = fullfile(cdir, ['faces_' pid '.json']); write_text(file, jsonencode(F));
    ok = true; if o.replay, ok = replay_faces(file); end
    for k = 1:2
        G = F.orientations{k}; parts = {};
        for f = G.faces
            if numel(f.S) == 1
                parts{end + 1} = sprintf('{%s}: A_%d%s, %s (%s)', strjoin(string(f.J), ','), f.S, tsym(G.orientation), ...
                    fmt_dec(parse_rationals(string(f.eta_lower)), 9, 'floor'), slowstr(f.slow)); %#ok<AGROW>
            else
                parts{end + 1} = sprintf('{%s}: modes %s', strjoin(string(f.J), ','), mat2str(f.S)); %#ok<AGROW>
            end
        end
        if isempty(parts), parts = {'none'}; end
        if G.all_slow, v = 'all invariant faces slow: lambda_min = lambda_check'; else, v = 'some face not certified slow'; end
        say('| %s | %s | %s | %s | %s | %s |', pid, G.orientation, fmt_dec(parse_rationals(string(G.c_upper)), 6, 'ceil'), ...
            strjoin(parts, '; '), v, tf(ok));
    end
    T.(pid) = F;
end
say('');
say('Paper: 1195 (as given) {2,3} under A_1, >= -0.289318, c <= -0.532919; 677 (as given) {2} under A_2, -0.237097,');
say('c <= -0.394492; 77 transposed {2} under A_2^T, 0.024892, c <= -0.148186; 1440 none; 77 as given -0.663165 and');
say('677 transposed -1.650614 (below c); 148: {2} under A_2, about -1.368010, transposed {3} under A_1^T, -3.669288,');
say('and {1,3} under A_2^T, >= 0.147543, c <= -1.236979.');
% ---------------------------------------------------------------------------------- case 148: the cut cone
say('');
say('## Case 148, transposed family: a cone that touches only slow faces (cor:eulerface)');
A = load_pair('test', 'case148'); kap = sym(1) / 20;
V = [1 0 kap 0; 0 1 0 kap; 0 0 1 - kap 1 - kap];
[C, info] = certify_cone(V, A, 'T', sym(1) / 5, struct('name', 'case148', 'kind', 'cut cone'));
file = fullfile(cdir, 'cone_case148_T.json'); write_text(file, jsonencode(C));
ok = true; msgs = {}; if o.replay, [ok, msgs] = replay_cone(file); end
say('- C = {x >= 0 : (19/20)(x1 + x2) >= x3/20}, h1 = 1/5: %d extreme rays, %d images outside, invariant: %s; replay %s %s', ...
    C.n_vertices, info.images_outside, yesno(C.ok), tf(ok), strjoin(msgs, '; '));
say('- support classes of C: %s (all nonempty subsets except {3})', strjoin(string(C.classes), ' | '));
F = T.case148.orientations{2};
inv = arrayfun(@(f) strjoin(string(f.J), ','), F.faces, 'UniformOutput', false);
cls = string(C.classes);
touched = inv(ismember(string(inv), cls));
say('- invariant faces of the transposed family: %s; touched by C: %s', strjoin(string(inv), ' | '), strjoin(string(touched), ' | '));
for f = F.faces
    if ismember(string(strjoin(string(f.J), ',')), cls)
        say('  face {%s} under A_%d^T: exponent >= %s, slow: %s (c <= %s)', strjoin(string(f.J), ','), f.S, ...
            fmt_dec(parse_rationals(string(f.eta_lower)), 9, 'floor'), yesno(f.slow), ...
            fmt_dec(parse_rationals(string(F.c_upper)), 6, 'ceil'));
    end
end
% ---------------------------------------------------------------------------------- cases 77 and 677: polygons
say('');
say('## Cases 77 (as given) and 677 (transposed): invariant polygons of prop:cone(ii), h1 = 0.2');
for c = {'case77', 'A', [1e-3 5e-4 2e-3]; 'case677', 'T', [1e-4 5e-5 5e-4 1e-3]}.'
    pid = c{1}; orient = c{2}; A = load_pair('test', pid);
    % the paper's inflation factor first; MATLAB's convhull treats nearly collinear points differently from Qhull,
    % and an image can end up on the boundary (exact check fails by one image): then the next factor is tried
    for ep = c{3}
        [Pv, pinfo] = orbit_polygon(A, orient, 0.2, ep);
        [C, info] = certify_cone(Pv, A, orient, '0.2', struct('name', pid, 'kind', 'polygon', 'require_positive', true));
        if C.ok, break, end
        say('- %s (%s), inflation %g: %d vertices, %d images outside (exact); next factor', pid, orient, ep, ...
            C.n_vertices, info.images_outside);
    end
    C.inflation = ep;
    file = fullfile(cdir, sprintf('polygon_%s_%s.json', pid, orient)); write_text(file, jsonencode(C));
    ok = true; msgs = {}; if o.replay, [ok, msgs] = replay_cone(file); end
    say('- %s (%s): %d vertices, min coordinate %.4f, inflation %g, %d rounds; exact invariance %s; replay %s %s', ...
        pid, orient, C.n_vertices, min(C.min_coordinates), ep, pinfo.inflation_rounds, yesno(C.ok), tf(ok), ...
        strjoin(msgs, '; '));
end
say('');
say('Paper: 21 and 47 vertices, coordinates at least 0.0796 and 0.1134 (the counts depend on the floating-point hull;');
say('only the exact checks matter).');
end


function s = tsym(o)
if o == 'T', s = '^T'; else, s = ''; end
end


function s = slowstr(b)
if b, s = 'slow'; else, s = 'not certified slow'; end
end


function s = tf(b)
if b, s = 'VALID'; else, s = 'INVALID'; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
