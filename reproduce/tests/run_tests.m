function results = run_tests(opts)
% RUN_TESTS  End-to-end tests of the package (extends the 44 checks of the earlier test_certified_bracket).
%   run_tests                 about 5 to 10 minutes; prints PASS/FAIL/SKIP per check and a summary, and raises an
%                             error (nonzero exit status under matlab -batch) if any check fails
%   run_tests(struct('python', false))   skip the optional Python cross-checks
%   run_tests(struct('python_executable', '/path/to/python3'))   select the Python with NumPy installed
%   Certificates written by the tests go to a fresh tempname folder; in the package only run logs are appended
%   (results/logs).
%   Searches are deterministic (pmax caps on processed vertices; the time caps are only a safety net).
%
%   (a)-(c) certified_bracket on the planar pair, case 1 (h = 0.1) and ex:sharp (eps = 0.2): the brackets contain
%           lambda_check and reproduce the paper's six decimals
%   (d)     negative controls: tampered lower and upper certificates are rejected
%   (e)     the MATLAB replays accept the paper's Python certificates (tests/fixtures) and reject a tampered one
%   (f)     the Python replays (independent_replay/python) accept the MATLAB certificates (optional)
%   (g)     exact-arithmetic helpers (parse_rationals, dec2rat, limit_denominator, fmt_dec, exact_min, exp/log bounds)
%   (h)     chattering and schedule certificates (certify_best_upper / replay_best_upper), with tampering
%   (i)     prop:polytope(ii): a polytope certified at h = 0.1 re-certified at h = 0.05; table certificates
%   (j)     Computation 3.9: invariant faces, the case-148 cone, the case-77 polygon, with tampering
%   (k)     Algorithm (L) at the paper's shift (planar, tau = 0.1) and the tangent bound
%   (l)     Table 1 row (ex:sharp, h = 0.1) and the Perron seeds of Obs. 4.8 (switching pair 4, h = 0.1)
%   (n)     a three-mode schedule (Ex. 4.6), a driver in search mode
%   (m)     no material of the companion paper: no identifier of its headline case in any file, no such number
%           outside digit strings in sources, data and documents, and every certificate refers to a pair of data/

if nargin < 1, opts = struct(); end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); startup_reproduce();
results = cell(0, 3); T0 = tic;
outdir = tempname; mkdir(outdir);
fix = fullfile(root, 'tests', 'fixtures', 'python_certificates');
base = struct('pmax', 2000, 'probe', 120, 'budget', 600, 'outdir', outdir, 'source', 'table2', 'verbose', false);
fprintf('run_tests: certificates written to %s\n\n', outdir);

    function check(label, cond, detail)
        if nargin < 3, detail = ''; end
        if cond, stat = 'PASS'; else, stat = 'FAIL'; end
        results(end + 1, :) = {label, stat, detail};
        if isempty(detail), fprintf('[%s] %s\n', stat, label); else, fprintf('[%s] %s  (%s)\n', stat, label, detail); end
    end

    function skip(label, detail)
        results(end + 1, :) = {label, 'SKIP', detail};
        fprintf('[SKIP] %s  (%s)\n', label, detail);
    end

    function rejected(label, ok)
        if ok, s = 'accepted'; else, s = 'rejected'; end
        check(label, ~ok, s);
    end

%% (a) planar pair
fprintf('---- (a) planar pair, h = 0.1\n');
Ap = load_pair('test', 'planar'); lam_p = -0.1949763083;
o = base; o.name = 'planar';
rp = certified_bracket(Ap, '0.1', o);
check('(a) replays of the new certificates VALID', replay_lower(rp.lower_cert, Ap) && replay_upper(rp.upper_cert, Ap));
check('(a) lower <= lambda_check <= upper', rp.lower_value <= lam_p && lam_p <= rp.upper_value, ...
      sprintf('%s <= %.10f <= %s', rp.lower, lam_p, rp.upper));
check('(a) lower beats -0.281832 (the estimate of Table S1)', rp.lower_value > -0.281832, rp.lower);
check('(a) lower within 2e-6 of the paper -0.272230', abs(rp.lower_value + 0.272230) <= 2e-6, ...
      sprintf('%s, |V| = %d', rp.lower, rp.search.nV));
check('(a) upper = paper -0.123619 with the same word (k = 31, up to rotation)', strcmp(rp.upper, '-0.123619') && ...
      is_rotation(rp.upper_cert.word, '2212212122121221221212212122121'), sprintf('%s, word %s', rp.upper, rp.upper_cert.word));

%% (b) case 1
fprintf('\n---- (b) case 1, h = 0.1\n');
A1 = load_pair('test', 'case1'); lam_1 = -1.0123958;
o = base; o.name = 'case1';
r1 = certified_bracket(A1, '0.1', o);
check('(b) replays of the new certificates VALID', replay_lower(r1.lower_cert, A1) && replay_upper(r1.upper_cert, A1));
check('(b) lower <= lambda_check <= upper', r1.lower_value <= lam_1 && lam_1 <= r1.upper_value, ...
      sprintf('%s <= %.7f <= %s', r1.lower, lam_1, r1.upper));
check('(b) lower = paper -1.044234', strcmp(r1.lower, '-1.044234'), sprintf('%s, |V| = %d', r1.lower, r1.search.nV));
check('(b) upper = paper -0.988047 with the same word (k = 17, up to rotation)', strcmp(r1.upper, '-0.988047') && ...
      is_rotation(r1.upper_cert.word, '12111112111121111'), sprintf('%s, word %s', r1.upper, r1.upper_cert.word));

%% (c) ex:sharp
fprintf('\n---- (c) ex:sharp, eps = 0.2, h = 0.1\n');
ep = 0.2; h = 0.1; As = load_pair('test', 'sharp');
dg = (1 - sqrt(1 - h^2)) / h; g = ep - dg; G = ep + dg;
o = base; o.name = 'sharp'; o.source = 'matlab';
rs = certified_bracket(As, '0.1', o);
check('(c) replays of the new certificates VALID', replay_lower(rs.lower_cert, As) && replay_upper(rs.upper_cert, As));
check('(c) lower <= g(h) < lambda_check = eps < G(h) <= upper', rs.lower_value <= g + 1e-12 && G <= rs.upper_value + 1e-12, ...
      sprintf('%s <= %.9f < %.1f < %.9f <= %s', rs.lower, g, ep, G, rs.upper));
check('(c) lower within 2e-6 of g(h)', g - rs.lower_value <= 2e-6, sprintf('g - lower = %.2e', g - rs.lower_value));
check('(c) upper within 2e-6 of G(h)', rs.upper_value - G <= 2e-6, sprintf('upper - G = %.2e', rs.upper_value - G));

%% (d) negative controls
fprintf('\n---- (d) negative controls (tampered certificates must be rejected)\n');
for cs = {{'planar', rp}, {'case1', r1}}
    nm = cs{1}{1}; r = cs{1}{2};
    C = r.lower_cert;
    T = C; T.t = char(parse_rationals(string(C.t_star)) * (1 + sym(1) / sym(10)^12));
    rejected(sprintf('(d) %s lower: t raised to t*(1 + 1e-12)', nm), replay_lower(T));
    T = C; T.lambda_lower = fmt_dec(parse_rationals(string(C.lambda_lower)) + sym(1) / sym(10)^6, 6);
    rejected(sprintf('(d) %s lower: claimed bound raised by 1e-6', nm), replay_lower(T));
    T = C; f = fieldnames(T.weights(1).lambda); x = parse_rationals(string(T.weights(1).lambda.(f{1})));
    T.weights(1).lambda.(f{1}) = char(x * (1 + sym(1) / sym(10)^9));
    rejected(sprintf('(d) %s lower: one weight changed by 1e-9 (sum ~= 1)', nm), replay_lower(T));
    T = C; jb = C.binding.vertex + 1;
    T.vertices{jb}{1} = char(parse_rationals(string(C.vertices{jb}{1})) * (1 + sym(1) / sym(10)^9));
    rejected(sprintf('(d) %s lower: binding vertex changed by 1e-9', nm), replay_lower(T));
    C = r.upper_cert; k = C.k; j = ceil(k / 2);
    T = C; T.word(j) = char('1' + '2' - T.word(j));
    rejected(sprintf('(d) %s upper: letter %d of the word flipped', nm, j), replay_upper(T));
    T = C; x = parse_rationals(string(T.versions.binary64.chain{j}{1}));
    T.versions.binary64.chain{j}{1} = char(x * (1 - sym(1) / sym(10)^9));
    rejected(sprintf('(d) %s upper: chain entry x_%d(1) shrunk by 1e-9', nm, j - 1), replay_upper(T));
    T = C; T.versions.decimal.r = char(parse_rationals(string(C.versions.decimal.r)) * (1 - sym(2)^-30));
    rejected(sprintf('(d) %s upper: r lowered by 2^-30', nm), replay_upper(T));
    T = C; P = parse_rationals(string(C.printed)) - sym(1) / sym(10)^6;
    T.printed = fmt_dec(P, 6); T.printed_exact = char(P);
    rejected(sprintf('(d) %s upper: printed bound lowered by 1e-6', nm), replay_upper(T));
end

%% (e) Python certificates replayed in MATLAB
fprintf('\n---- (e) the paper''s Python certificates, MATLAB replay\n');
for nm = {'planar_0.1_best', 'planar_0.1_table', 'case1_0.1_best', 'case148_0.1_rowsumN_trivial', 'case148_0.1_table'}
    C = jsondecode(fileread(fullfile(fix, ['lower_' nm{1} '.json']))); t1 = tic;
    [ok, msgs, neg] = replay_lower(C);
    check(sprintf('(e) lower %s VALID: lambda_check >= %s', nm{1}, C.lambda_lower), ok, ...
          sprintf('|V| = %d, %.1fs; %s', numel(C.vertices), toc(t1), strjoin([msgs, neg], '; ')));
end
C = jsondecode(fileread(fullfile(fix, 'lower_case1_0.1_best.json')));
C.t = char(parse_rationals(string(C.t_star)) * (1 + sym(1) / sym(10)^12));
rejected('(e) lower case1_0.1_best with t raised to t*(1 + 1e-12)', replay_lower(C));
for nm = {'planar_0.1', 'case1_0.1', 'case148_0.1', 'case1440_0.1'}
    C = jsondecode(fileread(fullfile(fix, ['upper_' nm{1} '.json']))); t1 = tic;
    [ok, msgs, neg, info] = replay_upper(C);
    check(sprintf('(e) upper %s VALID: lambda_check <= %s', nm{1}, C.printed), ok, ...
          sprintf('k = %d, slack %.1e, %.1fs; %s', C.k, info.slack, toc(t1), strjoin([msgs, neg], '; ')));
end
for nm = {'planar', 'case1440'}
    [ok, msgs] = replay_best_upper(fullfile(fix, ['best_upper_' nm{1} '.json']));
    check(sprintf('(e) chattering/schedule certificate %s (Python) VALID', nm{1}), ok, strjoin(msgs, '; '));
end

%% (g) exact arithmetic helpers
fprintf('\n---- (g) exact arithmetic helpers\n');
s = ["-0.272229892", "1e-05", "123456789012345678901234567890/7", "0.0", "-17", "3.5E+2"];
q = parse_rationals(s);
ref = [sym(-272229892) / 10^9, sym(1) / 10^5, sym('123456789012345678901234567890') / 7, 0, -17, 350];
check('(g) parse_rationals: decimals, exponents, big p/q', all(isAlways(q == ref)));
check('(g) dec2rat("0.1") is exactly 1/10 (str2sym("0.1") is a symbolic float)', isAlways(str2sym(dec2rat("0.1")) == sym(1) / 10));
bad = {};
for grp = {'test', 'scaling', 'switching', 'bcv'}
    [~, names] = load_pair(grp{1});
    for k = 1:numel(names)
        D = exact_data(load_pair(grp{1}, names{k}), '1');
        for i = 1:D.m
            if max(abs(double(D.ver(1).A{i} - D.ver(2).A{i})), [], 'all') > 1e-15 * max(1, max(abs(D.Af{i}), [], 'all'))
                bad{end + 1} = names{k}; %#ok<AGROW>
            end
        end
    end
end
check('(g) the decimal and binary64 readings of every pair in data/ agree to 1e-15 (parser check)', isempty(bad), ...
      strjoin(unique(bad), ', '));
check('(g) limit_denominator(pi, 1e6) = 3126535/995207 (as Python)', isAlways(limit_denominator(pi, 1e6) == sym(3126535) / 995207));
check('(g) fmt_dec floor/ceil', strcmp(fmt_dec(sym(-1) / 3, 6, 'floor'), '-0.333334') && strcmp(fmt_dec(sym(-1) / 3, 6, 'ceil'), '-0.333333'));
x = sym([3 1 2 1 5]) / 7 + [0 0 0 sym(1) / 10^40 0];
[m, i] = exact_min(x.');
check('(g) exact_min resolves differences of 1e-40', isAlways(m == sym(1) / 7) && i == 2);
E = exp_upper(sym([-1 2; 1 -3]), sym(3) / 2); Ef = expm(1.5 * [-1 2; 1 -3]);
check('(g) exp_upper encloses e^{tA} from above, tightly', all(double(E(:)) >= Ef(:) - 1e-15) && max(abs(double(E(:)) - Ef(:))) < 1e-14);
L2 = log_upper(sym(2));
check('(g) log_upper(2) >= log 2 within 1e-30', double(L2 - sym(log(sym(2)))) >= 0 && double(L2 - log(sym(2))) < 1e-30);

%% (h) chattering and schedule certificates
fprintf('\n---- (h) certified chattering value and schedules\n');
[R, json] = certify_best_upper(load_pair('test', 'planar'), struct('id', 'planar'));
fh = fullfile(outdir, 'planar_best_upper.json'); write_text(fh, json);
check('(h) planar: certified c <= -0.194976 (paper)', strcmp(R.printed, '-0.194976'), R.printed);
check('(h) planar: replay VALID', replay_best_upper(fh));
T = R; T.U_c = char(parse_rationals(string(R.U_c)) - sym(1) / sym(10)^9); T.U_best = T.U_c;
rejected('(h) planar: U_c lowered by 1e-9', replay_best_upper(T));
Lv = paper_levels('best_upper', 'case1440');
sch = struct('modes', [1 2], 'durations', {{Lv.t1, Lv.t2}});
[R, json] = certify_best_upper(load_pair('test', 'case1440'), struct('id', 'case1440', 'schedule', sch));
fh = fullfile(outdir, 'case1440_best_upper.json'); write_text(fh, json);
check('(h) case 1440: certified two-segment schedule 0.345675 (paper)', strcmp(R.printed, '0.345675') && ...
      strcmp(R.best_kind, 'schedule'), R.printed);
check('(h) case 1440: replay VALID (uniformization, no logarithms)', replay_best_upper(fh));
T = R; T.U_s = char(parse_rationals(string(R.U_s)) - sym(1) / sym(10)^7); T.U_best = T.U_s;
rejected('(h) case 1440: U_s lowered by 1e-7', replay_best_upper(T));

%% (i) transfer and table certificates
fprintf('\n---- (i) prop:polytope(ii) and table certificates\n');
Vq = cert_vertices(rp.lower_cert);
out = run_lower_row(Ap, '0.05', struct('name', 'planar', 'kind', 'transfer', 'vertices', Vq, ...
                    'file', fullfile(outdir, 'planar_0.05_transfer.json')));
check('(i) the h = 0.1 polytope certifies at h = 0.05 at least its h = 0.1 bound', out.replay_ok && ...
      str2double(out.lambda) >= rp.lower_value, sprintf('%s >= %s', out.lambda, rp.lower));
out = run_lower_row(Ap, '0.1', struct('name', 'planar', 'kind', 'table', 'level', (1 - 0.1 * 0.281832) * (1 + 1e-9), ...
                    'lambda_claim', '-0.281832', 'file', fullfile(outdir, 'planar_0.1_table.json')));
check('(i) table certificate of -0.281832 (the estimate of Table S1; planar, h = 0.1)', out.replay_ok && strcmp(out.lambda, '-0.281832'), ...
      sprintf('|V| = %d', out.nV));
threw = false;
try
    certify_lower_exact(double(Vq), exact_data(Ap, '0.1'), struct('lambda_claim', '-0.272000'));
catch
    threw = true;
end
check('(i) a claim above the exact level of the polytope is refused', threw);

%% (j) Computation 3.9
fprintf('\n---- (j) Computation 3.9: faces and cones\n');
F = struct('case', 'case1195', 'orientations', {{invariant_faces(load_pair('test', 'case1195'), 'A')}});
fh = fullfile(outdir, 'faces_case1195.json'); write_text(fh, jsonencode(F));
check('(j) case 1195: the only invariant proper face {2,3} is slow (>= -0.289318 > c)', F.orientations{1}.all_slow && ...
      numel(F.orientations{1}.faces) == 1 && strcmp(fmt_dec(parse_rationals(string(F.orientations{1}.faces(1).eta_lower)), 6, 'floor'), '-0.289318'));
check('(j) case 1195: replay_faces VALID', replay_faces(fh));
T = jsondecode(fileread(fh)); T.orientations.faces.eta_lower = char(parse_rationals(string(T.orientations.faces.eta_lower)) + sym(1) / 10^6);
rejected('(j) case 1195: face bound raised by 1e-6', replay_faces(T));
A148 = load_pair('test', 'case148'); kap = sym(1) / 20;
[C, info] = certify_cone([1 0 kap 0; 0 1 0 kap; 0 0 1 - kap 1 - kap], A148, 'T', sym(1) / 5, struct('name', 'case148', 'kind', 'cut cone'));
fh = fullfile(outdir, 'cone148.json'); write_text(fh, jsonencode(C));
check('(j) case 148: the cut cone is invariant (transposed family, h1 = 1/5)', C.ok && info.images_outside == 0);
check('(j) case 148: its support classes are all nonempty subsets but {3}; replay_cone VALID', ...
      isequal(sort(string(C.classes(:))), sort(["1"; "2"; "1,2"; "1,3"; "2,3"; "1,2,3"])) && replay_cone(fh));
k9 = sym(9) / 10;
[C2, info2] = certify_cone([1 0 k9 0; 0 1 0 k9; 0 0 1 - k9 1 - k9], A148, 'T', sym(1) / 5, struct('name', 'case148', 'kind', 'cut cone'));
check('(j) case 148: the cone with kappa = 9/10 is not invariant (negative control)', ~C2.ok && info2.images_outside > 0, ...
      sprintf('%d images outside', info2.images_outside));
T = jsondecode(fileread(fh)); T.vertices{3} = {'9/10', '0', '1/10'};
rejected('(j) case 148: replay_cone rejects the tampered cone (a generator moved to kappa = 9/10)', replay_cone(T));
[Pv, pinfo] = orbit_polygon(load_pair('test', 'case77'), 'A', 0.2, 1e-3);
[C, info] = certify_cone(Pv, load_pair('test', 'case77'), 'A', '0.2', struct('name', 'case77', 'kind', 'polygon', 'require_positive', true));
fh = fullfile(outdir, 'polygon77.json'); write_text(fh, jsonencode(C));
check('(j) case 77: invariant polygon with 21 vertices, coordinates >= 0.0796 (paper)', C.ok && C.n_vertices == 21 && ...
      min(C.min_coordinates) >= 0.0796, sprintf('%d vertices, min %.4f', C.n_vertices, min(C.min_coordinates)));
check('(j) case 77: replay_cone VALID', replay_cone(fh));

%% (k) Algorithm (L)
fprintf('\n---- (k) Algorithm (L) and the tangent bound\n');
Lg = paper_levels('glp', 'planar_0.1');
G = glp_algorithm_l(Ap, 0.1, struct('nu', Lg.nu));
check('(k) Algorithm (L), planar, tau = 0.1, at the paper''s shift: alpha = -0.215660 (paper)', strcmp(G.alpha_str, '-0.215660'), ...
      sprintf('%s, beta %.6f, |V| = %d', G.alpha_str, G.beta_l, G.nV));
[a, af] = tangent_bound(double(cert_vertices(rp.lower_cert)), Ap);
check('(k) the tangent bound on our planar polytope does not improve on its level (paper)', ...
      double(a) <= rp.search.lambda_star + 1e-6, sprintf('alpha %.9f, level %.9f', double(a), rp.search.lambda_star));

%% (l) Table 1 row and Perron seeds
fprintf('\n---- (l) Table 1 row and Obs. 4.8 seeds\n');
Af = cellfun(@str2double, As, 'UniformOutput', false);
[~, bv] = brute_words(cellfun(@(M) eye(2) + 0.1 * M, Af, 'UniformOutput', false), 16);
[~, rv] = brute_words(cellfun(@(M) inv(eye(2) - 0.1 * M), Af, 'UniformOutput', false), 16);
check('(l) Table 1, symmetric pair, h = 0.1: g_est 0.149874, G_up 0.250126', ...
      strcmp(sprintf('%.6f', (exp(min(bv)) - 1) / 0.1), '0.149874') && strcmp(sprintf('%.6f', (1 - exp(-min(rv))) / 0.1), '0.250126'));
Aw = load_pair('switching', 'sw5_4'); D = exact_data(Aw, '0.1');
B = cellfun(@(M) eye(5) + 0.1 * M, D.Af, 'UniformOutput', false);
w = [ones(1, 24), 2 * ones(1, 7)]; S = perron_seeds(B, w);
L = paper_levels('switching', 'sw5_4_0.1');
out = run_lower_row(Aw, '0.1', struct('name', 'sw5_4', 'seeds', S, 'level', L.t_build, 'nmax', 20000, ...
                    'file', fullfile(outdir, 'sw5_4_0.1_best.json')));
check('(l) switching pair 4, h = 0.1: Perron-seeded closure certifies -0.460394 (paper)', out.closed && out.replay_ok && ...
      strcmp(out.lambda, '-0.460394'), sprintf('%s, |V| = %d', out.lambda, out.nV));

%% (n) other components: three-mode schedules, a driver in search mode
fprintf('\n---- (n) three-mode schedule, search mode\n');
Ab = load_pair('bcv', 'bcv13');
[R, json] = certify_best_upper(Ab, struct('id', 'bcv13', 'chatter', false, ...
                                          'schedule', struct('modes', [3 2 1], 'durations', {{'30', '30', '30'}})));
fh = fullfile(outdir, 'bcv13_schedule.json'); write_text(fh, json);
check('(n) Ex. 4.6: three-mode schedule certified lambda_check <= -0.299986; replay VALID', strcmp(R.printed, '-0.299986') && ...
      replay_best_upper(fh), R.printed);
cs = fullfile(outdir, 'search'); 
table2_test_pairs('mode', 'search', 'rows', {'planar_0.1'}, 'part', 'lower', 'cdir', cs, 'pmax', 2000, 'maxprobes', 30);
C = jsondecode(fileread(fullfile(cs, 'planar_0.1_best.json')));
check('(n) table2_test_pairs in search mode (pmax 2000, 30 probes): planar h = 0.1 within 2e-6 of -0.272230', ...
      abs(str2double(C.lambda_lower) + 0.272230) <= 2e-6, sprintf('%s, %d probes, |V| = %d', C.lambda_lower, ...
      C.meta.probes, numel(C.vertices)));

%% (f) Python replays of the MATLAB certificates (optional)
fprintf('\n---- (f) MATLAB certificates, Python replay (optional)\n');
if isfield(opts, 'python') && ~opts.python
    skip('(f) Python replays', 'disabled');
else
    pyd = fullfile(root, 'independent_replay', 'python');
    groups = {'replay_lower.py', {'planar_0.1_best.json', 'case1_0.1_best.json', 'sharp_0.1_best.json', ...
                                  'planar_0.05_transfer.json', 'planar_0.1_table.json', 'sw5_4_0.1_best.json', ...
                                  fullfile('search', 'planar_0.1_best.json')}; ...
              'replay_upper.py', {'planar_0.1.json', 'case1_0.1.json'}; ...
              'replay_best_upper.py', {'planar_best_upper.json', 'case1440_best_upper.json', 'bcv13_schedule.json'}; ...
              'replay_cones.py', {'cone148.json', 'polygon77.json'}};
    for g = 1:size(groups, 1)
        files = cellfun(@(f) shq(fullfile(outdir, f)), groups{g, 2}, 'UniformOutput', false);
        pyexe = getopt(opts, 'python_executable', 'python3');
        [st, out] = system(sprintf('cd %s && %s -I -B %s %s 2>&1', shq(pyd), shq(pyexe), groups{g, 1}, strjoin(files, ' ')));
        if contains(out, 'No module named') || contains(out, 'command not found')
            skip(sprintf('(f) %s', groups{g, 1}), strtrim(out)); continue
        end
        check(sprintf('(f) %s accepts %d MATLAB certificates', groups{g, 1}, numel(files)), st == 0, ...
              strtrim(regexprep(out(max(1, end - 200):end), '\s+', ' ')));
    end
end

%% (m) no material of the companion paper
fprintf('\n---- (m) B only\n');
num = ['14' '65'];                       % built from pieces, so that this file does not match itself
files = dir(fullfile(root, '**', '*'));
files = files(~[files.isdir]);
txt = files(endsWith({files.name}, {'.m', '.py', '.json', '.md', '.cff', '.txt', '.sh', '.tex'}) | strcmp({files.name}, 'LICENSE'));
hits = {}; ids = {};
for f = txt.'
    p = fullfile(f.folder, f.name); s = fileread(p);
    if ~isempty(regexpi(s, ['case[_ ]?' num '|' num '_|_' num '|' num '[a-z]'], 'once')), ids{end + 1} = f.name; end %#ok<AGROW>
    incert = contains(p, [filesep 'certificates' filesep]) || contains(p, [filesep 'fixtures' filesep]);
    if ~incert && ~isempty(regexp(s, ['(?<![0-9])' num '(?![0-9])'], 'once')), hits{end + 1} = f.name; end %#ok<AGROW>
end
check(sprintf('(m) no identifier of case %s (case%s, %s_, ...) in any text file', num, num, num), isempty(ids), strjoin(ids, ', '));
check(sprintf('(m) no %s outside digit strings in sources, data and documents (certificates excluded: vertex indices)', num), ...
      isempty(hits), strjoin(hits, ', '));
[~, known] = cellfun(@(g) load_pair(g), {'test', 'scaling', 'switching', 'bcv', 'switching_journal', 'scaling_journal'}, ...
                     'UniformOutput', false);
known = [known{:}]; bad = {};
for f = dir(fullfile(root, 'certificates', '**', '*.json')).'
    C = jsondecode(fileread(fullfile(f.folder, f.name)));
    if isfield(C, 'case'), nm = C.('case'); elseif isfield(C, 'id'), nm = C.id; else, continue, end
    if ~ismember(nm, known), bad{end + 1} = f.name; end %#ok<AGROW>
end
check('(m) every certificate refers to a pair of data/ (the pairs of this paper)', isempty(bad), strjoin(bad, ', '));

%% proof prerequisites and malformed-certificate regression tests
nvalidation = run_validation_tests();
check('proof-prerequisite regression suite', nvalidation > 0, sprintf('%d checks passed', nvalidation));

%% summary
st = results(:, 2);
fprintf('\n==== summary: %d PASS, %d FAIL, %d SKIP, %.0f s\n', sum(strcmp(st, 'PASS')), sum(strcmp(st, 'FAIL')), ...
        sum(strcmp(st, 'SKIP')), toc(T0));
for k = find(strcmp(st, 'FAIL')).'
    fprintf('  FAILED: %s %s\n', results{k, 1}, results{k, 3});
end
if any(strcmp(st, 'FAIL'))
    error('run_tests:failed', '%d check(s) failed', sum(strcmp(st, 'FAIL')));
end
end


function tf = is_rotation(a, b)
tf = numel(a) == numel(b) && contains([b b], a);
end


function s = shq(p)
s = ['''' strrep(p, '''', '''\''''') ''''];
end
