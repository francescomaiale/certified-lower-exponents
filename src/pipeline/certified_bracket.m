function res = certified_bracket(A, h, opts)
% CERTIFIED_BRACKET  Certified bracket lower <= lambda_check <= upper for the lower Lyapunov exponent of the
% positive switched system x' = A_sigma x, from the explicit and implicit Euler families at step h (paper B, Sec. 4).
%
%   res = certified_bracket(A, h, opts)
%     A     cell array {A_1, ..., A_m} of d x d Metzler matrices: numeric, or decimal-string matrices
%           (paper_data, or load_pair in reproduce/).  Both readings (exact decimals and exact binary64 values) are certified.
%     h     positive step; its decimal is exact (0.1 -> 1/10). Requires h max|a_jj| < 1 and h eta(A_i) < 1.
%     opts  (all optional)
%       name        label stored in the certificates ('pair')
%       source      'source' field of the upper certificate ('matlab'; 'table2' lets the Python replay of
%                   reproduce/ find the data)
%       orientation 'N' (B_i = I + h A_i, default) or 'T' (transposes, same lsr) for the lower certificate
%       seeds       d x s starting vertices of the closure (default eye(d))
%       probe       time cap of one closure, s (10);  budget  time budget of the bisection, s (60)
%       nmax        vertex cap of one closure (6000); tol     bisection stops when hi - lo <= tol*h (2e-7)
%       pmax        cap on the vertices processed by one closure (Inf); set it (with generous probe/budget) for
%                   a search that does not depend on the speed of the machine
%       t_start     a level of lsr(B_h) to probe first, e.g. a lower estimate from another method (none)
%       word        certify this word instead of searching (row vector of modes)
%       search      options for word_search (L, n, Q, ntop, sweeps, kmax)
%       versions    data readings ({'decimal', 'binary64'})
%       outdir      if set, write the certificates there as <name>_<h>_best.json (lower) and <name>_<h>.json
%                   (upper), e.g. for the independent replays of reproduce/
%       verbose     print progress (true)
%   res fields
%     lower, upper              certified six-decimal bounds (char), lower_value, upper_value (double)
%     lower_cert, upper_cert    certificates (decoded JSON structs, in the format read by the replays of reproduce/)
%     lower_json, upper_json    their JSON text
%     search                    t_triv, t_up, chatter c, lsr product bound, bisection probes, word, U_float, ...
%     timing                    seconds per stage
%
%   Lower bound (prop:polytope): a finite V in R^d_+ with B_i v >= t V lambda_{i,v} (lambda >= 0, sum 1) for
%   every vertex and mode certifies lsr(B_h) >= t, hence lambda_check >= (t - 1)/h.  V is the closure of the
%   seeds under v -> B_i v / t, with t chosen by bisection between a level that closes and
%   t_up = min(1 + h c, rho(B_w)^(1/|w|) of the best explicit word, 1 + h * upper); certify_lower_exact checks
%   it exactly and prints floor_6((t* - 1)/h).
%   Upper bound (prop:upper): a word w of length k with a chain x_0 > 0, (I - h A_{w_j}) x_j >= x_{j-1}, and
%   r^k >= max(x_k ./ x_0) certifies lambda_check <= (1 - 1/r)/h, printed rounded up.

if nargin < 3, opts = struct(); end
T0 = tic;
name    = getopt(opts, 'name', 'pair');
source  = getopt(opts, 'source', 'matlab');
orient  = getopt(opts, 'orientation', 'N');
versions = getopt(opts, 'versions', {'decimal', 'binary64'});
verbose = getopt(opts, 'verbose', true);
sopts   = getopt(opts, 'search', struct());
outdir  = getopt(opts, 'outdir', '');
say = @(varargin) verbose && fprintf(varargin{:});

D = exact_data(A, h, versions);
d = D.d; m = D.m; hf = D.hf; I = eye(d);
check_euler_step(D, 'explicit');
check_euler_step(D, 'implicit');
if ~ismember(string(orient), ["N", "T"]), error('certified_bracket:orientation', 'orientation must be N or T'); end
say('== %s, d = %d, m = %d, h = %s (exact %s); data readings: %s\n', name, d, m, D.h_str, char(D.hq), ...
    strjoin(versions, ', '));
tim = struct();

% ---------------------------------------------------------------- upper bound (implicit Euler, prop:upper)
t1 = tic;
R = cell(1, m);
for i = 1:m, R{i} = inv(I - hf * D.Af{i}); end
w = getopt(opts, 'word', []);
if isempty(w)
    [w, v] = word_search(R, hf, sopts);
else
    w = w(:).';
    v = word_value(R, w);
end
U_float = (1 - exp(-v)) / hf;
tim.upper_search = toc(t1);
say('upper: best word length %d, U_float = %+.9f (%.1fs)\n', numel(w), U_float, tim.upper_search);
t1 = tic;
[CU, jsonU, infoU] = certify_upper_exact(w, D, struct('name', name, 'source', source, 'U_float_search', U_float));
tim.upper_exact = toc(t1);
say('upper: CERTIFIED lambda_check <= %s (k = %d; exact U = %+.9f; %.1fs)\n', CU.printed, numel(w), ...
    max(double(infoU.U)), tim.upper_exact);

% ---------------------------------------------------------------- lower bound (explicit Euler, prop:polytope)
t1 = tic;
F = cell(1, m);
for i = 1:m
    F{i} = I + hf * D.Af{i};
    if orient == 'T', F{i} = F{i}.'; end
end
c = chattering_value(D.Af);
eopts = sopts; eopts.ntop = getopt(sopts, 'ntop', 2);
[wB, vB] = word_search(F, hf, eopts);
lsr_up = exp(vB);
t_up = min([1 + hf * c, lsr_up, 1 + hf * str2double(CU.printed)]);
t_start = getopt(opts, 't_start', []);
tim.lower_bounds = toc(t1);
say('lower: mixture/mode estimate = %+.9f; best explicit word (k = %d): g(h) <= %+.9f; search up to %+.9f\n', ...
    c, numel(wB), (lsr_up - 1) / hf, (t_up - 1) / hf);
popts = struct('seeds', getopt(opts, 'seeds', eye(d)), 'probe', getopt(opts, 'probe', 10), ...
               'budget', getopt(opts, 'budget', 60), 'nmax', getopt(opts, 'nmax', 6000), ...
               'pmax', getopt(opts, 'pmax', Inf), ...
               'maxprobes', getopt(opts, 'maxprobes', Inf), ...
               'tol', getopt(opts, 'tol', 2e-7), 't_start', t_start, 'verbose', verbose);
S = polytope_search(F, hf, t_up, popts);
tim.lower_search = S.secs;
t1 = tic;
meta = struct('t_build', S.t_build, 't_fail', S.t_fail, 't_triv', S.t_triv, 't_up', S.t_up, ...
              'probes', numel(S.probes), 'search_secs', round(S.secs, 1));
[CL, jsonL, infoL] = certify_lower_exact(S.V, D, struct('name', name, 'kind', 'best', 'orientation', orient, ...
                                                       'meta', meta));
tim.lower_exact = toc(t1);
if S.t_fail < S.t_up, endnote = 'next probe failed at'; else, endnote = 'bisection reached the upper end'; end
say('lower: CERTIFIED lambda_check >= %s (exact (t*-1)/h = %+.9f; |V| = %d; %s %+.9f; %.1fs)\n', ...
    CL.lambda_lower, double((infoL.t_star - 1) / D.hq), infoL.nV, endnote, (S.t_fail - 1) / hf, tim.lower_exact);

% ---------------------------------------------------------------- output
res = struct('name', name, 'h', char(D.h_str));
res.lower = CL.lambda_lower; res.lower_value = str2double(CL.lambda_lower);
res.upper = CU.printed; res.upper_value = str2double(CU.printed);
res.lower_cert = CL; res.upper_cert = CU; res.lower_json = jsonL; res.upper_json = jsonU;
res.search = struct('t_triv', S.t_triv, 't_up', t_up, 'chatter', c, 'lsr_product_bound', lsr_up, ...
                    'explicit_word', sprintf('%d', wB), 't_start', t_start, 't_build', S.t_build, ...
                    't_fail', S.t_fail, 'probes', S.probes, 'nV', infoL.nV, 'word', sprintf('%d', w), ...
                    'k', numel(w), 'U_float', U_float, 't_star', double(infoL.t_star), ...
                    'lambda_star', double((infoL.t_star - 1) / D.hq));
if ~isempty(outdir)
    if ~isfolder(outdir), mkdir(outdir); end
    write_text(fullfile(outdir, sprintf('%s_%s_best.json', name, D.h_str)), jsonL);
    write_text(fullfile(outdir, sprintf('%s_%s.json', name, D.h_str)), jsonU);
end
tim.total = toc(T0);
res.timing = tim;
say('== %s, h = %s: %s <= lambda_check <= %s (width %.6f), total %.1fs\n', name, D.h_str, res.lower, res.upper, ...
    res.upper_value - res.lower_value, tim.total);
end


function v = word_value(F, w)
P = eye(size(F{1}, 1)); lg = 0;
for a = w
    P = F{a} * P; s = max(abs(P(:))); P = P / s; lg = lg + log(s);
end
v = (log(max(abs(eig(P)))) + lg) / numel(w);
end


function write_text(path, txt)
fid = fopen(path, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end


function v = getopt(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end
