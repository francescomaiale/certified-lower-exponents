function [R, json] = certify_best_upper(A, opts)
% CERTIFY_BEST_UPPER  Certified chattering and schedule upper bounds for lambda_check (Sec. 4, Construction).
%   [R, json] = certify_best_upper(A, opts)
%     A     cell array {A_1, ..., A_m} of decimal-string (or numeric) Metzler matrices
%     opts  id ('pair'); versions ({'decimal', 'binary64'});
%           chatter (true for m = 2): certify c = min_alpha eta(S_alpha) <= max_r (S_alpha x)_r / x_r at a
%                   rational alpha near the floating-point minimizer (Collatz-Wielandt, prop:basic(iii));
%           alpha   floating-point minimizer alpha_1 (default: chattering_value);
%           schedule struct with fields modes (row of mode indices 1..m, in the order they act) and durations
%                   (row of positive rationals: sym, 'p/q' strings or doubles, which are rounded to the closest
%                   rational with denominator <= 10^6), or [] for none.  It is certified as
%                   lambda_check <= log rho(Phi) / T, Phi = e^{t_k A_{i_k}} ... e^{t_1 A_{i_1}}, T = sum t_j, with
%                   entrywise upper enclosures of each exponential (exp_upper), Collatz-Wielandt for rho(Phi) and
%                   an artanh series for the logarithm (log_upper), rounded up to 1e-12.
%           meta    struct stored as is.
%   R is the certificate (decoded JSON, in the format read by the replays of reproduce/: readings.<name>.alpha/x/U_c
%   and, for a
%   schedule, t1/t2 (two segments) or modes/durations, x_sched, rho_upper, U_s; U_c, U_s and U_best are the
%   maxima over the readings).  replay_best_upper re-checks it independently (reproduce/replay; uniformization, no logarithms).

if nargin < 2, opts = struct(); end
T0 = tic;
id = getopt(opts, 'id', 'pair');
versions = getopt(opts, 'versions', {'decimal', 'binary64'});
D = exact_data(A, '1', versions);
m = D.m; d = D.d;
do_chat = getopt(opts, 'chatter', m == 2);
sch = getopt(opts, 'schedule', []);
if ~do_chat && isempty(sch), error('certify_best_upper:empty', 'request a mixture bound or a schedule'); end
R = struct('id', id, 'd', d);
if do_chat
    if m ~= 2, error('certify_best_upper:chatter', 'the chattering certificate is implemented for pairs'); end
    [cf, af] = chattering_value(D.Af);
    af = getopt(opts, 'alpha', af);
    if ~isnumeric(af) || ~isscalar(af) || ~isreal(af) || ~isfinite(af) || af < 0 || af > 1
        error('certify_best_upper:alpha', 'alpha must be a finite real number in [0,1]');
    end
    R.c_float = cf; R.alpha_float = af;
    a = limit_denominator(af, 1e9);
end
if ~isempty(sch)
    modes = sch.modes(:).'; dur = sch.durations;
    if isempty(modes) || ~isnumeric(modes) || ~isreal(modes) || ...
            any(~isfinite(modes) | modes ~= round(modes) | modes < 1 | modes > m) || numel(modes) ~= numel(dur)
        error('certify_best_upper:schedule', 'provide one positive duration per integer mode index 1..m');
    end
    if ~iscell(dur) && ~isa(dur, 'sym'), dur = num2cell(dur); end
    tq = sym(zeros(1, numel(modes)));
    for k = 1:numel(modes)
        if iscell(dur), x = dur{k}; else, x = dur(k); end
        if isa(x, 'sym'), tq(k) = x;
        elseif ischar(x) || isstring(x), tq(k) = parse_rationals(x);
        else, tq(k) = limit_denominator(x, 1e6);
        end
    end
    if ~all(isAlways(tq > 0)) || any(~isfinite(double(tq)))
        error('certify_best_upper:schedule', 'durations must be finite and positive after rationalization');
    end
    R.schedule_float = getopt(sch, 'rate_float', NaN);
    R.modes = modes - 1; R.durations = cellstr(string(tq));
    if numel(modes) == 2 && isequal(modes, [1 2])
        R.t1_float = double(tq(1)); R.t2_float = double(tq(2));
    end
end
R.readings = struct();
Ucs = sym([]); Uss = sym([]);
for v = 1:numel(D.ver)
    As = D.ver(v).A; rd = struct();
    if do_chat
        S = a * As{1} + (1 - a) * As{2};
        [EV, ev] = eig(double(S)); [~, k] = max(real(diag(ev)));
        pv = abs(real(EV(:, k))); pv = max(pv / max(pv), 1e-12);
        x = sym(pv, 'f');
        Uc = xmax_sym((S * x) ./ x);
        rd.alpha = char(a); rd.x = cellstr(string(x)).'; rd.U_c = char(Uc);
        Ucs(end + 1) = Uc; %#ok<AGROW>
    end
    if ~isempty(sch)
        Phi = sym(eye(d)); G = sym(2)^90;
        for k = 1:numel(modes)
            Phi = ceil(exp_upper(As{modes(k)}, tq(k)) * Phi * G) / G;
        end
        [EV, ev] = eig(double(Phi)); [~, k] = max(abs(diag(ev)));
        pv = abs(real(EV(:, k))); pv = max(pv / max(pv), 1e-12);
        xs = sym(pv, 'f');
        rho = xmax_sym((Phi * xs) ./ xs);
        U = log_upper(rho) / sum(tq);
        U12 = ceil(U * sym(10)^12) / sym(10)^12;          % rounded up to 1e-12 (slack for the replay)
        if numel(modes) == 2 && isequal(modes, [1 2])
            rd.t1 = char(tq(1)); rd.t2 = char(tq(2));
        end
        rd.modes = modes - 1; rd.durations = cellstr(string(tq));
        rd.x_sched = cellstr(string(xs)).'; rd.rho_upper = char(ceil(rho * G) / G); rd.U_s = char(U12);
        Uss(end + 1) = U12; %#ok<AGROW>
    end
    R.readings.(D.ver(v).name) = rd;
end
best = [];
if do_chat
    Uc = xmax_sym(Ucs); R.U_c = char(Uc); best = Uc; kind = 'chattering';
else
    R.U_c = [];
end
if ~isempty(sch)
    Us = xmax_sym(Uss); R.U_s = char(Us);
    if isempty(best) || double(sign(Us - best)) < 0, best = Us; kind = 'schedule'; end
else
    R.U_s = [];
end
R.U_best = char(best); R.U_best_float = double(best); R.best_kind = kind;
R.printed = fmt_dec(best, 6, 'ceil');
R.meta = getopt(opts, 'meta', struct()); R.meta.implementation = 'matlab';
R.secs = round(toc(T0), 1);
json = jsonencode(R);
R = jsondecode(json);
end


function m = xmax_sym(v)
v = v(:); m = v(1);
for k = 2:numel(v)
    if double(sign(v(k) - m)) > 0, m = v(k); end
end
end
