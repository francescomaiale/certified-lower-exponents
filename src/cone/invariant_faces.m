function F = invariant_faces(A, orient)
% INVARIANT_FACES  Invariant faces of K = R^d_+ with exact exponent bounds, against a certified bound for c
% (Computation 3.9).
%   F = invariant_faces(A, orient)
%     A       cell {A_1, A_2} of decimal-string or numeric matrices; the binary64 values are used (as in the paper)
%     orient  'A' (family as given) or 'T' (transposed family, same g and lambda_check by prop:cone(iv))
%   For every proper nonempty J, S(J) is the set of modes with (A_i)_{kj} = 0 for k outside J, j in J.  When
%   S(J) = {i}, lambda_J = eta(A_i[J,J]) is bounded below exactly by Collatz-Wielandt: with a rational x >= 0,
%   x ~= 0 (a rationalized Perron vector), M x >= L x entrywise gives eta(M) >= L, L = min_{x_r > 0} (Mx)_r/x_r.
%   The chattering value c = min_alpha eta(S_alpha) >= lambda_check is bounded above by Collatz-Wielandt at a rational
%   alpha near the minimizer (prop:basic(iii)).  A face is certified slow if its lower bound exceeds the bound for c.
%   F has fields orientation, alpha, x_c, c_upper (strings), c_upper_float, faces (struct array: J, S,
%   eta_lower, eta_lower_float, x, slow) and all_slow (true if every invariant proper face is certified slow; then
%   lambda_min = lambda_check by thm:faceformula).  replay_faces re-checks it independently.

if nargin < 2, orient = 'A'; end
if ~iscell(A) || numel(A) ~= 2 || ~ismember(string(orient), ["A", "T"])
    error('invariant_faces:input', 'provide a matrix pair and orientation A or T');
end
D = exact_data(A, '1', {'binary64'});
m = numel(A); d = size(A{1}, 1);
Aq = cell(1, m); Af = cell(1, m);
for i = 1:m
    Af{i} = D.Af{i}; Aq{i} = D.ver(1).A{i};
    if orient == 'T', Aq{i} = Aq{i}.'; Af{i} = Af{i}.'; end
end
% certified upper bound for c (pairs)
eta = @(a) max(real(eig(a * Af{1} + (1 - a) * Af{2})));
[af, ~] = fminbnd(eta, 0, 1, optimset('TolX', 1e-13));
if eta(0) < eta(af), af = 0; end
if eta(1) < eta(af), af = 1; end
a = limit_denominator(af, 1e6);
S = a * Aq{1} + (1 - a) * Aq{2};
[EV, ev] = eig(double(S)); [~, k] = max(real(diag(ev)));
v = abs(real(EV(:, k))); v = max(v / max(v), 1e-9);
x = sym(zeros(d, 1));
for r = 1:d, x(r) = limit_denominator(v(r), 1e12); end
UC = maxq((S * x) ./ x);
F = struct('orientation', orient, 'alpha', char(a), 'x_c', {cellstr(string(x)).'}, 'c_upper', char(UC), ...
           'c_upper_float', double(UC), 'faces', struct('J', {}, 'S', {}, 'eta_lower', {}, 'eta_lower_float', {}, ...
           'x', {}, 'slow', {}));
allslow = true;
for sz = 1:d - 1
    Js = nchoosek(1:d, sz);
    for q = 1:size(Js, 1)
        J = Js(q, :); out = setdiff(1:d, J);
        Si = [];
        for i = 1:m
            if all(all(Af{i}(out, J) == 0)), Si(end + 1) = i; end %#ok<AGROW>
        end
        if isempty(Si), continue, end
        f = struct('J', J, 'S', Si, 'eta_lower', '', 'eta_lower_float', NaN, 'x', {{}}, 'slow', false);
        if numel(Si) == 1
            M = Aq{Si}(J, J);
            [EV, ev] = eig(double(M)); [~, k] = max(real(diag(ev)));
            xv = abs(real(EV(:, k))); xv(xv < 1e-14 * max(xv)) = 0;
            xj = sym(zeros(numel(J), 1));
            for r = 1:numel(J), xj(r) = limit_denominator(xv(r), 1e12); end
            p = double(sign(xj)) > 0; Mx = M * xj;
            lo = minq(Mx(p) ./ xj(p));
            f.eta_lower = char(lo); f.eta_lower_float = double(lo); f.x = cellstr(string(xj)).';
            f.slow = double(sign(lo - UC)) > 0;
        end
        allslow = allslow && f.slow;
        F.faces(end + 1) = f;
    end
end
F.all_slow = allslow;
end


function m = maxq(v)
v = v(:); m = v(1);
for k = 2:numel(v), if double(sign(v(k) - m)) > 0, m = v(k); end, end
end


function m = minq(v)
v = v(:); m = v(1);
for k = 2:numel(v), if double(sign(v(k) - m)) < 0, m = v(k); end, end
end
