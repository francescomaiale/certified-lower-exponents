function out = run_upper_row(A, h, spec)
% RUN_UPPER_ROW  One certified implicit Euler upper bound (thm:implicit(i), prop:upper): word search, chain
% certificate, certificate file and independent replay.
%   out = run_upper_row(A, h, spec)
%     A, h  data (decimal strings) and step (decimal string)
%     spec  name (pair name), word (certify this word, row of modes; default: word_search on R_h),
%           search (word_search options), file (certificate path), replay (true), source ('matlab')
%   out: name, h, upper (char, rounded up), U_exact (double), U_float (float value of the word), word, k,
%   replay_ok, secs_search, secs_exact, secs_replay, file.

T0 = tic;
D = exact_data(A, h);
check_euler_step(D, 'implicit');
d = D.d; m = D.m;
R = cell(1, m);
for i = 1:m, R{i} = inv(eye(d) - D.hf * D.Af{i}); end
w = getopt(spec, 'word', []);
if ischar(w) || isstring(w), w = double(char(w)) - 48; end
if isempty(w)
    [w, v] = word_search(R, D.hf, getopt(spec, 'search', struct()));
else
    w = w(:).';
    v = word_value(R, w);
end
Uf = (1 - exp(-v)) / D.hf;
out = struct('name', spec.name, 'h', char(D.h_str), 'upper', '', 'U_exact', NaN, 'U_float', Uf, ...
             'word', sprintf('%d', w), 'k', numel(w), 'replay_ok', false, 'file', '', 'secs_search', toc(T0), ...
             'secs_exact', 0, 'secs_replay', 0);
T1 = tic;
[C, json, info] = certify_upper_exact(w, D, struct('name', spec.name, 'source', getopt(spec, 'source', 'matlab'), ...
                                                   'U_float_search', Uf));
out.secs_exact = toc(T1);
out.upper = C.printed; out.U_exact = max(double(info.U));
file = getopt(spec, 'file', '');
if ~isempty(file), write_text(file, json); out.file = file; end
if getopt(spec, 'replay', true)
    T2 = tic; [ok, msgs] = replay_upper(C);
    out.replay_ok = ok; out.replay_msgs = msgs; out.secs_replay = toc(T2);
end
end
