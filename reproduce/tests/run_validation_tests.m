function count = run_validation_tests()
% RUN_VALIDATION_TESTS  Proof prerequisites, malformed witnesses and scalar edge cases.
count = 0;
must_error(@() exact_data({}, '0.1'), 'exact_data:input');
must_error(@() exact_data({Inf}, '0.1'), 'exact_data:input');
must_error(@() exact_data({1 + 1i}, '0.1'), 'exact_data:input');
must_error(@() exact_data({0}, 1 + 1i), 'exact_data:input');
must_error(@() parse_rationals("1/0"), 'parse_rationals:parse');
must_error(@() exact_data({[0 -1; 0 0]}, '0.1'), 'exact_data:metzler');
must_error(@() exact_data({0}, '-0.1'), 'exact_data:input');
must_error(@() exact_data({0}, 'Inf'), 'exact_data:input');
must_error(@() certified_bracket({[-1 20; 20 -1]}, '0.1', struct('verbose', false)), 'check_euler_step:implicit');
must_error(@() certify_upper_exact([], exact_data({0}, '0.1')), 'certify_upper_exact:word');
must_error(@() certify_lower_exact(1, exact_data({-1}, '1')), 'check_euler_step:explicit');
must_error(@() certify_lower_exact(0, exact_data({0}, '0.1')), 'certify_lower_exact:vertices');
must_error(@() certify_lower_exact(-1, exact_data({0}, '0.1')), 'certify_lower_exact:vertices');
must_error(@() certify_best_upper({0, 1}, struct('alpha', 2)), 'certify_best_upper:alpha');
must_error(@() certify_best_upper({0, 1}, struct('chatter', false)), 'certify_best_upper:empty');
sch = struct('modes', 1, 'durations', [1 1]);
must_error(@() certify_best_upper({0, 1}, struct('schedule', sch)), 'certify_best_upper:schedule');
sch.durations = -1;
must_error(@() certify_best_upper({0, 1}, struct('schedule', sch)), 'certify_best_upper:schedule');
must_error(@() exp_upper(sym(0), sym(-1)), 'exp_upper:time');
must_error(@() log_upper(sym(0)), 'log_upper:input');
must_error(@() certify_cone(ones(3), {eye(3)}, 'A', '0.1'), 'certify_cone:generators');

% One mode and a column word must work; modes above nine use an integer array in JSON.
A = {zeros(2)}; D = exact_data(A, '0.1');
C = certify_upper_exact([1; 1], D);
check(replay_upper(C, A), 'single-mode column word');
Cs = certify_upper_exact(1, exact_data({0}, '0.1'));
check(replay_upper(Cs, {0}), 'scalar single-mode certificate');
A10 = repmat(A, 1, 10);
C10 = certify_upper_exact(10, exact_data(A10, '0.1'));
check(isnumeric(C10.word) && C10.word == 10 && replay_upper(C10, A10), 'ten-mode integer word');
for badh = {'0', '-0.1'}
    T = C; T.h = badh{1}; check(~replay_upper(T, A), 'nonpositive upper step rejected');
end
T = C; T.versions.decimal.z = {}; check(~replay_upper(T, A), 'missing mode witness rejected');
T = C; T.versions.decimal.chain = {}; check(~replay_upper(T, A), 'empty chain rejected');
[w, val] = word_search({eye(2)}, 0.1, struct('L', 1, 'kmax', 1));
check(isequal(w, 1) && val == 0, 'single-mode search with no local moves');

% A nonconvex mixture can give a false bound: 2*0 + (1-2)*I = -I, but lambda=0.
rd = struct('alpha', '2', 'x', {{'1', '1'}});
R = struct('U_c', '-1', 'U_s', [], 'U_best', '-1', 'printed', '-1', ...
           'readings', struct('decimal', rd, 'binary64', rd));
check(~replay_best_upper(R, {zeros(2), eye(2)}), 'false affine-mixture bound rejected');
% A truncated product with the full duration sum gives a false 0.75 bound for lambda=1.
rd = struct('modes', 0, 'durations', {{'1', '1'}}, 'x_sched', {{'1', '1'}});
R = struct('U_c', [], 'U_s', '3/4', 'U_best', '3/4', 'printed', '0.75', ...
           'readings', struct('decimal', rd, 'binary64', rd));
check(~replay_best_upper(R, {eye(2), eye(2)}), 'mismatched schedule rejected');
R = certify_best_upper({[0.3 0; 0 0.3], eye(2)}, struct('alpha', 1));
check(replay_best_upper(R, {[0.3 0; 0 0.3], eye(2)}), 'numeric decimal data replay consistently');
R.printed = '-1'; check(~replay_best_upper(R, {[0.3 0; 0 0.3], eye(2)}), 'printed mixture bound checked');
F = struct('orientations', []); check(~replay_faces(F, {eye(2), eye(2)}), 'empty face proof rejected');
check(~replay_all({'__missing_validation_fixture__'}), 'empty replay selection rejected');
fprintf('validation: %d PASS\n', count);

    function check(ok, label)
        assert(ok, 'run_validation_tests:failed', '%s', label);
        count = count + 1;
        fprintf('[PASS] validation: %s\n', label);
    end

    function must_error(f, id)
        caught = false;
        try
            f();
        catch err
            if ~strcmp(err.identifier, id), rethrow(err); end
            caught = true;
        end
        check(caught, [id ' rejected']);
    end
end
