%% This script can be used to rebuild Table 2 and the initial Table 3 closures of paper B.

% The stronger extended-search bounds marked + in the revised Table 3 are checked by the replays in reproduce/;
% rebuilding those searches needs table3_extended (and, for some seeds, the external GM run outputs).
% Each lower bound is rebuilt by one closure at the level recorded in paper_data.m and then checked in exact
% arithmetic; each upper bound is re-certified from the recorded word (Table 2) or schedule (Table 3).
% quick = true skips the long rows: case 1 at h = 0.025 (about an hour), and the Table 3 rows with d = 10 or
% h = 0.05. The quick run takes about 20 minutes on one core; the full run several hours.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));

quick = true;
P = paper_data();

%% Table 2: planar pair, cases 1, 1440 and 148
fprintf('\nTable 2: certified Euler brackets for the test pairs\n');
fprintf('%-9s %-7s %11s %11s   %11s %11s   %s\n', 'pair', 'h', 'lower', '(paper)', 'upper', '(paper)', 'match');
polytopes = struct();                       % certified vertex sets, reused for the transfers to smaller steps
ok2 = true(1, 0);
for r = P.table2.'
    row = matlab.lang.makeValidName([r.pair '_' r.h]);
    if quick && strcmp(row, 'case1_0_025'), continue, end
    A = P.pairs.(r.pair); D = exact_data(A, r.h); d = D.d;
    if isempty(r.from)
        F = cellfun(@(M) eye(d) + D.hf * M, D.Af, 'UniformOutput', false);
        if strcmp(r.seeds, 'e1e3'), seeds = [1 0; 0 0; 0 1]; else, seeds = eye(d); end
        V = polytope_closure(F, r.level, seeds, Inf, 6000, Inf);
    else                                    % Prop. 4.3(ii): the polytope of the larger step is certified again
        V = polytopes.(matlab.lang.makeValidName([r.pair '_' r.from]));
    end
    CL = certify_lower_exact(V, D, struct('name', r.pair));
    polytopes.(row) = cert_vertices(CL);
    CU = certify_upper_exact(r.upper_word - '0', D, struct('name', r.pair));
    if isempty(r.from)
        match = strcmp(CL.lambda_lower, r.lower);
    else                                    % a transfer may certify more than the bound of the larger step,
        match = str2double(CL.lambda_lower) >= str2double(r.lower);     % which is what the paper prints
    end
    match = match && strcmp(CU.printed, r.upper);
    ok2(end + 1) = match; %#ok<SAGROW>
    fprintf('%-9s %-7s %11s %11s   %11s %11s   %s\n', r.pair, r.h, CL.lambda_lower, r.lower, CU.printed, r.upper, ...
            yesno(match));
end

%% Table 3: random pairs selected by numerical schedule-versus-mixture comparisons
fprintf('\nTable 3 initial closures (before extended refinements): lower bounds for g(h) and schedule upper bounds\n');
fprintf('%-8s %-5s %-3s %11s %11s   %s\n', 'pair', 'h', 'or.', 'g(h) >=', '(recorded)', 'match');
ok3 = true(1, 0);
for r = P.table3.'
    if isempty(r.orientation)
        fprintf('%-8s %-5s: no initial closure recorded; see reproduce/ for extended results.\n', r.pair, r.h);
        continue
    end
    d = size(P.pairs.(r.pair){1}, 1);
    if quick && (d == 10 || strcmp(r.h, '0.05')), continue, end
    A = P.pairs.(r.pair); D = exact_data(A, r.h);
    B = cellfun(@(M) eye(d) + D.hf * M, D.Af, 'UniformOutput', false);
    w = [ones(1, r.word(1)), 2 * ones(1, r.word(2))];
    if r.orientation == 'T'                 % Prop. 4.3(iii): the transposed family has the same lower spectral radius
        F = cellfun(@(M) M.', B, 'UniformOutput', false); seeds = perron_seeds(F, fliplr(w));
    else
        F = B; seeds = perron_seeds(B, w);
    end
    V = polytope_closure(F, r.level, seeds, Inf, 20000, Inf);
    CL = certify_lower_exact(V, D, struct('name', r.pair, 'orientation', r.orientation));
    match = strcmp(CL.lambda_lower, r.lower);
    ok3(end + 1) = match; %#ok<SAGROW>
    fprintf('%-8s %-5s %-3s %11s %11s   %s\n', r.pair, r.h, r.orientation, CL.lambda_lower, r.lower, yesno(match));
end
fprintf('\n%-8s %11s %11s   %s\n', 'pair', 'upper', '(paper)', 'match');
for s = P.schedules.'
    if quick && startsWith(s.pair, 'sw10'), continue, end
    R = certify_best_upper(P.pairs.(s.pair), struct('chatter', false, ...
        'schedule', struct('modes', [1 2], 'durations', {{s.t1, s.t2}})));
    match = strcmp(R.printed, s.upper);
    ok3(end + 1) = match; %#ok<SAGROW>
    fprintf('%-8s %11s %11s   %s\n', s.pair, R.printed, s.upper, yesno(match));
end

fprintf('\n%d of %d bounds match the paper.\n', nnz([ok2 ok3]), numel([ok2 ok3]));


function s = yesno(b)
if b, s = 'yes'; else, s = 'NO'; end
end
