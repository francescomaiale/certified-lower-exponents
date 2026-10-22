function run_gm_tests()
% RUN_GM_TESTS  Timed-level fallback, portable job outputs and historical seeds.
root = paperB_root(); oldpath = path; oldpwd = pwd;
restore = onCleanup(@() path(oldpath)); %#ok<NASGU>
addpath(fullfile(root, 'gm'), '-begin');
% A gap after level one makes the zero-time cap exercise the fallback branch.
A = [0.8 0.4; 0.1 0.9]; B = [0.9 0.1; 0.4 0.8]; X = [A B];
for adaptive = [false true]
    if adaptive
        [r0, m0] = adaptive_capped(X, [], 1e-6, 2, eye(2), [], false, true, Inf);
        [r1, m1, ~, ~, ~, status] = adaptive_capped(X, [], 1e-6, 20000, eye(2), [], false, true, 0);
    else
        [r0, m0] = lsr_comp_capped(X, [], 1e-6, 2, eye(2), [], false, Inf);
        [r1, m1, ~, status] = lsr_comp_capped(X, [], 1e-6, 20000, eye(2), [], false, 0);
    end
    assert(strcmp(status, 'time') && max(abs(r1-r0)) < 1e-10 && isequal(m1, m0), ...
        'run_gm_tests:timeout', 'A partial level changed the completed-level bounds or counters.');
end
work = tempname; mkdir(work); cleanup = onCleanup(@() rmdir(work, 's')); %#ok<NASGU>
for variant = {'adaptive', 'fixed'}
    for side = 'EI'
        R = gm_fixed_job('test', 'case1', '0.1', side, variant{1}, 'rdir', work, ...
            'cap', 0, 'products', 2, 'certify', false);
        assert(~R.paper_protocol && ~R.certification_requested && R.depth == 1);
        f = fullfile(work, variant{1}, 'test', ['case1_0.1_' side '.json']);
        saved = jsondecode(fileread(f)); assert(saved.products == R.products);
        assert(strcmp(pwd, oldpwd), 'run_gm_tests:cwd', 'Job changed working directory.');
    end
end
% Exercise the exact-certification path, using an easy coordinate lower level.
R = gm_fixed_job('fixed', 'case1', '0.1', 'E', 'fixed', 'rdir', work, ...
    'cap', 0, 'products', 2);
assert(R.certification_requested && R.seed.closed && R.seed.replay_ok && ~R.paper_protocol, ...
    'run_gm_tests:certificate', 'The short GM run failed exact a posteriori certification.');
manifest = jsondecode(fileread(fullfile(root, 'results', 'gm_original', 'MANIFEST.json')));
count = 0;
for k = 1:numel(manifest.files)
    if iscell(manifest.files), F = manifest.files{k}; else, F = manifest.files(k); end
    if ~endsWith(F.path, '.mat'), continue, end
    S = load(fullfile(root, 'results', 'gm_original', F.path), 'V');
    assert(isequal(size(S.V), F.shape(:).') && all(isfinite(S.V(:))) && ...
        all(S.V(:) >= 0) && all(sum(S.V, 1) > 0), 'run_gm_tests:seed', 'Invalid saved seed %s', F.path);
    record = jsondecode(fileread(fullfile(root, 'results', 'gm_original', replace(F.path, '_V.mat', '.json'))));
    assert(record.nV == size(S.V, 2)); count = count + 1;
end
assert(count == 142);
fprintf('GM checks PASS: both timeout fallbacks, four portable jobs, one exact certificate, %d original vertex sets.\n', count);
end
