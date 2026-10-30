function run_all(profile)
% RUN_ALL  Reproduce the computations of the paper, in order.
%   run_all('quick')   (default) everything that takes minutes: Table 1, Computation 3.9, Table 2 without its
%                      longest row (case 1 at h = 0.025), the Table 4 certificates, Ex. 4.6, and the original
%                      d = 5 scaling and switching pairs of Obs. 4.8 at h = 0.1
%   run_all('full')    the original suite with the paper's settings (hours; see README.md).
%                      Run table3_transposed, table3_extended and larger_samples separately.
%   The optional Gripenberg runs of Appendix B and Tables S1-S3 use gm/gm_fixed_job; see gm/README.md.
%   They are bundled with this repository but run separately from this driver.
%   Most certificate drivers skip files that already exist; the rate searches and the BCV driver run again. The summaries
%   are written to results/*.md.

if nargin < 1, profile = 'quick'; end
if ~ismember(string(profile), ["quick", "full"]), error('run_all:profile', 'profile must be quick or full'); end
full = strcmp(profile, 'full');
t0 = tic;
fprintf('== Table 1\n'); table1_sandwich();
fprintf('== Computation 3.9\n'); comp38_reducible_pairs();
fprintf('== Table 2\n');
rows = {'planar_0.1', 'planar_0.05', 'planar_0.025', 'planar_0.0125', 'case1_0.1', 'case1_0.05', ...
        'case1440_0.1', 'case1440_0.05', 'case1440_0.025', 'case148_0.1', 'case148_0.05', 'case148_0.025'};
if full, rows{end + 1} = 'case1_0.025'; end
table2_test_pairs('rows', rows, 'part', 'lower'); table2_test_pairs('rows', rows, 'part', 'transfer');
table2_test_pairs('part', 'upper'); table2_test_pairs('part', 'summary');
fprintf('== Ex. 4.6\n'); ex410_bcv();
fprintf('== Obs. 4.8 / Table S4\n');
if full
    table7_scaling();
else
    table7_scaling('pairs', {'d5_dense1', 'd5_dense2', 'd5_dense3', 'd5_sparse1', 'd5_sparse2', 'd5_sparse3'}, ...
                   'hs', {'0.1'}, 'part', 'lower');
    table7_scaling('part', 'upper'); table7_scaling('part', 'summary');
end
fprintf('== Table 4\n'); table4_best_brackets('part', 'certify');
if full, table4_algorithm_l(); end
table4_best_brackets('part', 'summary');
fprintf('== Obs. 4.8 / Table 3\n');
if full
    table3_switching();
else
    table3_switching('pairs', arrayfun(@(k) sprintf('sw5_%d', k), 1:10, 'UniformOutput', false), 'hs', {'0.1'}, ...
                     'part', 'all');
end
fprintf('== Appendix C / Table 5\n'); table8_rates();
if full
    fprintf('== Fig. 2\n'); fig1_work_precision();
end
fprintf('== done in %.0f s; summaries in results/*.md\n', toc(t0));
end
