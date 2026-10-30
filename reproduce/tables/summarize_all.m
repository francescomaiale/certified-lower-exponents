function summarize_all()
% SUMMARIZE_ALL  Rebuild every summary in results/ from the certificates and run outputs that exist (no new
% computation, except the tangent bounds of Obs. 4.9 if missing): results/table2.md, scaling_sample.md, table3.md, table4.md,
% algorithm_l.md, work_precision.md and work_precision.pdf.  Rows not yet computed show as missing.
table2_test_pairs('part', 'summary');
scaling_sample('part', 'summary');
table3_switching('part', 'summary');
table4_algorithm_l('part', 'tangent');
table4_best_brackets('part', 'summary');
root = paperB_root();
if ~isempty(dir(fullfile(root, 'results', 'wp', '*.json'))), work_precision('part', 'plot'); end
end
