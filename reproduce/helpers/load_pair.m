function [A, names, file] = load_pair(file, name)
% LOAD_PAIR  The matrices of one pair (or family) of a data file, as decimal strings exactly as written.
%   A = load_pair(file, name) returns a cell array {A_1, ..., A_m} of d x d string matrices for the entry "name"
%   of a JSON file of the form {"name": [A_1, A_2, ...], ...}.  file is a path or one of the shortcuts
%     'test'       data/test_pairs.json       planar, sharp, case1, case77, case148, case677, case1195, case1440
%     'scaling'    data/scaling_pairs.json    d5_dense1, ..., d10_sparse3 (Obs. 4.8)
%     'switching'  data/switching_pairs.json  sw5_1, ..., sw10_10 (Obs. 4.8)
%     'bcv'        data/bcv13.json            bcv13 (Ex. 4.6)
%     'switching_journal'  data/switching_pairs_journal.json  sw5_11, ..., sw10_30, sw20_1, ..., sw20_5 (larger sample, Obs. 4.8)
%     'scaling_journal'    data/scaling_pairs_journal.json    d5_dense4, ..., d20_sparse2 (larger sample, Obs. 4.8)
%   jsondecode would turn the entries into doubles and lose their decimal reading, so the numbers are taken from
%   the raw text.  [~, names] = load_pair(file) lists the names; the third output is the resolved path.
%   With name omitted and a single output, A is a struct mapping every name to its matrices.

short = struct('test', 'test_pairs.json', 'scaling', 'scaling_pairs.json', 'switching', 'switching_pairs.json', ...
               'bcv', 'bcv13.json', 'switching_journal', 'switching_pairs_journal.json', ...
               'scaling_journal', 'scaling_pairs_journal.json');
file = char(file);
if isfield(short, file), file = fullfile(paperB_root(), 'data', short.(file)); end
txt = fileread(file);
names = regexp(txt, '"([^"]+)"\s*:\s*\[', 'tokens'); names = cellfun(@(c) c{1}, names, 'UniformOutput', false);
if nargin < 2
    if nargout <= 1
        A = struct();
        for k = 1:numel(names), A.(matlab.lang.makeValidName(names{k})) = load_pair(file, names{k}); end
    else
        A = {};
    end
    return
end
key = regexp(txt, ['"' regexptranslate('escape', char(name)) '"\s*:\s*\['], 'end', 'once');
if isempty(key), error('load_pair:name', 'pair "%s" not found in %s', name, file); end
depth = 0; e = key;                       % the balanced [...] that starts at txt(key)
for p = key:numel(txt)
    if txt(p) == '[', depth = depth + 1; elseif txt(p) == ']', depth = depth - 1; end
    if depth == 0, e = p; break, end
end
seg = txt(key:e);
tok = regexp(seg, '-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', 'match');
J = jsondecode(seg);                      % only for the shape: m x d x d array or cell of matrices
if iscell(J), m = numel(J); d = size(J{1}, 1); else, m = size(J, 1); d = size(J, 2); end
if numel(tok) ~= m * d * d, error('load_pair:shape', 'pair "%s" is not a list of square matrices', name); end
A = cell(1, m);
for i = 1:m
    A{i} = reshape(string(tok((i - 1) * d * d + (1:d * d))), d, d).';   % rows as in the file
end
end
