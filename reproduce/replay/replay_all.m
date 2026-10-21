function [ok, T] = replay_all(folders)
% REPLAY_ALL  Re-check every certificate under certificates/ with the independent MATLAB replays.
%   [ok, T] = replay_all              all subfolders of certificates/
%   [ok, T] = replay_all({'table2'})  some of them
%   Each file is classified by its fields: polytope ('weights') -> replay_lower, implicit Euler chain ('versions') ->
%   replay_upper, chattering/schedule ('readings') -> replay_best_upper, invariant cone ('h1') -> replay_cone,
%   invariant faces ('orientations') -> replay_faces.  The data are re-read from data/ by the replays themselves.
%   T lists file, kind, ok and seconds; ok is true iff every certificate is valid.
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(folders)
    d = dir(fullfile(root, 'certificates')); folders = {d([d.isdir] & ~startsWith({d.name}, '.')).name};
end
T = struct('file', {}, 'kind', {}, 'ok', {}, 'secs', {});
selection_ok = true;
for f = folders(:).'
    files = dir(fullfile(root, 'certificates', f{1}, '*.json'));
    if isempty(files), selection_ok = false; fprintf('INVALID selection: no certificates in %s\n', f{1}); end
    for k = 1:numel(files)
        p = fullfile(files(k).folder, files(k).name); C = jsondecode(fileread(p)); t0 = tic;
        if isfield(C, 'weights'), kind = 'lower'; good = replay_lower(C);
        elseif isfield(C, 'versions'), kind = 'upper'; good = replay_upper(C);
        elseif isfield(C, 'readings'), kind = 'best_upper'; good = replay_best_upper(C);
        elseif isfield(C, 'h1'), kind = 'cone'; good = replay_cone(C);
        elseif isfield(C, 'orientations'), kind = 'faces'; good = replay_faces(C);
        else, kind = 'unknown'; good = false;
        end
        T(end + 1) = struct('file', fullfile(f{1}, files(k).name), 'kind', kind, 'ok', good, 'secs', toc(t0)); %#ok<AGROW>
        if good, s = 'VALID'; else, s = 'INVALID'; end
        fprintf('%-45s %-10s %s (%.1fs)\n', T(end).file, kind, s, T(end).secs);
    end
end
ok = selection_ok && ~isempty(T) && all([T.ok]);
kinds = unique({T.kind});
parts = cellfun(@(k) sprintf('%s %d', k, nnz(strcmp({T.kind}, k))), kinds, 'UniformOutput', false);
if ~selection_ok || isempty(T), status = 'INVALID selection';
else, status = ternary(ok, 'all VALID', sprintf('%d INVALID', nnz(~[T.ok]))); end
fprintf('replay_all: %d certificates (%s): %s\n', numel(T), strjoin(parts, ', '), status);
if nargout == 0 && ~ok, error('replay_all:invalid', 'certificate replay failed or selection was empty'); end
end


function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
