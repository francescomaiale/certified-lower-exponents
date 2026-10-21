function say = logger(path)
% LOGGER  say = logger(path) returns a function say(fmt, ...) that prints a line and appends it to the file
% path.  The file is truncated now, except for files in a folder named logs (run logs accumulate, with a time
% stamp per session).  logger('') only prints.
if ~isempty(path) && contains(path, [filesep 'logs' filesep])
    folder = fileparts(path);
    if ~isfolder(folder), mkdir(folder); end
    fid = fopen(path, 'a'); fprintf(fid, '---- %s\n', char(datetime('now'))); fclose(fid);
elseif ~isempty(path)
    folder = fileparts(path);
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    fid = fopen(path, 'w'); fclose(fid);
end
say = @(varargin) emit(path, varargin{:});
end

function emit(path, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
if ~isempty(path)
    fid = fopen(path, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
