function write_text(path, txt)
% WRITE_TEXT  Write a char/string to a file (UTF-8), creating the folder if needed.
folder = fileparts(path);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid < 0, error('write_text:open', 'cannot open %s for writing', path); end
fprintf(fid, '%s\n', txt); fclose(fid);
end
