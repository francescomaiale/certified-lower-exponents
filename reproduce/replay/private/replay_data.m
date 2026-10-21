function A = replay_data(name)
% REPLAY_DATA  (replay side) The matrices of the pair "name" as decimal strings, read from data/*.json.
%   Independent of the certify path (src/core/load_pair): every data file of the package is scanned for the key,
%   and the numbers of its value are taken verbatim from the raw text.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
files = dir(fullfile(root, 'data', '*.json'));
for f = files.'
    txt = fileread(fullfile(f.folder, f.name));
    [s, e] = regexp(txt, ['"' regexptranslate('escape', char(name)) '"\s*:\s*\['], 'start', 'end', 'once');
    if isempty(s), continue, end
    lvl = 0; stop = e; open3 = []; 
    for p = e:numel(txt)
        if txt(p) == '[', lvl = lvl + 1; if lvl == 2, open3(end + 1) = p; end, end %#ok<AGROW>
        if txt(p) == ']', lvl = lvl - 1; end
        if lvl == 0, stop = p; break, end
    end
    m = numel(open3); A = cell(1, m);
    for i = 1:m
        if i < m, seg = txt(open3(i):open3(i + 1) - 1); else, seg = txt(open3(i):stop - 1); end
        rows = regexp(seg, '\[([^\[\]]*)\]', 'tokens');
        d = numel(rows); M = strings(d, d);
        for r = 1:d
            v = strtrim(split(string(rows{r}{1}), ','));
            if numel(v) ~= d, error('replay:data', '%s: matrix %d is not square', name, i); end
            M(r, :) = v.';
        end
        A{i} = M;
    end
    return
end
error('replay:data', 'pair "%s" not found in data/*.json', name);
end
