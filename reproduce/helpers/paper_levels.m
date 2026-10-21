function v = paper_levels(group, key)
% PAPER_LEVELS  Recorded outputs of the paper's floating-point searches (data/paper_search_results.json).
%   v = paper_levels(group, key), e.g. paper_levels('table2', 'planar_0.1').t_build: the level at which the paper's
%   bisection last closed.  Groups: table2, scaling, switching, best_upper, glp, bcv, wp.  These numbers are only
%   starting points: every certificate built from them is checked in exact arithmetic.
persistent J
if isempty(J), J = jsondecode(fileread(fullfile(paperB_root(), 'data', 'paper_search_results.json'))); end
v = J.(group);
if nargin > 1, v = v.(matlab.lang.makeValidName(key)); end
end
