function startup_reproduce()
% STARTUP_REPRODUCE  Add the code (../src) and the reproduction folders to the MATLAB path (run once per session).
r = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(r), 'src')), fullfile(r, 'helpers'), fullfile(r, 'replay'), ...
        fullfile(r, 'tables'), fullfile(r, 'tests'));
end
