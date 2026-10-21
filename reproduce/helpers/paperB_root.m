function r = paperB_root()
% PAPERB_ROOT  The reproduce/ folder (it contains data/, certificates/, results/, tables/, ...).
r = fileparts(fileparts(mfilename('fullpath')));
end
