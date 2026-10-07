function setup()
%SETUP Add MNAD-MRTM source directories to the MATLAB path.

rootDir = fileparts(mfilename('fullpath'));
srcDir  = fullfile(rootDir, 'src');

addpath(genpath(srcDir));

fprintf('MNAD-MRTM source code added to MATLAB path.\n');
end