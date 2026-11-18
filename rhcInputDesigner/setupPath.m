function setupPath()
%SETUPPATH  Put the project's source, test and example folders on the MATLAB path.
%
%   Adds ONLY src/, tests/ and examples/. Never use genpath from the repository
%   root: rhcInputDesignerOld/ contains legacy functions of similar name that
%   would shadow the new implementation in src/ (see rhcInputDesignerOld/README.md).
%
%   Usage:
%       setupPath
%
%   Inputs:  none
%   Outputs: none

here = fileparts(mfilename("fullpath"));

addpath(fullfile(here, "src"));                    % namespace +rid resolves from here
addpath(genpath(fullfile(here, "tests")));
addpath(fullfile(here, "examples"));

end
