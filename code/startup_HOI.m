function root = startup_HOI()
%STARTUP_HOI Add the active 1D HOI project entry points to the MATLAB path.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'runs'));
addpath(fullfile(root, 'post'));
addpath(fullfile(root, 'tests'));

resultFolders = {
    fullfile(root, 'results', 'single')
    fullfile(root, 'results', 'mesh_refinement')
    fullfile(root, 'results', 'regularization')
    fullfile(root, 'results', 'diagnostics')
};
for j = 1:numel(resultFolders)
    if ~isfolder(resultFolders{j})
        mkdir(resultFolders{j});
    end
end

fprintf('HOI density-formulation project ready: %s\n', root);
end
