%RUN_GRADIENT_CONSISTENCY Check all active regularization gradients.
clearvars; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'tests'));
startup_HOI();

stats = test_gradient_consistency();
fprintf('maximum best relative error: %.3e\n', ...
    max(stats.best_relative_error));
