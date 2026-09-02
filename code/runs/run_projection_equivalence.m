%RUN_PROJECTION_EQUIVALENCE Compare simplex and semismooth projections.
clearvars; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root, 'tests'));
startup_HOI();

stats = test_projection_equivalence();
fprintf('maximum projection difference: %.3e\n', stats.max_difference);
