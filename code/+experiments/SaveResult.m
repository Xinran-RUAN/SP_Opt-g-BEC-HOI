function filename = SaveResult(result, category, baseName, overwriteExisting)
%SAVERESULT Save a result using a uniform, post-processing-friendly schema.

if nargin < 2 || isempty(category)
    category = 'single';
end
if nargin < 3 || isempty(baseName)
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    baseName = ['result_' timestamp '.mat'];
end
if nargin < 4 || isempty(overwriteExisting)
    overwriteExisting = false;
end
if ~endsWith(baseName, '.mat')
    baseName = [baseName '.mat'];
end

root = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(root, 'results', category);
if ~isfolder(folder)
    mkdir(folder);
end
filename = fullfile(folder, baseName);
if isfile(filename) && ~overwriteExisting
    warning('experiments:SaveResult:Exists', ...
        'Result exists and overwrite is disabled: %s', filename);
    filename = '';
    return;
end

parameters = result.parameters;
regularization = result.regularization;
grid = result.grid;
solver = result.solver;
rho = result.rho;
energy = result.energy;
if isfield(result, 'target_energy')
    target_energy = result.target_energy;
else
    target_energy = energy;
end
if isfield(result, 'baseline_energy')
    baseline_energy = result.baseline_energy;
else
    baseline_energy = [];
end
if isfield(result, 'physical_energy')
    physical_energy = result.physical_energy;
    entropy_value = result.entropy_value;
    augmented_energy = result.augmented_energy;
else
    physical_energy = energy;
    entropy_value = [];
    augmented_energy = energy;
end
diagnostics = result.diagnostics;
history = result.history;
if isfield(result, 'fisher_regularization')
    fisher_regularization = result.fisher_regularization;
else
    fisher_regularization = [];
end
if isfield(result, 'potential_regularization')
    potential_regularization = result.potential_regularization;
else
    potential_regularization = [];
end
if isfield(result, 'trapping_potential')
    trapping_potential = result.trapping_potential;
else
    trapping_potential = [];
end
save(filename, 'parameters', 'regularization', 'grid', 'solver', ...
    'rho', 'energy', 'target_energy', 'baseline_energy', ...
    'physical_energy', 'entropy_value', ...
    'augmented_energy', 'fisher_regularization', ...
    'potential_regularization', 'trapping_potential', ...
    'diagnostics', 'history', '-v7.3');
end
