%POST_REGULARIZATION_COMPARISON Plot a saved fixed-grid epsilon study.
root = fileparts(fileparts(mfilename('fullpath')));
if ~exist('result_file', 'var') || isempty(result_file)
    files = dir(fullfile(root, 'results', 'regularization', '*.mat'));
    if isempty(files)
        error('No saved regularization-comparison result was found.');
    end
    [~, newest] = max([files.datenum]);
    result_file = fullfile(files(newest).folder, files(newest).name);
end
data = load(result_file);
names = data.regularization.names;
epsilon = data.regularization.epsilon(:);
referenceIndex = numel(epsilon);

figure('Name', 'HOI regularization diagnostics');
tiles = tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for j = 1:numel(names)
    errors = abs(data.energy(j, 1:referenceIndex-1) ...
        - data.energy(j, referenceIndex));
    loglog(epsilon(1:end-1), errors, '-o', 'DisplayName', names{j});
end
xlabel('\epsilon'); ylabel('|E_\epsilon-E_{\epsilon_{min}}|');
grid on; legend('Location', 'best');

nexttile; hold on;
for j = 1:numel(names)
    errors = zeros(referenceIndex-1, 1);
    referenceDensity = data.rho{j, referenceIndex};
    for ell = 1:referenceIndex-1
        difference = data.rho{j, ell} - referenceDensity;
        errors(ell) = sqrt(data.grid.h * sum(difference .^ 2));
    end
    loglog(epsilon(1:end-1), errors, '-o', 'DisplayName', names{j});
end
xlabel('\epsilon'); ylabel('density L^2 difference');
grid on; legend('Location', 'best');

nexttile; hold on;
for j = 1:numel(names)
    loglog(epsilon, data.diagnostics.pg_residual(j, :), '-o', ...
        'DisplayName', names{j});
end
xlabel('\epsilon'); ylabel('projected-gradient residual'); grid on;
legend('Location', 'best');

nexttile; hold on;
for j = 1:numel(names)
    semilogx(epsilon, data.diagnostics.iterations(j, :), '-o', ...
        'DisplayName', names{j});
end
xlabel('\epsilon'); ylabel('iterations'); grid on; legend('Location', 'best');
title(tiles, sprintf('Fixed-grid regularization comparison, N=%d', data.grid.N));
