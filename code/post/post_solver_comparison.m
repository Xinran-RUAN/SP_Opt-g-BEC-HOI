%POST_SOLVER_COMPARISON Plot saved first-order histories and common polish.
root = fileparts(fileparts(mfilename('fullpath')));
if ~exist('result_file', 'var') || isempty(result_file)
    files = dir(fullfile(root, 'results', 'diagnostics', ...
        'solver_comparison_*.mat'));
    if isempty(files)
        error('No saved solver-comparison result was found.');
    end
    [~, newest] = max([files.datenum]);
    result_file = fullfile(files(newest).folder, files(newest).name);
end
data = load(result_file);
referenceEnergy = min(data.energy);

figure('Name', 'HOI solver comparison');
tiles = tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on;
for j = 1:numel(data.history)
    main = data.history{j}.main;
    polish = data.history{j}.polish;
    time = main.elapsed_time;
    energyTrace = main.energy;
    if isfield(polish, 'elapsed_time') && ~isempty(polish.elapsed_time)
        time = [time; data.diagnostics.records(j).first_stage_time ...
            + polish.elapsed_time]; %#ok<AGROW>
        energyTrace = [energyTrace; polish.energy]; %#ok<AGROW>
    end
    gap = abs(energyTrace - referenceEnergy);
    gap(gap == 0) = NaN;
    semilogy(time, gap, 'LineWidth', 1.2, ...
        'DisplayName', upper(strrep(data.diagnostics.solver_names{j}, '_', '-')));
end
xlabel('elapsed time'); ylabel('|E-E_{best final}|'); grid on;
legend('Location', 'best');

nexttile; hold on;
for j = 1:numel(data.history)
    main = data.history{j}.main;
    polish = data.history{j}.polish;
    time = main.elapsed_time;
    residual = main.pg_residual;
    if isfield(polish, 'elapsed_time') && ~isempty(polish.elapsed_time)
        time = [time; data.diagnostics.records(j).first_stage_time ...
            + polish.elapsed_time]; %#ok<AGROW>
        residual = [residual; polish.pg_residual]; %#ok<AGROW>
    end
    semilogy(time, residual, 'LineWidth', 1.2, ...
        'DisplayName', upper(strrep(data.diagnostics.solver_names{j}, '_', '-')));
end
xlabel('elapsed time'); ylabel('projected-gradient residual'); grid on;
legend('Location', 'best');
title(tiles, sprintf('Shared problem and PDAS polish, N=%d', data.grid.N));
