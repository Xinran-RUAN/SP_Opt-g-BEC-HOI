%POST_MESH_REFINEMENT Plot the new spectral-comparison result schema.
root = fileparts(fileparts(mfilename('fullpath')));
if ~exist('result_file', 'var') || isempty(result_file)
    files = dir(fullfile(root, 'results', 'mesh_refinement', '*.mat'));
    if isempty(files)
        error('No saved mesh-refinement result was found.');
    end
    [~, newest] = max([files.datenum]);
    result_file = fullfile(files(newest).folder, files(newest).name);
end
data = load(result_file);
names = data.regularization.names;

figure('Name', 'Spectral state convergence'); hold on;
for j = 1:numel(names)
    r = data.diagnostics.records{j};
    loglog([r.N], [r.resolved_L2_error], 'o-', ...
        'DisplayName', [names{j} ' resolved']);
    loglog([r.N], [r.reference_tail_L2], '--', ...
        'DisplayName', [names{j} ' tail']);
end
xlabel('N'); ylabel('L2 error'); grid on; legend('Location', 'best');

figure('Name', 'Energy diagnostics'); hold on;
for j = 1:numel(names)
    r = data.diagnostics.records{j};
    loglog([r.N], [r.native_augmented_energy_error], 'o-', ...
        'DisplayName', [names{j} ' native']);
    loglog([r.N], [r.common_physical_energy_error], '--', ...
        'DisplayName', [names{j} ' common']);
end
xlabel('N'); ylabel('energy diagnostic'); grid on;
legend('Location', 'best');

figure('Name', 'Fourier coefficient decay'); hold on;
for j = 1:numel(names)
    referenceState = data.rho{j, end};
    [coefficients, modes] = ...
        src.discretization.ps.FourierCoefficients(referenceState);
    semilogy(abs(modes), abs(coefficients), '.', ...
        'DisplayName', names{j});
end
set(gca, 'YScale', 'log'); xlabel('|k|'); ylabel('abs(rho hat)');
grid on; legend('Location', 'best');

figure('Name', 'Optimization versus spatial error'); hold on;
for j = 1:numel(names)
    r = data.diagnostics.records{j};
    loglog([r.N], [r.total_spectral_L2_error], 'o-', ...
        'DisplayName', [names{j} ' total L2']);
    loglog([r.N], [r.pg_residual], '--', ...
        'DisplayName', [names{j} ' PG']);
end
xlabel('N'); ylabel('residual or state error'); grid on;
legend('Location', 'best');
