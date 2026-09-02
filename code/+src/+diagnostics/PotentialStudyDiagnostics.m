function record = PotentialStudyDiagnostics(result, baselineResult)
%POTENTIALSTUDYDIAGNOSTICS Target/baseline, support, tail, and solver cost.

zeroOptions.zero_tol_list = [0, 1e-14, 1e-12, 1e-10, 1e-8];
active = src.diagnostics.ActiveSetDiagnostics( ...
    result.rho, result.grid, zeroOptions);
relativeSupport = src.diagnostics.RelativeSupportDiagnostics( ...
    result.rho, result.grid);
farField = src.diagnostics.FarFieldDiagnostics( ...
    result.rho, result.grid, result.problem.plan);
fourier = src.diagnostics.FourierTailDiagnostics(result.rho);

record.variant = result.potential_regularization.name;
record.power = result.potential_regularization.power;
record.sigma = result.potential_regularization.sigma;
record.target_energy = result.target_energy;
record.baseline_energy = result.baseline_energy;
record.baseline_energy_bias = ...
    result.baseline_energy - baselineResult.baseline_energy;
record.density_L2_bias = sqrt(result.grid.h ...
    * sum((result.rho - baselineResult.rho) .^ 2));
record.density_Linf_bias = max(abs(result.rho - baselineResult.rho));
record.pg_residual = result.diagnostics.pg_residual;
record.final_pg_residual = result.diagnostics.final_pg_residual;
record.projected_active_count = ...
    result.diagnostics.projected_active_count_handoff;
record.projected_active_count_handoff = ...
    result.diagnostics.projected_active_count_handoff;
record.projected_active_count_final = ...
    result.diagnostics.projected_active_count_final;
record.main_stop_reason = result.diagnostics.main_stop_reason;
record.polish_entered = result.diagnostics.polish_entered;
record.polish_status = result.diagnostics.polish_status;
record.polish_solver = result.diagnostics.polish_solver;
record.polish_failure_message = result.diagnostics.polish_failure_message;
record.polish_fallback_reason = result.diagnostics.polish_fallback_reason;
record.iterations = result.diagnostics.main_iterations;
record.main_iterations = result.diagnostics.main_iterations;
record.polish_iterations = result.diagnostics.polish_iterations;
record.elapsed_time = result.diagnostics.total_elapsed_time;
record.main_elapsed_time = result.diagnostics.main_elapsed_time;
record.polish_elapsed_time = result.diagnostics.polish_elapsed_time;
record.accepted_L_final = result.diagnostics.accepted_L_final;
record.accepted_L_max = result.diagnostics.accepted_L_max;
record.mean_backtracks = result.diagnostics.mean_backtracks;
record.max_backtracks = result.diagnostics.max_backtracks;
record.restart_count = result.diagnostics.restart_count;
record.mean_prox_lambda_iterations = ...
    result.diagnostics.mean_prox_lambda_iterations;
record.max_prox_lambda_iterations = ...
    result.diagnostics.max_prox_lambda_iterations;
record.mean_prox_inner_iterations = ...
    result.diagnostics.mean_prox_inner_iterations;
record.max_prox_inner_iterations = ...
    result.diagnostics.max_prox_inner_iterations;
record.total_pcg_z_iterations = ...
    result.diagnostics.total_pcg_z_iterations;
record.total_pcg_w_iterations = ...
    result.diagnostics.total_pcg_w_iterations;
record.max_pcg_iterations = result.diagnostics.max_pcg_iterations;

record.exact_zero_count = active.exact_zero_count;
record.strict_positive_count = nnz(result.rho > 0);
record.exact_zero_fraction = active.zero_fraction_by_tol(1);
record.zero_fraction_1e14 = active.zero_fraction_by_tol(2);
record.zero_fraction_1e12 = active.zero_fraction_by_tol(3);
record.zero_fraction_1e10 = active.zero_fraction_by_tol(4);
record.zero_fraction_1e8 = active.zero_fraction_by_tol(5);
record.min_positive_density = active.min_positive_density;
record.radius_rel_1e3 = relativeSupport.radius_rel_1e3;
record.radius_rel_1e6 = relativeSupport.radius_rel_1e6;
record.radius_rel_1e10 = relativeSupport.radius_rel_1e10;
record.radius_rel_1e14 = relativeSupport.radius_rel_1e14;

record.tail_mass = farField.tail_mass;
record.tail_max = farField.tail_max;
record.boundary_strip_mass = farField.boundary_strip_mass;
record.boundary_strip_max = farField.boundary_strip_max;
record.min_density = farField.min_density;
record.edge_density_left = farField.edge_density_left;
record.edge_density_right = farField.edge_density_right;
record.edge_density = farField.edge_density;
record.edge_abs_drho_left = farField.edge_abs_drho_left;
record.edge_abs_drho_right = farField.edge_abs_drho_right;
record.edge_abs_drho = farField.edge_abs_drho;
record.boundary_strip_max_abs_drho = ...
    farField.boundary_strip_max_abs_drho;
record.rho_hat_abs = fourier.rho_hat_abs;
record.fourier_modes = fourier.modes;
record.fft_tail_ratio_quarter = fourier.tail_ratio_quarter;
record.fft_tail_ratio_third = fourier.tail_ratio_third;
record.fourier_decay_slope = fourier.fourier_decay_slope;
if record.sigma > 0
    record.tail_mass_over_sigma = record.tail_mass / record.sigma;
    record.tail_max_over_sigma = record.tail_max / record.sigma;
    record.edge_density_over_sigma = record.edge_density / record.sigma;
else
    record.tail_mass_over_sigma = NaN;
    record.tail_max_over_sigma = NaN;
    record.edge_density_over_sigma = NaN;
end
end
