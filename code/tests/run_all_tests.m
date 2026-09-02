%RUN_ALL_TESTS Low-cost verification suite for the active package API.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();

fprintf('\nRunning HOI density-formulation tests...\n');
test_results.fourier_derivative = test_fourier_derivative();
test_results.fourier_derivative_2d = test_fourier_derivative_2d();
test_results.energy_gradient_hessian_2d = ...
    test_energy_gradient_hessian_2d();
test_results.inline_trapping_potential_interface = ...
    test_inline_trapping_potential_interface();
test_results.smooth_step_cinf = test_smooth_step_cinf();
test_results.harmonic_cinf_periodic_potential = ...
    test_harmonic_cinf_periodic_potential();
test_results.evaluate_fourier_series_1d = ...
    test_evaluate_fourier_series_1d();
test_results.fourier_prolongation = test_fourier_prolongation();
test_results.spectral_reference_projection = ...
    test_spectral_reference_projection();
test_results.spectral_error_decomposition = ...
    test_spectral_error_decomposition();
test_results.fourier_envelope_fit = test_fourier_envelope_fit();
test_results.projection_equivalence = test_projection_equivalence();
test_results.potential_regularization_derivative = ...
    test_potential_regularization_derivative();
test_results.potential_small_sigma = test_potential_small_sigma();
test_results.inline_regularization_interface = ...
    test_inline_regularization_interface();
test_results.fisher_regularization_family = ...
    test_fisher_regularization_family();
test_results.fisher_regularization_consistency = ...
    test_fisher_regularization_consistency();
test_results.general_fisher_gradient_consistency = ...
    test_general_fisher_gradient_consistency();
test_results.general_fisher_hessian_consistency = ...
    test_general_fisher_hessian_consistency();
test_results.general_potential_prox = test_general_potential_prox();
test_results.gradient_component_sum = test_gradient_component_sum();
test_results.transition_layer_diagnostics = ...
    test_transition_layer_diagnostics();
test_results.actual_endpoint_compatibility = ...
    test_actual_endpoint_compatibility();
test_results.small_sigma_interior_dispatch = ...
    test_small_sigma_interior_dispatch();
test_results.potential_positive_conservative_prox = ...
    test_potential_positive_conservative_prox();
test_results.potential_splitting_consistency = ...
    test_potential_splitting_consistency();
test_results.gradient_consistency = test_gradient_consistency();
test_results.gradient_consistency_all_potentials = ...
    test_gradient_consistency_all_potentials();
test_results.hessian_consistency = test_hessian_consistency();
test_results.interior_newton_pcg_schur = ...
    test_interior_newton_pcg_schur();
test_results.fd_hessian_preconditioner = ...
    test_fd_hessian_preconditioner();
test_results.interior_vs_pdas = test_interior_vs_pdas();
test_results.fista_to_polish_handoff = ...
    test_fista_to_polish_handoff();
test_results.discrete_convexity = test_discrete_convexity();
test_results.convexity_all_potentials = test_convexity_all_potentials();
test_results.spg_feasibility = test_spg_feasibility();
test_results.spg_descent = test_spg_descent();
test_results.pdas_small_problem = test_pdas_small_problem();
test_results.entropy_prox_mass = test_entropy_prox_mass();
test_results.entropy_prox_optimality = test_entropy_prox_optimality();
test_results.entropy_prox_eta_zero_limit = test_entropy_prox_eta_zero_limit();
test_results.entropy_gradient_hessian = test_entropy_gradient_hessian();
test_results.entropy_newton_small_problem = test_entropy_newton_small_problem();
fprintf('All HOI tests passed.\n');
