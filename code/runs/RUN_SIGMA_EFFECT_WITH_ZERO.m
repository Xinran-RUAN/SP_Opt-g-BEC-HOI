%RUN_SIGMA_EFFECT_WITH_ZERO Reproduce Figure 5.3 including sigma=0.
%
% The implementation and compatibility checks live in the existing
% fixed-sigma mesh-refinement driver.  Its positive-sigma case archives
% are reused without alteration; the new linear case is stored separately.
run(fullfile(fileparts(mfilename('fullpath')), ...
    'run_plot_sigma_spectral_accuracy.m'));
