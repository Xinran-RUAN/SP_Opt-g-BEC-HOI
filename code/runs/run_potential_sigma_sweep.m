%RUN_POTENTIAL_SIGMA_SWEEP Compatibility entry point for the final study.
%
% The former L=8, N=512 experiment is retired. The paper data are now
% generated only by the L=32 high-resolution regularization-effects run.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root, 'runs', ...
    'run_potential_regularization_effects_L32.m'));
