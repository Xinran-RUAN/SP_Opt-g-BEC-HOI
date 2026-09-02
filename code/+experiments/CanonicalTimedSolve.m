function [result, timing] = CanonicalTimedSolve(config, rho0, repeats)
%CANONICALTIMEDSOLVE Run production solves under one timing protocol.
%
% This wrapper deliberately reads the timers maintained inside the
% production FISTA and KKT-refinement routines.  Consequently, experiment
% setup and postprocessing performed by experiments.SolveGroundState do not
% enter the reported paper time.

if nargin < 3 || isempty(repeats)
    repeats = 1;
end
if ~isscalar(repeats) || repeats < 1 || repeats ~= round(repeats)
    error('experiments:CanonicalTimedSolve:InvalidRepeats', ...
        'repeats must be a positive integer.');
end
states = cell(repeats, 1);
stageTiming = cell(repeats, 1);
times = nan(repeats, 1);
for repeat = 1:repeats
    states{repeat} = experiments.SolveGroundState(config, rho0);
    stageTiming{repeat} = experiments.CanonicalSolverTiming(states{repeat});
    times(repeat) = stageTiming{repeat}.total_solver_time;
end
medianTime = median(times);
[~, selected] = min(abs(times - medianTime));
result = states{selected};
timing = stageTiming{selected};
timing.repeats = repeats;
timing.solver_times = times;
timing.median_solver_time = medianTime;
timing.selected_repeat = selected;
result.diagnostics.canonical_solver_time = timing.total_solver_time;
result.diagnostics.canonical_timing_protocol = timing.protocol;
end
