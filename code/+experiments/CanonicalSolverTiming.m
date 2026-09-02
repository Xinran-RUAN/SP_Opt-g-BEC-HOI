function timing = CanonicalSolverTiming(result)
%CANONICALSOLVERTIMING Extract the paper's canonical solver wall time.
%
% The canonical interval is the time spent in the production FISTA stage
% plus the production KKT-refinement stage.  It includes the standard
% backtracking, proximal solves, switch checks, Newton globalization,
% preconditioner construction, PCG solves, and convergence diagnostics
% performed inside those stages.  Problem construction, initial-density
% generation, file I/O, plotting, and experiment-only postprocessing are
% outside this interval.

if ~isstruct(result) || ~isfield(result, 'diagnostics')
    error('experiments:CanonicalSolverTiming:InvalidResult', ...
        'A production result with diagnostics is required.');
end
d = result.diagnostics;
required = {'main_elapsed_time', 'polish_elapsed_time'};
if ~all(isfield(d, required))
    error('experiments:CanonicalSolverTiming:MissingStageTimes', ...
        'The result does not contain production FISTA/Newton stage times.');
end
timing.protocol = 'canonical_solver_time_v1';
timing.definition = [ ...
    'production FISTA stage plus production KKT-refinement stage; ' ...
    'problem setup, initialization, I/O, plotting, and experiment-only ' ...
    'postprocessing excluded'];
timing.fista_time = d.main_elapsed_time;
timing.newton_time = d.polish_elapsed_time;
timing.total_solver_time = timing.fista_time + timing.newton_time;
timing.repeats = 1;
timing.solver_times = timing.total_solver_time;
timing.median_solver_time = timing.total_solver_time;
timing.selected_repeat = 1;
if isfield(d, 'total_elapsed_time')
    scale = max(1, abs(timing.total_solver_time));
    tolerance = 100 * eps(scale);
    if abs(d.total_elapsed_time - timing.total_solver_time) > tolerance
        error('experiments:CanonicalSolverTiming:InconsistentTotal', ...
            ['diagnostics.total_elapsed_time is inconsistent with the ' ...
            'sum of the production stage timers.']);
    end
end
end
