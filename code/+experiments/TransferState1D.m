function rhoFine = TransferState1D(rhoCoarse, gridFine, mass, problem, solver)
%TRANSFERSTATE1D Fourier-prolong and restore nodal positivity/fixed mass.

rhoFine = src.discretization.ps.Prolong(rhoCoarse, gridFine.N);
if nargin >= 5
    problem.grid = gridFine;
    if isfield(problem, 'entropy') && problem.entropy.enabled ...
            && problem.entropy.eta > 0
        transferTau = 1;
        if isfield(problem.entropy, 'transfer_tau') ...
                && ~isempty(problem.entropy.transfer_tau)
            transferTau = problem.entropy.transfer_tau;
        end
    else
        transferTau = 1;
    end
    rhoFine = src.solvers.ApplyProx(rhoFine, transferTau, problem, solver);
else
    rhoFine = src.constraints.ProjectSimplex(rhoFine, mass, gridFine.h);
end
end
