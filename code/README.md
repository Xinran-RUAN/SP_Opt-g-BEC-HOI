# Convex optimization of the HOI density formulation

This directory is the self-contained, package-based MATLAB implementation
used by the current paper. The main model and experiments are one-dimensional;
the same vector optimization path is also exercised by two representative 2D
examples.

- Start with this file for the mathematical model and solver choices.
- See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the code layers,
  data structures, numerical invariants, and safe modification boundaries.
- See [`runs/README.md`](runs/README.md) for the complete, purpose-based run
  catalog.

The dated snapshot in `legacy/` is provenance only and is not placed on the
MATLAB path. Saved numerical data and figures live in `results/` and `figs/`;
they are outputs, not alternative source trees.

## Model

On the periodic collocation grid

```matlab
h = 2*L/N;
x = -L + (0:N-1)'*h;
```

the code minimizes

\[
E_h(\rho)=h\sum_j\left[
\frac{(D\rho)_j^2}{8s_\varepsilon(\rho_j)}+V_j\rho_j
+\frac{\beta}{2}\rho_j^2+\frac{\delta}{2}(D\rho)_j^2
\right]
\]

subject to nodal (collocation-point) positivity and fixed mass,
\(\rho_j\ge 0\) and \(h\sum_j\rho_j=1\). The main line is a convex
regularized density formulation with positivity- and mass-preserving
Fourier pseudospectral optimization. The linear-potential model remains the
baseline; it has no biharmonic penalty or spectral filter.

The paper-aligned core interface is:

```text
s_epsilon, ds_epsilon, d2s_epsilon
p_sigma,   dp_sigma,   d2p_sigma
```

`Energy`, `Gradient`, the matrix-free Hessian, the FD preconditioner, and
the potential prox call these handles directly. Simple formulas are defined
at the top of the formal run file. The active convexity-preserving legacy
Fisher factories remain:

- `shift_smooth`: \(s=\rho+\varepsilon\);
- `piecewise_c1`, `piecewise_c2`, `piecewise_c3`: the concave transition
  formulas implemented by `src.regularization.EvaluateDenominator`.

Legacy names are converted once to the same handle interface before the
solver starts; they are not mathematical switches inside the core solver.

## Layout

```text
+model/                    parameters, 1D/2D grids, potentials, initial density
+src/+regularization/      active denominator definitions and validation
+src/+potential/           linear and convex smoothed potential maps
+src/+discretization/+ps/  Fourier plans, D, D^T, energy, gradient, transfer
+src/+constraints/         mass/feasibility and two independent projections
+src/+entropy/             entropy value/derivatives, simplex prox, KKT
+src/+diagnostics/         active-set/support and Fourier-tail evidence
+src/+solvers/             ISTA, SPG, FISTA-CD, residuals, two polishers
+experiments/              configuration, assembly, transfer, timing, save schema
runs/                      paper, comparison, validation, and diagnostic entries
post/                      read-only plotting of saved studies
tests/                     low-cost numerical verification
results/, figs/            saved data and generated figures
docs/                      architecture and maintenance notes
legacy/                    dated source snapshot; never added to the path
```

## Running

Start MATLAB in this directory and initialize paths once:

```matlab
startup_HOI
```

The shortest useful workflow is:

```matlab
run('runs/run_single_ground_state.m')
run('tests/run_all_tests.m')
```

Each run keeps its editable parameters together at the top. For a single
change, edit only the corresponding value there (for example `N`,
`regularization_name`, `epsilon`, `solver_name`, or `projection_name`).

The current final potential-smoothing entry is
`run_potential_regularization_effects_L32.m`. The old
`run_potential_sigma_sweep.m` name is retained only as a compatibility
forwarder. Large paper, 2D, entropy, and diagnostic entries are listed with
their prerequisites in [`runs/README.md`](runs/README.md).

## Optimization architecture

- `ISTA` is the reference projected-gradient method for correctness and
  small-problem comparisons.
- `SPG` is the default first-stage convex minimizer. It uses safeguarded BB
  step lengths and a nonmonotone projected line search.
- `FISTA-CD` is an optional accelerated first-order comparison, not the
  production default. Its feasibility and monotone restarts mean that no
  complete unrestarted-FISTA rate claim is made here.
- `PDAS / semismooth Newton` is the optional high-accuracy matrix-free KKT
  polish used after a first-order method enters the local region.

The word **spectral** in SPG refers to the Barzilai-Borwein spectral step
length and is unrelated to the Fourier pseudospectral spatial
discretization.

All methods stop and compare using one fixed-step projected-gradient
residual (`residual_step = 1`). Energy or state plateaus are diagnostics,
not stationarity certificates. The default pipeline is:

```text
SPG (PG about 1e-8) -> optional PDAS polish (PG about 1e-12)
```

### Potential regularization study

The baseline potential contribution \(V\rho\) can permit an artificial
vacuum active region in the fixed-\(\varepsilon\) obstacle problem. The
potential study therefore compares it with

\[
V p_\sigma(\rho),\qquad
p_\sigma(\rho)=\sqrt{\rho^2+\sigma^2}-\sigma,
\]

For example, a formal run defines the formulas directly:

```matlab
epsilon = 1e-3;
s_epsilon   = @(rho) rho + epsilon;
ds_epsilon  = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));

sigma = epsilon^4;
p_sigma = @(rho) rho.^2 ./ (hypot(rho,sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho,sigma);
d2p_sigma = @(rho) ...
    (sigma./hypot(rho,sigma)).^2 ./ hypot(rho,sigma);
```

The fixed-grid scale study supports powers 1--4. The legacy names
`sqrt_same_scale` and `sqrt_squared_scale` remain aliases for powers 1 and
2. Since
\(p_\sigma'(0)=0\) and \(p_\sigma\) is convex, the modified potential term
preserves convexity when \(V\ge0\); the solver enforces this condition.
Energy, gradient, Hessian, and prox therefore share the exact same function
handles. The stable value formula avoids cancellation while its label may
retain the paper form
\(p_\sigma(\rho)=\sqrt{\rho^2+\sigma^2}-\sigma\). Both handles and labels are
saved in MAT results. `src.potential.Expression` and
`src.regularization.KineticExpression` remain documentation/legacy helpers,
not the mathematical source used by formal solves.

The KKT sign convention is uniform throughout:

\[
G(\rho)+\lambda\mathbf 1-\mu=0.
\]

For strictly positive states, projected-active counts are diagnostics only.
Interior Newton-PCG eligibility uses exact positivity, the vacuum slope
\(p_\sigma'(0)\), and the strong-convexity guard; a tiny positive tail is not
sent to PDAS merely because a projected point contains zeros.

The fixed-grid A/B/C experiment uses the same positive initial state,
FISTA-CD settings, and Bo Lin semismooth positive conservative projection
for all three variants. Entropy is disabled and no polish is applied. This
study is designed to test whether smoothing the vacuum slope removes the
artificial active set and improves Fourier accuracy; it does not assume
that spectral accuracy is restored. Target energy, baseline linear-potential
energy bias, far-field behavior, Fourier tails, and first-order solver cost
are reported separately.

For the sqrt variants, the formal FISTA-CD splitting places the potential
term together with positivity and mass conservation in a separable
positive-conservative prox. Smooth backtracking then majorizes only the
Fisher, beta, and delta terms. The complete target `Energy` and `Gradient`
remain unchanged and a separate Bo Lin full-gradient mapping certifies the
final state. `legacy_full_gradient` remains available only for the local
splitting A/B run:

```matlab
run('runs/run_potential_splitting_comparison.m')
```

### Vanishing entropy option (experimental diagnostic branch)

The optional term

\[
\eta H_h(\rho)=\eta h\sum_j\rho_j(\log\rho_j-1)
\]

is a **computational regularization**, not a new physical term. Setting
`entropy.enabled = false` or `entropy.eta = 0` recovers the original convex
density problem and its existing simplex-projection path. For `eta>0`, the
first-order method uses FISTA-CD (or reference ISTA) with an
entropy-simplex composite prox. Its optional interior equality-constrained
Newton polish does not use the PDAS active-set classifier. Convexity is
preserved.

This branch is designed to test whether removal of the artificial active
set improves Fourier spatial regularity. It does not assume that entropy
restores spectral accuracy. Entropy bias and same-`eta` spatial error must
be studied separately, and finite-box tails must be checked before a
Fourier-convergence interpretation. Classic SPG and PDAS remain the
`eta=0` baseline.

The default mesh study uses `N_list = [32 64 128 256 512]` and
`N_ref = 1024`. Its primary spectral state comparison is

```text
rho_N versus P_N rho_ref,
```

where `P_N` is an orthogonal restriction of the reference Fourier
coefficients to the nested coarse spectral space. It reports the
resolved-mode discrepancy, the omitted reference tail directly through
Parseval, and their orthogonal total. No positivity or mass projection is
applied to `P_N rho_ref`. Spectral prolongation of `rho_N` to the finest
grid remains only a secondary visualization/comparison diagnostic.

Native-grid nonlinear energy differences are kept separate from guarded
common-grid energy diagnostics. A common-grid entropy energy is invalid
if the unchanged oversampled trigonometric interpolant is nonpositive; it
is never clipped. Tail mass and tail maximum are saved separately so a
finite-box plateau is not mislabeled as Fourier error. The finest state is
called a reference only when its optimization residual and Fourier tail
pass the configured adequacy tolerances. The full study is intentionally
started only by:

```matlab
run('runs/run_mesh_refinement.m')
```

Before selecting a production entropy value, run the fixed-grid tradeoff:

```matlab
run('runs/run_entropy_bias_resolution_tradeoff.m')
```

It searches for an eta with both controlled physical-energy bias and at
least the requested number of nodal cells across a configured relative
density layer. Only after an overlap is found and an eta is selected should
`run_entropy_mesh_refinement.m` be launched.

The regularization study fixes the spatial grid and decreases `epsilon`
from large to small with a simple warm start:

```matlab
run('runs/run_regularization_comparison.m')
```

Mesh and regularization results are saved separately. Files in `post/`
only load saved MAT files and plot; they never invoke a solver. The observed
spectral behavior of smooth versus piecewise regularization is a numerical
study here, not a stated theorem.

## Tests

```matlab
run('tests/run_all_tests.m')
```

The suite checks Fourier differentiation/skew-adjointness, zero-padding
prolongation, exact reference-spectrum restriction, Parseval resolved/tail
decomposition, projection equivalence, centered energy-gradient consistency,
matrix-free Hessian consistency, random discrete convexity, SPG feasibility
and descent, FISTA-to-polish handoff, interior Newton-PCG/PDAS agreement,
entropy prox behavior, and 2D energy-gradient-Hessian consistency.

## Maintenance rule

Core numerical routines are deliberately not deduplicated merely because two
implementations look similar: independent projection paths and compatibility
adapters are part of the verification strategy. In particular, changes to
energy, gradient, Hessian actions, Fourier normalization, projection/prox,
solver iteration order, stopping criteria, or defaults require a mathematical
equivalence argument and the full test suite. The detailed risk classification
and compatibility-entry list are in
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
