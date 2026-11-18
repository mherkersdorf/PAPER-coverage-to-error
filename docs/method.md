# The method: criterion and receding-horizon loop

How the method of the paper works as implemented in `rhcInputDesigner/`. Symbols and code names are in
[`notation.md`](notation.md). The paper gives the theory (Secs. 2 and 3); this page adds what the
code does around it.

## 1 What it does

The excitation signal is not planned in advance. At every step the inputs of the whole prediction
horizon are optimized: each future sample is a decision variable. The optimizer simulates a candidate
through a surrogate model, collects the predicted points in the regressor space, and scores them with
the space-filling criterion `J`. Only the **first sample** of the best solution is applied. The rest of
the horizon is look-ahead, never a plan. Then the horizon moves by one sample and the loop repeats.

```
while not terminated
    feasible set        <- hard constraints on the decision variables
    start values        <- warm start from the previous step + target sequences
    best solution       <- multi-start local search on the objective (J + penalties)
    apply the first sample of the best solution
    measure the process (online); nothing is measured offline
    append to the signal history
    add the realized points to the coverage state
    poll for a surrogate update
end
```

Two points are easy to get wrong:

- **What is committed is not what was scored.** The scored points are a *prediction*. What enters the
  coverage state are the *realized* points, built from the signal history. Committing the predicted
  points instead would let the design believe its own forecast.
- **Scoring has no side effects.** Evaluating a candidate never changes the coverage state; only adding
  realized points does. This is what makes it safe to evaluate candidates in parallel.

Offline and online are the same loop. Offline, no process is applied and no output is measured; the
loop designs a whole signal in advance, as in the demo
(`rhcInputDesigner/examples/demo_StateSpace_Offline.m`) and in the paper's studies. Online, the
process is measured after every sample, and the surrogate can be updated from the data.

## 2 Regressor space, region of interest, supporting points

The criterion works in a **regressor space**: a vector $\vec z(k)$ built from the signal history. Two
representations are implemented. A *state-space* representation uses $[u(k), \vec x(k)]$, the joint
input--state space. A *NARX* representation uses past inputs and outputs. The space of the design need
not be that of the process or of the identified model.

The **region of interest** (ROI) is a box in that space. All distances are computed in normalized
coordinates; physical units appear only where signals are generated and plotted.

The ROI is discretized by **supporting points** $\mathcal{Q}$ on a cell-centered grid. The criterion is
evaluated only at these points. Their fill distance $h_{\mathcal{Q}}$ is the discretization slack of
the certificate in section 4.

## 3 The criterion

For the samples $\vec Z$ and the supporting points $\vec q_j$, with the distance
$d_j = \min_i \lVert \hat{\vec z}(i) - \vec q_j \rVert_2$ to the nearest sample:

$$\psi_j = \min\bigl\{ (d_j / h_0 - 1)_+ ,\ \psi_{\max} \bigr\} , \qquad J = \Bigl( \tfrac{1}{N_{\mathrm{sup}}} \sum_j \psi_j^{\,p} \Bigr)^{1/p} .$$

`J` is minimized. In the paper $p = 1$, so `J` is the mean of the deviations $\psi_j$. The reference
implementation is four lines:

```matlab
dRaw   = min(pdist2(Q, Znorm), [], 2);       % (nSup x 1), cached
dTilde = dRaw ./ epsReq;                     % distance in units of the demanded resolution
psi    = min(max(dTilde - 1, 0), psiMax);
J      = mean(psi.^pNorm)^(1/pNorm);
```

How to read it:

- **Only positive deviations count.** A supporting point that is resolved, $d_j \le h_0$, contributes
  nothing; resolving it more finely earns no further gain. This separates the criterion from a plain
  max-min space-filling score.
- **`psiMax` clips from above.** It sets an *outer radius* $r = (1 + \psi_{\max})\,h_0$: a supporting
  point whose nearest sample is farther than $r$ counts as much as one at exactly $r$. With a small
  `psiMax`, `J` is a *volumetric* measure, indifferent to how empty a region is. With `psiMax` large,
  `J` approaches the mean nearest-sample distance, biggest hole first. The default `psiMax = 1`
  is a starting point; a study sets its own value as part of its settings.
- **The criterion is flat outside the shell.** Between $h_0$ and $r$ the gradient is nonzero. Inside
  $h_0$ and beyond $r$ it is exactly zero. So only supporting points in the shell give the optimizer a
  direction. A start whose predicted points have no supporting point in the shell cannot move under
  the criterion (see "dead starts", section 7).
- **`pNorm`** interpolates between the mean ($p = 1$, the default) and the worst case
  ($p \to \infty$).

### Why not minimize the fill distance directly

The paper's design problem minimizes the fill distance $h$. As an objective it is a poor choice:

1. It is set by the largest empty ball alone. The few predicted samples of a horizon usually cannot
   reach that ball, so $h$ does not change with the horizon inputs at all: its gradient is exactly
   zero.
2. Closing any gap other than the largest leaves $h$ unchanged, so there is no reward for partial
   progress and no direction for a local search.
3. It is a max of a min, hence not smooth.
4. Evaluating it needs a global search for the largest empty ball inside the innermost loop.

`J` grades every gap at the supporting points, so any sample that comes within the outer radius of an
unresolved supporting point lowers it, and its zero still certifies $h$ (section 4).

### Demanded resolution

The baseline `epsReq0` $= h_0$ is the coarsest resolution demanded anywhere in the ROI. **The paper uses
this uniform case.** The code can refine the demand locally with three fields on the supporting points:
an application relevance $a_j$, a process complexity $c_j$ and a noise level $v_j$:

$$\varphi_j = \min\bigl\{ (a_j c_j)^{-\beta} v_j^{-\beta_v},\ 1 \bigr\} , \qquad h^{\mathrm{req}}_j = h_0\,\varphi_j .$$

Each field is calibrated so that **1 is nominal**; only values above nominal refine. The factor is
clipped at 1, so no region is ever demanded coarser than $h_0$: adaptation adds resolution where it is
needed and never takes it away elsewhere. The criterion applies no normalization of its own to the
fields, because dividing by a mean or maximum would make the demand at one point depend on every other
point. These fields shape *where* samples are placed; they are not used in the paper's studies.

## 4 What $J = 0$ certifies

- **Certificate.** If $J = 0$ for the samples, then $h \le h_0 + h_{\mathcal{Q}}$ (paper, Sec. 3.1). With
  the coverage bound of paper Sec. 2, $J = 0$ thus bounds the worst-case error of any model identified
  from these samples.
- **Surrogate.** If $J = 0$ holds for the points *predicted* by a surrogate, then
  $h \le h_0 + h_{\mathcal{Q}} + \delta_{\mathrm{S}}$, where $\delta_{\mathrm{S}}$ is the surrogate's largest
  distance between a predicted and a realized point. $\delta_{\mathrm{S}}$ is not set by the design;
  offline it is the surrogate's worst-case simulation error on the designed input.
- **Converse.** Under uniform demand, $h \le h_0$ implies $J = 0$. The zero of `J` is pinned between
  $h \le h_0$ (sufficient) and $h \le h_0 + h_{\mathcal{Q}}$ (implied).
- **Once zero, always zero.** Adding samples never increases a nearest-sample distance, so once `J`
  reaches zero at fixed demand and supporting points, the certificate holds for every later sample.
- **Bias only.** These statements concern coverage and hence bias. They are no total-error guarantee.
  Model bias, the approximation error of the chosen model class, is outside what an excitation signal
  can reduce.
- **Density of the supporting points.** The certificate balls cover the ROI only if the supporting
  points are dense relative to the demanded resolution. A coarser set would make the criterion partly
  blind: samples between supporting points would register no gain.

## 5 Objective, constraints and penalties

**`J` never contains a constraint term.** It measures distribution quality and nothing else, so it
stays comparable across runs. What the optimizer minimizes is a different function, the **objective**:
`J` plus soft-constraint penalties. Both are logged separately.

- **Hard constraints** on the decision variables are structural: bounds and linear inequalities, so a violating point is not in the feasible set at all.
- **Soft constraints** on simulated quantities, such as the state bounds $\mathcal{X}$, cannot be
  structural, because they are nonlinear in the decision variables through the surrogate. They enter
  as smooth penalties in the objective. A violation is measured as a fraction of the admissible range,
  and a violation of `violationScale` costs as much as a saturated supporting point.
- **"Penalty" means only this soft-constraint term.** The deviation $\psi_j$ is part of `J`; it is not
  a penalty.

The default parametrization optimizes **every future sample** of the horizon. There is no candidate set
to enumerate: the problem is continuous with `nHorizon * nInputs` variables. Reduced parametrizations
(splines, steps) are not implemented.

## 6 Horizon

The *prediction horizon* `nHorizon` is the number of simulated future samples that are scored. The
*applied horizon* is what is committed per iteration: one sample. With the default parametrization the
two coincide in dimension, so `nHorizon` alone sets the size of the optimization problem.

`nHorizon` is chosen by the user and stays fixed during a run. The surrogate's dominant time constant
is a guide: the horizon should be long enough for the transient that a step triggers to become
visible. It is not computed from it, so a surrogate update cannot silently redefine the horizon and
make scores incomparable across the update.

The horizon has a floor. A NARX representation with an input delay `nk` contains the applied sample
first in row `N + 1 + nk`, so a horizon of `nk` samples or fewer would score nothing that depends on
the applied input. The code states the floor (`nk + 1` for NARX, 1 otherwise) and refuses a shorter
horizon.

## 7 The optimizer

The optimizer is a replaceable component with one contract: it is handed the generator of candidates
and a score function, and returns the best decision vector with diagnostics.

- **`rid.opt.MultiStart`** (the default): local searches from several start vectors, one of them the
  warm start from the previous step, the others sequences that steer the surrogate to target input
  levels spanning the ROI. The number of starts is the option `opts.opt.nStarts`.
- **`rid.opt.GridSearch`**: exhaustive search over a discrete parametrization. It is the correctness
  oracle: on a small discrete problem any optimizer must find the grid-search optimum.

The score function returns the objective. It can return a gradient as a second output; the code uses
the analytic gradient when the score function provides one and differences it otherwise. A third
output holds `J`, the penalty and the number of supporting points in the shell, and is requested once,
at the optimum.

**Dead starts and plateaus.** In the clipped regime the criterion has exactly zero gradient outside
the shell, so a gradient-based search can stall on a flat region and report convergence. The
diagnostics therefore count the **dead starts**: starts whose predicted points have no supporting
point in the shell, so that the criterion gradient is zero at the start. (A dead start may still move
under the penalty gradient, toward feasibility rather than coverage.) For this reason `psiMax` is an
optimization parameter as much as a statistical one.

## 8 Termination

A run ends on exactly one of the following, and never stops silently:

- **Certification** (`mode = "certification"`): `J == 0`. The run reports the achieved bound
  $h_0 + h_{\mathcal{Q}}$. In this mode the demanded resolution is floored at $h_{\mathcal{Q}}$, since
  no claim finer than the supporting-point grid can be made.
- **Budget** (`mode = "budget"`): the prescribed experiment length is reached. The length is the
  commitment, so `J` reaching zero does not stop the run. Because a run with every supporting point
  resolved has no gradient left, the coverage state is reset at that moment, so that the remaining
  budget keeps producing data that covers.
- **A guard**: maximum steps, maximum wall time, or no improvement of `J` over a configured number
  of steps.

The result carries a status: `certified`, `budgetReached`, `stalled` or `guardTripped`. A stall is a
legitimate outcome: it usually means the remaining supporting points cannot be reached dynamically.

In budget mode the demanded resolution only shapes the design and no certificate is claimed.

## 9 Runtime

The criterion sits in the innermost loop, so it is implemented for speed:

1. **The raw nearest-sample distances `dRaw` are cached**, not the relative ones. The demanded
   resolution is a constant factor per point, so the nearest neighbour does not depend on it, and a
   resolution update costs only $O(N_{\mathrm{sup}})$.
2. **Candidate scoring is incremental.** A candidate costs $O(N_{\mathrm{sup}} \cdot \texttt{nHorizon})$:
   the cached distances are combined with the distances to the `nHorizon` predicted points. The full
   history is never recomputed.
3. **Active set.** A supporting point with $\psi_j = 0$ cannot contribute again while the demanded
   resolution and the surrogate hold, so scoring skips it. Adding realized points still updates every
   point, so a later change of the demand cannot expose a stale distance.
4. **A surrogate update triggers a rebuild** of the regressor points from the full history only where
   the points actually depend on the surrogate (simulated past values, or a representation built from
   it). It is rare, so a full recomputation is acceptable.

## 10 Terms in the code comments

The code comments use short names for rules of the implementation. The ones that matter:

| Name | Meaning |
|---|---|
| the criterion/objective split | `J` is the criterion; the optimizer minimizes the objective, `J` plus penalties (section 5) |
| the hard/soft constraint split | hard constraints are structural, soft ones are penalties (section 5) |
| the coverage reset | in budget mode, when `J` reaches zero the coverage state is cleared so the run keeps covering (section 8) |
| the certifiability floor | in certification mode the demanded resolution is not set finer than $h_{\mathcal{Q}}$ (section 8) |
| the active-set rule | resolved supporting points are skipped in scoring (section 9) |
| the starts-provider rule | the optimizer asks a function for the number of start vectors it wants, so it never needs to know what a target level is (section 7) |
| the scoreFcn arity rule | the score function has up to three outputs; the gradient path is detected by calling it (section 7) |
| the warm-start initialization | one start of every step is built from the previous step's solution (section 7) |
| the history struct | the signal history `hist` (`.u`, `.y`, `.x`, `.tS`, `.x0`) is the single source of truth for what has been applied and measured |
| the pastSource option | whether past output coordinates in the regressors are `"measured"` or `"simulated"` |
| the direct-parametrization default | every future sample is a decision variable (section 5) |
| the psiMax = 1 default | the default clipping level (section 3) |
| the refinement clipping at 1 | no region is demanded coarser than $h_0$ (section 3) |
| the calibration of the fields | 1 is nominal for $a_j$, $c_j$, $v_j$ (section 3) |
