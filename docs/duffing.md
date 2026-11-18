# The Duffing study (paper, Sec. 4)

The full setup of the first study. The paper's Sec. 4 gives the essentials; this page adds the
settings behind them. The files that hold the settings are in `experiments/` (see
[`studies.md`](studies.md)); the plotted values are in `data/`.

## 1 What the study asks

The process has measured states, so coverage can be assessed in the true joint input--state space
$\vec z(k) = [u(k), x_1(k), x_2(k)]$. The study asks whether the realized model error follows the
coverage, which the bound of paper Sec. 2 only limits. It also isolates the effects of the surrogate
and of the NARX space through one-factor comparisons against a reference arm, and compares the RHC
design with baseline methods.

## 2 Process and region of interest

A Duffing oscillator with a hardening spring,

$$\ddot x + 2\zeta\omega_0\dot x + \omega_0^2\,(x + \alpha x^3) = \omega_0^2\,u , \qquad y = x ,$$

with states $x_1 = x$, $x_2 = \dot x$, $\omega_0 = 1$ rad/s, $\zeta = 0.3$ and $\alpha = 1$. It is integrated
with fourth-order Runge--Kutta and a zero-order-hold input at the sample time $t_S = 0.2$ s. For every
constant input it has a unique, globally asymptotically stable equilibrium. The process is noise-free.

The **ROI** is the box $u \in [-1, 1]$, $x_1 \in [-0.9, 0.9]$, $x_2 \in [-1.2, 1.2]$. It lies inside the
set of states reachable from rest with $|u| \le 1$; this was verified by a grid-based forward
reachability computation, in which every cell of the box is reached. The box extends beyond the
static range of $x_1$ ($\pm 0.68$), so its outer part is visited only in transit. All distances in the
study are computed in the ROI scaled to the unit cube.

The nonlinearity matters for model accuracy. Over the ROI, the best affine model of the one-step
increment $\Delta x_2$ has an RMSE of 10.9 % of its standard deviation (1.7 % for $\Delta x_1$), measured
on a $41^3$ grid.

**Admissible region of the design.** The input is limited to the ROI's own range, $[-1, 1]$, for every
method. The states are a soft box at the limits that the process can reach, $x_1 \in [-1.42, 1.42]$ and
$x_2 \in [-2.0, 2.0]$, so the soft-constraint penalty never acts on a transit outside the ROI. The
input rate is limited only by the amplitude width, so the input axis is trivially coverable. What the
study measures is therefore coverage of the **states**.

## 3 Surrogates

- **S0** is the exact process, used as a surrogate. It is the upper end of the comparison.
- **S1** is a linear state-space model, fitted by equation-error least squares with an intercept on
  the *measured states* of a pilot signal. Fitting on measured states puts S1's states in the process's
  coordinates, which the joint-space design needs. The same S1 serves arms (ii) and (iii). Its
  matrices are in `experiments/p2/surrogates/S1.json`.
  - The pilot is a 10 000-sample APRBS from rest with amplitude at the bounds of $u$, mean hold 16
    samples (see section 5).
  - S1 has a dominant time constant of 3.17 s. It is stiffer than the process at the origin, because
    the fit averages the hardening spring.
  - Its one-step error is at most 1 % RMS of the ROI width. Its free-run error on a validation APRBS
    is up to 29 % of the width of $x_1$, largest in the big excursions near the edge of $x_1$. That
    free-run error is the premise on which arms (ii) and (iii) run, because they re-simulate the whole
    history with S1.

States are never measured during the design. Measured states enter in two places only: the S1 fit
and the evaluation. The design is offline.

## 4 The three design arms

All arms are offline, in budget mode, with past states simulated by the surrogate, starting at rest.
They share $N = 8000$, the horizon $L = 17$, the outer radius rule below and the number of seeds
(1--100, each seed one signal).

| Arm | Space | $n_z$ | Surrogate | $h_0$ | $h_{\mathcal{Q}}$ | $N_{\mathrm{sup}}$ | Isolates |
|---|---|---|---|---|---|---|---|
| (i) | joint $[u(k), x_1(k), x_2(k)]$ | 3 | S0 | 0.043 | 0.048 | 5 832 (18 per axis) | the surrogate |
| (ii) | joint | 3 | S1 | 0.043 | 0.048 | 5 832 | the reference |
| (iii) | NARX $[u(k{-}1), u(k{-}2), y(k{-}1), y(k{-}2)]$ | 4 | S1 | 0.050 | 0.050 | 160 000 (20 per axis) | the representation |

- **Arms (i) against (ii)** differ only in the surrogate. **Arms (ii) against (iii)** differ only in
  the space, with S1 in both.
- **Demanded resolution.** $h_0 = 0.043\sqrt{n_z/3}$: the fill distance of a grid with 20 points per
  axis, derived from a grid with $N = 8000$ points in 3-D, carried over to every space with the same
  share of the ROI per axis. It therefore grows with $\sqrt{n_z}$.
- **Outer radius.** $r = 0.2\sqrt{n_z/3}$, so $\psi_{\max} = r/h_0 - 1 \approx 3.65$ in every space.
- **Supporting points.** A cell-centered grid, with the number of points per axis chosen as the smallest
  that reaches the absolute fill distance 0.05. This keeps the slack $h_{\mathcal{Q}}$ approximately
  equal across dimensions. Its fill distance is $\sqrt{n_z}/(2n)$ for $n$ points per axis.
- **NARX representation of arm (iii).** Two input lags and two output lags, input delay 1, no input
  difference.
- **Optimizer.** Multi-start with a gradient-based local solver; `rid.opt.MultiStart` is the default.
- **Seeds.** The seed varies the multi-start draws. The supporting points are a deterministic grid.
- **Arm (iii) is run in a different regime.** In 4-D the demand (0.050) is finer than what 8 000 samples
  can reach (about 0.106), so hardly any supporting point is resolved. That gap is the cost of the
  extra dimension, and it shows in the measured fill distance.

## 5 Baseline methods

All methods use $u \in [-1, 1]$ and $N$ up to 8 000.

- **APRBS.** Exponential holds, mean hold 16 samples (3.2 s, about one time constant). The levels are
  drawn without repetition from equally spaced values across $[-1, 1]$. A switch occurs at each sample
  with probability 1/16, so holds are geometric with minimum 1 and no maximum.
- **Multisine.** Schroeder phases, flat band from the lowest line to just below the Nyquist frequency,
  peak 1. One complete period per dataset size.
- **Chirp.** Linear, from MATLAB's `chirp`, one complete sweep per dataset size over the same band.

Every method is evaluated at **80 dataset sizes**, log-spaced from 50 to 8 000, with 100 seeds for the
methods with randomness. For the RHC arms and the APRBS, a shorter dataset is the first samples of a
run. The multisine and the chirp are generated anew for each size.

## 6 Evaluation

Two model architectures with different inductive biases translate coverage into realized error, so
that the findings do not rest on one. Both are trained on noise-free measured states and predict the
one-step increment $\Delta\vec x(k)$ from $[u(k), \vec x(k)]$, whatever the space in which the data were
designed. Their hyperparameters are identical across methods and sizes, so that performance
differences stem from the data alone.

- **Gaussian process (GP).** Matérn 5/2 kernel, one fixed length scale of 0.25 in the normalized ROI
  (0.5 in the $[-1, 1]$ scaling of the ROI), affine mean, a small relative jitter. It is deliberately
  local, as fortuitous extrapolation cannot be relied upon.
- **Local model network (LMN).** Trained with LOLIMOT (LMNtool 1.5.3), affine local models, one network
  with both outputs, inputs $[u, x_1, x_2]$. The number of local models grows with the dataset size,
  $M = \max\bigl(1, \mathrm{round}(\min(2.25\sqrt N,\, N/10))\bigr)$, identically for every method.
- An **affine model** is evaluated as a reference at every size.

**Measures.**

- **Fill distance $h$**, measured with $10^5$ probe points in the true input--state space.
- **Worst-case one-step error $E_{\max}$.** The maximum over both components of $\Delta\vec x$, each
  normalized by the ROI half-width of its state (0.9 and 1.2). It is evaluated on a $41^3$ grid over the
  ROI plus the center of the largest empty ball. As maxima over finite point sets, the measured $h$
  and $E_{\max}$ never exceed the suprema they approximate.
- **NRMSE of the simulated output** (RMSE divided by the standard deviation of $y$). The test signal is
  an APRBS, a multisine and a chirp of 20 000 samples each, applied in one free run from rest. The
  score is taken inside the ROI, over a grid of $5^3$ cubes; a cube with fewer than 30 scored test
  points is not scored. The single run favors the baseline methods as a group, but not any individual
  method.
- **Gap map.** For every grid point, the distance to its nearest training sample, in units of $h$.
  $d^{\star}$ is this distance at the location of the worst one-step error; the hole analysis compares
  it with the largest gap $h$ per run and size.

## 7 A check that is not in the paper: the GP length scale

At $N = 8000$ the evaluation was repeated with the GP length scale halved and doubled (0.125 and 0.5
in the normalized ROI, around the paper's 0.25). The median $E_{\max}$ of the three RHC arms
(i)/(ii)/(iii) against the three baselines APRBS/multisine/chirp:

| Length scale | RHC arms (i)/(ii)/(iii) | Baselines APRBS/multisine/chirp | Best baseline / arm (ii) |
|---|---|---|---|
| 0.125 | 0.0442 / 0.0486 / 0.0480 | 0.0812 / 0.0571 / 0.0519 | 1.07 |
| 0.25 (paper) | 0.0220 / 0.0243 / 0.0226 | 0.0687 / 0.0606 / 0.0563 | 2.3 |
| 0.5 | 0.0088 / 0.0094 / 0.0085 | 0.0446 / 0.0568 / 0.0559 | 4.7 |

The ranking of the methods holds at all three length scales, but the lead of the RHC arms shrinks
with a shorter length scale. The settings of the study were fixed in advance and not tuned.
