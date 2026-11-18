# Notation

The symbols of the paper, with the names the code uses for them. A quantity has exactly one code
name.

## Conventions

- **Samples in rows.** Every point set is `(nPoints x nDim)`.
- **Indices.** $i$ is a sample, $j$ a supporting point, $k$ a time step.
- **Normalized coordinates.** All criterion and distance computations happen in the normalized
  regressor space. Physical units appear only at the boundaries (signal generation, plotting).
- **Transpose** is written $^{\top}$. Vectors are columns; a matrix whose rows are samples is
  $\vec Z^{\top} = [\vec z(1), \ldots, \vec z(N)]$.
- **Hats mark model outputs.** A quantity that comes out of a model carries a hat; a quantity that
  comes out of the process carries none. An upright superscript names the model: $\mathrm{S}$ is
  the surrogate used inside the design. Example: the process state is $\vec x(k)$, the surrogate's
  simulated state is $\hat{\vec x}^{\mathrm{S}}(k)$.
- **Ideal and measured values.** $h$, $E_{\max}$ and the moduli of continuity are defined as
  suprema over the region of interest. When they are evaluated numerically, the symbol is kept and
  the value is called *measured*: $h$ from $10^5$ probe points, $E_{\max}$ on a $41^3$ grid plus the
  center of the largest empty ball. A measured value is a lower estimate of the ideal one.

## Spaces and sets

| Symbol | Meaning | Code |
|---|---|---|
| $\mathcal{R}$ | region of interest (ROI): a box in the joint input--state space $\mathcal{U} \times \mathcal{X}$ (in the NARX space, a box of past inputs and outputs) | `roi` |
| $\mathcal{U}$ | admissible inputs; hard bounds in the design | -- |
| $\mathcal{X}$ | admissible states; soft bounds, penalized in the objective | -- |
| $\mathcal{D}$ | dataset $\{(\vec z(k), \Delta\vec x(k)) : k = 1, \ldots, N\}$ | -- |
| $\mathcal{Q}$ | supporting points that discretize the ROI | `Q` `(nSup x nDim)` |
| $\vec Z$ | regressor points whose output is measured | -- |
| $\hat{\vec Z}^{\mathrm{S}}$ | regressor points whose output comes from the surrogate | -- |
| -- | the matrix the criterion is evaluated on, in physical / normalized units | `Zphys` / `Znorm` `(nSamples x nDim)` |
| $n_z$ | dimension of the regressor space | `nDim` |
| $N$, $N_{\mathrm{sup}}$ | number of samples / supporting points | `nSamples`, `nSup` |

## Points and signals

| Symbol | Meaning | Code |
|---|---|---|
| $\vec z(k)$, $\hat{\vec z}^{\mathrm{S}}(k)$ | one point in regressor space, from the process / predicted by the surrogate | `z` `(1 x nDim)` |
| $\vec q_j$ | one supporting point | `q` `(1 x nDim)` |
| $u(k)$ | process input (the excitation signal) | `u` `(nSteps x nInputs)` |
| $\vec u$ | input sequence $[u(1), \ldots, u(N)]$, the decision variable of the design | `u` |
| $y(k)$ | measured process output | `y` |
| $\hat y^{\mathrm{S}}(k)$ | surrogate output | `yHat` |
| $\vec x(k)$ | process state measured at the instant $u(k)$ is applied | `hist.x` |
| $\hat{\vec x}^{\mathrm{S}}(k)$ | surrogate state trajectory | `xTraj` |
| $\Delta\vec x(k)$ | measured state increment $\vec x(k+1) - \vec x(k)$, the identification target | -- |
| $m_u$, $m_y$ | input / output lag orders of the NARX representation | `nLagU`, `nLagY` |

## Process, model and error

| Symbol | Meaning |
|---|---|
| $f$ | the map to be identified, as a function on the regressor space. For the Duffing study: the one-step state increment as a function of $[u(k), \vec x(k)]$ |
| $\mathcal{P}$ | the process, $\vec x(k+1) = \vec x(k) + f(\vec z(k))$, $y(k) = g(\vec x(k))$, assumed noise-free |
| $\mathcal{M}$ | an identified model of $f$, of any architecture; not the surrogate |
| S | the surrogate used inside the design to predict samples |
| $\omega_f$, $\omega_{\mathcal{M}}$ | moduli of continuity of $f$ and $\mathcal{M}$, on the ROI enlarged by $h$ |
| $e(\vec z)$ | one-step error of $\mathcal{M}$ at $\vec z$, $\lVert f(\vec z) - \mathcal{M}(\vec z) \rVert$ |
| $E_{\mathrm{train}}$ | largest error of $\mathcal{M}$ at its training samples, $\max_k e(\vec z(k))$; 0 for interpolating models |
| $E_{\max}$ | worst-case one-step error of $\mathcal{M}$ over the ROI, $\sup_{\vec z \in \mathcal{R}} e(\vec z)$; the norm is the maximum over the state components, each in units of its ROI half-width. The coverage bound is $E_{\max} \le E_{\mathrm{train}} + \omega_f(h) + \omega_{\mathcal{M}}(h)$ |
| $E_{\max}^{\mathrm{sim}}$ | maximum absolute error of a model's simulated output on a test set, in the units of $y$ (Bouc--Wen study). Not the one-step $E_{\max}$ over the ROI |

## Distances, resolution and the criterion

| Symbol | Meaning | Code |
|---|---|---|
| $d(\vec z)$ | distance from $\vec z$ to its nearest sample, $\min_k \lVert \vec z - \vec z(k) \rVert_2$ in normalized coordinates | -- |
| $d_j$ | distance of the supporting point $\vec q_j$ to its nearest sample, $d_j = d(\vec q_j)$ | `dRaw` `(nSup x 1)` |
| $h$ | fill distance of the samples over the ROI, $\sup_{\vec z \in \mathcal{R}} d(\vec z)$ | `fillDistance` |
| $h_{\mathcal{Q}}$ | fill distance of the supporting-point set itself, the discretization slack | `fillDistanceQ` |
| $h_0$ | demanded resolution under uniform demand | `epsReq0` |
| $d^{\star}$ | distance from the location of the worst one-step error to its nearest training sample, in units of $h$ | `dStar` |
| $\psi_j$ | clipped deviation of $\vec q_j$ beyond the demanded resolution, $\psi_j = \min\{(d_j / h_0 - 1)_+, \psi_{\max}\}$ | `psi` `(nSup x 1)` |
| $\psi_{\max}$ | clipping level of $\psi_j$ | `psiMax` |
| $r$ | outer radius, $r = (1 + \psi_{\max})\,h_0$; beyond it $\psi_j$ saturates | -- |
| $J$ | criterion value, the mean of the deviations, $J = \frac{1}{N_{\mathrm{sup}}} \sum_j \psi_j$; minimized | `J` |
| $\delta_{\mathrm{S}}$ | largest distance between a predicted sample and its realized counterpart, $\max_k \lVert \hat{\vec z}^{\mathrm{S}}(k) - \vec z(k) \rVert_2$ | -- |

"Penalty" always means the soft-constraint term of the optimization objective. The deviation
$\psi_j$ is part of $J$; the penalty is not.

## Loop and optimization

| Symbol | Meaning | Code |
|---|---|---|
| $L$ | prediction horizon: the number of simulated future samples | `nHorizon` |
| -- | decision vector of the parametrization (`nHorizon * nInputs` entries) | `theta` |
| -- | objective minimized by the optimizer: $J$ plus the soft-constraint penalties, never $J$ itself | `objective` |

## Terms

| Term | Meaning |
|---|---|
| receding-horizon control (RHC) design | the proposed method; its signals are RHC signals |
| input design | the general activity; "space-filling design" and "optimal experiment design" keep their literature meaning |
| excitation method, method | any way of generating an input: the RHC design, APRBS, multisine, chirp |
| baseline methods | APRBS, multisine and chirp, the methods the RHC design is compared with |
| signal | one concrete input sequence of a method, of length $N$ |
| realization | one signal of a seeded method (the RHC design, APRBS) |
| dataset | the measured trajectories a signal produces on the process |
| sample | one point $\vec z(k)$ of a dataset in a regressor space |
| measured (value) | the numerical evaluation of an ideally defined quantity (see Conventions) |
| arms (i)--(iii) | the three variants of the RHC design in the Duffing study; in the Bouc--Wen study the variants are named by their lag order $m$ |
