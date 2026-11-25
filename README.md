# From Coverage to Worst-Case Error Bounds: Space-Filling Excitation Design

Code, settings, plotted values and a demo for the paper

> M. Herkersdorf, O. Nelles and K. A. Prajapati,
> *From Coverage to Worst-Case Error Bounds: Space-Filling Excitation Design*,
> submitted to the IFAC Symposium on System Identification (SYSID).

University of Siegen, Automatic Control and Mechatronics.

> **Status:** release **v1.0-submission** (2026-11-23), the state of the submitted paper.

The repository contains:

1. **The method.** `rhcInputDesigner/` designs an excitation signal step by step by receding-horizon
   control (RHC). Candidate inputs are simulated through a surrogate model, and the predicted points
   in the regressor space are scored by the space-filling criterion `J`. The first step of the best
   candidate is applied, and the loop repeats. You can run it on your own system by supplying a
   surrogate and a region of interest.
2. **A demo** that runs in about a minute: `rhcInputDesigner/examples/demo_StateSpace_Offline.m` designs
   a 100-sample excitation offline on a linear state-space surrogate (no real process) and plots
   the input, the surrogate states and the samples in the joint input-state space.
3. **Documentation:** the help text of every function (`help rid.rhc.run`) and the pages in `docs/`:
   - `docs/notation.md`: the notation of the paper and of the code;
   - `docs/method.md`: the method;
   - `docs/duffing.md`: the Duffing study;
   - `docs/boucwen.md`: the Bouc--Wen study;
   - `docs/studies.md`: where the studies' values and settings are, prerequisites.
4. **The plotted values** of Figs. 3, 5 and 6 and Table 2 in `data/`.
5. **A study reference** in `experiments/`: the process models, the settings of the two studies, the
   error measures and the baseline signal generators. It is a reference, **not runnable end to end**
   (`docs/studies.md`).

The method's settings in both studies were fixed in advance and are documented as such. They were
not tuned or validated for robustness; a sensitivity study is not part of this paper.

## Requirements

- MATLAB **R2025b**, the release this repository is tested on.
- Optimization Toolbox and Statistics and Machine Learning Toolbox: the method, the demo and the tests.
- No third-party code is needed. The third-party code that the paper's studies used is **not
  included**; `THIRD_PARTY.md` says what it is and where to get it.

## Quick start

```matlab
cd rhcInputDesigner
setupPath                                   % adds src/, tests/, examples/; never genpath from the root
r = runtests("tests", "IncludeSubfolders", true); disp(r)
run examples/demo_StateSpace_Offline.m      % a small offline design on a state-space surrogate
```

To use the method on your own system, copy the demo and replace its surrogate, regressor
representation, region of interest and (optionally) process; every function has a help text.

## Structure

| Path | Content |
|---|---|
| `rhcInputDesigner/src/` | the method, namespace `+rid` (criterion, RHC loop, surrogates, regressor spaces) |
| `rhcInputDesigner/tests/`, `examples/` | its unit tests (with golden vectors in `tests/data/`) and the demo |
| `experiments/` | the study reference (processes, settings, metrics, baseline signals) |
| `data/` | the plotted values (`data/README.md`) |
| `docs/` | notation, method, the two studies (`studies.md`: where their values and settings are, prerequisites) |
| `THIRD_PARTY.md` | third-party code and data the studies used; none is included |

## Citation

Cite the paper above. `CITATION.cff` describes this software release.

## Related work

The design builds on the RHC framework of

> M. Herkersdorf and O. Nelles, *Online and offline space-filling input design for nonlinear system
> identification: A receding horizon control-based approach*, American Control Conference (ACC), 2026.

## Contact

max.herkersdorf@uni-siegen.de

## License

© 2026 University of Siegen, Institute of Mechanics and Automatic Control (M. Herkersdorf,
O. Nelles). Code, documentation and data: Creative Commons Attribution-NonCommercial-NoDerivatives
4.0 International (CC BY-NC-ND 4.0), see `LICENSE`.
