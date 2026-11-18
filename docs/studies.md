# The two studies of the paper

This repository gives you the method and a demo. The paper's two studies were run on a cluster and
are not meant to be rerun from here. This page says what they were, where their plotted values are
and where their settings are.

## What the studies were

- **Duffing study (Sec. 4).** RHC designs of an excitation signal for a Duffing oscillator
  (the RHC arms (i), (ii) and (iii)) compared with an APRBS, a linear chirp and a Schroeder
  multisine. Seeds 1-100 for every seeded family; the RHC arms and the APRBS are paired by seed. The
  designs are evaluated with a local GP and a local model network (LMN); the measures are the fill
  distance and the worst one-step error. Details: [duffing.md](duffing.md).
- **Bouc-Wen study (Sec. 5).** RHC designs in NARX regressor spaces of lag order *m* = 1, 2, 3
  (*N* = 500 to 5000) compared with an APRBS, a multisine and a chirp, on the Bouc-Wen benchmark. A
  model (LMSSN) is trained on each signal and tested on the benchmark's multisine and sine-sweep
  test sets; the paper evaluates 30 realizations per method. Details: [boucwen.md](boucwen.md).

## Plotted values

`data/` holds the values behind Figs. 3, 5 and 6 and Table 2, as stored, with the names of the code.
`data/README.md` describes each file and its columns. Fig. 4 (the median gap and error maps) is not
included, and the generation records of arm (iii) were not retained.

## Settings

The settings of both studies were fixed in advance and are documented as such; they were not tuned
or validated for robustness. They are in `experiments/`, which is a **reference, not runnable end
to end**: the drivers, figure scripts and run records are not part of this repository.

| What | Where |
|---|---|
| Duffing process, its parameters and the linear surrogate S1 | `experiments/processes/` (the `p2Duffing` files), `experiments/p2/surrogates/S1.json` |
| Duffing RHC generation settings | `experiments/p2/scripts/private/p2GenConfig.m` |
| Duffing evaluation settings (test signal, measures) | `experiments/p2/scripts/private/p2EvalConfig.m`, `lolimotPilotConfig.m`, `p2EvaluationSizes.m` |
| Local GP of the evaluation (Matérn kernel) | `experiments/models/gpMaternFit.m`, `gpMaternPredict.m`, `private/gpMaternKernel.m` |
| Size of the LMN as a function of *N* | `experiments/p2/scripts/private/p2EvalLmnSize.m` |
| Bouc-Wen process (simulation) and its BLA surrogate | `experiments/processes/p3BoucWenSimulate.m`, `p3BoucWenBla.m`, `experiments/p3/surrogates/bla.json` |
| Bouc-Wen settings (spaces, signals, limits, seeds, training) | `experiments/p3/scripts/private/p3Config.m` |
| Error measures | `experiments/metrics/` |
| Baseline signals (APRBS, chirp, multisine) | `experiments/signals/` |

Some of these files call code that is not included (the third-party code in `THIRD_PARTY.md`, and the
study drivers). Read them for the settings; do not expect them to run.

## Prerequisites of the demo and the tests

- MATLAB **R2025b**, the release this repository is tested on.
- Optimization Toolbox and Statistics and Machine Learning Toolbox. No third-party code is needed.
- The demo: `rhcInputDesigner/examples/demo_StateSpace_Offline.m`, after `setupPath` in
  `rhcInputDesigner/`. It writes no file.
- The tests: `matlab -batch "cd rhcInputDesigner; setupPath; r = runtests('tests','IncludeSubfolders',true); disp(r)"`.
