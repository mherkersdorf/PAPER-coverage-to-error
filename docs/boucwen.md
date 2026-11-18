# The Bouc--Wen study (paper, Sec. 5)

The full setup of the second study, and the comparison over the full record that the paper's figure
summarizes with medians. The settings are in `experiments/p3/scripts/private/p3Config.m` and the
process files in `experiments/processes/` (see [`studies.md`](studies.md)); the plotted values are in
`data/`. The benchmark data and its integrator are not included; see `THIRD_PARTY.md`.

## 1 What the study asks

The RHC design is tested under practical conditions on the Bouc--Wen benchmark (Noël and Schoukens,
2016): only the input and a noise-distorted output are measured, only a linear model is known, and the
test data are fixed by the benchmark. The process is a hysteretic oscillator whose internal variable
cannot be measured. Designing in NARX spaces with lag orders $m = 1, 2, 3$ shows whether the outcome
depends on the choice of space.

- **Excitation is varied, not the identification method.** One model structure, one training setting,
  for every excitation method, length and realization. The models are local model state-space networks
  (LMSSN).
- **Output-only.** There are no states in which to measure the fill distance of the true process, and
  the hysteretic variable is not determined by the NARX regressor, so the one-step error $E_{\max}$ over
  the ROI is not defined. The models are graded by $E_{\max}^{\mathrm{sim}}$, the maximum absolute error
  of the simulated output on each test set.
- The study is empirical. It is not evidence for the coverage bound of Sec. 2, and the lag-order
  comparison is heuristic insight, not a general result.

## 2 Process and benchmark conditions

The conditions are those the benchmark prescribes.

- **Sampling and simulation.** Working sampling frequency 750 Hz. The integrator runs at 15 000 Hz
  (the input is upsampled by 20); the output is low-pass filtered and decimated back to 750 Hz.
  Initial conditions are zero. A signal that exists only at 750 Hz reaches 15 kHz by band-limited
  interpolation (`interp(u, 20)`), as the benchmark's own test data do; a zero-order hold would add a
  phase lag of about half a sample. The 750-Hz samples stay within the input range; the interpolated
  input may overshoot between samples.
- **Noise.** Training and validation data carry the benchmark's measurement noise on the output:
  Gaussian, band-limited to 0--375 Hz, RMS $8 \cdot 10^{-3}$ mm (about 40 dB SNR at 50 N). The input is
  noiseless. Noiseless simulation is not used for training or validation. The test data are noiseless.
- **Test data.** Two official test sets: a random-phase multisine (one steady-state period of
  8192 samples, 5--150 Hz, RMS 50 N) and a sine sweep (153 000 samples, 20--50 Hz, amplitude 40 N). They
  are used for **evaluation only**: not for training, validation, model selection, tuning, the
  surrogate, the ROI or anything else in the design, as the benchmark requires.

## 3 What each method knows

A fair comparison fixes what each method may know about the process:

1. **Shared prior knowledge:** the benchmark document's description of its test classes (random-phase
   multisines, 5--150 Hz, RMS 50 N, 8192 samples; the sweep's settings).
2. **The ROI**, from a pilot experiment on the process: 100 multisines of that class, simulated on the
   Bouc--Wen simulator. The input bounds are the medians of the realizations' minima and maxima, widened
   by 5 % of the span; the output bounds likewise. The 100 realizations only make the bounds
   independent of a seed. The ROI stands for the operating range a practitioner knows (input range and
   admissible output range). The pilot data reach the methods only as these four bounds.
3. **The BLA**, from a second pilot experiment: 4 new random-phase multisines (5--150 Hz), each scaled
   to fill the ROI's input range.
4. **Design and scaling use only the ROI and the BLA.** The RHC design fills the ROI with the BLA as its
   surrogate. Every baseline signal is scaled so that its input and its BLA-predicted output stay inside
   the ROI. The APRBS hold comes from the BLA.

No method uses process data beyond the two pilots and its own validation signal. One asymmetry is
inherent to model-based design: the BLA pilot (4 realizations of 6 periods of 8192 samples) is the RHC
design's surrogate, while the baselines use the BLA only for the APRBS hold and the output check. The
pilots and the validation signals are not counted in $N$.

## 4 Design spaces

Contiguous input and output lag orders without feedthrough $u(k)$, input delay 1:

| $m$ | $n_z$ | Regressor |
|---|---|---|
| 1 | 2 | $[u(k{-}1),\, y(k{-}1)]$ |
| 2 | 4 | $[u(k{-}1),\, u(k{-}2),\, y(k{-}1),\, y(k{-}2)]$ |
| 3 | 6 | $[u(k{-}1), \ldots, u(k{-}3),\, y(k{-}1), \ldots, y(k{-}3)]$ |

Spaces with feedthrough ($n_z = 3, 5, 7$) are not run.

## 5 ROI and surrogate

**ROI.** From the 100 pilot multisines (section 3): $u \in [-204.1,\, 207.1]$ N and
$y \in [-2.357,\, 2.342]$ mm. Every lag coordinate of $u$ (of $y$) has the same range. All distances are
computed with every axis scaled to width 1.

**Surrogate.** The best linear approximation (BLA) of the Bouc--Wen process, of order 3, estimated by the
robust method (Pintelon and Schoukens, 2012):

- $M = 4$ random-phase multisines, 5--150 Hz, one period of 8192 samples, zero mean and scaled to the
  largest amplitude inside the ROI's input range (an RMS of about 45--55 N). They are simulated for six
  periods with the benchmark noise, and the first period is dropped.
- The nonparametric BLA is the mean over the realizations of each realization's period-averaged
  frequency response. The noise variance comes from the variation between periods, the total variance
  (noise plus stochastic nonlinear distortion) from the variation between realizations.
- A parametric model of order 3 is fitted to it. Of the two fitted candidates (the plain subspace model
  and its refinement), the refined one was kept, because it is stable and its cost is lower.
- Its poles are a real pole at $|z| = 0.973$ (about 3.3 Hz, below the excited band) and a pair at 40 Hz,
  $|z| = 0.971$. The stochastic nonlinear distortion exceeds the noise by 11--40 dB. A direct term of
  $7.9 \cdot 10^{-8}$ m/N comes from the sampling chain, not from the process, and the test data carry
  it too. The model is in `experiments/p3/surrogates/bla.json`.
- The BLA is never estimated on the test data: the surrogate shapes the RHC signals, so test
  information would leak into the design and into no baseline.

## 6 Excitation signals

| Method | Realizations | Generation |
|---|---|---|
| RHC, $m = 1, 2, 3$ | 30 per lag order | $N = 5000$ per realization; shorter lengths are nested prefixes |
| APRBS | 30 | levels drawn uniformly from the input range within the output limit; mean hold = one time constant of the BLA's dominant pole, 36.75 samples |
| Multisine | 1 per length | one full period of $N$ samples, Schroeder phases, band 5--150 Hz, zero mean, scaled to the largest amplitude within the limits below |
| Chirp | 1 per length | one linear sweep from 5 to 150 Hz over $N$ samples, zero mean, scaled like the multisine |

- **Limits.** Every signal, validation signals included, keeps its input inside the ROI's input range and
  its BLA-predicted output inside the ROI's output range. Each baseline gets the largest amplitude that
  meets both. For the RHC design, the output limit is a soft penalty on its BLA-predicted output. There
  is **no rate limit**.
- **Band.** 5--150 Hz, the band the benchmark states for its multisine test set. The sweep's 20--50 Hz
  lies inside.
- **Realizations.** 40 realizations per seeded method were generated, and realizations 1--30 are the
  sample, fixed in advance; 31--40 enter no reported number. A training that diverges, keeps splitting or
  hits a time limit stays in the sample and is reported. Seeds are distinct for every random signal, and
  validation seeds are disjoint from every training seed.

## 7 Space-filling settings

| $n_z$ | $h_0$ | Points per axis | $N_{\mathrm{sup}}$ | $h_{\mathcal{Q}}$ | $r$ | $\psi_{\max}$ |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.010 | 148 | 21 904 | 0.0048 | 0.42 | 41.4 |
| 4 | 0.119 | 12 | 20 736 | 0.083 | 0.60 | 4.05 |
| 6 | 0.296 | 5 | 15 625 | 0.245 | 0.73 | 1.48 |

- **Demanded resolution.** Matched to the budget: $h_0 = \sqrt{n_z} / (2 N^{1/n_z})$ at $N = 5000$, the
  largest gap that a perfect grid of $N$ samples leaves in that space. It is a target scale, not a
  guarantee, and it shows the curse of dimensionality.
- **Supporting points.** The largest cell-centered full grid with at most 22 000 points. Every grid is
  finer than $h_0$. Short smoke runs showed that the time per step and supporting point does not grow
  with $n_z$ (0.8, 1.0 and 0.9 ms). In 6-D the grid is fixed at 5 points per axis by the demand and the
  cost: 4 per axis would give $h_{\mathcal{Q}} = 0.306 > h_0$, and 6 would give 46 656 points. This
  differs from the Duffing study, which holds the fill distance equal across spaces.
- **Outer radius.** $r = 0.3\sqrt{n_z}$, which is 60 % of the ROI half-diagonal, and
  $\psi_{\max} = r/h_0 - 1$. Fixed in advance. $r$ stays well inside the box, so regions the process
  cannot reach saturate at $\psi_{\max}$ instead of dominating $J$. In 2-D the criterion weights hole
  size strongly; in 6-D it behaves more like counting unresolved supporting points.
- **Horizon.** $L = 40$. A rule declared in advance gave $L = 20$ (median fill-distance ratio of
  $L = 20$ to $L = 40$: 0.962, in a pilot at $N = 1000$ with 3 seeds). $L = 40$ was chosen because the
  pilot cannot test whether a longer horizon reaches hard-to-reach regions better once the easy ones
  are covered. It costs about twice as much per step.
- **Mode.** Budget mode, $N = 5000$.

## 8 Training lengths and validation

- **Lengths.** 20, log-spaced from 500 to 5000. RHC and APRBS use the first $N$ samples of the same
  realization. Multisine and chirp have their own signal per length.
- **Validation.** One signal of 1200 samples per realization, of that method's own kind, with its own
  seed. The same signal is used for all lengths of that realization. Multisine and chirp keep one
  validation signal each: a random-phase multisine (5--150 Hz) and a downward chirp (150 to 5 Hz). Same
  input range, output limit and noise as the training data. A validation signal never repeats a
  training signal. It decides each training's termination and model size. Giving every realization of
  a method its own validation signal avoids a bias that the spread over realizations would not show.

## 9 LMSSN training

- **One structure and one training setting** for every method, length and realization, and **one
  training per signal**. The training is deterministic for a given data set: the starting model is
  computed from the data, so only signal randomness enters. Repeated initializations would be better
  and were not run.
- **Termination.** Automatic stop by validation (3 splits without improvement), with the model selected
  by validation. In a pilot (APRBS, $N = 5000$) the automatic stop selected the same model as a fixed
  cap of 30 local models, at 29 % of its cost.
- **Bounds.** A safety cap of 30 local models; a hard deadline of 30 h per training, after which the
  model is selected from the last saved history (the result counts, flagged); and at most 760 000
  function evaluations per child optimization. In the first runs some trainings hung inside one
  optimizer iteration for days, which is why the deadline exists.
- **Start of every simulation.** A pre-roll of 1500 zero-input samples, with weight 0, in front of
  every training set, validation signal and the sine-sweep test input, so that the model settles to its
  own rest before the scored data begin, as the data start from rest. Against an estimated start state,
  the pre-roll halved the validation error in a pilot, at 5.6 times the training time.

## 10 Measures, flags and failures

- **$E_{\max}^{\mathrm{sim}}$**, the maximum absolute error of the simulated output, on the multisine
  test set and on the sine sweep: the measure of the paper's Fig. 6. The multisine test set is simulated
  over two periods and the second is scored; the sine sweep is scored after the pre-roll.
- **RMSE of the simulated output** (the benchmark's figure of merit) on both test sets, and their mean:
  computed and stored, not in the paper.
- **Counted, flagged and failed runs.** A training counts with its test error if the selected model has
  a finite error on both test sets, whatever happened on the way. These are reported per method and
  length as **flags**, not failures: the child optimization failed, the training error is NaN or
  large, the model has not settled to rest at the end of the pre-roll, the selected model's
  simulation diverges on a test set (it enters with its error), or the deadline was hit. A run
  **fails** (infinite error in the medians) if the training or the test simulation errors out, the
  test output is not finite, no model exists after the deadline, or the run was cancelled.
- Methods are compared by their test errors, never by the criterion value $J$, which is not comparable
  across spaces.

## 11 The full record behind Fig. 6

Fig. 6 shows the RHC medians with IQR bands over 30 realizations (nested prefixes of one run, so
smooth over $N$) next to single multisine and chirp signals generated per length. A comparison over
all counted trainings shows how reliably each method gives a small error: 595--600 trainings per RHC
lag order, 599 for the APRBS and 20 each for multisine and chirp. The share of trainings with
$E_{\max}^{\mathrm{sim}} > 3 \cdot 10^{-3}$ m, more than 10 times the RHC median, is:

| Method | Multisine test set | Sine sweep |
|---|---:|---:|
| RHC, $m = 1$ | 0 % | 0.2 % |
| RHC, $m = 2$ | 1.2 % | 1.8 % |
| RHC, $m = 3$ | 0.5 % | 0.8 % |
| APRBS | 1.2 % | 0.5 % |
| Multisine | 20 % | 5 % |
| Chirp | 10 % | 10 % |

- **Diverged simulations.** Every $E_{\max}^{\mathrm{sim}} > 10^{-2}$ m is a diverged test simulation,
  in every method: 1--4 per RHC lag order, 2 for the APRBS, 1 for the multisine (at $N = 2728$) and 1 for
  the chirp (at $N = 1318$). The baselines' jumps on the sine sweep in Fig. 6 are these single events.
- **Not counted** in these shares: 5 trainings for $m = 2$, 1 for $m = 3$, 1 for the APRBS.
- **Multisine test set.** The RHC designs lie below every baseline in 58 of 60 pairs of lag order and
  length. The exceptions are $m = 2$ at $N = 719$ and $N = 916$.
- **Sine sweep.** The RHC designs lie below the APRBS at 18, 14 and 20 of 20 lengths for $m = 1, 2, 3$.
  The multisine or the chirp lies below the best RHC median at 9 of 20 lengths. With one signal per
  length, 20 baseline signals cannot resolve a difference on this test set.

## 12 What the BLA does to the signals

The BLA is the only model of the process that the design and the scaling of the signals may use
(section 3). It is a poor model of the Bouc--Wen process, and two consequences follow for every method.

- **The RHC design with a poor surrogate.** The design drives the BLA-predicted output to the edges of
  the ROI, where the BLA's error is largest. It keeps the BLA-predicted output inside the ROI (within
  0.13 %), but the true output of the generated data sets exceeds the ROI's output range by up to 44 %
  for $m = 3$ (33 % in a pilot), against up to 5 % for the baselines. The data are nevertheless
  better: the RHC models have the smaller errors in Fig. 6 and in section 11, even though the surrogate
  is poor.
- **Staying within the ROI is something a practitioner has to consider.** The BLA overestimates the
  response near the resonance by about 30 %. Every baseline is scaled with the same BLA and the same
  knowledge as the RHC design, so that its predicted output stays inside the ROI. The multisine and the
  chirp therefore run at an RMS of about 37 N ($N = 5000$), below the test class's 50 N, and use only
  32--34 % of the input range. Their true outputs use about 76 % of the output range, the APRBS's
  91--105 %.

## 13 Known limitations

- **One training per signal.** The spread over realizations is signal randomness only. Multisine and
  chirp have no spread.
- **The pilots.** The BLA pilot, the ROI pilot and the validation signals are not counted in $N$. The
  BLA is used by the RHC design as its surrogate and by the baselines to set the APRBS hold and to
  check that their signals stay inside the ROI.
- **No certificate in 6-D.** The supporting-point fill distance (0.245) is close to the demand (0.296),
  so the criterion's zero says little there. The paper's claims for this study are empirical.
- **Training time is bounded at 30 h.** A training stopped there keeps the best model of its completed
  splits and may be less trained than without the bound.
