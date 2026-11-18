# Golden vectors for the criterion

Fifteen frozen cases: inputs, the expected `J`, and the intermediates that produced it. They exist
to catch a class of error that no property test can reach -- an implementation that is *consistently
wrong*.

## Why an independent implementation

Every property of the criterion's test contract is satisfied by a criterion with a systematic error in it.
Drop the `- 1` from the penalty and permutation invariance, monotonicity, the scaling equivalence,
boundedness and locality all still hold; only the value is wrong, everywhere. The only thing that
catches that is a value computed from a **second, independent reading of the specification**.

`generate_golden.py` is that second reading. It was written, in Python, using no numerical library, **before any MATLAB existed in `src/+rid/+crit`**. It
is not a translation of the implementation and must never be "corrected" by consulting it. If the
two disagree, one of them is wrong, and finding out which is the point.

## Rules

- **Never hand-edit a value to make a test pass.** If a golden value looks wrong, say so and stop.
- Regenerate only when the *specification* changes, by re-running the generator -- and commit the
  regenerated files **alone**, with a one-line reason, never together with a code change
 .
- The generator is not run by the MATLAB test suite. It is kept for provenance and regeneration.

## File format

One JSON file per case:

```
name, description
source     : which document, which generator, its hash prefix, the date
inputs     : Znorm, Q, relevance, complexity, noiseLevel, epsReq0, beta, betaNoise,
             psiMax, pNorm, mode, fillDistanceQ
expected   : refinement, epsReq, dRaw, dTilde, psi, J
```

Arrays are always numeric. Non-finite values appear only as the **scalar string** `"Inf"` (used for
`psiMax` in `psiMaxInf`), because MATLAB's `jsondecode` rejects a bare `Infinity` literal; the reader
converts it with `str2double`. Where an intermediate would be infinite -- `emptySamples`, in which
every `dRaw` is `Inf` -- the intermediate arrays are **omitted** and only `J` is given.

## Coverage

| Case | What it pins down |
|---|---|
| `singleSampleSingleSupport`, `grid1dOneSample` | values computable by hand |
| `allResolved` | `J = 0` exactly when everything is resolved (P4) |
| `allClipped`, `emptySamples` | saturation at `psiMax` (P5) |
| `mixedRegimes` | resolved, graded and clipped points in one evaluation |
| `fieldsNonNominal` | refinement clipped at 1 for below-nominal fields |
| `betaNoiseDiffers` | the noise field's own exponent |
| `pNormTwo`, `pNormLarge` | the aggregation exponent |
| `psiMaxInf` | the unclipped, distance-based limit (section 10) |
| `floorInactive`, `floorActive` | the certification floor |
| `random2d`, `random3d` | the vectorized path in 2-D and 3-D, with non-nominal fields |

Coordinates in the random cases are stored explicitly, so no case depends on any random number
generator.

## Reproducibility: what the stored values actually are

**Re-running the generator does not reproduce these files bit-for-bit, and that is not a defect in
the files.** Measured on 2026-09-04 (Python 3.12.3, macOS arm64) against the committed
`generate_golden.py`, whose SHA-256 prefix `71d25d2f85357706` matches the `generatorSha256Prefix`
stored in every case: 13 of 15 files reproduce exactly, and 9 of 2450 numeric values differ across
the other two, by one or two ULP (worst relative difference `5.05e-16`). The test tolerance is
`1e-12`, four orders larger, so every case passes either way.

The two files differ for **different reasons**, and only one of them is drift:

- **`grid1dOneSample.json` stores the exact answer, and the generator is the weaker route.** This is
  one of the hand-computable cases: `epsReq = 0.25` throughout, `psi = [1, .6, .2, 0, 0, 0, 0, 0,
  .2, .6, 1]`, so `J = 3.6/11 = 18/55`. The stored `J` is `float(18/55)` -- the correctly rounded
  exact value, confirmed with `fractions.Fraction`. The generator's `sum(psi_floats)/11` lands one
  ULP away. **Do not "fix" this by regenerating**: a golden that reproduces the implementation's
  accumulation order is a weaker golden, and independence is the entire point of this directory.

- **`random3d.json` is genuine platform drift, and is not reproducible here by any route.** Its four
  differing `dRaw` values are neither this machine's `math.sqrt(sum((zi-qi)**2))` nor the correctly
  rounded exact value: compared against exact rational arithmetic, three of the four are *further*
  from the true value than the naive route and one is closer. That asymmetry is the signature of the
  same naive expression evaluated with different rounding, not of a more accurate computation --
  `(zi-qi)**2` goes through libm `pow`, which IEEE-754 does not require to be correctly rounded.
  `math.dist` and `math.hypot`, which are correctly rounded, differ from the stored values in ten
  places rather than four, so switching to them would change more of this file, not less.

**`generatorSha256Prefix` therefore identifies the generator version; it does not promise that
re-running that generator returns these bytes.** Read it as provenance, not as a reproduction
guarantee.

## Regenerating

```bash
python3 generate_golden.py
```

Writes the JSON files beside the script and prints one line per case.

**Regenerating overwrites `grid1dOneSample.json`'s exact `J` with a less accurate one and rewrites
four values in `random3d.json` to this machine's rounding.** Neither changes any test outcome, and
both lose information. Regenerate only when a case is genuinely added or a specification reading
changes -- never as a side effect of an unrelated commit -- and say in the commit message which
values moved and why.
