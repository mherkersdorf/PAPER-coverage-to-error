# Plotted values

The values behind Figs. 3, 5 and 6 and Table 2 of the paper, extracted from the run records named
below. Values are copied as stored, at full precision, and never recomputed. The records themselves
are kept by the authors and are available on request. Licensed CC BY-NC-ND 4.0 like the rest of this
repository (`LICENSE`).

**Names.** The CSVs use the code's names. `eps` is the paper's worst-case distance (fill distance) *h*.
`worst` (`worstGp`, `worstLmn`) is the worst one-step error *E*<sub>max</sub>. `dstarGp` is *d*<sup>⋆</sup>,
the distance from the location of the GP's worst one-step error to its nearest training sample.
`eMaxMultisine` and `eMaxSweep` are *E*<sub>max</sub><sup>sim</sup> on the two Bouc-Wen test sets.
Families: `i`, `ii`, `iii` are the RHC arms (Duffing), `rhcM1`-`rhcM3` the RHC designs with lag
order *m* = 1-3 (Bouc-Wen), and `aprbs`, `multisine`, `chirp` the baselines.

| File | Figure / table | Rows | Columns | Source record |
|---|---|---|---|---|
| `fig3_epsVsWorst.csv` | Fig. 3 (a) GP, (b) LMN: one marker per run | 1 600 | `family` (ii, aprbs), `seed` (1-100), `N` (123, 505, 1947, 8000), `model` (gp, lmn), `eps`, `worst` | `p2Eval_20261019T214542` (`measures.csv`) |
| `fig5_studyB.csv` | Fig. 5 (a)-(f): one row per panel, family and *N* | 2 880 | `panel`, `family`, `N`, `median`, `q25`, `q75`, `min`, `max` | `figP2StudyB_20261020T114708` (`series.csv`) |
| `table2_reach.csv` | Table 2 | 15 | `measure`, `statistic`, `reference` (the best baseline at *N* = 8000), `referenceValue`, `arm`, `firstN` (the samples the arm needs), `factor` (8000 / `firstN`) | `figP2StudyB_20261020T114708` (`reachTable.csv`) |
| `fig6_boucwen.csv` | Fig. 6 (a) multisine test set, (b) sine sweep | 240 | `panel`, `method`, `N`, `median`, `q25`, `q75` | `figP3Own_20261109T000048` (`series.csv`) |
| `boucwen_measures_summary.csv` | every Bouc-Wen measure of the seeded methods (`rhcM1`-`rhcM3`, `aprbs`) and length behind Fig. 6 | 400 | `method`, `N`, `measure`, `median`, `q25`, `q75`, `nRealizations`, `nValues`, `nInf` | `p3OwnMeasures_20261109T000002` (`summary.csv`) |

Notes:

- **Seeded methods** (the RHC arms and designs, APRBS): median and quartiles over seeds (Duffing:
  100) or realizations (Bouc-Wen: 30). **Multisine and chirp**: one independent signal per size, so
  `median` is that signal's value and `q25`, `q75` are `NaN`.
- **Fig. 5 (a)**: the regular-grid reference line in the panel is the formula √3 / (2 *N*<sup>1/3</sup>),
  not data, so it is not in the CSV.
- **Fig. 6**: a training that does not count as a result under the study's collection rule, or
  whose model cannot be re-simulated on a test set, scores `Inf` for *E*<sub>max</sub><sup>sim</sup>. The medians and quartiles include these
  values; `nInf` in `boucwen_measures_summary.csv` counts them. The paper's upper y limits (1e-2,
  1e-3) clip the plotted lines but do not change the values.
- **Fig. 4** (median gap and error maps) is not included in this release. The maps are not stored
  numerically in the figure record.
