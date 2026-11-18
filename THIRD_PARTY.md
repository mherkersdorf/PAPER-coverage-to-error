# Third-party code and data used by the studies

**None of it is included in this repository.** The method (`rhcInputDesigner/`), the demo and the
tests need none of it. The study reference in `experiments/` mentions these items because the
paper's studies used them.

| Item | What it is | Why the studies used it | Not included; where to get it |
|---|---|---|---|
| LMNtool 1.5.3 | local model networks (LOLIMOT), O. Nelles and co-workers, University of Siegen | the LMN model that evaluates the Duffing designs | not included. The institute's copy of 1.5.3 was used; the public lmn-tool, <https://github.com/mrtsiegen/lmn-tool>, is a different version |
| LMSSN | local model state-space networks, own code (O. Nelles, M. Schüssler) | the model trained on each Bouc-Wen signal | not included, no link |
| APRBS.m | APRBS / PRBS generator, C. Illg, University of Siegen | the APRBS baseline signals | not included, no link |
| `multisine.m` | Schroeder-phase multisine generator, B. Holmes (adapted from M. van Walstijn) | the multisine baseline signals | not included (all rights reserved); <https://bholmesqub.github.io/DAFx19/big-muff-pi-tonestack/multisine-func/> |
| matlab2tikz | TikZ export of MATLAB figures, N. Schlömer and contributors | the TikZ files of the paper's figures | not included; <https://github.com/matlab2tikz/matlab2tikz> |
| Bouc-Wen benchmark | hysteretic benchmark: integrator and test signals, J. P. Noël and M. Schoukens; CC BY-SA 4.0 | the process of the Bouc-Wen study and its test sets | not included; nonlinearbenchmark.org, doi:10.4121/12967592 |
