# HereditaryUnit

Analytical figure code for *When Can a Hereditary Unit Function as an Evolutionary Gene?*

## Run

Requires R (tested with R 4.6.1). Only base R is used; no packages need to be installed.

```sh
Rscript genes8PT_analytical_figures.R figures
```

The script writes PDF and 600-dpi PNG versions of Figures 1 and 2 (recombination only, physical span in Mb) and S1 and S2 (mutation extension, physical span in Mb). It also writes a CSV of illustrative human-scale bounds and R session information.

## Model and parameters

The main figures omit mutation. Recombination follows the Haldane endpoint map, c = (1 - exp(-2 L))/2, with map length L in Morgans. Physical spans use the illustrative uniform calibration of 1 cM/Mb; this is not a universal recombination map. 100 cM = 1 Morgan. For short intervals, 1 cM corresponds approximately to 1% endpoint recombination; the equivalence is not exact at larger distances.

The half-sweep retention model uses closed-form continuous two-class dynamics and explicit analytical boundaries, without stochastic simulations or integer-generation rounding. Mutation figures assume a disruptive per-base mutation probability of 1.2e-8 and that every base in the physical span is required; their scalar root calculations evaluate analytical functions rather than simulate trajectories. These are conditional illustrative models, not estimates of measured gene lengths or a full multilocus fixation model.

Parameters are near the beginning of the script. Changing the recombination map or mutation assumptions changes physical-span results. PDF figures retain vector curves. Legends have enlarged text and sit immediately below the panels.

The script checks the defining selection and retention thresholds and verifies integrated exposure by independent numerical quadrature before drawing figures. Additional numerical checks can be run with `Rscript validation.R` from this directory.
