# Momentum-resolution tuning

Copy `TestKinfit.cc` and `TestKinfit.hh` into the existing fitter directory,
alongside its `src/` and `include/` folders. In ROOT:

```cpp
.L TestKinfit.cc+
auto phi = TunePhiPhiResolutions();
auto xi  = TuneXiResolutions();
auto eta = TuneMissingEtaResolutions();
```

Each function defaults to 10,000 events, at most **four fit rounds**, and an
initial momentum-resolution estimate of **0.5 times the true resolution**.
Iteration 1 uses that initial guess; subsequent iterations use the updated
estimates. The corresponding initial momentum variances are one quarter of
their true values. Angular resolutions and their covariance entries stay fixed.

The true resolutions are copied from the existing `ResP`, `ResPi1`, etc.
settings at the start of each function. Events are generated and smeared once,
then the exact same measurements are fitted in every iteration. Momentum
variances use the event's true momentum, matching the existing test framework.
The calibration functions tune diagonal momentum variances.

| Function | Estimated parameters |
| --- | --- |
| `TunePhiPhiResolutions()` | `ResPbarBeamEstm`; one shared `ResPEstm` for the four kaons |
| `TuneXiResolutions()` | `ResPEstm`, `ResPi1Estm`, `ResPi2Estm` |
| `TuneMissingEtaResolutions()` | `ResPEstm`, `ResPi1Estm`, `ResKBeamEstm` |

For each estimated parameter, a Gaussian is fitted to the momentum correction
pulls. The four kaon momentum-pull distributions are pooled for their shared
parameter. A parameter is stable when **0.9 < fitted width < 1.1**. An unstable
parameter with a valid fit is updated as:

```cpp
ResEstm *= fittedPullWidth;
```

Thus a width of 2 doubles the estimated resolution. A stable parameter keeps
its estimate for that round and is checked again if other parameters change.
The function stops when every parameter is stable in the same round, or the
iteration limit is reached. Failed fits are counted; a failed Gaussian fit
does not update its parameter or count as stable. Forbidden kinematics raise
an exception.

The returned result contains `converged`, `groups`, the final tested
`estimates`, and the per-round `iterations`. The corresponding global
`Res*Estm` settings are set to the final tested values. Each reaction has its
own calibration; running another one can update shared estimated parameters.

```cpp
// Events, maximum fit rounds, initial scale, random seed:
auto result = TuneXiResolutions(20000, 4, 0.5, 12345);
std::cout << result.converged << std::endl;

// PhiPhi additionally accepts beam momentum in GeV/c; -1 uses pPbar:
auto resultPhi = TunePhiPhiResolutions(10000, 4, 0.5, 12345, 1.0);
```

The optional final argument is the output directory. Defaults are
`figs/TunePhiPhiResolutions`, `figs/TuneXiResolutions`, and
`figs/TuneMissingEtaResolutions`.

Each directory contains:

- `ResolutionIterations.pdf` and `.png`: estimated resolutions by fit round,
  with dashed true-resolution lines; correction-pull widths with the stable band.
- `ResolutionIterations.tsv`: numerical estimates, fitted means and widths,
  uncertainties, stability flags, and event counts for every round.
- `MomentumPulls_IterationN.pdf`: Gaussian pull fits with means and widths.
- `tuning.root`: pull histograms with fitted functions, canvases, and history graphs.

`examples/` contains the three 10,000-event runs using seed 12345 and the
default settings. All three reached the width band after two fit rounds with
zero failed fits. The stopping rule uses pull widths; the true-resolution
lines provide a reference for judging the resulting estimates.
