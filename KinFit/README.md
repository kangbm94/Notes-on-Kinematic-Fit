# KinFit

KinFit is a ROOT-based framework for kinematic fitting. It adjusts measured
parameters within their uncertainties to satisfy constraints defined for a
particular reaction or decay. Each concrete fitter supplies the kinematics;
the common base class handles the iterative fit and its diagnostics.

## Project structure

| Files or directory | Purpose |
| --- | --- |
| `include/KinFit.hh`, `src/KinFit.cc` | The common base class, `KinematicFitter`: covariance handling, iteration, chi-square, pulls, and result storage. |
| `include/<Fitter>.hh`, `src/<Fitter>.cc` | The declaration and implementation of each concrete fitter. |
| `TestKinfit.hh` | Generic tools for smearing measurements and plotting residuals, pulls, constraints, and fit-quality distributions. |
| `TestKinfit.cc` | Reaction-specific generation and performance-test examples. |
| `figs/` | Output plots and ROOT histogram files from the tests. |

The concrete classes inherit from `KinematicFitter`, the abstract/interface
layer implemented in the `KinFit` files:

| Concrete class | Supported kinematics |
| --- | --- |
| `FourVectorFitter` | A two-body decay, `R -> P + Q`, with a specified parent mass. |
| `CascadeFitter` | A cascade decay, `X -> L + R`, followed by `L -> P + Q`, with specified intermediate and parent masses. |
| `MissingMassFitter` | A reaction `Beam + Target -> R + Missing`, followed by `R -> P + Q`, with specified missing and parent masses. |
| `BalanceFitter` | Four-momentum conservation for `Beam + Target -> P1 + P2 + ...`, with a stationary target. All beam and final-particle parameters are measured; there are no unknown parameters. |

The current base class has default bodies for its virtual hooks, so it is not
formally abstract in C++. Its intended role is the common interface: use a
concrete fitter for an actual event.

## Using a fitter

For each event:

1. Construct the concrete fitter with the measured particles or other inputs.
2. Set any reaction-specific quantities, such as parent or missing masses.
3. Supply the diagonal covariance entries with `SetVariance()`.
4. If measurements are correlated, add the off-diagonal entries with
   `AddOffdiagonals()`.
5. Run `DoKinematicFit()`.
6. Retrieve the fitted parameters and construct the observables needed by your
   analysis, or use the concrete class's result getters.

For example, after loading the balance fitter source in a ROOT macro:

```cpp
#include "src/BalanceFitter.cc"

// User inputs: beam, target, finalParticles, and their resolutions in sigma.
BalanceFitter fitter(beam, target, finalParticles);

std::vector<double> variance(sigma.size());
for (std::size_t i = 0; i < sigma.size(); ++i)
    variance[i] = sigma[i] * sigma[i];

fitter.SetVariance(variance.data());

// Optional: covariance is a symmetric TMatrixD in the same parameter order.
// fitter.AddOffdiagonals(covariance);

double chi2 = fitter.DoKinematicFit();
double pvalue = fitter.GetPValue();
auto pulls = fitter.GetPull();
auto fittedMeasurements = fitter.GetBestMeasurements();
auto fittedParticles = fitter.GetFittedLV();
```

The `sigma` values and covariance rows must follow the fitter's measured
parameter order. In `BalanceFitter`, that order is:

```text
p_beam, theta_beam, phi_beam, p_P1, theta_P1, phi_P1, ...
```

There are `3 * (1 + number_of_final_particles)` measured parameters. Angles
are in radians; momenta and masses must use consistent units.
`GetFittedLV()` returns the fitted beam followed by the fitted final particles.
Other concrete fitters define their own output order in their headers.

`SetVariance()` accepts **variances, not standard deviations**. Pass one
positive variance per measured parameter. `AddOffdiagonals()` accepts a
covariance matrix and ignores its diagonal entries. Supply covariances in
physical parameter units; the base class handles its internal scaling.
The complete covariance matrix should be symmetric and positive definite.

For fitters with unknown parameters, the base class also provides
`GetBestUnknowns()` and `GetUnmeasuredCovariance()`. Use them only when the
concrete fitter defines unknowns. A fresh fitter object per event keeps the
inputs and fit history separate.

## Adding unsupported kinematics

If the existing fitters do not describe your reaction or constraints, implement
a new concrete class inherited from `KinematicFitter`. Place its declaration
in `include/` and its implementation in `src/`. Use an existing fitter as a
template and document the order and units of every parameter.

The concrete class should provide the following:

| Part | Responsibility |
| --- | --- |
| **Class definition and constructor** | Declare the input particles, fixed quantities, output observables, and getters. Assign the measured inputs in the constructor and call `Initialize()`. This is the class declaration/constructor, not a required method named `ClassDefinition()`. |
| **`Initialize()`** | Clear the previous fit state, set `nMeas`, `nUnkn`, and `nConst`, and build the initial parameter columns. Store the initial measurements in `Measurements`; initialize `Unknowns` when `nUnkn > 0`. Initialize the required pull and chi-square history entries and mark the fitter initialized. |
| **`SetConstraints()`** | Define the constraints and their derivatives for the current step. Store the constraint column in `FMats`, derivatives with respect to measured parameters in `dFdMs`, and derivatives with respect to unknown parameters in `dFdUs` when unknowns are present. The reaction-specific mathematics belongs here. |
| **`SampleStepPoint(int steps)`** | Convert the parameters at the requested step into particles or other observables. Read `Measurements.at(steps)` and, when applicable, `Unknowns.at(steps)`, then update the results exposed by the concrete class. |

Keep the parameter order consistent across initialization, covariance assignment,
constraint derivatives, observable reconstruction, and tests. For a fitter
without unknowns, set `nUnkn = 0` and omit unknown-parameter handling.

## Checking fit performance

Use representative signal events, a calibration sample, or simulated events
to examine:

- **Chi-square:** returned by `DoKinematicFit()`. Compare its distribution with
  the expected number of degrees of freedom, available through `GetNDF()`.
- **P-value:** obtained with `GetPValue()`. For a well-modeled signal sample
  with appropriate Gaussian uncertainties, its distribution should be
  approximately uniform before fit-quality selections.
- **Pulls:** obtained with `GetPull()`. Their distributions should be centered
  near zero and have a Gaussian width near one when the error model is appropriate.
- **Constraint residuals:** obtained with `GetKinematicConstraints()`. Check
  them separately to verify that the selected fit result satisfies the constraints.

The measured pulls describe the measurement corrections normalized by the
uncertainties of those corrections. They are not simply raw measurement
residuals divided by the input resolutions. When generated truth is available,
also compare initial and fitted parameters with that truth. For a resonance
with a finite width, define an invariant-mass residual relative to the mass
generated for that event.

## Tuning the covariance for practical use

Start with an initial estimate of each measurement resolution, then calibrate
the covariance using the resulting pull distributions:

1. Choose a representative calibration sample and initial resolutions.
2. Assign their squared values through `SetVariance()`, add any correlations,
   and fit the events.
3. Fit the central part of each pull distribution with a Gaussian. Record its
   mean and width, `w`.
4. Rescale the corresponding resolution by that width:
   **`sigma_new = w * sigma_old`**. Consequently, its diagonal covariance entry
   becomes **`variance_new = w * w * variance_old`**.
5. Rebuild the covariance, refit the same calibration sample, and repeat until
   the pull widths are consistent with one within their statistical uncertainties.

For example, a pull width of `1.2` calls for increasing the resolution estimate
by 20%, which increases the corresponding variance by a factor of `1.44`.
A width below one calls for a smaller resolution estimate.

When retaining the same correlation coefficients, rescale an off-diagonal
entry consistently: `covariance_new(i,j) = w_i * w_j * covariance_old(i,j)`.
Check the complete covariance again after updating it. Because parameters
are coupled by the constraints, tuning generally requires several iterations.

Keep the calibration events and their underlying measurement errors or simulated
smearing fixed while tuning the covariance supplied to the fitter. Check pull
means and tails as well as widths: a displaced mean can indicate a measurement
bias or an incorrect constraint model. Validate the tuned covariance on an
independent sample and recheck chi-square and p-value distributions.

## Running the examples

ROOT 6 is required. From the project directory, load the test macro:

```cpp
.L TestKinfit.cc+
TestKinfit(); // Print the available examples and apply the example style.
```

Examples include:

```cpp
TestFourVectorFit(10000);
TestMissingEtaFit(10000);
TestPbarElasticFit(10000);
TestPbarPhiPhiFit(10000);
// TestBalanceFit(10000); // Run both antiproton balance tests.
```

The antiproton examples default to a 1 GeV/c beam. The phi-phi example samples
finite-width phi masses and plots invariant-mass residuals relative to event
truth. `TestKinfit.hh` allows variable names, units, histogram ranges, and
smearing resolutions to be configured for each test. Output plots and ROOT
histograms are saved below `figs/`.
