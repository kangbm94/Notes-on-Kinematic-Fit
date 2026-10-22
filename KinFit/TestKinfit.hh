#ifndef TEST_KINFIT_HH
#define TEST_KINFIT_HH

// Generic ROOT-macro utilities for kinematic-fit performance tests.
// All matrix getters must return UN-SCALED physical parameters as column vectors.
// Variable order must match the rows of the corresponding fitter getter.
// This header has no dependency on a particular fitter or reaction.
//
// Typical loop:
//   KFTestDiagnostics test("MyFit", measuredVariables, unknownVariables);
//   // ... generate trueMeasured and trueUnknown ...
//   auto sample = test.Smear(trueMeasured);
//   // ... construct a fitter from sample.measurements ...
//   sample.SetVariance(fitter);
//   double chi2 = fitter.DoKinematicFit();
//   test.FillHist(fitter, chi2, trueMeasured, trueUnknown);
//   // After the loop: test.DrawHist("figs/MyFit");

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

struct KFTestRange {
  int bins;
  double min, max;
  KFTestRange(int n = 100, double lo = -1., double hi = 1.)
      : bins(n), min(lo), max(hi) {}
};

typedef std::function<double(double)> KFTestResolution;

inline KFTestResolution KFAbsoluteResolution(double sigma) {
  return [sigma](double) { return sigma; };
}

inline KFTestResolution KFRelativeResolution(double fraction) {
  return [fraction](double value) { return fraction * std::abs(value); };
}

struct KFTestVariable {
  TString name;                    // Unique within the measured/unknown list.
  TString axisTitle;               // ROOT axis title, including units.
  KFTestRange residualRange;       // Range of initial/fitted minus true.
  KFTestRange pullRange;
  KFTestResolution resolution;     // True sigma used only for smearing.
  // Estimated sigma for the fit covariance, evaluated at the same true value.
  // Leave empty to use the smearing sigma (backward-compatible default).
  KFTestResolution resolutionEstm;
  double meanScale = 1.;           // Smearing mean = scale * truth + shift.
  double meanShift = 0.;
  double period = 0.;              // E.g. 2*pi for phi; zero = not periodic.
  // Optional non-Gaussian generator: sampler(mean, sigma, random).
  std::function<double(double, double, TRandom&)> sampler;

  KFTestVariable(const TString& key, const TString& title,
                 KFTestRange range,
                 KFTestResolution sigma = KFTestResolution(),
                 double wrapPeriod = 0.)
      : name(key), axisTitle(title), residualRange(range),
        pullRange(100, -5., 5.), resolution(sigma), period(wrapPeriod) {}

  double Difference(double a, double b) const {
    return period > 0. ? std::remainder(a - b, period) : a - b;
  }
};

inline void KFCheckRange(const KFTestRange& r) {
  if (r.bins <= 0 || !std::isfinite(r.min) || !std::isfinite(r.max) ||
      r.min >= r.max)
    throw std::invalid_argument("Invalid histogram binning");
}

inline void KFCheckColumn(const TMatrixD& v, std::size_t size,
                          const char* description) {
  if (size == 0 && v.GetNrows() == 0) return; // ROOT empty matrix: 0 x 0.
  if (v.GetNrows() != static_cast<int>(size) || v.GetNcols() != 1)
    throw std::invalid_argument(TString::Format(
        "%s: expected %zu x 1, got %d x %d", description, size,
        v.GetNrows(), v.GetNcols()).Data());
}

inline TMatrixD KFColumn(const std::vector<double>& values) {
  if (values.empty()) return TMatrixD();
  TMatrixD out(static_cast<int>(values.size()), 1);
  for (std::size_t i = 0; i < values.size(); ++i) out(i, 0) = values[i];
  return out;
}

// Convenience adapters for the spherical (p, theta, phi) fitters in the examples.
// Generic smearing and diagnostics below do not assume groups of three.
inline TMatrixD KFSphericalParameters(
    std::initializer_list<TLorentzVector> particles) {
  std::vector<double> values;
  for (const auto& p : particles) {
    values.push_back(p.Vect().Mag());
    values.push_back(p.Theta());
    values.push_back(p.Phi());
  }
  return KFColumn(values);
}

inline TLorentzVector KFFourVector(const TMatrixD& v, int first, double mass) {
  if (first < 0 || first + 2 >= v.GetNrows() || v.GetNcols() != 1)
    throw std::invalid_argument("KFFourVector: invalid parameter rows");
  const double p = v(first, 0), th = v(first + 1, 0), ph = v(first + 2, 0);
  return TLorentzVector(p * std::sin(th) * std::cos(ph),
                        p * std::sin(th) * std::sin(ph), p * std::cos(th),
                        std::hypot(p, mass));
}

struct KFTestSample {
  TMatrixD measurements;
  std::vector<double> variances; // Estimated sigmas squared, for SetVariance().

  template<class Fitter> void SetVariance(Fitter& fitter) const {
    std::vector<double> diagonal = variances; // Existing API accepts double*.
    const TMatrixD initial = fitter.GetInitialMeasurements();
    if (diagonal.empty() || diagonal.size() != static_cast<std::size_t>(initial.GetNrows()))
      throw std::invalid_argument("Variance count does not match fitter measurements");
    for (double v : diagonal)
      if (!std::isfinite(v) || v <= 0.)
        throw std::invalid_argument("Variances must be finite and positive");
    fitter.SetVariance(diagonal.data());
  }
};

inline KFTestSample KFSmear(const TMatrixD& truth,
                            const std::vector<KFTestVariable>& variables,
                            TRandom& random) {
  KFCheckColumn(truth, variables.size(), "Smearing truth");
  KFTestSample out{TMatrixD(truth), {}};
  for (std::size_t i = 0; i < variables.size(); ++i) {
    const auto& v = variables[i];
    if (!v.resolution)
      throw std::invalid_argument(("Missing resolution for " + v.name).Data());
    const double value = truth(i, 0);
    const double sigma = v.resolution(value);
    const double sigmaEstm = v.resolutionEstm ? v.resolutionEstm(value) : sigma;
    const double mean = v.meanScale * value + v.meanShift;
    if (!std::isfinite(value) || !std::isfinite(mean) ||
        !std::isfinite(sigma) || sigma <= 0.)
      throw std::invalid_argument(("Invalid smearing parameters for " + v.name).Data());
    if (!std::isfinite(sigmaEstm) || sigmaEstm <= 0.)
      throw std::invalid_argument(("Invalid covariance resolution for " + v.name).Data());
    double smeared = v.sampler ? v.sampler(mean, sigma, random)
                               : random.Gaus(mean, sigma);
    if (!std::isfinite(smeared))
      throw std::runtime_error(("Nonfinite smeared value for " + v.name).Data());
    if (v.period > 0.) smeared = std::remainder(smeared, v.period);
    out.measurements(i, 0) = smeared;
    out.variances.push_back(sigmaEstm * sigmaEstm);
  }
  return out;
}

// For sequential fits or externally generated correlated measurements.
// Covariance is in physical units, not the fitter's scaled coordinates.
template<class Fitter>
inline void KFSetCovariance(Fitter& fitter, const TMatrixD& covariance) {
  const int n = covariance.GetNrows();
  if (n == 0 || covariance.GetNcols() != n ||
      fitter.GetInitialMeasurements().GetNrows() != n)
    throw std::invalid_argument("Covariance dimensions do not match fitter measurements");
  std::vector<double> diagonal(n);
  for (int i = 0; i < n; ++i) {
    diagonal[i] = covariance(i, i);
    if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.)
      throw std::invalid_argument("Invalid covariance diagonal");
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      if (!std::isfinite(covariance(i, j)) ||
          std::abs(covariance(i, j) - covariance(j, i)) >
              1.e-10 * std::sqrt(diagonal[i] * diagonal[j]))
        throw std::invalid_argument("Covariance must be finite and symmetric");
  }
  fitter.SetVariance(diagonal.data());
  fitter.AddOffdiagonals(covariance);
}

class KFTestDiagnostics {
public:
  struct VariableHistograms {
    std::unique_ptr<TH1D> initialResidual, fittedResidual, pull;
    unsigned long skippedPulls = 0;
  };
  struct HistogramPair {
    std::unique_ptr<TH1D> initial, fitted;
  };
  // Measured-variable indices in ROOT pad order (left to right, top to bottom).
  struct CanvasPage {
    int columns, rows;
    std::vector<std::size_t> variables;
  };

  KFTestDiagnostics(const TString& name,
                    std::vector<KFTestVariable> measured,
                    std::vector<KFTestVariable> unknown = {},
                    std::vector<KFTestVariable> constraints = {},
                    KFTestRange chi2Range = KFTestRange(200, 0., 50.))
      : measured_(measured), unknown_(unknown), constraintSpecs_(constraints) {
    static unsigned long serial = 0;
    prefix_ = name + TString::Format("_%lu", ++serial);
    CheckVariables(measured_);
    CheckVariables(unknown_);
    CheckVariables(constraintSpecs_);
    MakeHist(chi2Range);
  }

  const std::vector<KFTestVariable>& MeasuredVariables() const { return measured_; }
  const std::vector<KFTestVariable>& UnknownVariables() const { return unknown_; }

  KFTestSample Smear(const TMatrixD& truth, TRandom& random) const {
    return KFSmear(truth, measured_, random);
  }
  KFTestSample Smear(const TMatrixD& truth) const { return Smear(truth, *gRandom); }

  // Measured pull: fitter.GetPull(), normally (initial-best)/sigma_correction.
  // Unknown pull: (best-true)/sqrt(fitted unknown covariance diagonal).
  // GetUPull() is deliberately not used: initial unknowns are starting estimates,
  // not independent measurements with known errors.
  // This records the fitter's selected result; it imposes no chi2/p-value or
  // convergence cut. "Finite result" does not imply satisfied constraints.
  template<class Fitter>
  bool FillHist(Fitter& fitter, double chi2, const TMatrixD& trueMeasured,
                const TMatrixD& trueUnknown = TMatrixD(), double weight = 1.) {
    const TMatrixD initial = fitter.GetInitialMeasurements();
    const TMatrixD best = fitter.GetBestMeasurements();
    // A measured-only fit (e.g. BalanceFitter) must not query unknowns.
    const TMatrixD initialU = unknown_.empty() ? TMatrixD() : fitter.GetInitialUnknowns();
    const TMatrixD bestU = unknown_.empty() ? TMatrixD() : fitter.GetBestUnknowns();
    KFCheckColumn(trueMeasured, measured_.size(), "Measured truth");
    KFCheckColumn(initial, measured_.size(), "Initial measurements");
    KFCheckColumn(best, measured_.size(), "Best measurements");
    KFCheckColumn(trueUnknown, unknown_.size(), "Unknown truth");
    KFCheckColumn(initialU, unknown_.size(), "Initial unknowns");
    KFCheckColumn(bestU, unknown_.size(), "Best unknowns");
    if (!std::isfinite(weight)) throw std::invalid_argument("Nonfinite event weight");

    if (!std::isfinite(chi2) || chi2 < 0. || !Finite(initial) || !Finite(best) ||
        !Finite(initialU) || !Finite(bestU) ||
        !Finite(trueMeasured) || !Finite(trueUnknown)) {
      status_->Fill(1., weight);
      return false;
    }
    const auto pulls = fitter.GetPull();
    if (pulls.size() != measured_.size())
      throw std::invalid_argument("GetPull size does not match measured configuration");
    const TMatrixD unknownCovariance = unknown_.empty() ? TMatrixD()
                                                       : fitter.GetUnmeasuredCovariance();
    if (!unknown_.empty()) {
      if (unknownCovariance.GetNrows() != static_cast<int>(unknown_.size()) ||
          unknownCovariance.GetNcols() != static_cast<int>(unknown_.size()))
        throw std::invalid_argument("Unknown covariance size does not match configuration");
    }
    const auto initialConstraints = fitter.GetInitialConstraints();
    const auto bestConstraints = fitter.GetKinematicConstraints();
    if (initialConstraints.size() != bestConstraints.size())
      throw std::invalid_argument("Initial/best constraint sizes differ");
    if (!constraintsBooked_) BookConstraints(initialConstraints.size());
    if (constraints_.size() != initialConstraints.size())
      throw std::invalid_argument("Constraint count changed or configuration is incorrect");

    status_->Fill(0., weight);
    chi2_->Fill(chi2, weight);
    const double pval = fitter.GetPValue();
    // ROOT treats exactly 1 as overflow on a [0,1) axis.
    if (std::isfinite(pval) && pval >= 0. && pval <= 1.)
      pvalue_->Fill(pval == 1. ? std::nextafter(1., 0.) : pval, weight);
    steps_->Fill(fitter.GetNStep(), weight);
    bestStep_->Fill(fitter.GetBestStep(), weight);

    for (std::size_t i = 0; i < measured_.size(); ++i) {
      FillResidual(measuredHist_[i], measured_[i], initial(i, 0), best(i, 0),
                   trueMeasured(i, 0), weight);
      FillPull(measuredHist_[i], pulls[i], weight);
    }
    for (std::size_t i = 0; i < unknown_.size(); ++i) {
      FillResidual(unknownHist_[i], unknown_[i], initialU(i, 0), bestU(i, 0),
                   trueUnknown(i, 0), weight);
      const double variance = unknownCovariance(i, i);
      if (std::isfinite(variance) && variance > 0.)
        FillPull(unknownHist_[i], unknown_[i].Difference(bestU(i, 0), trueUnknown(i, 0)) /
                                      std::sqrt(variance), weight);
      else
        ++unknownHist_[i].skippedPulls;
    }
    for (std::size_t i = 0; i < constraints_.size(); ++i)
      FillPair(constraints_[i], initialConstraints[i], bestConstraints[i], weight);
    return true;
  }

  // Extra observables (invariant mass, missing mass, derived residuals, ...).
  // Their range is the value range passed to FillObservable, not an automatic
  // truth subtraction. Add them before the event loop.
  void AddObservable(const KFTestVariable& variable) {
    for (const auto& v : observableSpecs_)
      if (v.name == variable.name) throw std::invalid_argument("Duplicate observable");
    CheckVariables({variable});
    observableSpecs_.push_back(variable);
    observables_.push_back(MakePair("Observable_" + variable.name,
                                     variable.axisTitle, variable.residualRange));
  }

  void FillObservable(const TString& name, double initial, double fitted,
                      double weight = 1.) {
    FillPair(observables_.at(ObservableIndex(name)), initial, fitted, weight);
  }

  TH1D* Chi2Histogram() const { return chi2_.get(); }
  TH1D* PValueHistogram() const { return pvalue_.get(); }
  TH1D* StatusHistogram() const { return status_.get(); }
  TH1D* ObservableHistogram(const TString& name, bool fitted = false) const {
    const auto& pair = observables_.at(ObservableIndex(name));
    return fitted ? pair.fitted.get() : pair.initial.get();
  }
  const VariableHistograms& MeasuredHistograms(std::size_t i) const {
    return measuredHist_.at(i);
  }
  const VariableHistograms& UnknownHistograms(std::size_t i) const {
    return unknownHist_.at(i);
  }

  // refit=false reuses the stored Gaussian for additional monitors.
  static void DrawCorrectionPull(TH1D* histogram, bool refit = true) {
    if (!histogram) throw std::invalid_argument("Missing correction-pull histogram");
    DrawPull(histogram, true, refit);
  }

  // Optionally prepend derived-observable overlays to the first unknown
  // residual canvas. Unknown parameter lists and pull canvases are unchanged.
  // Explicit measured pages override the default nine-variable pagination.
  void DrawHist(const TString& directory = "",
                const std::vector<TString>& unknownResidualTopRow = {},
                const std::vector<CanvasPage>& measuredPages = {}) const {
    if (!directory.IsNull()) gSystem->mkdir(directory, true);
    if (measuredPages.empty()) DrawVariables(measuredHist_, "Measured", directory);
    else DrawMeasuredPages(measuredPages, directory);
    DrawVariables(unknownHist_, "Unknown", directory, unknownResidualTopRow);
    DrawPairs(constraints_, "Constraints", directory);
    DrawPairs(observables_, "Observables", directory);
    auto* canvas = new TCanvas(prefix_ + "_Summary", prefix_ + " summary", 1200, 800);
    canvas->Divide(3, 2);
    const TH1D* h[] = {chi2_.get(), pvalue_.get(), steps_.get(), bestStep_.get(), status_.get()};
    for (int i = 0; i < 5; ++i) { canvas->cd(i + 1); h[i]->DrawCopy("hist"); }
    SaveCanvas(canvas, directory);
    PrintSummary();
  }

  void WriteHist(const TString& filename) const {
    TFile output(filename, "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create histogram ROOT file");
    chi2_->Write(); pvalue_->Write(); steps_->Write(); bestStep_->Write(); status_->Write();
    for (const auto* group : {&measuredHist_, &unknownHist_})
      for (const auto& h : *group) {
        h.initialResidual->Write(); h.fittedResidual->Write(); h.pull->Write();
      }
    for (const auto* group : {&constraints_, &observables_})
      for (const auto& h : *group) { h.initial->Write(); h.fitted->Write(); }
  }

  void PrintSummary() const {
    std::cout << prefix_ << ": " << status_->GetBinContent(1) << " finite results, "
              << status_->GetBinContent(2) << " invalid results" << std::endl;
    for (const auto* group : {&measuredHist_, &unknownHist_})
      for (const auto& h : *group)
        if (h.skippedPulls)
          std::cout << "  " << h.pull->GetName() << ": skipped " << h.skippedPulls
                    << " undefined pulls" << std::endl;
  }

private:
  TString prefix_;
  std::vector<KFTestVariable> measured_, unknown_, constraintSpecs_, observableSpecs_;
  std::vector<VariableHistograms> measuredHist_, unknownHist_;
  std::vector<HistogramPair> constraints_, observables_;
  std::unique_ptr<TH1D> chi2_, pvalue_, steps_, bestStep_, status_;
  bool constraintsBooked_ = false;

  static void CheckVariables(const std::vector<KFTestVariable>& variables) {
    std::set<TString> names;
    for (const auto& v : variables) {
      if (v.name.IsNull() || !names.insert(v.name).second)
        throw std::invalid_argument("Empty or duplicate variable name");
      KFCheckRange(v.residualRange);
      KFCheckRange(v.pullRange);
      if (!std::isfinite(v.period) || v.period < 0.)
        throw std::invalid_argument("Invalid variable period");
    }
  }
  static bool Finite(const TMatrixD& v) {
    for (int i = 0; i < v.GetNrows(); ++i)
      for (int j = 0; j < v.GetNcols(); ++j)
        if (!std::isfinite(v(i, j))) return false;
    return true;
  }
  std::unique_ptr<TH1D> MakeOne(const TString& key, const TString& title,
                               KFTestRange range, int color = kBlack) const {
    KFCheckRange(range);
    std::unique_ptr<TH1D> h(new TH1D(prefix_ + "_" + key,
                                     ";" + title + ";Events", range.bins, range.min, range.max));
    h->SetDirectory(nullptr);
    h->Sumw2();
    h->SetLineColor(color);
    return h;
  }
  HistogramPair MakePair(const TString& key, const TString& title, KFTestRange range) const {
    HistogramPair h;
    h.initial = MakeOne(key + "_Initial", title, range, kRed);
    h.fitted = MakeOne(key + "_Fitted", title, range, kBlue);
    return h;
  }
  void MakeHist(KFTestRange chi2Range) {
    for (int group = 0; group < 2; ++group) {
      const auto& vars = group == 0 ? measured_ : unknown_;
      auto& hist = group == 0 ? measuredHist_ : unknownHist_;
      const TString label = group == 0 ? "Measured_" : "Unknown_";
      for (const auto& v : vars) {
        VariableHistograms h;
        h.initialResidual = MakeOne(label + v.name + "_InitialResidual",
                                     "#Delta " + v.axisTitle, v.residualRange, kRed);
        h.fittedResidual = MakeOne(label + v.name + "_FittedResidual",
                                    "#Delta " + v.axisTitle, v.residualRange, kBlue);
        h.pull = MakeOne(label + v.name + "_Pull",
                         (group == 0 ? "Correction pull: " : "Truth pull: ") + v.axisTitle,
                         v.pullRange);
        hist.push_back(std::move(h));
      }
    }
    chi2_ = MakeOne("Chi2", "#chi^{2}", chi2Range);
    pvalue_ = MakeOne("PValue", "P-value", KFTestRange(1000, 0., 1.));
    steps_ = MakeOne("Steps", "Number of iterations", KFTestRange(301, -.5, 300.5));
    bestStep_ = MakeOne("BestStep", "Selected iteration", KFTestRange(301, -.5, 300.5));
    status_ = MakeOne("Status", "Result status", KFTestRange(2, -.5, 1.5));
    status_->GetXaxis()->SetBinLabel(1, "Finite result");
    status_->GetXaxis()->SetBinLabel(2, "Invalid result");
    if (!constraintSpecs_.empty()) BookConstraints(constraintSpecs_.size());
  }
  void BookConstraints(std::size_t count) {
    if (constraintSpecs_.empty())
      for (std::size_t i = 0; i < count; ++i)
        constraintSpecs_.push_back(KFTestVariable(
            TString::Format("F%zu", i + 1), TString::Format("F_{%zu}", i + 1),
            KFTestRange(200, -1., 1.)));
    if (constraintSpecs_.size() != count)
      throw std::invalid_argument("Constraint configuration size mismatch");
    for (const auto& v : constraintSpecs_)
      constraints_.push_back(MakePair("Constraint_" + v.name, v.axisTitle, v.residualRange));
    constraintsBooked_ = true;
  }
  std::size_t ObservableIndex(const TString& name) const {
    for (std::size_t i = 0; i < observableSpecs_.size(); ++i)
      if (observableSpecs_[i].name == name) return i;
    throw std::invalid_argument(("Unknown observable: " + name).Data());
  }
  static void FillPair(HistogramPair& h, double initial, double best, double weight) {
    if (!std::isfinite(weight)) throw std::invalid_argument("Nonfinite event weight");
    if (std::isfinite(initial)) h.initial->Fill(initial, weight);
    if (std::isfinite(best)) h.fitted->Fill(best, weight);
  }
  static void FillResidual(VariableHistograms& h, const KFTestVariable& v,
                            double initial, double best, double truth, double weight) {
    h.initialResidual->Fill(v.Difference(initial, truth), weight);
    h.fittedResidual->Fill(v.Difference(best, truth), weight);
  }
  static void FillPull(VariableHistograms& h, double pull, double weight) {
    if (std::isfinite(pull)) h.pull->Fill(pull, weight);
    else ++h.skippedPulls;
  }
  static void DrawOverlay(const TH1D* initial, const TH1D* fitted) {
    TH1* first = initial->DrawCopy("hist");
    first->SetMaximum(1.2 * std::max(1., std::max(initial->GetMaximum(), fitted->GetMaximum())));
    TH1* second = fitted->DrawCopy("hist same");
    auto* legend = new TLegend(.60, .74, .90, .90);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(first, "Initial", "l");
    legend->AddEntry(second, "Fitted", "l");
    legend->Draw();
  }
  static void DrawPull(TH1D* histogram, bool fitGaussian, bool refit = true) {
    if (!fitGaussian) {
      histogram->DrawCopy("hist");
      return;
    }
    // Reserve a band above the axes for the fit result. Transparent text alone
    // would still overlap a wide distribution if placed inside the plot.
    if (gPad) gPad->SetTopMargin(std::max(.26, double(gPad->GetTopMargin())));
    // Fit the complete displayed pull range. Use the stored bin errors, also
    // for weighted samples; do not replace a failed fit by histogram moments.
    const TString fitName = TString(histogram->GetName()) + "_Gaussian";
    if (refit)
      if (auto* previous = histogram->GetFunction(fitName)) {
        histogram->GetListOfFunctions()->Remove(previous);
        delete previous;
      }
    const double lo = histogram->GetXaxis()->GetXmin();
    const double hi = histogram->GetXaxis()->GetXmax();
    const double rms = histogram->GetRMS();
    const bool enoughData = histogram->GetEffectiveEntries() >= 20. &&
                            std::isfinite(rms) && rms > 0.;
    TF1 gaussian(fitName, "gaus", lo, hi);
    bool valid = false;
    if (!refit) {
      if (const auto* stored = histogram->GetFunction(fitName)) {
        for (int parameter = 0; parameter < 3; ++parameter) {
          gaussian.SetParameter(parameter, stored->GetParameter(parameter));
          gaussian.SetParError(parameter, stored->GetParError(parameter));
        }
        valid = true;
      }
    } else if (enoughData) {
      gaussian.SetParameters(histogram->GetMaximum(), histogram->GetMean(), rms);
      gaussian.SetParLimits(2, 1.e-6 * (hi - lo), 10. * (hi - lo));
      // ROOT owns the stored function; WriteHist() retains its fit parameters.
      const TFitResultPtr result = histogram->Fit(&gaussian, "QRS0");
      valid = result.Get() && int(result) == 0 && result->IsValid() &&
              std::isfinite(gaussian.GetParameter(1)) &&
              std::isfinite(gaussian.GetParameter(2)) && gaussian.GetParameter(2) > 0. &&
              std::isfinite(gaussian.GetParError(1)) &&
              std::isfinite(gaussian.GetParError(2));
    }
    if (auto* stored = histogram->GetFunction(fitName)) {
      if (valid) {
        stored->SetLineColor(kRed + 1);
        stored->SetLineWidth(2);
        stored->ResetBit(TF1::kNotDraw);
      } else {
        histogram->GetListOfFunctions()->Remove(stored);
        delete stored;
      }
    }
    TH1* displayed = histogram->DrawCopy("hist");
    if (valid) displayed->GetFunction(fitName)->Draw("same");

    auto* text = new TPaveText(.17, .77, .94, .98, "NDC");
    text->SetBorderSize(0); text->SetFillStyle(0);
    text->SetTextFont(42); text->SetTextAlign(12); text->SetTextSize(.035);
    if (valid) {
      text->AddText("Gaussian fit");
      text->AddText(TString::Format("Mean = %.3f #pm %.3f",
                                    gaussian.GetParameter(1), gaussian.GetParError(1)));
      text->AddText(TString::Format("Width = %.3f #pm %.3f",
                                    gaussian.GetParameter(2), gaussian.GetParError(2)));
    } else {
      text->AddText(enoughData ? "Gaussian fit failed" : "Too few entries for Gaussian fit");
    }
    text->Draw();
  }
  void DrawMeasuredPages(const std::vector<CanvasPage>& pages,
                         const TString& directory) const {
    // Validate all pages before allocating canvases.
    for (const auto& page : pages) {
      if (page.columns <= 0 || page.rows <= 0 || page.variables.empty() ||
          page.variables.size() > static_cast<std::size_t>(page.columns) * page.rows)
        throw std::invalid_argument("Invalid measured canvas layout");
      for (auto index : page.variables)
        if (index >= measuredHist_.size())
          throw std::out_of_range("Measured canvas variable index");
    }
    for (std::size_t pageIndex = 0; pageIndex < pages.size(); ++pageIndex) {
      const auto& page = pages[pageIndex];
      const TString suffix = TString::Format("_%zu", pageIndex + 1);
      auto* residual = new TCanvas(prefix_ + "_MeasuredResidual" + suffix,
                                   "Measured residuals", 400 * page.columns, 300 * page.rows);
      auto* pull = new TCanvas(prefix_ + "_MeasuredPull" + suffix,
                               "Measured pulls", 400 * page.columns, 300 * page.rows);
      residual->Divide(page.columns, page.rows);
      pull->Divide(page.columns, page.rows);
      for (std::size_t pad = 0; pad < page.variables.size(); ++pad) {
        const auto& hist = measuredHist_[page.variables[pad]];
        residual->cd(pad + 1);
        DrawOverlay(hist.initialResidual.get(), hist.fittedResidual.get());
        pull->cd(pad + 1); DrawPull(hist.pull.get(), true);
      }
      SaveCanvas(residual, directory); SaveCanvas(pull, directory);
    }
  }
  void DrawVariables(const std::vector<VariableHistograms>& hist,
                      const TString& label, const TString& directory,
                      const std::vector<TString>& topRow = {}) const {
    for (std::size_t begin = 0; begin < hist.size(); begin += 9) {
      const int count = std::min<std::size_t>(9, hist.size() - begin);
      const TString page = TString::Format("_%zu", begin / 9 + 1);
      auto* residual = new TCanvas(prefix_ + "_" + label + "Residual" + page,
                                   label + " residuals", 1200, 800);
      auto* pull = new TCanvas(prefix_ + "_" + label + "Pull" + page,
                               label + " pulls", 1200, 800);
      const int cols = std::min(3, count), rows = (count + cols - 1) / cols;
      const int extra = begin == 0 ? static_cast<int>(topRow.size()) : 0;
      const int residualCols = extra ? 3 : cols;
      const int extraRows = (extra + residualCols - 1) / residualCols;
      const int offset = extraRows * residualCols;
      residual->Divide(residualCols, extraRows + (count + residualCols - 1) / residualCols);
      pull->Divide(cols, rows);
      for (int i = 0; i < extra; ++i) {
        residual->cd(i + 1);
        const auto& pair = observables_.at(ObservableIndex(topRow[i]));
        DrawOverlay(pair.initial.get(), pair.fitted.get());
      }
      for (int i = 0; i < count; ++i) {
        residual->cd(offset + i + 1);
        DrawOverlay(hist[begin + i].initialResidual.get(), hist[begin + i].fittedResidual.get());
        pull->cd(i + 1); DrawPull(hist[begin + i].pull.get(), true);
      }
      SaveCanvas(residual, directory); SaveCanvas(pull, directory);
    }
  }
  void DrawPairs(const std::vector<HistogramPair>& hist, const TString& label,
                   const TString& directory) const {
    for (std::size_t begin = 0; begin < hist.size(); begin += 9) {
      const int count = std::min<std::size_t>(9, hist.size() - begin);
      auto* canvas = new TCanvas(prefix_ + "_" + label + TString::Format("_%zu", begin / 9 + 1),
                                 label, 1200, 800);
      const int cols = std::min(3, count);
      canvas->Divide(cols, (count + cols - 1) / cols);
      for (int i = 0; i < count; ++i) {
        canvas->cd(i + 1);
        DrawOverlay(hist[begin + i].initial.get(), hist[begin + i].fitted.get());
      }
      SaveCanvas(canvas, directory);
    }
  }
  static void SaveCanvas(TCanvas* canvas, const TString& directory) {
    if (!directory.IsNull())
      canvas->SaveAs(directory + "/" + canvas->GetName() + ".pdf");
  }
};

// These runtime helpers contain functions and unique histogram owners; ROOT
// I/O collection dictionaries are unnecessary and cannot copy these owners.
#ifdef __ROOTCLING__
#pragma link off class KFTestDiagnostics;
#pragma link off class KFTestDiagnostics::VariableHistograms;
#pragma link off class KFTestDiagnostics::HistogramPair;
#pragma link off class std::vector<KFTestVariable>;
#pragma link off class std::vector<KFTestDiagnostics::VariableHistograms>;
#pragma link off class std::vector<KFTestDiagnostics::HistogramPair>;
#endif

#endif
