#include "src/FourVectorFitter.cc"
#include "src/MissingMassFitter.cc"
#include "src/CascadeFitter.cc"
#include "src/BalanceFitter.cc"
#include "TestKinfit.hh"

#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <limits>
#include <utility>

// Load with: .L TestKinfit.cc+
// Existing fitter sources are used unchanged. Besides their existing API, they
// must provide GetBestMeasurements(), GetInitialMeasurements(), GetBestUnknowns()
// and GetInitialUnknowns(), all returning physical, unscaled TMatrixD columns.
// Histogram configuration, smearing, filling and drawing are reusable helpers.
// Every Test... function shows generation, fitter input, SetVariance(),
// DoKinematicFit(), and fitted-parameter extraction directly in its event loop.
// All generated decays are two-body, so their Generate() weights are constant.

int nev = 100000;
double mXi = 1.321, mXi0 = 1.315, mXi1530 = 1.535;
double mL = 1.115, mP = .938, mK = .493, mPi = .139, mPi0 = .135, mEta = .5478;
double pK18 = 1.8;
double ResP = .05, ResPi1 = .05, ResPi2 = .05;
double ResThP = .005, ResTh = .005, ResThPi2 = .005;
double ResPhP = .01, ResPh = .02, ResPhPi2 = .02;
double ScalePP = 1., ScalePPi1 = 1., ScalePPi2 = 1.;

// Balance examples: a fixed proton target and a measured antiproton beam.
double pPbar = 1.;
double mPhi = 1.019461;
double WidthPhi = .004249; // Full Breit-Wigner width in GeV (PDG: 4.249 MeV).
double ResPbarBeam = .001, ResThPbarBeam = .001, ResPhPbarBeam = .005;

// A particle helper is convenient for these examples. The generic tester also
// accepts individual Cartesian components, energies, masses, etc.; no grouping
// or fixed number of rows is required by KFTestDiagnostics.
std::vector<KFTestVariable> KFParticleVariables(
    TString name, TString title, double momentumResidualRange,
    double relativeMomentumSigma = 0., double thetaSigma = 0., double phiSigma = 0.,
    double momentumScale = 1.) {
  const double pi = std::acos(-1.);
  std::vector<KFTestVariable> vars = {
      {name + "_p", "p_{" + title + "} [GeV/c]",
       {100, -momentumResidualRange, momentumResidualRange},
       relativeMomentumSigma > 0. ? KFRelativeResolution(relativeMomentumSigma) : KFTestResolution()},
      {name + "_theta", "#theta_{" + title + "} [rad]", {100, -.1, .1},
       thetaSigma > 0. ? KFAbsoluteResolution(thetaSigma) : KFTestResolution()},
      {name + "_phi", "#phi_{" + title + "} [rad]", {100, -.3, .3},
       phiSigma > 0. ? KFAbsoluteResolution(phiSigma) : KFTestResolution(), 2. * pi}};
  vars[0].meanScale = momentumScale;
  return vars;
}

std::vector<KFTestVariable> KFJoinVariables(
    std::initializer_list<std::vector<KFTestVariable>> groups) {
  std::vector<KFTestVariable> result;
  for (const auto& group : groups) result.insert(result.end(), group.begin(), group.end());
  return result;
}

std::vector<KFTestVariable> KFCascadeMeasurements() {
  // Must match CascadeFitter's measured rows: p, pi1, pi2; (p,theta,phi) each.
  return KFJoinVariables({
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP),
      KFParticleVariables("pi1", "#pi_{1}", .1, ResPi1, ResTh, ResPh, ScalePPi1),
      KFParticleVariables("pi2", "#pi_{2}", .1, ResPi2, ResThPi2, ResPhPi2, ScalePPi2)});
}

std::vector<KFTestVariable> KFCascadeUnknowns() {
  // CascadeFitter has only the Xi (p, theta, phi) as unknowns.
  // Lambda is reconstructed from the proton and pi1 measurements.
  return KFParticleVariables("xi", "#Xi", .5);
}

void SetStyle() {
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gStyle->SetTitleSize(.05, "XY"); gStyle->SetLabelSize(.045, "XY");
  gStyle->SetPadLeftMargin(.15); gStyle->SetPadBottomMargin(.15);
  gROOT->ForceStyle();
}

// Lambda is derived from measured daughters, never added to the unknowns.
std::vector<KFTestVariable> KFCascadeLambdaResiduals() {
  auto variables = KFParticleVariables("lambda", "#Lambda", .5);
  for (auto& variable : variables) {
    variable.name += "Residual";
    variable.axisTitle = "#Delta " + variable.axisTitle;
  }
  return variables;
}

void KFAddCascadeObservables(KFTestDiagnostics& test) {
  test.AddObservable({"LambdaMass", "M(p#pi_{1}) [GeV/c^{2}]", {500, 1.04, 1.2}});
  test.AddObservable({"XiMass", "M(p#pi_{1}#pi_{2}) [GeV/c^{2}]", {500, mXi - .075, mXi + .075}});
  for (const auto& variable : KFCascadeLambdaResiduals()) test.AddObservable(variable);
}

void KFDrawCascadeHist(KFTestDiagnostics& test, const TString& directory) {
  test.DrawHist(directory, {"lambda_pResidual", "lambda_thetaResidual", "lambda_phiResidual"});
}

struct KFCascadeFitResult {
  double initialMass = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double pvalue = std::numeric_limits<double>::quiet_NaN();
  double piXiPullScore = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
};

// pi2 is the charged pion assigned to Xi -> Lambda pi (pi_Xi).
// Measured rows are proton[0..2], pi1[3..5], pi2[6..8].
double KFPiXiPullScore(const std::vector<double>& pulls) {
  if (pulls.size() != 9)
    throw std::invalid_argument("Cascade pull score requires nine measured pulls");
  double score = 0.;
  for (std::size_t i = 6; i < 9; ++i) {
    if (!std::isfinite(pulls[i])) return std::numeric_limits<double>::quiet_NaN();
    score += pulls[i] * pulls[i];
  }
  return std::isfinite(score) ? score : std::numeric_limits<double>::quiet_NaN();
}

// Histogram-only helper: generation, fitter setup and parameter extraction
// remain visible in the examples below. No fitter API is called here.
void KFFillCascadeObservables(KFTestDiagnostics& test,
                              const TLorentzVector& trueLambda,
                              const TLorentzVector& initialLambda,
                              const TLorentzVector& fittedLambda,
                              const TLorentzVector& initialXi,
                              const TLorentzVector& fittedXi) {
  const TMatrixD truth = KFSphericalParameters({trueLambda});
  const TMatrixD initial = KFSphericalParameters({initialLambda});
  const TMatrixD fitted = KFSphericalParameters({fittedLambda});
  const auto variables = KFCascadeLambdaResiduals();
  for (std::size_t i = 0; i < variables.size(); ++i) {
    const auto& variable = variables[i];
    test.FillObservable(variable.name, variable.Difference(initial(i, 0), truth(i, 0)),
                        variable.Difference(fitted(i, 0), truth(i, 0)));
  }
  test.FillObservable("LambdaMass", initialLambda.M(), fittedLambda.M());
  test.FillObservable("XiMass", initialXi.M(), fittedXi.M());
}


// Exact event-level cut scans. Efficiencies include all generated events after
// the K+ acceptance selection; invalid KF results fail KF-based selections.
struct KFSeparationScan {
  std::vector<double> cuts, efficiency, purity;
};

KFSeparationScan KFScanSeparation(std::vector<double> signal,
                                  std::vector<double> background,
                                  std::size_t generatedSignal,
                                  std::size_t generatedBackground,
                                  double signalFraction, bool keepAbove,
                                  double maximumCut, int steps = 1000,
                                  bool logarithmicCuts = false) {
  if (!generatedSignal || !generatedBackground ||
      !std::isfinite(signalFraction) || signalFraction <= 0. || signalFraction >= 1. ||
      !std::isfinite(maximumCut) || maximumCut <= 0. || steps <= 0)
    throw std::invalid_argument("Invalid separation scan configuration");
  std::sort(signal.begin(), signal.end());
  std::sort(background.begin(), background.end());
  auto count = [keepAbove](const std::vector<double>& values, double cut) {
    return keepAbove ? values.end() - std::lower_bound(values.begin(), values.end(), cut)
                     : std::upper_bound(values.begin(), values.end(), cut) - values.begin();
  };
  KFSeparationScan scan;
  for (int i = 0; i <= steps; ++i) {
    // A score can have a long tail: log1p spacing resolves small cuts while
    // retaining the largest event. Force the last cut to include the endpoint.
    const double cut = i == steps ? maximumCut :
        (logarithmicCuts ? std::expm1(std::log1p(maximumCut) * i / steps)
                         : maximumCut * i / steps);
    const double es = static_cast<double>(count(signal, cut)) / generatedSignal;
    const double eb = static_cast<double>(count(background, cut)) / generatedBackground;
    const double selected = signalFraction * es + (1. - signalFraction) * eb;
    if (selected <= 0.) continue; // Purity is undefined for an empty selection.
    scan.cuts.push_back(cut);
    scan.efficiency.push_back(es);
    scan.purity.push_back(signalFraction * es / selected);
  }
  return scan;
}

// Xi1530- -> Xi- pi0 is signal; Xi1530- -> Xi0 pi- is the competing channel.
// Both end in p pi- pi- pi0 and are reconstructed as Xi- -> Lambda pi-.
void KFCompareSamples(const std::vector<KFCascadeFitResult>& signal,
                       const std::vector<KFCascadeFitResult>& background,
                       const TString& directory, double signalFraction = 1. / 3.) {
  if (signal.empty() || background.empty() || !std::isfinite(signalFraction) ||
      signalFraction <= 0. || signalFraction >= 1.)
    throw std::invalid_argument("Separation requires both channels and 0 < signalFraction < 1");
  gSystem->mkdir(directory, true);
  TH1D sChi2("Xi1530SignalChi2", ";#chi^{2};Events", 500, 0., 50.);
  TH1D bChi2("Xi1530BackgroundChi2", ";#chi^{2};Events", 500, 0., 50.);
  TH1D sPval("Xi1530SignalPValue", ";P-value;Events", 1000, 0., 1.);
  TH1D bPval("Xi1530BackgroundPValue", ";P-value;Events", 1000, 0., 1.);
  TH1D sMass("Xi1530SignalMass", ";M(p#pi^{-}#pi^{-}) [GeV/c^{2}];Events", 600, mXi - .15, mXi + .15);
  TH1D bMass("Xi1530BackgroundMass", ";M(p#pi^{-}#pi^{-}) [GeV/c^{2}];Events", 600, mXi - .15, mXi + .15);
  const char* scoreAxis = ";S = u_{p,#pi_{#Xi}}^{2} + u_{#theta,#pi_{#Xi}}^{2} + u_{#phi,#pi_{#Xi}}^{2};Events";
  TH1D sScore("Xi1530SignalPiXiPullScore", scoreAxis, 600, 0., 30.);
  TH1D bScore("Xi1530BackgroundPiXiPullScore", scoreAxis, 600, 0., 30.);
  for (auto* h : {&sChi2, &bChi2, &sPval, &bPval, &sMass, &bMass, &sScore, &bScore}) {
    h->SetDirectory(nullptr); h->Sumw2();
  }
  std::vector<double> sp, bp, sm, bm, ss, bs;
  auto fill = [&](const std::vector<KFCascadeFitResult>& events,
                  TH1D& chi2, TH1D& pval, TH1D& mass, TH1D& score,
                  std::vector<double>& probabilities, std::vector<double>& errors,
                  std::vector<double>& scores) {
    for (const auto& event : events) {
      if (std::isfinite(event.initialMass)) {
        mass.Fill(event.initialMass);
        errors.push_back(std::abs(event.initialMass - mXi));
      }
      if (!event.valid) continue;
      if (std::isfinite(event.piXiPullScore) && event.piXiPullScore >= 0.) {
        score.Fill(event.piXiPullScore);
        scores.push_back(event.piXiPullScore);
      }
      if (std::isfinite(event.chi2) && event.chi2 >= 0.) chi2.Fill(event.chi2);
      if (std::isfinite(event.pvalue) && event.pvalue >= 0. && event.pvalue <= 1.) {
        pval.Fill(event.pvalue == 1. ? std::nextafter(1., 0.) : event.pvalue);
        probabilities.push_back(event.pvalue);
      }
    }
  };
  fill(signal, sChi2, sPval, sMass, sScore, sp, sm, ss);
  fill(background, bChi2, bPval, bMass, bScore, bp, bm, bs);
  const auto pscan = KFScanSeparation(sp, bp, signal.size(), background.size(), signalFraction, true, 1.);
  const auto mscan = KFScanSeparation(sm, bm, signal.size(), background.size(), signalFraction, false, .15);
  double maximumScoreCut = 30.;
  for (const auto* scores : {&ss, &bs})
    for (double score : *scores) maximumScoreCut = std::max(maximumScoreCut, score);
  const auto sscan = KFScanSeparation(ss, bs, signal.size(), background.size(),
                                     signalFraction, false, maximumScoreCut, 2000, true);

  // Scale channel overlays to the same assumed mixture used for purity.
  const double total = signal.size() + background.size();
  for (auto* h : {&sChi2, &sPval, &sMass, &sScore}) h->Scale(total * signalFraction / signal.size());
  for (auto* h : {&bChi2, &bPval, &bMass, &bScore}) h->Scale(total * (1. - signalFraction) / background.size());
  auto overlay = [&](const TH1D& s, const TH1D& b, bool logY = false, bool zoomChi2 = false) {
    TH1D sum(s); sum.SetDirectory(nullptr); sum.Add(&b); sum.SetLineColor(kBlack);
    sum.SetMaximum((logY ? 2. : 1.2) * std::max(1., sum.GetMaximum()));
    sum.SetMinimum(logY ? .5 : 0.);
    if (zoomChi2) sum.GetXaxis()->SetRangeUser(0., 10.);
    if (logY) gPad->SetLogy();
    auto* all = sum.DrawCopy("hist");
    auto* sh = s.DrawCopy("hist same"); sh->SetLineColor(kBlue);
    auto* bh = b.DrawCopy("hist same"); bh->SetLineColor(kRed);
    auto* legend = new TLegend(.47, .69, .91, .91);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(all, "All events", "l");
    legend->AddEntry(sh, "#Xi^{*-} #rightarrow #Xi^{-}#pi^{0}", "l");
    legend->AddEntry(bh, "#Xi^{*-} #rightarrow #Xi^{0}#pi^{-}", "l");
    legend->Draw();
  };
  auto* sum = new TCanvas("c_sum", "Xi1530 channel separation", 1200, 600);
  sum->Divide(2, 1);
  sum->cd(1); overlay(sChi2, bChi2, false, true);
  sum->cd(2); overlay(sPval, bPval, true);
  sum->SaveAs(directory + "/c_sum.pdf");
  auto* mass = new TCanvas("c_imxi", "Unfit Xi candidate mass", 800, 600);
  mass->cd(); overlay(sMass, bMass);
  mass->SaveAs(directory + "/c_imxi.pdf");
  auto* scoreDistribution = new TCanvas("c_pullscore", "Xi-pion squared-pull score", 800, 600);
  scoreDistribution->cd(); overlay(sScore, bScore, true);
  scoreDistribution->SaveAs(directory + "/c_pullscore.pdf");

  auto graph = [](const std::vector<double>& x, const std::vector<double>& y,
                  const TString& name, int color, bool positiveX = false) {
    auto* g = new TGraph();
    for (std::size_t i = 0; i < x.size(); ++i) {
      // Zero efficiency has no position on a logarithmic axis.
      if (positiveX && x[i] <= 0.) continue;
      g->SetPoint(g->GetN(), x[i], y[i]);
    }
    g->SetName(name); g->SetLineColor(color); g->SetLineWidth(2);
    return g;
  };
  auto* pEfficiency = graph(pscan.cuts, pscan.efficiency, "PValueEfficiency", kBlue);
  auto* pPurity = graph(pscan.cuts, pscan.purity, "PValuePurity", kRed);
  auto* mEfficiency = graph(mscan.cuts, mscan.efficiency, "MassEfficiency", kBlue);
  auto* mPurity = graph(mscan.cuts, mscan.purity, "MassPurity", kRed);
  auto* sEfficiency = graph(sscan.cuts, sscan.efficiency, "PiXiPullScoreEfficiency", kBlue);
  auto* sPurity = graph(sscan.cuts, sscan.purity, "PiXiPullScorePurity", kRed);
  auto drawCuts = [&](TGraph* efficiency, TGraph* purity, const char* xTitle, double xmax) {
    // DrawFrame also handles a scan with no selected events.
    auto* frame = gPad->DrawFrame(0., 0., xmax, 1.05);
    frame->SetTitle(TString(";") + xTitle + ";Efficiency / purity");
    efficiency->Draw("L same"); purity->Draw("L same");
    auto* legend = new TLegend(.52, .68, .90, .90);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(efficiency, "Signal efficiency", "l");
    legend->AddEntry(purity, TString::Format("Purity (f_{S}=%.3g)", signalFraction), "l");
    legend->Draw();
  };
  auto* pval = new TCanvas("c_pval", "P-value selection", 800, 600);
  pval->cd(); drawCuts(pEfficiency, pPurity, "P-value cut", 1.);
  pval->SaveAs(directory + "/c_pval.pdf");
  auto* mcut = new TCanvas("c_mass", "Mass-window selection", 800, 600);
  mcut->cd(); drawCuts(mEfficiency, mPurity, "Mass half-window [GeV/c^{2}]", .15);
  mcut->SaveAs(directory + "/c_mass.pdf");
  auto* scut = new TCanvas("c_score", "Xi-pion squared-pull selection", 800, 600);
  scut->cd(); drawCuts(sEfficiency, sPurity, "Squared-pull score upper cut", 30.);
  scut->SaveAs(directory + "/c_score.pdf");
  auto* pTradeoff = graph(pscan.efficiency, pscan.purity, "PValueEfficiencyPurity", kBlue, true);
  auto* mTradeoff = graph(mscan.efficiency, mscan.purity, "MassEfficiencyPurity", kRed, true);
  auto* sTradeoff = graph(sscan.efficiency, sscan.purity, "PiXiPullScoreEfficiencyPurity", kGreen + 2, true);
  double minimumEfficiency = 0.9; // Retain the supplied high-efficiency zoom.
  #if 0
  for (const auto* scan : {&pscan, &mscan})
    for (double efficiency : scan->efficiency)
      if (efficiency > 0.) minimumEfficiency = std::min(minimumEfficiency, .5 * efficiency);
  #endif
  auto drawTradeoff = [&]() {
    gPad->SetLogx();
    auto* frame = gPad->DrawFrame(minimumEfficiency, 0., 1., 1.05);
    frame->SetTitle(";Signal efficiency;Purity");
    pTradeoff->Draw("L same"); mTradeoff->Draw("L same"); sTradeoff->Draw("L same");
    auto* legend = new TLegend(.20, .20, .70, .45);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(pTradeoff, "P-value cut", "l");
    legend->AddEntry(mTradeoff, "Mass-window cut", "l");
    legend->AddEntry(sTradeoff, "#pi_{#Xi} squared-pull score cut", "l"); legend->Draw();
  };
  auto* tradeoff = new TCanvas("c_eff_pur", "Selection comparison", 800, 600);
  tradeoff->cd(); drawTradeoff();
  tradeoff->SaveAs(directory + "/c_eff_pur.pdf");

  // Retain the existing views and add the score selection to the comparison.
  auto* comparison = new TCanvas("Xi1530Comparison", "Xi1530 comparison", 1500, 600);
  comparison->Divide(3, 1);
  comparison->cd(1); overlay(sScore, bScore, true);
  comparison->cd(2); overlay(sPval, bPval, true);
  comparison->cd(3); overlay(sMass, bMass);
  comparison->SaveAs(directory + "/comparison.pdf");
  auto* selection = new TCanvas("Xi1530Cuts", "Selection studies", 1400, 1000);
  selection->Divide(2, 2);
  selection->cd(1); drawCuts(pEfficiency, pPurity, "P-value cut", 1.);
  selection->cd(2); drawCuts(mEfficiency, mPurity, "Mass half-window [GeV/c^{2}]", .15);
  selection->cd(3); drawCuts(sEfficiency, sPurity, "Squared-pull score upper cut", 30.);
  selection->cd(4); drawTradeoff();
  selection->SaveAs(directory + "/selection.pdf");
  TFile output(directory + "/separation.root", "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create Xi1530 separation ROOT file");
  for (auto* h : {&sChi2, &bChi2, &sPval, &bPval, &sMass, &bMass, &sScore, &bScore}) h->Write();
  for (auto* g : {pEfficiency, pPurity, mEfficiency, mPurity, sEfficiency, sPurity,
                  pTradeoff, mTradeoff, sTradeoff}) g->Write();
  for (auto* c : {sum, mass, scoreDistribution, pval, mcut, scut, tradeoff, comparison, selection}) c->Write();
  std::cout << "Xi1530 separation: signal = Xi- pi0, competing channel = Xi0 pi-\n"
            << "  Assumed signal fraction before selection: " << signalFraction
            << "; efficiency denominator: all generated events after K+ acceptance.\n"
            << "  P-value rule: P >= cut; mass rule: |M_initial - mXi| <= half-window.\n"
            << "  Score rule: sum of pi_Xi (p,theta,phi) squared correction pulls <= cut.\n"
            << "  Invalid fits/pulls fail the score selection; all events remain in efficiency denominators.\n"
            << "  Score display: 0..30; tradeoff scan includes the full score range.\n"
            << "  Invalid KF results fail the P-value selection; mass selection is independent of KF."
            << std::endl;
}

// K- p -> K+ Xi-, Xi- -> Lambda pi2, Lambda -> p pi1.
void TestXiFit(int nEvents = -1) {
  KFTestDiagnostics test("XiCascade", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFAddCascadeObservables(test);
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector primary = beam + target;
  double primaryMasses[2] = {mXi, mK};
  double xiMasses[2] = {mL, mPi};
  double lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace evVert, evXi, evLambda;
  if (!evVert.SetDecay(primary, 2, primaryMasses))
    throw std::runtime_error("Kinematically forbidden Xi production");
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    // 1. Generate the entire decay chain explicitly.
    evVert.Generate();
    TLorentzVector xi = *evVert.GetDecay(0);
    if (!evXi.SetDecay(xi, 2, xiMasses))
      throw std::runtime_error("Kinematically forbidden Xi decay");
    evXi.Generate();
    TLorentzVector lambda = *evXi.GetDecay(0);
    const TLorentzVector truePi2 = *evXi.GetDecay(1);
    if (!evLambda.SetDecay(lambda, 2, lambdaMasses))
      throw std::runtime_error("Kinematically forbidden Lambda decay");
    evLambda.Generate();
    const TLorentzVector trueP = *evLambda.GetDecay(0);
    const TLorentzVector truePi1 = *evLambda.GetDecay(1);

    // 2. Smear the measured rows: p[0..2], pi1[3..5], pi2[6..8].
    const TMatrixD truth = KFSphericalParameters({trueP, truePi1, truePi2});
    const TMatrixD trueUnknown = KFSphericalParameters({trueP + truePi1 + truePi2});
    const KFTestSample sample = test.Smear(truth);
    const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pi1 = KFFourVector(sample.measurements, 3, mPi);
    const TLorentzVector pi2 = KFFourVector(sample.measurements, 6, mPi);

    // 3. Declare the fitter, set constraints and covariance, then fit.
    CascadeFitter fitter(p, pi1, pi2);
    fitter.SetInvMass(mL, mXi);
    fitter.SetMaximumStep(5);
    fitter.UpdateVariance(true);
    std::vector<double> variance = sample.variances; // sigma^2 in measured-row order
    fitter.SetVariance(variance.data());
    // For correlated measurements, also call fitter.AddOffdiagonals(covariance).
    const double chi2 = fitter.DoKinematicFit();
    if (!test.FillHist(fitter, chi2, truth, trueUnknown)) continue;

    // 4. Extract physical fitted parameters and reconstruct the daughters.
    const TMatrixD best = fitter.GetBestMeasurements();
    const TMatrixD bestXi = fitter.GetBestUnknowns(); // Xi (p,theta,phi) only
    KFCheckColumn(bestXi, 3, "Fitted Xi unknowns");
    const TLorentzVector fittedP = KFFourVector(best, 0, mP);
    const TLorentzVector fittedPi1 = KFFourVector(best, 3, mPi);
    const TLorentzVector fittedPi2 = KFFourVector(best, 6, mPi);
    const TLorentzVector fittedLambda = fittedP + fittedPi1;
    const TLorentzVector fittedXi = fittedLambda + fittedPi2;
    KFFillCascadeObservables(test, trueP + truePi1, p + pi1, fittedLambda,
                            p + pi1 + pi2, fittedXi);
  }
  KFDrawCascadeHist(test, "figs/XiFit");
  test.WriteHist("figs/XiFit/histograms.root");
}

// Isospin limit: Xi1530- -> Xi0 pi- : Xi- pi0 = 2:1, so signalFraction = 1/3.
// signalFraction is the assumed Xi- pi0 fraction AFTER K+ acceptance, before cuts.
void TestXi1530Fit(int nEvents = -1, double signalFraction = 1. / 3.,
                   bool updateVariance = false) {
  const int count = nEvents < 0 ? nev : nEvents;
  if (count <= 0 || !std::isfinite(signalFraction) || signalFraction <= 0. || signalFraction >= 1.)
    throw std::invalid_argument("TestXi1530Fit needs positive events and 0 < signalFraction < 1");
  std::vector<KFCascadeFitResult> signalResults, backgroundResults;
  signalResults.reserve(count); backgroundResults.reserve(count);
  KFTestDiagnostics signal("Xi1530Signal", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFTestDiagnostics background("Xi1530Background", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFAddCascadeObservables(signal); KFAddCascadeObservables(background);
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector primary = beam + target;
  double primaryMasses[2] = {mXi1530, mK};
  double lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace evVert, evXi1530, evXi, evLambda;
  if (!evVert.SetDecay(primary, 2, primaryMasses))
    throw std::runtime_error("Kinematically forbidden Xi1530 production");
  for (int i = 0; i < count; ++i) {
    int attempts = 0;
    do {
      evVert.Generate();
      if (++attempts > 1000000) throw std::runtime_error("K+ selection has no generated acceptance");
    } while (evVert.GetDecay(1)->CosTheta() < .97);
    TLorentzVector xi1530 = *evVert.GetDecay(0);

    // Generate and fit one event per channel using the same Xi- hypothesis.
    for (int channel = 0; channel < 2; ++channel) {
      const bool isSignal = channel == 0;
      KFTestDiagnostics& test = isSignal ? signal : background;
      auto& results = isSignal ? signalResults : backgroundResults;
      TLorentzVector lambda, truePi2;
      if (isSignal) {
        // Signal: Xi1530- -> Xi- pi0, Xi- -> Lambda pi-.
        double resonanceMasses[2] = {mXi, mPi0};
        if (!evXi1530.SetDecay(xi1530, 2, resonanceMasses))
          throw std::runtime_error("Kinematically forbidden Xi- pi0 decay");
        evXi1530.Generate();
        TLorentzVector xi = *evXi1530.GetDecay(0);
        double xiMasses[2] = {mL, mPi};
        if (!evXi.SetDecay(xi, 2, xiMasses))
          throw std::runtime_error("Kinematically forbidden Xi- decay");
        evXi.Generate();
        lambda = *evXi.GetDecay(0);
        truePi2 = *evXi.GetDecay(1);
      } else {
        // Background: Xi1530- -> Xi0 pi-, Xi0 -> Lambda pi0.
        // Assign the charged pion from Xi1530 to the Xi- decay hypothesis.
        double resonanceMasses[2] = {mXi0, mPi};
        if (!evXi1530.SetDecay(xi1530, 2, resonanceMasses))
          throw std::runtime_error("Kinematically forbidden Xi0 pi- decay");
        evXi1530.Generate();
        TLorentzVector xi0 = *evXi1530.GetDecay(0);
        truePi2 = *evXi1530.GetDecay(1);
        double xi0Masses[2] = {mL, mPi0};
        if (!evXi.SetDecay(xi0, 2, xi0Masses))
          throw std::runtime_error("Kinematically forbidden Xi0 decay");
        evXi.Generate();
        lambda = *evXi.GetDecay(0);
      }
      if (!evLambda.SetDecay(lambda, 2, lambdaMasses))
        throw std::runtime_error("Kinematically forbidden Lambda decay");
      evLambda.Generate();
      const TLorentzVector trueP = *evLambda.GetDecay(0);
      const TLorentzVector truePi1 = *evLambda.GetDecay(1);
      const TMatrixD truth = KFSphericalParameters({trueP, truePi1, truePi2});
      // For background this is the wrongly assigned Lambda+pi- system.
      const TMatrixD trueUnknown = KFSphericalParameters({trueP + truePi1 + truePi2});
      const KFTestSample sample = test.Smear(truth);
      const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
      const TLorentzVector pi1 = KFFourVector(sample.measurements, 3, mPi);
      const TLorentzVector pi2 = KFFourVector(sample.measurements, 6, mPi);
      KFCascadeFitResult result;
      result.initialMass = (p + pi1 + pi2).M(); // Preserve even if the fit fails.

      // Fitter input and settings are identical for both generated channels.
      CascadeFitter fitter(p, pi1, pi2);
      fitter.SetInvMass(mL, mXi);
      fitter.SetMaximumStep(5);
      fitter.UpdateVariance(updateVariance);
      std::vector<double> variance = sample.variances;
      fitter.SetVariance(variance.data());
      const double chi2 = fitter.DoKinematicFit();
      result.chi2 = chi2;
      result.pvalue = fitter.GetPValue();
      result.valid = test.FillHist(fitter, chi2, truth, trueUnknown);
      if (result.valid) {
        const TMatrixD best = fitter.GetBestMeasurements();
        const TMatrixD bestXi = fitter.GetBestUnknowns();
        KFCheckColumn(bestXi, 3, "Fitted Xi unknowns");
        const auto pulls = fitter.GetPull(); // Same measured-row order as truth.
        result.piXiPullScore = KFPiXiPullScore(pulls); // rows 6,7,8: pi_Xi
        const TLorentzVector fittedP = KFFourVector(best, 0, mP);
        const TLorentzVector fittedPi1 = KFFourVector(best, 3, mPi);
        const TLorentzVector fittedPi2 = KFFourVector(best, 6, mPi);
        const TLorentzVector fittedLambda = fittedP + fittedPi1;
        const TLorentzVector fittedXi = fittedLambda + fittedPi2;
        KFFillCascadeObservables(test, trueP + truePi1, p + pi1, fittedLambda,
                                p + pi1 + pi2, fittedXi);
      }
      results.push_back(result); // Failed fits remain in efficiency denominators.
    }
  }
  KFDrawCascadeHist(signal, "figs/Xi1530Fit");
  KFDrawCascadeHist(background, "figs/Xi1530Fit");
  signal.WriteHist("figs/Xi1530Fit/signal.root");
  background.WriteHist("figs/Xi1530Fit/background.root");
  KFCompareSamples(signalResults, backgroundResults, "figs/Xi1530Fit", signalFraction);
}

// Fit Lambda -> p pi1 only; the generated Xi decay supplies Lambda kinematics.
void TestFourVectorFit(int nEvents = -1) {
  const auto lambdaMeasured = KFJoinVariables({
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP),
      KFParticleVariables("pi1", "#pi_{1}", .1, ResPi1, ResTh, ResPh, ScalePPi1)});
  KFTestDiagnostics lambdaTest("LambdaFourVector", lambdaMeasured,
                                 KFParticleVariables("lambda", "#Lambda", .5));
  lambdaTest.AddObservable({"Mass", "M(p#pi_{1}) [GeV/c^{2}]", {500, 1.04, 1.2}});
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector primary = beam + target;
  double primaryMasses[2] = {mXi, mK};
  double xiMasses[2] = {mL, mPi};
  double lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace evVert, evXi, evLambda;
  if (!evVert.SetDecay(primary, 2, primaryMasses))
    throw std::runtime_error("Kinematically forbidden Xi production");
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    evVert.Generate();
    TLorentzVector xi = *evVert.GetDecay(0);
    if (!evXi.SetDecay(xi, 2, xiMasses))
      throw std::runtime_error("Kinematically forbidden Xi decay");
    evXi.Generate();
    TLorentzVector lambda = *evXi.GetDecay(0);
    if (!evLambda.SetDecay(lambda, 2, lambdaMasses))
      throw std::runtime_error("Kinematically forbidden Lambda decay");
    evLambda.Generate();
    const TLorentzVector trueP = *evLambda.GetDecay(0);
    const TLorentzVector truePi1 = *evLambda.GetDecay(1);
    const TLorentzVector trueLambda = trueP + truePi1;
    const TMatrixD truth = KFSphericalParameters({trueP, truePi1});
    const KFTestSample sample = lambdaTest.Smear(truth);
    const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pi1 = KFFourVector(sample.measurements, 3, mPi);

    FourVectorFitter fitter(p, pi1, p + pi1); // Last argument: initial Lambda guess.
    fitter.SetInvMass(mL);
    std::vector<double> variance = sample.variances; // p then pi1; p,theta,phi each
    fitter.SetVariance(variance.data());
    const double chi2 = fitter.DoKinematicFit();
    if (!lambdaTest.FillHist(fitter, chi2, truth, KFSphericalParameters({trueLambda}))) continue;
    const TMatrixD best = fitter.GetBestMeasurements();
    const TMatrixD bestLambda = fitter.GetBestUnknowns();
    KFCheckColumn(bestLambda, 3, "Fitted Lambda unknowns");
    const TLorentzVector fittedP = KFFourVector(best, 0, mP);
    const TLorentzVector fittedPi1 = KFFourVector(best, 3, mPi);
    lambdaTest.FillObservable("Mass", (p + pi1).M(), (fittedP + fittedPi1).M());
  }
  lambdaTest.DrawHist("figs/FourVectorFit");
  lambdaTest.WriteHist("figs/FourVectorFit/lambda.root");
}

// K- p -> Lambda eta, Lambda -> p pi-; eta is not measured.
// Measured rows: p, pion, beam; unknown rows: Lambda (p,theta,phi).
double pE72 = 0.78;
void TestMissingEtaFit(int nEvents = -1) {
  auto variables = KFJoinVariables({
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP),
      KFParticleVariables("pion", "#pi", .1, ResPi1, ResTh, ResPh, ScalePPi1),
      KFParticleVariables("beam", "K^{-}", .02, .001, .001, 1.)});
  KFTestDiagnostics test("MissingEta", variables, KFParticleVariables("lambda", "#Lambda", .5));
  test.AddObservable({"LambdaMass", "M(p#pi) [GeV/c^{2}]", {500, 1.04, 1.2}});
  test.AddObservable({"MissingMass", "M_{X} [GeV/c^{2}]", {500, mEta - .2, mEta + .2}});
  // A small coordinate tilt keeps the spherical beam angles away from the pole.
  TVector3 beamE72(0.1,0,1);
  beamE72.SetMag(pE72);
  const TLorentzVector beam(beamE72, std::hypot(mK, pE72));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector primary = beam + target;
  double primaryMasses[2] = {mL, mEta};
  double lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace evVert, evLambda;
  if (!evVert.SetDecay(primary, 2, primaryMasses))
    throw std::runtime_error("Kinematically forbidden Lambda eta production");
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    evVert.Generate();
    TLorentzVector lambda = *evVert.GetDecay(0);
    if (!evLambda.SetDecay(lambda, 2, lambdaMasses))
      throw std::runtime_error("Kinematically forbidden Lambda decay");
    evLambda.Generate();
    const TMatrixD truth = KFSphericalParameters({*evLambda.GetDecay(0), *evLambda.GetDecay(1), beam});
    const KFTestSample sample = test.Smear(truth);
    const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pion = KFFourVector(sample.measurements, 3, mPi);
    const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 6, mK);
    MissingMassFitter fitter(measuredBeam, target, mEta, p, pion, p + pion);
    fitter.SetInvMass(mL);
    std::vector<double> variance = sample.variances; // p, pion, beam; p,theta,phi each
    fitter.SetVariance(variance.data());
    const double chi2 = fitter.DoKinematicFit();
    if (!test.FillHist(fitter, chi2, truth, KFSphericalParameters({lambda}))) continue;
    const TMatrixD best = fitter.GetBestMeasurements();
    const TMatrixD bestUnknown = fitter.GetBestUnknowns();
    const TLorentzVector fittedP = KFFourVector(best, 0, mP);
    const TLorentzVector fittedPion = KFFourVector(best, 3, mPi);
    const TLorentzVector fittedBeam = KFFourVector(best, 6, mK);
    const TLorentzVector fittedLambda = KFFourVector(bestUnknown, 0, mL);
    test.FillObservable("LambdaMass", (p + pion).M(),
                         (fittedP + fittedPion).M());
    test.FillObservable("MissingMass", (measuredBeam + target - (p + pion)).M(),
                         (fittedBeam + target - fittedLambda).M());
  }
  test.DrawHist("figs/MissingEtaFit");
  test.WriteHist("figs/MissingEtaFit/histograms.root");
}

// All balance-fit rows are measured: beam first, then each final particle,
// with (p,theta,phi) for each. There are no unknowns or intermediate phi rows.
std::vector<KFTestVariable> KFPbarBalanceMeasurements(bool phiPhi) {
  auto variables = KFParticleVariables("beam", "#bar{p}_{beam}", .01,
                                        ResPbarBeam, ResThPbarBeam, ResPhPbarBeam);
  const auto finals = phiPhi
      ? KFJoinVariables({
          KFParticleVariables("km1", "K^{-}_{1}", .15, ResP, ResThP, ResPhP),
          KFParticleVariables("kp1", "K^{+}_{1}", .15, ResP, ResThP, ResPhP),
          KFParticleVariables("km2", "K^{-}_{2}", .15, ResP, ResThP, ResPhP),
          KFParticleVariables("kp2", "K^{+}_{2}", .15, ResP, ResThP, ResPhP)})
      : KFJoinVariables({
          KFParticleVariables("pbar", "#bar{p}", .2, ResP, ResThP, ResPhP),
          KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP)});
  variables.insert(variables.end(), finals.begin(), finals.end());
  return variables;
}

std::vector<KFTestVariable> KFPbarBalanceConstraints() {
  return {
    {"Px", "#Sigma p_{x,final} - p_{x,beam} [GeV/c]", {200, -.25, .25}},
    {"Py", "#Sigma p_{y,final} - p_{y,beam} [GeV/c]", {200, -.25, .25}},
    {"Pz", "#Sigma p_{z,final} - p_{z,beam} [GeV/c]", {200, -.25, .25}},
    {"Energy", "#Sigma E_{final} - E_{beam} - m_{p} [GeV]", {200, -.25, .25}}
  };
}

struct KFPbarBalanceSummary {
  unsigned long invalid = 0, nonconserving = 0;
  double maxResidual = 0.;
};

// Diagnostic only: inputs are already extracted fitted four-vectors.
void KFCheckBalanceConservation(const TLorentzVector& fittedBeam,
                                const TLorentzVector& target,
                                const std::vector<TLorentzVector>& fittedFinal,
                                KFPbarBalanceSummary& summary) {
  TLorentzVector residual = -fittedBeam - target;
  for (const auto& particle : fittedFinal) residual += particle;
  const double components[4] = {residual.Px(), residual.Py(), residual.Pz(), residual.E()};
  bool conserved = true;
  for (double component : components) {
    if (!std::isfinite(component) || std::abs(component) > 1.e-8) conserved = false;
    if (std::isfinite(component))
      summary.maxResidual = std::max(summary.maxResidual, std::abs(component));
  }
  if (!conserved) ++summary.nonconserving;
}

void KFFinishPbarBalance(KFTestDiagnostics& test, const TString& directory,
                          const KFPbarBalanceSummary& summary,
                          const std::vector<KFTestDiagnostics::CanvasPage>& measuredPages = {}) {
  test.DrawHist(directory, {}, measuredPages);
  test.WriteHist(directory + "/histograms.root");
  std::cout << "Balance check: mean chi2 = " << test.Chi2Histogram()->GetMean()
            << " (expect about 4), mean p-value = " << test.PValueHistogram()->GetMean()
            << " (expect about 0.5)\n"
            << "  Invalid results: " << summary.invalid
            << ", nonconserving finite results: " << summary.nonconserving
            << ", max |F_i| = " << summary.maxResidual << " GeV" << std::endl;
  // Finite and unconverged fits remain in the histograms. Invalid results and
  // failed conservation are counted and reported without aborting the run.
}

void TestPbarElasticFit(int nEvents = -1, double beamMomentum = -1.) {
  const int count = nEvents < 0 ? nev : nEvents;
  if (count <= 0) throw std::invalid_argument("TestPbarElasticFit: event count must be positive");
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (!std::isfinite(momentum) || momentum <= 0.)
    throw std::invalid_argument("Antiproton beam momentum must be finite and positive");
  // Slightly tilt the coordinate axes to avoid the spherical-angle pole.
  TVector3 beamDirection(.1, .05, 1.);
  beamDirection.SetMag(momentum);
  const TLorentzVector beam(beamDirection, std::hypot(momentum, mP));
  const TLorentzVector target(0., 0., 0., mP);
  const auto variables = KFPbarBalanceMeasurements(false);
  KFTestDiagnostics test("PbarElastic", variables, {}, KFPbarBalanceConstraints());
  const double sqrtS = (beam + target).M();
  test.AddObservable({"FinalMassResidual", "M(#bar{p}p)-M^{true}(#bar{p}p) [GeV/c^{2}]", {200, -.15, .15}});
  test.AddObservable({"MissingMassResidual", "M(B+T-#bar{p})-m_{p} [GeV/c^{2}]", {200, -.15, .15}});
  TLorentzVector initialSystem = beam + target;
  double finalMasses[2] = {mP, mP};
  TGenPhaseSpace primary;
  if (!primary.SetDecay(initialSystem, 2, finalMasses))
    throw std::runtime_error("Kinematically forbidden elastic production");
  KFPbarBalanceSummary summary;
  // Isotropic two-body phase space is used as a kinematic performance test,
  // not as a model of the elastic differential cross section.
  for (int i = 0; i < count; ++i) {
    primary.Generate();
    const TLorentzVector truePbar = *primary.GetDecay(0);
    const TLorentzVector trueProton = *primary.GetDecay(1);
    const TMatrixD truth = KFSphericalParameters({beam, truePbar, trueProton});
    const KFTestSample sample = test.Smear(truth);
    // BalanceFitter row order: beam first, then the two final particles.
    const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pbar = KFFourVector(sample.measurements, 3, mP);
    const TLorentzVector proton = KFFourVector(sample.measurements, 6, mP);
    const std::vector<TLorentzVector> measuredFinal = {pbar, proton};

    BalanceFitter fitter(measuredBeam, target, measuredFinal);
    if (fitter.GetNDF() != 4)
      throw std::runtime_error("BalanceFitter must have four constraints and no unknowns");
    fitter.SetChi2DifCut(1.e-10);
    std::vector<double> variance = sample.variances;
    fitter.SetVariance(variance.data());
    const double chi2 = fitter.DoKinematicFit();
    // Measured-only fit: no unknown getters or unknown truth column.
    if (!test.FillHist(fitter, chi2, truth)) { ++summary.invalid; continue; }

    const TMatrixD best = fitter.GetBestMeasurements(); // (p,theta,phi) per particle
    const auto fittedLV = fitter.GetFittedLV(); // beam, pbar, proton
    if (fittedLV.size() != 3)
      throw std::runtime_error("Unexpected BalanceFitter four-vector count");
    KFCheckColumn(best, 9, "Fitted balance measurements");
    const TLorentzVector fittedBeam = fittedLV[0];
    const TLorentzVector fittedPbar = fittedLV[1];
    const TLorentzVector fittedProton = fittedLV[2];
    KFCheckBalanceConservation(fittedBeam, target, {fittedPbar, fittedProton}, summary);
    test.FillObservable("FinalMassResidual", (pbar + proton).M() - sqrtS,
                         (fittedPbar + fittedProton).M() - sqrtS);
    test.FillObservable("MissingMassResidual", (measuredBeam + target - pbar).M() - mP,
                         (fittedBeam + target - fittedPbar).M() - mP);
  }
  KFFinishPbarBalance(test, "figs/PbarElasticFit", summary);
}

void TestPbarPhiPhiFit(int nEvents = -1, double beamMomentum = -1.) {
  const int count = nEvents < 0 ? nev : nEvents;
  if (count <= 0) throw std::invalid_argument("TestPbarPhiPhiFit: event count must be positive");
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (!std::isfinite(momentum) || momentum <= 0.)
    throw std::invalid_argument("Antiproton beam momentum must be finite and positive");
  // Slightly tilt the coordinate axes to avoid the spherical-angle pole.
  TVector3 beamDirection(.1, .05, 1.);
  beamDirection.SetMag(momentum);
  const TLorentzVector beam(beamDirection, std::hypot(momentum, mP));
  const TLorentzVector target(0., 0., 0., mP);
  const double sqrtS = (beam + target).M();
  if (!std::isfinite(mPhi) || !std::isfinite(mK) || mK <= 0. ||
      !std::isfinite(WidthPhi) || WidthPhi < 0. || mPhi <= 2. * mK ||
      sqrtS <= 4. * mK || (WidthPhi == 0. && sqrtS <= 2. * mPhi))
    throw std::invalid_argument("TestPbarPhiPhiFit: invalid phi width/mass or production is forbidden");
  std::cout << "pbar p -> phi phi: p_beam = " << beam.P()
            << " GeV/c, sqrt(s) = " << sqrtS
            << " GeV, excess above two pole masses = " << sqrtS - 2. * mPhi
            << " GeV, phi width = " << WidthPhi << " GeV" << std::endl;
  const auto variables = KFPbarBalanceMeasurements(true);
  KFTestDiagnostics test("PbarPhiPhi", variables, {}, KFPbarBalanceConstraints());
  test.AddObservable({"Phi1MassResidual", "M(K^{-}_{1}K^{+}_{1})-M^{true}_{#phi_{1}} [GeV/c^{2}]",
                       {200, -.04, .04}});
  test.AddObservable({"Phi2MassResidual", "M(K^{-}_{2}K^{+}_{2})-M^{true}_{#phi_{2}} [GeV/c^{2}]",
                       {200, -.04, .04}});
  test.AddObservable({"FourKaonMassResidual", "M(4K)-M^{true}(4K) [GeV/c^{2}]",
                       {200, -.15, .15}});
  TGenPhaseSpace primary, decay1, decay2;
  KFPbarBalanceSummary summary;
  // Draw two masses for each event, then generate the production and decays.
  // No phi mass constraints are imposed in the balance fit: it is still 4C.
  for (int i = 0; i < count; ++i) {
    // 1. Generate two Breit-Wigner masses with allowed production phase space.
    // ROOT's width is FWHM. Angular distributions remain isotropic.
    double phiMass1 = mPhi, phiMass2 = mPhi;
    if (WidthPhi > 0.) {
      const double minimum = 2. * mK;
      const double s = sqrtS * sqrtS;
      const double maxMomentum = .5 * std::sqrt(s - 4. * minimum * minimum);
      bool accepted = false;
      for (int attempt = 0; attempt < 1000000; ++attempt) {
        phiMass1 = gRandom->BreitWigner(mPhi, WidthPhi);
        phiMass2 = gRandom->BreitWigner(mPhi, WidthPhi);
        if (!std::isfinite(phiMass1) || !std::isfinite(phiMass2) ||
            phiMass1 <= minimum || phiMass2 <= minimum || phiMass1 + phiMass2 >= sqrtS)
          continue;
        const double sum = phiMass1 + phiMass2, difference = phiMass1 - phiMass2;
        const double pStar = std::sqrt((s - sum * sum) * (s - difference * difference)) / (2. * sqrtS);
        if (gRandom->Rndm() * maxMomentum < pStar) { accepted = true; break; }
      }
      if (!accepted) throw std::runtime_error("Unable to generate an allowed phi mass pair");
    }
    TLorentzVector initialSystem = beam + target;
    double phiMasses[2] = {phiMass1, phiMass2};
    double kaonMasses[2] = {mK, mK};
    if (!primary.SetDecay(initialSystem, 2, phiMasses))
      throw std::runtime_error("Kinematically forbidden phi phi production");
    primary.Generate();
    TLorentzVector phi1 = *primary.GetDecay(0), phi2 = *primary.GetDecay(1);
    if (!decay1.SetDecay(phi1, 2, kaonMasses) || !decay2.SetDecay(phi2, 2, kaonMasses))
      throw std::runtime_error("Kinematically forbidden phi decay");
    decay1.Generate();
    decay2.Generate();
    const TLorentzVector trueKm1 = *decay1.GetDecay(0);
    const TLorentzVector trueKp1 = *decay1.GetDecay(1);
    const TLorentzVector trueKm2 = *decay2.GetDecay(0);
    const TLorentzVector trueKp2 = *decay2.GetDecay(1);

    // 2. Smear beam + four kaons, in exactly the fitter's measured-row order.
    const TMatrixD truth = KFSphericalParameters({beam, trueKm1, trueKp1, trueKm2, trueKp2});
    const KFTestSample sample = test.Smear(truth);
    const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector km1 = KFFourVector(sample.measurements, 3, mK);
    const TLorentzVector kp1 = KFFourVector(sample.measurements, 6, mK);
    const TLorentzVector km2 = KFFourVector(sample.measurements, 9, mK);
    const TLorentzVector kp2 = KFFourVector(sample.measurements, 12, mK);
    const std::vector<TLorentzVector> measuredFinal = {km1, kp1, km2, kp2};

    // 3. Fit total four-momentum conservation, with no phi mass constraints.
    BalanceFitter fitter(measuredBeam, target, measuredFinal);
    if (fitter.GetNDF() != 4)
      throw std::runtime_error("BalanceFitter must have four constraints and no unknowns");
    fitter.SetChi2DifCut(1.e-10);
    std::vector<double> variance = sample.variances; // 15 diagonal variances
    fitter.SetVariance(variance.data());
    // fitter.AddOffdiagonals(covariance); // Optional for correlated inputs.
    const double chi2 = fitter.DoKinematicFit();
    if (!test.FillHist(fitter, chi2, truth)) { ++summary.invalid; continue; }

    // 4. Extract fitted measurements and four-vectors explicitly.
    const TMatrixD best = fitter.GetBestMeasurements();
    const auto fittedLV = fitter.GetFittedLV(); // beam, Km1, Kp1, Km2, Kp2
    if (fittedLV.size() != 5)
      throw std::runtime_error("Unexpected BalanceFitter four-vector count");
    KFCheckColumn(best, 15, "Fitted balance measurements");
    const TLorentzVector fittedBeam = fittedLV[0];
    const TLorentzVector fittedKm1 = fittedLV[1];
    const TLorentzVector fittedKp1 = fittedLV[2];
    const TLorentzVector fittedKm2 = fittedLV[3];
    const TLorentzVector fittedKp2 = fittedLV[4];
    KFCheckBalanceConservation(fittedBeam, target,
                               {fittedKm1, fittedKp1, fittedKm2, fittedKp2}, summary);
    // Subtract each event's actual generated mass, not the phi pole mass.
    test.FillObservable("Phi1MassResidual", (km1 + kp1).M() - phi1.M(),
                         (fittedKm1 + fittedKp1).M() - phi1.M());
    test.FillObservable("Phi2MassResidual", (km2 + kp2).M() - phi2.M(),
                         (fittedKm2 + fittedKp2).M() - phi2.M());
    test.FillObservable("FourKaonMassResidual", (km1 + kp1 + km2 + kp2).M() - sqrtS,
                         (fittedKm1 + fittedKp1 + fittedKm2 + fittedKp2).M() - sqrtS);
  }
  // Beam on its own page; four kaons in columns, with p/theta/phi in rows.
  const std::vector<KFTestDiagnostics::CanvasPage> measuredPages = {
      {3, 1, {0, 1, 2}},
      {4, 3, {3, 6, 9, 12,
              4, 7, 10, 13,
              5, 8, 11, 14}}};
  KFFinishPbarBalance(test, "figs/PbarPhiPhiFit", summary, measuredPages);
}

// Both tests default to the 1 GeV/c pbar beam configured by pPbar.
void TestBalanceFit(int nEvents = -1, double beamMomentum = -1.) {
  TestPbarElasticFit(nEvents, beamMomentum);
  TestPbarPhiPhiFit(nEvents, beamMomentum);
}

void TestKinfit() {
  SetStyle();
  std::cout << "TestXiFit([nEvents])\n"
               "TestXi1530Fit([nEvents], [signalFraction=1./3.], [updateVariance=false])\n"
               "TestFourVectorFit([nEvents])\n"
               "TestMissingEtaFit([nEvents])\n"
               "TestPbarElasticFit([nEvents], [beamMomentum])\n"
               "TestPbarPhiPhiFit([nEvents], [beamMomentum])\n"
               "TestBalanceFit([nEvents], [beamMomentum]) runs both pbar tests (default 1 GeV/c)\n"
               "Variable setup: KFTestVariable(name, axisTitle, {bins,min,max}, sigma, period)\n"
               "Residual convention: initial/fitted minus truth; matrix row order must match.\n";
}
