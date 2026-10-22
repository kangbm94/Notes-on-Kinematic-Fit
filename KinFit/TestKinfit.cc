#include "src/FourVectorFitter.cc"
#include "src/MissingMassFitter.cc"
#include "src/CascadeFitter.cc"
#include "src/BalanceFitter.cc"
#include "TestKinfit.hh"

#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <limits>
#include <fstream>
#include <iomanip>
#include <utility>

double margin = 0.02;
// Load with: .L TestKinfit.cc+
// Existing fitter sources are used unchanged. Besides their existing API, they
// must provide GetBestMeasurements(), GetInitialMeasurements(), GetBestUnknowns()
// and GetInitialUnknowns(), all returning physical, unscaled TMatrixD columns.
// Histogram configuration, smearing, filling and drawing are reusable helpers.
// Every Test... function shows generation, fitter input, SetVariance(),
// DoKinematicFit(), and fitted-parameter extraction directly in its event loop.
// Two-body angular weights are constant. The sequential four-pion BG example
// samples intermediate masses with the appropriate four-body phase-space weight.

int nev = 100000;
double mXi = 1.321, mXi0 = 1.315, mXi1530 = 1.535;
double mL = 1.115, mP = .938, mK = .493, mPi = .139, mPi0 = .135, mEta = .5478;
double pK18 = 1.8;
// True relative momentum resolutions used to smear generated momenta.
double ResP = .05, ResPi1 = .05, ResPi2 = .05;
// Independent relative resolutions used only for the fitter covariance.
// These start equal to the true resolutions; changing one does not change the other.
double ResPEstm = ResP, ResPi1Estm = ResPi1, ResPi2Estm = ResPi2;
double ResThP = .005, ResTh = .005, ResThPi2 = .005;
double ResPhP = .01, ResPh = .02, ResPhPi2 = .02;
double ScalePP = 1., ScalePPi1 = 1., ScalePPi2 = 1.;

// Balance examples: a fixed proton target and a measured antiproton beam.
double pPbar = 1.;
double mPhi = 1.019461;
double WidthPhi = .004249; // Full Breit-Wigner width in GeV (PDG: 4.249 MeV).
double ResPbarBeam = .001, ResThPbarBeam = .001, ResPhPbarBeam = .005;
double ResPbarBeamEstm = ResPbarBeam;
double ResKBeam = .02, ResKBeamEstm = ResKBeam; // Missing-eta K- beam.

// A particle helper is convenient for these examples. The generic tester also
// accepts individual Cartesian components, energies, masses, etc.; no grouping
// or fixed number of rows is required by KFTestDiagnostics.
std::vector<KFTestVariable> KFParticleVariables(
    TString name, TString title, double momentumResidualRange,
    double relativeMomentumSigma = 0., double thetaSigma = 0., double phiSigma = 0.,
    double momentumScale = 1., double relativeMomentumSigmaEstm = -1.) {
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
  // The optional final argument sets the covariance resolution independently.
  // Negative (omitted) means use the same sigma as smearing, as before.
  if (relativeMomentumSigmaEstm >= 0.)
    vars[0].resolutionEstm = KFRelativeResolution(relativeMomentumSigmaEstm);
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
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP, ResPEstm),
      KFParticleVariables("pi1", "#pi_{1}", .1, ResPi1, ResTh, ResPh, ScalePPi1, ResPi1Estm),
      KFParticleVariables("pi2", "#pi_{2}", .1, ResPi2, ResThPi2, ResPhPi2, ScalePPi2, ResPi2Estm)});
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
    std::vector<double> variance = sample.variances; // estimated sigma^2 in measured-row order
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
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP, ResPEstm),
      KFParticleVariables("pi1", "#pi_{1}", .1, ResPi1, ResTh, ResPh, ScalePPi1, ResPi1Estm)});
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
      KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP, ResPEstm),
      KFParticleVariables("pion", "#pi", .1, ResPi1, ResTh, ResPh, ScalePPi1, ResPi1Estm),
      KFParticleVariables("beam", "K^{-}", .05, ResKBeam, .01, .01, 1., ResKBeamEstm)});
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
                                        ResPbarBeam, ResThPbarBeam, ResPhPbarBeam, 1., ResPbarBeamEstm);
  const auto finals = phiPhi
      ? KFJoinVariables({
          KFParticleVariables("km1", "K^{-}_{1}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
          KFParticleVariables("kp1", "K^{+}_{1}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
          KFParticleVariables("km2", "K^{-}_{2}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
          KFParticleVariables("kp2", "K^{+}_{2}", .15, ResP, ResThP, ResPhP, 1., ResPEstm)})
      : KFJoinVariables({
          KFParticleVariables("pbar", "#bar{p}", .2, ResP, ResThP, ResPhP, 1., ResPEstm),
          KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, 1., ResPEstm)});
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

// One event record for the no-PID separation plots. Invalid/unconverged fits
// remain in the efficiency denominator, but fail every KF selection.
struct KFPbarBGResult {
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  double pvalue = std::numeric_limits<double>::quiet_NaN();
  double energyResidual = std::numeric_limits<double>::quiet_NaN();
  double pairMass1 = std::numeric_limits<double>::quiet_NaN();
  double pairMass2 = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
};

struct KFPbarBGScan {
  std::vector<double> cuts, signalEfficiency, backgroundEfficiency, purity;
};

// Plot/counting helper only; event generation and fitter usage are below.
KFPbarBGScan KFScanPbarBG(const std::vector<KFPbarBGResult>& signal,
                         const std::vector<KFPbarBGResult>& background,
                         double signalFraction) {
  if (signal.empty() || background.empty() || !std::isfinite(signalFraction) ||
      signalFraction <= 0. || signalFraction >= 1.)
    throw std::invalid_argument("Invalid pbar background separation samples/mixture");
  std::vector<double> sp, bp;
  auto collect = [](const std::vector<KFPbarBGResult>& events, std::vector<double>& probabilities) {
    for (const auto& event : events)
      if (event.valid && std::isfinite(event.pvalue) && event.pvalue >= 0. && event.pvalue <= 1.)
        probabilities.push_back(event.pvalue);
    std::sort(probabilities.begin(), probabilities.end());
  };
  collect(signal, sp); collect(background, bp);
  // Include P=0 and P=1, ordinary cuts, and tiny probabilities from the wrong
  // mass hypothesis. Counts use event values, never histogram bin integrals.
  std::vector<double> cuts = {0.};
  for (int i = 1; i <= 1000; ++i) cuts.push_back(i / 1000.);
  for (int i = 0; i <= 1280; ++i) cuts.push_back(std::pow(10., -320. + .25 * i));
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  KFPbarBGScan scan;
  for (double cut : cuts) {
    const double es = double(sp.end() - std::lower_bound(sp.begin(), sp.end(), cut)) / signal.size();
    const double eb = double(bp.end() - std::lower_bound(bp.begin(), bp.end(), cut)) / background.size();
    const double selected = signalFraction * es + (1. - signalFraction) * eb;
    scan.cuts.push_back(cut);
    scan.signalEfficiency.push_back(es);
    scan.backgroundEfficiency.push_back(eb);
    scan.purity.push_back(selected > 0. ? signalFraction * es / selected
                                       : std::numeric_limits<double>::quiet_NaN());
  }
  return scan;
}

void KFDrawPbarBGSeparation(const std::vector<KFPbarBGResult>& signal,
                            const std::vector<KFPbarBGResult>& background,
                            const TString& directory, double signalFraction) {
  const auto scan = KFScanPbarBG(signal, background, signalFraction);
  gSystem->mkdir(directory, true);
  TH1D sChi2("PbarBGSignalChi2", ";#chi^{2} (4K hypothesis);Events", 1000, 0., 2000.);
  TH1D bChi2("PbarBGBackgroundChi2", ";#chi^{2} (4K hypothesis);Events", 1000, 0., 2000.);
  TH1D sPval("PbarBGSignalPValue", ";P-value (4K hypothesis);Events", 200, 0., 1.);
  TH1D bPval("PbarBGBackgroundPValue", ";P-value (4K hypothesis);Events", 200, 0., 1.);
  TH1D sEnergy("PbarBGSignalEnergy", ";#Sigma E_{K} - E_{beam} - m_{p} [GeV];Events", 500, -.5, 2.5);
  TH1D bEnergy("PbarBGBackgroundEnergy", ";#Sigma E_{K} - E_{beam} - m_{p} [GeV];Events", 500, -.5, 2.5);
  TH1D sMass1("PbarBGSignalPair1", ";M_{1} (K mass assignment) [GeV/c^{2}];Events", 500, .7, 2.2);
  TH1D bMass1("PbarBGBackgroundPair1", ";M_{1} (K mass assignment) [GeV/c^{2}];Events", 500, .7, 2.2);
  TH1D sMass2("PbarBGSignalPair2", ";M_{2} (K mass assignment) [GeV/c^{2}];Events", 500, .7, 2.2);
  TH1D bMass2("PbarBGBackgroundPair2", ";M_{2} (K mass assignment) [GeV/c^{2}];Events", 500, .7, 2.2);
  for (auto* h : {&sChi2, &bChi2, &sPval, &bPval, &sEnergy, &bEnergy, &sMass1, &bMass1, &sMass2, &bMass2}) {
    h->SetDirectory(nullptr); h->Sumw2();
  }
  auto fill = [](const std::vector<KFPbarBGResult>& events, TH1D& chi2, TH1D& pval,
                 TH1D& energy, TH1D& mass1, TH1D& mass2) {
    for (const auto& event : events) {
      if (std::isfinite(event.energyResidual)) energy.Fill(event.energyResidual);
      if (std::isfinite(event.pairMass1)) mass1.Fill(event.pairMass1);
      if (std::isfinite(event.pairMass2)) mass2.Fill(event.pairMass2);
      if (!event.valid) continue;
      if (std::isfinite(event.chi2) && event.chi2 >= 0.) chi2.Fill(event.chi2);
      if (std::isfinite(event.pvalue) && event.pvalue >= 0. && event.pvalue <= 1.)
        pval.Fill(event.pvalue == 1. ? std::nextafter(1., 0.) : event.pvalue);
    }
  };
  fill(signal, sChi2, sPval, sEnergy, sMass1, sMass2);
  fill(background, bChi2, bPval, bEnergy, bMass1, bMass2);
  const double total = signal.size() + background.size();
  for (auto* h : {&sChi2, &sPval, &sEnergy, &sMass1, &sMass2}) h->Scale(total * signalFraction / signal.size());
  for (auto* h : {&bChi2, &bPval, &bEnergy, &bMass1, &bMass2}) h->Scale(total * (1. - signalFraction) / background.size());
  auto overlay = [](const TH1D& signalHist, const TH1D& backgroundHist) {
    gPad->SetLogy();
    TH1D sum(signalHist); sum.SetDirectory(nullptr); sum.Add(&backgroundHist);
    sum.SetLineColor(kBlack); sum.SetMinimum(.1);
    sum.SetMaximum(2. * std::max(1., sum.GetMaximum()));
    auto* all = sum.DrawCopy("hist");
    auto* sh = signalHist.DrawCopy("hist same"); sh->SetLineColor(kBlue);
    auto* bh = backgroundHist.DrawCopy("hist same"); bh->SetLineColor(kRed);
    auto* legend = new TLegend(.48, .70, .91, .91);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(all, "Assumed mixture", "l");
    legend->AddEntry(sh, "#phi#phi #rightarrow 4K", "l");
    legend->AddEntry(bh, "Selected 4#pi #rightarrow 4K hypothesis", "l"); legend->Draw();
  };
  auto* comparison = new TCanvas("PbarBGComparison", "No-PID separation", 1500, 600);
  comparison->Divide(3, 1);
  comparison->cd(1); overlay(sChi2, bChi2);
  comparison->cd(2); overlay(sPval, bPval);
  comparison->cd(3); overlay(sEnergy, bEnergy);
  comparison->SaveAs(directory + "/comparison.pdf");
  auto* masses = new TCanvas("PbarBGPairMasses", "Pair masses before fitting", 1200, 600);
  masses->Divide(2, 1);
  masses->cd(1); overlay(sMass1, bMass1);
  masses->cd(2); overlay(sMass2, bMass2);
  masses->SaveAs(directory + "/pair_masses.pdf");
  auto graph = [](const std::vector<double>& x, const std::vector<double>& y,
                  const char* name, int color, bool positiveX) {
    auto* g = new TGraph(); g->SetName(name); g->SetLineColor(color); g->SetLineWidth(2);
    for (std::size_t i = 0; i < x.size(); ++i)
      if (std::isfinite(x[i]) && std::isfinite(y[i]) && (!positiveX || x[i] > 0.))
        g->SetPoint(g->GetN(), x[i], y[i]);
    return g;
  };
  auto* efficiency = graph(scan.cuts, scan.signalEfficiency, "PbarBGSignalEfficiency", kBlue, true);
  auto* backgroundSurvival = graph(scan.cuts, scan.backgroundEfficiency, "PbarBGBackgroundEfficiency", kGreen + 2, true);
  auto* purity = graph(scan.cuts, scan.purity, "PbarBGPurity", kRed, true);
  auto* tradeoff = graph(scan.signalEfficiency, scan.purity, "PbarBGEfficiencyPurity", kBlue, true);
  std::vector<double> rejection;
  for (double eb : scan.backgroundEfficiency) rejection.push_back(1. - eb);
  auto* roc = graph(scan.signalEfficiency, rejection, "PbarBGRejection", kRed, true);
  auto* cuts = new TCanvas("PbarBGCuts", "P-value selection", 900, 650);
  cuts->cd(); gPad->SetLogx();
  gPad->DrawFrame(1.e-8, 0., 1., 1.05)->SetTitle(";P-value lower cut;Efficiency / purity");
  efficiency->Draw("L same"); backgroundSurvival->Draw("L same"); purity->Draw("L same");
  auto* legend = new TLegend(.20, .18, .65, .42);
  legend->SetBorderSize(0); legend->SetFillStyle(0);
  legend->AddEntry(efficiency, "Signal efficiency", "l");
  legend->AddEntry(backgroundSurvival, "Background survival", "l");
  legend->AddEntry(purity, TString::Format("Purity (f_{S}=%.3g)", signalFraction), "l"); legend->Draw();
  cuts->SaveAs(directory + "/c_pval.pdf");
  const double xmin = .5 / signal.size();
  auto* efficiencyPurity = new TCanvas("PbarBGEffPur", "Efficiency versus purity", 900, 650);
  efficiencyPurity->cd(); gPad->SetLogx();
  gPad->DrawFrame(xmin, 0., 1., 1.05)->SetTitle(";Signal efficiency;Purity");
  tradeoff->Draw("L same"); efficiencyPurity->SaveAs(directory + "/c_eff_pur.pdf");
  auto* rejectionCanvas = new TCanvas("PbarBGROC", "Background rejection", 900, 650);
  rejectionCanvas->cd(); gPad->SetLogx();
  gPad->DrawFrame(xmin, 0., 1., 1.05)->SetTitle(";Signal efficiency;Background rejection (1-#epsilon_{B})");
  roc->Draw("L same"); rejectionCanvas->SaveAs(directory + "/c_rejection.pdf");
  TFile output(directory + "/separation.root", "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create pbar background separation ROOT file");
  for (auto* h : {&sChi2, &bChi2, &sPval, &bPval, &sEnergy, &bEnergy, &sMass1, &bMass1, &sMass2, &bMass2}) h->Write();
  for (auto* g : {efficiency, backgroundSurvival, purity, tradeoff, roc}) g->Write();
  for (auto* c : {comparison, masses, cuts, efficiencyPurity, rejectionCanvas}) c->Write();
  std::cout << "No-PID separation: efficiencies relative to all signal events and all PRESELECTED background events.\n"
            << "Assumed signal fraction in these samples: " << signalFraction
            << " (a configurable mixture, not a branching-ratio prediction).\n"
            << "Invalid/nonconserving fits fail all P-value cuts.\n"
            << "P cut    signal efficiency    background survival    purity\n";
  for (double cut : {0., .001, .01, .05, .1, .5}) {
    const auto i = std::lower_bound(scan.cuts.begin(), scan.cuts.end(), cut) - scan.cuts.begin();
    std::cout << scan.cuts[i] << "    " << scan.signalEfficiency[i] << "    "
              << scan.backgroundEfficiency[i] << "    " << scan.purity[i] << '\n';
  }
}

// pbar p -> phi phi -> 4K versus preselected, uniform four-body 4pi phase space.
// No external PID: measured momenta in BOTH channels receive kaon masses.
// BG selection is before smearing, using kaon-assigned masses by default.
// Set useKaonMassWindow=false for the previous true-pion-mass selection.
// nEvents is the number of signal events AND the number of accepted BG events.
// signalFraction refers to these samples AFTER the BG-only preselection.
void TestPbarBGSeparation(int nEvents = -1, double beamMomentum = -1.,
                          double signalFraction = .5, double massWindow = .100,
                          bool useKaonMassWindow = true) {
  const int count = nEvents < 0 ? nev : nEvents;
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (count <= 0 || !std::isfinite(momentum) || momentum <= 0. ||
      !std::isfinite(signalFraction) || signalFraction <= 0. || signalFraction >= 1. ||
      !std::isfinite(massWindow) || massWindow <= 0.)
    throw std::invalid_argument("Invalid TestPbarBGSeparation arguments");
  TVector3 direction(.1, .05, 1.); direction.SetMag(momentum);
  const TLorentzVector beam(direction, std::hypot(momentum, mP));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector initialSystem = beam + target;
  const double sqrtS = initialSystem.M();
  if (!std::isfinite(mPhi) || !std::isfinite(mK) || mK <= 0. ||
      !std::isfinite(mPi) || mPi <= 0. || !std::isfinite(WidthPhi) || WidthPhi < 0. ||
      mPhi <= 2. * mK || sqrtS <= 4. * mK ||
      (WidthPhi == 0. && sqrtS <= 2. * mPhi) || mK < mPi)
    throw std::invalid_argument("Signal production or selected four-pion background is kinematically forbidden");
  const TString directory = "figs/PbarBGSeparation";
  // Generic track labels: the BG selection permits all disjoint pairings.
  const auto variables = KFJoinVariables({
      KFParticleVariables("beam", "#bar{p}_{beam}", .01, ResPbarBeam, ResThPbarBeam, ResPhPbarBeam, 1., ResPbarBeamEstm),
      KFParticleVariables("track1", "K_{1}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("track2", "K_{2}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("track3", "K_{3}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("track4", "K_{4}", .15, ResP, ResThP, ResPhP, 1., ResPEstm)});
  KFTestDiagnostics signal("PbarBGSignal", variables, {}, KFPbarBalanceConstraints());
  KFTestDiagnostics background("PbarBGBackground", variables, {}, KFPbarBalanceConstraints());
  for (auto* test : {&signal, &background}) {
    test->AddObservable({"Pair1Mass", "M_{1} (4K hypothesis) [GeV/c^{2}]", {500, .7, 2.2}});
    test->AddObservable({"Pair2Mass", "M_{2} (4K hypothesis) [GeV/c^{2}]", {500, .7, 2.2}});
    test->AddObservable({"EnergyResidual", "#Sigma E_{K} - E_{beam} - m_{p} [GeV]", {500, -.5, 2.5}});
  }
  std::vector<KFPbarBGResult> signalResults, backgroundResults;
  signalResults.reserve(count); backgroundResults.reserve(count);
  KFPbarBalanceSummary signalSummary, backgroundSummary;
  // These are ALL disjoint pairings of four tracks; no generated-order or
  // charge-pairing restriction is imposed on the requested BG preselection.
  const int pairings[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
  // Restrict the TRUE intermediate pion-pair masses before generating angles.
  // For identical daughters with their momenta held fixed,
  // M_KK^2 >= M_pipi^2 + 4*(mK^2-mPi^2). Thus any KK pair below the upper
  // window edge MUST have M_pipi below the bound used here, in any frame.
  const double pairUpper = mPhi + massWindow;
  double bgMassLow = std::max(2. * mPi, mPhi - massWindow);
  double bgMassHigh = pairUpper;
  if (useKaonMassWindow) {
    bgMassLow = 2. * mPi;
    if (pairUpper <= 2. * mK)
      throw std::invalid_argument("The kaon-mass pair window lies below the KK threshold");
    bgMassHigh = std::sqrt(pairUpper * pairUpper - 4. * (mK * mK - mPi * mPi));
  }
  bgMassHigh = std::min(bgMassHigh, sqrtS - bgMassLow);
  if (bgMassHigh <= bgMassLow || sqrtS <= 2. * bgMassLow)
    throw std::invalid_argument("No allowed intermediate masses for the BG selection");
  // For uniform dM1 dM2 and isotropic decays, dPhi4 is proportional to
  // p*(sqrtS,M1,M2) * p*(M1,mpi,mpi) * p*(M2,mpi,mpi).
  // Bound each factor separately to obtain a safe rejection envelope.
  // Reference: PDG Kinematics, recursive phase-space factorization.
  const double bgProductionMax = .5 * std::sqrt(sqrtS * sqrtS - 4. * bgMassLow * bgMassLow);
  const double bgDaughterMax = .5 * std::sqrt(std::max(0., bgMassHigh * bgMassHigh - 4. * mPi * mPi));
  const double bgWeightMax = bgProductionMax * bgDaughterMax * bgDaughterMax;
  if (!std::isfinite(bgWeightMax) || bgWeightMax <= 0.)
    throw std::invalid_argument("Invalid sequential BG phase-space envelope");
  double pionPairMasses[2] = {mPi, mPi};
  TGenPhaseSpace evSignal, evPhi1, evPhi2, evBackground, evBG1, evBG2;
  unsigned long long bgProposals = 0, bgPhaseSpaceEvents = 0, bgWindowCandidates = 0;
  for (int event = 0; event < count; ++event) {
    for (int channel = 0; channel < 2; ++channel) {
      const bool isSignal = channel == 0;
      KFTestDiagnostics& test = isSignal ? signal : background;
      auto& results = isSignal ? signalResults : backgroundResults;
      auto& summary = isSignal ? signalSummary : backgroundSummary;
      TLorentzVector tracks[4];
      if (isSignal) {
        // 1a. Signal: finite-width phi phi production followed by two 2K decays.
        double phiMass1 = mPhi, phiMass2 = mPhi;
        if (WidthPhi > 0.) {
          const double minimum = 2. * mK, s = sqrtS * sqrtS;
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
        double phiMasses[2] = {phiMass1, phiMass2};
        double kaonMasses[2] = {mK, mK};
        if (!evSignal.SetDecay(initialSystem, 2, phiMasses))
          throw std::runtime_error("Kinematically forbidden phi phi production");
        evSignal.Generate();
        TLorentzVector phi1 = *evSignal.GetDecay(0), phi2 = *evSignal.GetDecay(1);
        if (!evPhi1.SetDecay(phi1, 2, kaonMasses) || !evPhi2.SetDecay(phi2, 2, kaonMasses))
          throw std::runtime_error("Kinematically forbidden phi decay");
        evPhi1.Generate(); evPhi2.Generate();
        tracks[0] = *evPhi1.GetDecay(0); tracks[1] = *evPhi1.GetDecay(1);
        tracks[2] = *evPhi2.GetDecay(0); tracks[3] = *evPhi2.GetDecay(1);
      } else {
        // 1b. BG: pbar p -> M1 M2, M1 -> pi pi, M2 -> pi pi.
        // This is a phase-space factorization, not a new resonance model.
        bool selected = false;
        for (int attempt = 0; attempt < 1000000; ++attempt) {
          ++bgProposals;
          const double mass1 = bgMassLow + (bgMassHigh - bgMassLow) * gRandom->Rndm();
          const double mass2 = bgMassLow + (bgMassHigh - bgMassLow) * gRandom->Rndm();
          if (mass1 + mass2 >= sqrtS) continue;
          const double sum = mass1 + mass2, difference = mass1 - mass2;
          const double pProduction = std::sqrt((sqrtS * sqrtS - sum * sum) *
                                               (sqrtS * sqrtS - difference * difference)) / (2. * sqrtS);
          const double pDecay1 = .5 * std::sqrt(std::max(0., mass1 * mass1 - 4. * mPi * mPi));
          const double pDecay2 = .5 * std::sqrt(std::max(0., mass2 * mass2 - 4. * mPi * mPi));
          const double weight = pProduction * pDecay1 * pDecay2;
          if (gRandom->Rndm() * bgWeightMax >= weight) continue;

          double intermediateMasses[2] = {mass1, mass2};
          if (!evBackground.SetDecay(initialSystem, 2, intermediateMasses))
            throw std::runtime_error("Kinematically forbidden BG intermediate production");
          evBackground.Generate();
          TLorentzVector parent1 = *evBackground.GetDecay(0), parent2 = *evBackground.GetDecay(1);
          if (!evBG1.SetDecay(parent1, 2, pionPairMasses) ||
              !evBG2.SetDecay(parent2, 2, pionPairMasses))
            throw std::runtime_error("Kinematically forbidden BG pion-pair decay");
          evBG1.Generate(); evBG2.Generate();
          ++bgPhaseSpaceEvents;
          // Choose the generating partition uniformly among all three choices.
          const auto& generatedPairing = pairings[int(3. * gRandom->Rndm())];
          tracks[generatedPairing[0]] = *evBG1.GetDecay(0);
          tracks[generatedPairing[1]] = *evBG1.GetDecay(1);
          tracks[generatedPairing[2]] = *evBG2.GetDecay(0);
          tracks[generatedPairing[3]] = *evBG2.GetDecay(1);

          // Test the requested window under the selected mass hypothesis.
          TLorentzVector massHypothesis[4];
          for (int track = 0; track < 4; ++track) {
            const auto& pion = tracks[track];
            massHypothesis[track] = useKaonMassWindow
                ? TLorentzVector(pion.Px(), pion.Py(), pion.Pz(), std::hypot(pion.P(), mK))
                : pion;
          }
          int proposalMultiplicity = 0, passingPairings = 0;
          const double tolerance = 1.e-10 * std::max(1., sqrtS);
          for (const auto& pairing : pairings) {
            const double pionMass1 = (tracks[pairing[0]] + tracks[pairing[1]]).M();
            const double pionMass2 = (tracks[pairing[2]] + tracks[pairing[3]]).M();
            if (pionMass1 >= bgMassLow - tolerance && pionMass1 <= bgMassHigh + tolerance &&
                pionMass2 >= bgMassLow - tolerance && pionMass2 <= bgMassHigh + tolerance)
              ++proposalMultiplicity;
            const double assignedMass1 = (massHypothesis[pairing[0]] + massHypothesis[pairing[1]]).M();
            const double assignedMass2 = (massHypothesis[pairing[2]] + massHypothesis[pairing[3]]).M();
            if (std::abs(assignedMass1 - mPhi) < massWindow &&
                std::abs(assignedMass2 - mPhi) < massWindow)
              ++passingPairings;
          }
          if (proposalMultiplicity == 0)
            throw std::runtime_error("Sequential BG event outside its proposal domain");
          if (passingPairings == 0) continue;
          ++bgWindowCandidates;
          // The three generating partitions overlap. Divide by the number of
          // proposal domains containing this event, NOT by passingPairings:
          // the narrow pion-mass domain and the KK-mass cut differ. This keeps
          // the original flat 4-body phase space conditioned on the final cut.
          if (gRandom->Rndm() * proposalMultiplicity >= 1.) continue;
          selected = true; break;
        }
        if (!selected) throw std::runtime_error("Unable to generate a sequential BG event passing both pair-mass windows");
      }

      // 2. Smear MOMENTUM measurements, independent of particle identification.
      // Order: beam, track1, track2, track3, track4; (p,theta,phi) for each.
      const TMatrixD truth = KFSphericalParameters({beam, tracks[0], tracks[1], tracks[2], tracks[3]});
      const KFTestSample sample = test.Smear(truth);
      const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 0, mP);
      // Assign KAON masses even when the true generated particles are pions.
      const TLorentzVector k1 = KFFourVector(sample.measurements, 3, mK);
      const TLorentzVector k2 = KFFourVector(sample.measurements, 6, mK);
      const TLorentzVector k3 = KFFourVector(sample.measurements, 9, mK);
      const TLorentzVector k4 = KFFourVector(sample.measurements, 12, mK);
      const std::vector<TLorentzVector> measuredFinal = {k1, k2, k3, k4};
      KFPbarBGResult result;
      result.energyResidual = (k1 + k2 + k3 + k4 - measuredBeam - target).E();
      // For the mass plots, choose the pairing closest to two phi masses under
      // the KAON hypothesis. This does not impose another event selection.
      int bestPairing = 0;
      double bestDistance = std::numeric_limits<double>::infinity();
      for (int p = 0; p < 3; ++p) {
        const auto& pairing = pairings[p];
        const double mass1 = (measuredFinal[pairing[0]] + measuredFinal[pairing[1]]).M();
        const double mass2 = (measuredFinal[pairing[2]] + measuredFinal[pairing[3]]).M();
        const double distance = std::max(std::abs(mass1 - mPhi), std::abs(mass2 - mPhi));
        if (distance < bestDistance) {
          bestDistance = distance; bestPairing = p;
          result.pairMass1 = mass1; result.pairMass2 = mass2;
        }
      }

      // 3. Same four-kaon 4C fit in both channels, with no phi mass constraints.
      BalanceFitter fitter(measuredBeam, target, measuredFinal);
      if (fitter.GetNDF() != 4)
        throw std::runtime_error("BalanceFitter must have four constraints and no unknowns");
      fitter.SetChi2DifCut(1.e-10);
      std::vector<double> variance = sample.variances;
      fitter.SetVariance(variance.data());
      const double chi2 = fitter.DoKinematicFit();
      result.chi2 = chi2; result.pvalue = fitter.GetPValue();
      const bool finite = test.FillHist(fitter, chi2, truth);
      if (finite) {
        // 4. Extract the fitted measurements/four-vectors explicitly.
        const TMatrixD best = fitter.GetBestMeasurements();
        const auto fittedLV = fitter.GetFittedLV(); // beam, track1, track2, track3, track4
        KFCheckColumn(best, 15, "Fitted balance measurements");
        if (fittedLV.size() != 5) throw std::runtime_error("Unexpected fitted four-vector count");
        const TLorentzVector fittedBeam = fittedLV[0];
        const TLorentzVector fittedK1 = fittedLV[1], fittedK2 = fittedLV[2];
        const TLorentzVector fittedK3 = fittedLV[3], fittedK4 = fittedLV[4];
        const std::vector<TLorentzVector> fittedFinal = {fittedK1, fittedK2, fittedK3, fittedK4};
        const auto previousNonconserving = summary.nonconserving;
        KFCheckBalanceConservation(fittedBeam, target, fittedFinal, summary);
        result.valid = summary.nonconserving == previousNonconserving;
        const auto& pairing = pairings[bestPairing];
        test.FillObservable("Pair1Mass", result.pairMass1,
                             (fittedFinal[pairing[0]] + fittedFinal[pairing[1]]).M());
        test.FillObservable("Pair2Mass", result.pairMass2,
                             (fittedFinal[pairing[2]] + fittedFinal[pairing[3]]).M());
        test.FillObservable("EnergyResidual", result.energyResidual,
                             (fittedK1 + fittedK2 + fittedK3 + fittedK4 - fittedBeam - target).E());
      } else ++summary.invalid;
      results.push_back(result); // Preserve failed fits in the cut denominators.
    }
  }
  const std::vector<KFTestDiagnostics::CanvasPage> measuredPages = {
      {3, 1, {0, 1, 2}}, {4, 3, {3, 6, 9, 12, 4, 7, 10, 13, 5, 8, 11, 14}}};
  gSystem->mkdir(directory, true);
  signal.DrawHist(directory, {}, measuredPages);
  background.DrawHist(directory, {}, measuredPages);
  signal.WriteHist(directory + "/signal.root");
  background.WriteHist(directory + "/background.root");
  std::cout << "Sequential BG mass range: [" << bgMassLow << ", " << bgMassHigh << "] GeV/c^2\n"
            << "BG mass proposals: " << bgProposals << ", generated sequential events: " << bgPhaseSpaceEvents
            << ", window candidates: " << bgWindowCandidates << ", accepted events: " << count
            << ", proposal efficiency: " << double(count) / bgProposals << '\n'
            << "Selection: exists a disjoint pairing with both |M_pair - mPhi| < " << massWindow
            << " GeV/c^2, " << (useKaonMassWindow ? "kaon-assigned" : "true pion")
            << " masses before smearing.\n"
            << "Proposal efficiency is NOT the selection acceptance of unrestricted 4pi phase space.\n"
            << "Signal invalid/nonconserving fits: " << signalSummary.invalid << "/" << signalSummary.nonconserving
            << "; BG invalid/nonconserving fits: " << backgroundSummary.invalid << "/" << backgroundSummary.nonconserving << '\n';
  KFDrawPbarBGSeparation(signalResults, backgroundResults, directory, signalFraction);
}

// Histogram-only comparison. Reactions, generated four-vectors, smearing and
// BalanceFitter calls are explicit in the two Test... functions below.
class KFPbarMomentumComparison {
public:
  KFPbarMomentumComparison(const TString& prefix, const TString& exclusiveLabel,
                           const TString& missingLabel, double momentumRange)
      : prefix_(prefix), labels_{exclusiveLabel, missingLabel} {
    const char* axes[7] = {
        ";#Delta p_{x} [GeV/c];Fraction / bin", ";#Delta p_{y} [GeV/c];Fraction / bin",
        ";#Delta p_{z} [GeV/c];Fraction / bin", ";|#Delta #vec{p}| [GeV/c];Fraction / bin",
        ";|#Delta #vec{p}_{T}| [GeV/c];Fraction / bin", ";#Delta E [GeV];Fraction / bin",
        ";M_{miss}^{2} [GeV^{2}/c^{4}];Fraction / bin"};
    const double mass2Range = std::max(.2, .2 * momentumRange * momentumRange);
    for (int stage = 0; stage < 3; ++stage)
      for (int channel = 0; channel < 2; ++channel)
        for (int variable = 0; variable < 7; ++variable) {
          const double high = variable == 6 ? mass2Range : momentumRange;
          const double low = variable == 3 || variable == 4 ? 0. : -high;
          const TString name = TString::Format("%s_S%d_C%d_V%d", prefix_.Data(), stage, channel, variable);
          std::unique_ptr<TH1D> h(new TH1D(name, axes[variable], 400, low, high));
          h->SetDirectory(nullptr); h->Sumw2();
          hist_[stage][channel].push_back(std::move(h));
        }
    for (int channel = 0; channel < 2; ++channel) {
      chi2_[channel].reset(new TH1D(TString::Format("%s_Chi2_C%d", prefix_.Data(), channel),
                                    ";#chi^{2} (exclusive 4C hypothesis);Fraction / bin", 1000, 0., 1000.));
      pvalue_[channel].reset(new TH1D(TString::Format("%s_PValue_C%d", prefix_.Data(), channel),
                                      ";P-value (exclusive 4C hypothesis);Fraction / bin", 200, 0., 1.));
      chi2_[channel]->SetDirectory(nullptr); pvalue_[channel]->SetDirectory(nullptr);
      chi2_[channel]->Sumw2(); pvalue_[channel]->Sumw2();
    }
  }

  void FillBeforeFit(int channel, const TLorentzVector& generatedBalance,
                     const TLorentzVector& measuredBalance, double weight) {
    sumWeights_[channel] += weight;
    FillBalance(0, channel, generatedBalance, weight);
    FillBalance(1, channel, measuredBalance, weight);
  }
  void FillAfterFit(int channel, const TLorentzVector& fittedBalance,
                    double chi2, double pvalue, bool conserved, double weight) {
    FillBalance(2, channel, fittedBalance, weight);
    finiteWeights_[channel] += weight;
    if (!conserved || !std::isfinite(pvalue) || pvalue < 0. || pvalue > 1.) return;
    usableWeights_[channel] += weight;
    chi2_[channel]->Fill(chi2, weight);
    pvalue_[channel]->Fill(pvalue == 1. ? std::nextafter(1., 0.) : pvalue, weight);
    if (pvalue >= .05) retainedWeights_[channel] += weight;
  }
  void DrawAndWrite(const TString& directory) const {
    gSystem->mkdir(directory, true);
    const char* names[3] = {"truth_balance", "momentum_balance", "fitted_balance"};
    const char* titles[3] = {"Generated charged-track balance", "Smeared charged-track balance",
                             "Fitted balance under the exclusive 4C hypothesis"};
    std::vector<TCanvas*> canvases;
    for (int stage = 0; stage < 3; ++stage) {
      auto* canvas = new TCanvas(prefix_ + "_" + names[stage], titles[stage], 1600, 850);
      canvas->Divide(4, 2);
      for (int variable = 0; variable < 7; ++variable) {
        canvas->cd(variable + 1);
        Overlay(*hist_[stage][0][variable], *hist_[stage][1][variable]);
      }
      canvas->SaveAs(directory + "/" + names[stage] + ".pdf");
      canvases.push_back(canvas);
    }
    auto* fitCanvas = new TCanvas(prefix_ + "_FitComparison", "Exclusive 4C fit comparison", 1200, 600);
    fitCanvas->Divide(2, 1);
    fitCanvas->cd(1); Overlay(*chi2_[0], *chi2_[1]);
    fitCanvas->cd(2); Overlay(*pvalue_[0], *pvalue_[1]);
    fitCanvas->SaveAs(directory + "/fit_comparison.pdf"); canvases.push_back(fitCanvas);
    TFile output(directory + "/momentum_balance.root", "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create momentum-balance ROOT file");
    for (int stage = 0; stage < 3; ++stage)
      for (int channel = 0; channel < 2; ++channel)
        for (const auto& hist : hist_[stage][channel]) hist->Write();
    for (int channel = 0; channel < 2; ++channel) { chi2_[channel]->Write(); pvalue_[channel]->Write(); }
    for (auto* canvas : canvases) canvas->Write();
    std::cout << "Balance convention: charged final state minus measured beam and fixed target.\n"
              << "The pi0 is omitted from reconstruction. Generated balance is zero for the exclusive channel\n"
              << "and -p_pi0 for the missing-pi0 channel; true missing-mass squared is 0 or m_pi0^2.\n"
              << "Plots are normalized per generated sample using phase-space weights; no production ratio is assumed.\n";
    for (int channel = 0; channel < 2; ++channel) {
      const double denominator = sumWeights_[channel];
      if (denominator <= 0.) continue;
      std::cout << labels_[channel] << ": finite fits = " << finiteWeights_[channel] / denominator
                << ", usable conserving fits = " << usableWeights_[channel] / denominator
                << ", retained at P >= 0.05 = " << retainedWeights_[channel] / denominator << '\n';
    }
  }
private:
  TString prefix_, labels_[2];
  std::vector<std::unique_ptr<TH1D>> hist_[3][2];
  std::unique_ptr<TH1D> chi2_[2], pvalue_[2];
  double sumWeights_[2] = {0., 0.}, finiteWeights_[2] = {0., 0.};
  double usableWeights_[2] = {0., 0.}, retainedWeights_[2] = {0., 0.};
  void FillBalance(int stage, int channel, const TLorentzVector& delta, double weight) {
    const double values[7] = {delta.Px(), delta.Py(), delta.Pz(), delta.P(),
                              std::hypot(delta.Px(), delta.Py()), delta.E(), delta.M2()};
    for (int i = 0; i < 7; ++i) hist_[stage][channel][i]->Fill(values[i], weight);
  }
  void Overlay(const TH1D& exclusive, const TH1D& missing) const {
    TH1D first(exclusive), second(missing);
    first.SetDirectory(nullptr); second.SetDirectory(nullptr);
    if (sumWeights_[0] > 0.) first.Scale(1. / sumWeights_[0]);
    if (sumWeights_[1] > 0.) second.Scale(1. / sumWeights_[1]);
    first.SetLineColor(kBlue); second.SetLineColor(kRed);
    first.SetMinimum(1.e-7);
    first.SetMaximum(2. * std::max(1.e-6, std::max(first.GetMaximum(), second.GetMaximum())));
    gPad->SetLogy();
    auto* a = first.DrawCopy("hist"); auto* b = second.DrawCopy("hist same");
    auto* legend = new TLegend(.42, .73, .92, .91);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(a, labels_[0], "l"); legend->AddEntry(b, labels_[1], "l"); legend->Draw();
  }
};

// Runtime histogram owner only; no ROOT I/O dictionary is needed.
#ifdef __ROOTCLING__
#pragma link off class KFPbarMomentumComparison;
#endif

// Compare pbar p -> pi- pi+ pi- pi+ with pbar p -> pi- pi+ pi- pi+ pi0.
// Pion masses throughout: no phi production, mass-window selection, or kaon hypothesis.
void TestPbar4PiPi0Fit(int nEvents = -1, double beamMomentum = -1.) {
  const int count = nEvents < 0 ? nev : nEvents;
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (count <= 0 || !std::isfinite(momentum) || momentum <= 0.)
    throw std::invalid_argument("Invalid TestPbar4PiPi0Fit arguments");
  TVector3 direction(.1, .05, 1.); direction.SetMag(momentum);
  const TLorentzVector beam(direction, std::hypot(momentum, mP));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector initialSystem = beam + target;
  const double sqrtS = initialSystem.M();
  const auto variables = KFJoinVariables({
      KFParticleVariables("beam", "#bar{p}_{beam}", .01, ResPbarBeam, ResThPbarBeam, ResPhPbarBeam, 1., ResPbarBeamEstm),
      KFParticleVariables("pi1", "#pi^{-}_{1}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("pi2", "#pi^{+}_{1}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("pi3", "#pi^{-}_{2}", .15, ResP, ResThP, ResPhP, 1., ResPEstm),
      KFParticleVariables("pi4", "#pi^{+}_{2}", .15, ResP, ResThP, ResPhP, 1., ResPEstm)});
  double finalMasses[5] = {mPi, mPi, mPi, mPi, mPi0};
  const int nCharged = 4;
  KFTestDiagnostics exclusive("FourPiExclusive", variables, {}, KFPbarBalanceConstraints());
  KFTestDiagnostics missing("FourPiMissingPi0", variables, {}, KFPbarBalanceConstraints());
  KFPbarMomentumComparison plots("Pbar4PiPi0", "4#pi", "4#pi + #pi^{0}", std::max(1., momentum + .5 * sqrtS));
  TGenPhaseSpace generators[2];
  // Check both thresholds before processing any events. No phi-related cuts.
  if (!generators[0].SetDecay(initialSystem, nCharged, finalMasses) ||
      !generators[1].SetDecay(initialSystem, nCharged + 1, finalMasses))
    throw std::invalid_argument("Exclusive or missing-pi0 reaction is kinematically forbidden");
  KFPbarBalanceSummary summaries[2];
  for (int event = 0; event < count; ++event) {
    for (int channel = 0; channel < 2; ++channel) {
      KFTestDiagnostics& test = channel == 0 ? exclusive : missing;
      auto& generator = generators[channel];
      auto& summary = summaries[channel];
      // 1. Generate directly and KEEP the phase-space weight. This avoids an
      // unweighting/rejection loop for the multi-body comparison samples.
      const double weight = generator.Generate();
      if (!std::isfinite(weight) || weight < 0.)
        throw std::runtime_error("Invalid phase-space weight");
      if (weight == 0.) continue;
      const TLorentzVector truePi1 = *generator.GetDecay(0);
      const TLorentzVector truePi2 = *generator.GetDecay(1);
      const TLorentzVector truePi3 = *generator.GetDecay(2);
      const TLorentzVector truePi4 = *generator.GetDecay(3);
      const TLorentzVector generatedBalance = truePi1 + truePi2 + truePi3 + truePi4 - beam - target;
      const TMatrixD truth = KFSphericalParameters({beam, truePi1, truePi2, truePi3, truePi4});
      const KFTestSample sample = test.Smear(truth);
      const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 0, mP);
      const TLorentzVector pi1 = KFFourVector(sample.measurements, 3, mPi);
      const TLorentzVector pi2 = KFFourVector(sample.measurements, 6, mPi);
      const TLorentzVector pi3 = KFFourVector(sample.measurements, 9, mPi);
      const TLorentzVector pi4 = KFFourVector(sample.measurements, 12, mPi);
      const std::vector<TLorentzVector> measuredFinal = {pi1, pi2, pi3, pi4};
      const TLorentzVector measuredBalance = pi1 + pi2 + pi3 + pi4 - measuredBeam - target;
      // generator.GetDecay(nCharged) is the undetected pi0 in channel 1.
      // It is deliberately absent from the measurement vector and fitter.
      if (channel == 1) {
        const TLorentzVector truePi0 = *generator.GetDecay(nCharged);
        const TLorentzVector closure = generatedBalance + truePi0;
        if (closure.P() > 1.e-8 || std::abs(closure.E()) > 1.e-8)
          throw std::runtime_error("Generated missing-pi0 event violates four-momentum conservation");
      } else if (generatedBalance.P() > 1.e-8 || std::abs(generatedBalance.E()) > 1.e-8)
        throw std::runtime_error("Generated exclusive event violates four-momentum conservation");
      plots.FillBeforeFit(channel, generatedBalance, measuredBalance, weight);

      // 2. Fit the SAME exclusive charged-particle 4C hypothesis in both samples.
      BalanceFitter fitter(measuredBeam, target, measuredFinal);
      if (fitter.GetNDF() != 4) throw std::runtime_error("Expected a four-constraint balance fit");
      fitter.SetChi2DifCut(1.e-10);
      std::vector<double> variance = sample.variances;
      fitter.SetVariance(variance.data());
      const double chi2 = fitter.DoKinematicFit();
      const double pvalue = fitter.GetPValue();
      // Weighted histograms; no unknown parameters are introduced.
      if (!test.FillHist(fitter, chi2, truth, TMatrixD(), weight)) { ++summary.invalid; continue; }

      // 3. Extract fitted measurements and four-vectors explicitly.
      const TMatrixD best = fitter.GetBestMeasurements();
      const auto fittedLV = fitter.GetFittedLV();
      KFCheckColumn(best, 15, "Fitted four-pion measurements");
      if (fittedLV.size() != 5) throw std::runtime_error("Unexpected fitted four-vector count");
      const TLorentzVector fittedBeam = fittedLV[0];
      const TLorentzVector fittedPi1 = fittedLV[1], fittedPi2 = fittedLV[2];
      const TLorentzVector fittedPi3 = fittedLV[3], fittedPi4 = fittedLV[4];
      const std::vector<TLorentzVector> fittedFinal = {fittedPi1, fittedPi2, fittedPi3, fittedPi4};
      const TLorentzVector fittedBalance = fittedPi1 + fittedPi2 + fittedPi3 + fittedPi4 - fittedBeam - target;
      const auto before = summary.nonconserving;
      KFCheckBalanceConservation(fittedBeam, target, fittedFinal, summary);
      const bool conserved = summary.nonconserving == before;
      plots.FillAfterFit(channel, fittedBalance, chi2, pvalue, conserved, weight);
    }
  }
  const TString directory = "figs/Pbar4PiPi0Fit";
  const std::vector<KFTestDiagnostics::CanvasPage> measuredPages = {
      {3, 1, {0, 1, 2}}, {4, 3, {3, 6, 9, 12, 4, 7, 10, 13, 5, 8, 11, 14}}};
  exclusive.DrawHist(directory, {}, measuredPages);
  missing.DrawHist(directory, {}, measuredPages);
  exclusive.WriteHist(directory + "/exclusive.root");
  missing.WriteHist(directory + "/missing_pi0.root");
  plots.DrawAndWrite(directory);
}

// Compare pbar-p elastic scattering with nonresonant pbar p -> pbar p pi0.
// The pi0 is present in generation only; both channels are fitted as exclusive pbar-p.
void TestPbarElasticPi0Fit(int nEvents = -1, double beamMomentum = -1.) {
  const int count = nEvents < 0 ? nev : nEvents;
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (count <= 0 || !std::isfinite(momentum) || momentum <= 0.)
    throw std::invalid_argument("Invalid TestPbarElasticPi0Fit arguments");
  TVector3 direction(.1, .05, 1.); direction.SetMag(momentum);
  const TLorentzVector beam(direction, std::hypot(momentum, mP));
  const TLorentzVector target(0., 0., 0., mP);
  TLorentzVector initialSystem = beam + target;
  const double sqrtS = initialSystem.M();
  const auto variables = KFPbarBalanceMeasurements(false);
  double finalMasses[3] = {mP, mP, mPi0};
  const int nCharged = 2;
  KFTestDiagnostics exclusive("PbarElasticExclusive", variables, {}, KFPbarBalanceConstraints());
  KFTestDiagnostics missing("PbarElasticMissingPi0", variables, {}, KFPbarBalanceConstraints());
  KFPbarMomentumComparison plots("PbarElasticPi0", "#bar{p}p elastic", "#bar{p}p + #pi^{0}",
                                 std::max(1., momentum + .5 * sqrtS));
  TGenPhaseSpace generators[2];
  // Check both thresholds before processing any events. No phi-related cuts.
  if (!generators[0].SetDecay(initialSystem, nCharged, finalMasses) ||
      !generators[1].SetDecay(initialSystem, nCharged + 1, finalMasses))
    throw std::invalid_argument("Exclusive or missing-pi0 reaction is kinematically forbidden");
  KFPbarBalanceSummary summaries[2];
  for (int event = 0; event < count; ++event) {
    for (int channel = 0; channel < 2; ++channel) {
      KFTestDiagnostics& test = channel == 0 ? exclusive : missing;
      auto& generator = generators[channel];
      auto& summary = summaries[channel];
      // 1. Generate directly and KEEP the phase-space weight. This avoids an
      // unweighting/rejection loop for the multi-body comparison samples.
      const double weight = generator.Generate();
      if (!std::isfinite(weight) || weight < 0.)
        throw std::runtime_error("Invalid phase-space weight");
      if (weight == 0.) continue;
      const TLorentzVector truePbar = *generator.GetDecay(0);
      const TLorentzVector trueProton = *generator.GetDecay(1);
      const TLorentzVector generatedBalance = truePbar + trueProton - beam - target;
      const TMatrixD truth = KFSphericalParameters({beam, truePbar, trueProton});
      const KFTestSample sample = test.Smear(truth);
      const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 0, mP);
      const TLorentzVector pbar = KFFourVector(sample.measurements, 3, mP);
      const TLorentzVector proton = KFFourVector(sample.measurements, 6, mP);
      const std::vector<TLorentzVector> measuredFinal = {pbar, proton};
      const TLorentzVector measuredBalance = pbar + proton - measuredBeam - target;
      // generator.GetDecay(nCharged) is the undetected pi0 in channel 1.
      // It is deliberately absent from the measurement vector and fitter.
      if (channel == 1) {
        const TLorentzVector truePi0 = *generator.GetDecay(nCharged);
        const TLorentzVector closure = generatedBalance + truePi0;
        if (closure.P() > 1.e-8 || std::abs(closure.E()) > 1.e-8)
          throw std::runtime_error("Generated missing-pi0 event violates four-momentum conservation");
      } else if (generatedBalance.P() > 1.e-8 || std::abs(generatedBalance.E()) > 1.e-8)
        throw std::runtime_error("Generated exclusive event violates four-momentum conservation");
      plots.FillBeforeFit(channel, generatedBalance, measuredBalance, weight);

      // 2. Fit the SAME exclusive charged-particle 4C hypothesis in both samples.
      BalanceFitter fitter(measuredBeam, target, measuredFinal);
      if (fitter.GetNDF() != 4) throw std::runtime_error("Expected a four-constraint balance fit");
      fitter.SetChi2DifCut(1.e-10);
      std::vector<double> variance = sample.variances;
      fitter.SetVariance(variance.data());
      const double chi2 = fitter.DoKinematicFit();
      const double pvalue = fitter.GetPValue();
      // Weighted histograms; no unknown parameters are introduced.
      if (!test.FillHist(fitter, chi2, truth, TMatrixD(), weight)) { ++summary.invalid; continue; }

      // 3. Extract fitted measurements and four-vectors explicitly.
      const TMatrixD best = fitter.GetBestMeasurements();
      const auto fittedLV = fitter.GetFittedLV();
      KFCheckColumn(best, 9, "Fitted pbar-p measurements");
      if (fittedLV.size() != 3) throw std::runtime_error("Unexpected fitted four-vector count");
      const TLorentzVector fittedBeam = fittedLV[0];
      const TLorentzVector fittedPbar = fittedLV[1], fittedProton = fittedLV[2];
      const std::vector<TLorentzVector> fittedFinal = {fittedPbar, fittedProton};
      const TLorentzVector fittedBalance = fittedPbar + fittedProton - fittedBeam - target;
      const auto before = summary.nonconserving;
      KFCheckBalanceConservation(fittedBeam, target, fittedFinal, summary);
      const bool conserved = summary.nonconserving == before;
      plots.FillAfterFit(channel, fittedBalance, chi2, pvalue, conserved, weight);
    }
  }
  const TString directory = "figs/PbarElasticPi0Fit";
  const std::vector<KFTestDiagnostics::CanvasPage> measuredPages = {
      {3, 3, {0, 1, 2, 3, 4, 5, 6, 7, 8}}};
  exclusive.DrawHist(directory, {}, measuredPages);
  missing.DrawHist(directory, {}, measuredPages);
  exclusive.WriteHist(directory + "/exclusive.root");
  missing.WriteHist(directory + "/missing_pi0.root");
  plots.DrawAndWrite(directory);
}

// Both tests default to the 1 GeV/c pbar beam configured by pPbar.
void TestBalanceFit(int nEvents = -1, double beamMomentum = -1.) {
  TestPbarElasticFit(nEvents, beamMomentum);
  TestPbarPhiPhiFit(nEvents, beamMomentum);
}

// Iterative momentum-resolution calibration. The true variables are copied
// once, and generated/smeared events are cached before any covariance update.
// The first fit uses initialScale * true sigma, hence one quarter of the true
// momentum variance for the default initialScale=0.5. Angular errors stay fixed.
struct KFResolutionGroup {
  TString name, title;
  double trueResolution;
  std::vector<int> rows;
};

struct KFResolutionEvent {
  TMatrixD truth;
  KFTestSample sample;
};

struct KFResolutionIteration {
  int step = 0;
  unsigned long finiteEvents = 0, failedEvents = 0;
  std::vector<double> estimates, means, widths, widthErrors;
  std::vector<bool> valid, stable;
};

struct KFResolutionTuningResult {
  bool converged = false;
  std::vector<KFResolutionGroup> groups;
  std::vector<double> estimates;
  std::vector<KFResolutionIteration> iterations;
};

// TGenPhaseSpace uses gRandom. Use a private generator and restore the previous
// generator even if a forbidden decay or an I/O error raises an exception.
class KFResolutionRandomScope {
  TRandom* previous_;
public:
  explicit KFResolutionRandomScope(TRandom& random) : previous_(gRandom) { gRandom = &random; }
  ~KFResolutionRandomScope() { gRandom = previous_; }
  KFResolutionRandomScope(const KFResolutionRandomScope&) = delete;
  KFResolutionRandomScope& operator=(const KFResolutionRandomScope&) = delete;
};
#ifdef __ROOTCLING__
#pragma link off class KFResolutionRandomScope;
#endif

int KFCheckResolutionTuning(int nEvents, int maxIterations, double initialScale) {
  const int count = nEvents < 0 ? nev : nEvents;
  if (count <= 0 || maxIterations <= 0 ||
      !std::isfinite(initialScale) || initialScale <= 0.)
    throw std::invalid_argument("Resolution tuning requires positive events, iterations and initial scale");
  return count;
}

// Reset estimates to the requested initial guess without altering true sigmas.
void KFInitialTuningVariables(std::vector<KFTestVariable>& variables,
                              const std::vector<KFResolutionGroup>& groups,
                              double initialScale) {
  for (const auto& group : groups) {
    if (!std::isfinite(group.trueResolution) || group.trueResolution <= 0.)
      throw std::invalid_argument("True momentum resolutions must be positive");
    for (int row : group.rows) {
      if (row < 0 || row >= static_cast<int>(variables.size()))
        throw std::invalid_argument("Invalid momentum row for resolution tuning");
      variables[row].resolutionEstm = KFRelativeResolution(initialScale * group.trueResolution);
    }
  }
}

// The leading entries are measured correction pulls and control momentum tuning.
// Optional appended unknown truth pulls are diagnostic only. BalanceFitter uses
// an empty unknown list and never queries unknown getters.
template<class Fitter>
std::vector<double> KFFitResolutionSample(Fitter& fitter, std::vector<double>& variances,
    const TMatrixD& trueUnknown = TMatrixD(),
    const std::vector<KFTestVariable>& unknownVariables = {}) {
  fitter.SetVariance(variances.data());
  fitter.UpdateVariance(false);
  fitter.SetChi2DifCut(1.e-10);
  const double chi2 = fitter.DoKinematicFit();
  if (!std::isfinite(chi2) || chi2 < 0.) return {};
  const TMatrixD best = fitter.GetBestMeasurements();
  for (int row = 0; row < best.GetNrows(); ++row)
    if (!std::isfinite(best(row, 0))) return {};
  auto pulls = fitter.GetPull();
  if (pulls.size() != variances.size()) return {};
  if (!unknownVariables.empty()) {
    const std::size_t first = pulls.size();
    pulls.resize(first + unknownVariables.size(), std::numeric_limits<double>::quiet_NaN());
    // A diagnostic getter must not remove an otherwise valid measured fit
    // from the calibration sample if an unknown covariance is unavailable.
    try {
      const TMatrixD bestUnknown = fitter.GetBestUnknowns();
      const TMatrixD covariance = fitter.GetUnmeasuredCovariance();
      KFCheckColumn(bestUnknown, unknownVariables.size(), "Fitted unknowns for diagnostic pulls");
      KFCheckColumn(trueUnknown, unknownVariables.size(), "True unknowns for diagnostic pulls");
      if (covariance.GetNrows() == int(unknownVariables.size()) &&
          covariance.GetNcols() == int(unknownVariables.size()))
        for (std::size_t row = 0; row < unknownVariables.size(); ++row)
          if (std::isfinite(covariance(row, row)) && covariance(row, row) > 0.)
            pulls[first + row] = unknownVariables[row].Difference(bestUnknown(row, 0), trueUnknown(row, 0)) /
                                 std::sqrt(covariance(row, row));
    } catch (const std::exception&) { /* Keep undefined diagnostic pulls as NaN. */ }
  }
  return pulls;
}

TCanvas* KFDrawResolutionHistory(const TString& name, const KFResolutionTuningResult& result,
                                  const TString& directory) {
  const int count = static_cast<int>(result.groups.size());
  const int steps = static_cast<int>(result.iterations.size());
  auto* canvas = new TCanvas(name + "_ResolutionIterations", "Momentum resolution calibration",
                             480 * count, 760);
  canvas->Divide(count, 2);
  const double xmin = .5, xmax = steps + .5;
  auto drawFrame = [&](const TString& key, const TString& ytitle, double ymax) {
    gPad->SetLeftMargin(.20); gPad->SetBottomMargin(.16);
    TH1D frame(name + key, ";Fit iteration;" + ytitle, steps, xmin, xmax);
    frame.SetDirectory(nullptr); frame.SetStats(false);
    frame.SetMinimum(0.); frame.SetMaximum(ymax);
    frame.GetXaxis()->SetTitleSize(.042); frame.GetYaxis()->SetTitleSize(.042);
    frame.GetXaxis()->SetLabelSize(.038); frame.GetYaxis()->SetLabelSize(.038);
    frame.GetYaxis()->SetTitleOffset(1.55);
    for (int step = 1; step <= steps; ++step)
      frame.GetXaxis()->SetBinLabel(step, TString::Format("%d", step));
    frame.DrawCopy("axis");
  };
  for (int group = 0; group < count; ++group) {
    const auto& spec = result.groups[group];
    double maximum = spec.trueResolution;
    for (const auto& step : result.iterations) maximum = std::max(maximum, step.estimates[group]);
    canvas->cd(group + 1);
    drawFrame(spec.name + "_ResolutionFrame", "ResEstm (#sigma_{p}/p)", 1.35 * maximum);
    auto* estimate = new TGraph();
    estimate->SetName(name + "_" + spec.name);
    estimate->SetLineColor(kBlue + 1); estimate->SetMarkerColor(kBlue + 1);
    estimate->SetMarkerStyle(20); estimate->SetLineWidth(2);
    for (const auto& step : result.iterations)
      estimate->SetPoint(estimate->GetN(), step.step, step.estimates[group]);
    estimate->Draw("LP same");
    auto* truth = new TLine(xmin, spec.trueResolution, xmax, spec.trueResolution);
    truth->SetLineColor(kRed + 1); truth->SetLineStyle(2); truth->SetLineWidth(2); truth->Draw();
    auto* legend = new TLegend(.32, .74, .93, .92);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(.034);
    legend->AddEntry(estimate, spec.name + " (" + spec.title + ")", "lp");
    legend->AddEntry(truth, TString::Format("True = %.5g", spec.trueResolution), "l");
    legend->Draw();

    canvas->cd(count + group + 1);
    maximum = 1.3;
    for (const auto& step : result.iterations)
      if (step.valid[group]) maximum = std::max(maximum, step.widths[group] + step.widthErrors[group]);
    drawFrame(spec.name + "_PullFrame", "Correction-pull width", 1.3 * maximum);
    auto* widths = new TGraphErrors();
    widths->SetName(name + "_" + spec.name + "_PullWidth");
    widths->SetLineColor(kBlue + 1); widths->SetMarkerColor(kBlue + 1); widths->SetMarkerStyle(20);
    for (const auto& step : result.iterations) {
      if (!step.valid[group]) continue;
      const int point = widths->GetN();
      widths->SetPoint(point, step.step, step.widths[group]);
      widths->SetPointError(point, 0., step.widthErrors[group]);
    }
    widths->Draw("LP same");
    for (double value : {.9, 1., 1.1}) {
      auto* reference = new TLine(xmin, value, xmax, value);
      reference->SetLineColor(value == 1. ? kRed + 1 : kGreen + 2);
      reference->SetLineStyle(value == 1. ? 2 : 3); reference->Draw();
    }
    auto* label = new TPaveText(.34, .79, .93, .92, "NDC");
    label->SetBorderSize(0); label->SetFillStyle(0); label->SetTextSize(.034);
    label->AddText(result.converged ? "All widths stable: 0.9 < w < 1.1" : "Iteration limit: not all widths stable");
    label->Draw();
  }
  canvas->SaveAs(directory + "/ResolutionIterations.pdf");
  canvas->SaveAs(directory + "/ResolutionIterations.png");
  return canvas;
}

struct KFResolutionMonitorVariable {
  TString name;
  std::size_t group;
  std::vector<TH1D*> pulls; // One already-fitted histogram per completed round.
};

// One 2 x 1 canvas per momentum variable: all iteration pulls on the left,
// the matching tested resolution history on the right. Fits are reused.
std::vector<TCanvas*> KFDrawCombinedResolutionMonitor(
    const TString& name, const KFResolutionTuningResult& result,
    const std::vector<KFResolutionMonitorVariable>& variables, const TString& directory) {
  if (result.iterations.empty() || variables.empty())
    throw std::invalid_argument("Combined monitor requires pull histograms and resolution history");
  const int steps = static_cast<int>(result.iterations.size());
  for (const auto& variable : variables) {
    if (variable.group >= result.groups.size() || variable.pulls.size() != result.iterations.size())
      throw std::invalid_argument("Invalid variable history in combined resolution monitor");
    for (const auto* pull : variable.pulls)
      if (!pull) throw std::invalid_argument("Missing iteration pull in combined resolution monitor");
  }
  const int iteration = result.iterations.back().step;
  const std::size_t pages = variables.size();
  const double xmin = -.5, xmax = steps - .5;
  const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                        kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2};
  const std::size_t colorCount = sizeof(colors) / sizeof(colors[0]);
  std::vector<TCanvas*> canvases;
  for (std::size_t page = 0; page < pages; ++page) {
    const auto& variable = variables[page];
    const auto& spec = result.groups[variable.group];
    const TString key = name + TString::Format("_CombinedMonitor_Iteration%d_Page%zu", iteration, page + 1);
    auto* canvas = new TCanvas(key, variable.name + TString::Format(" - iter0 through iter%d", steps - 1), 1500, 680);
    canvas->Divide(2, 1, .005, .005);
    canvas->cd(1);
    {
      const double topMargin = std::min(.60, .12 + .04 * (steps + 1));
      gPad->SetLeftMargin(.16); gPad->SetRightMargin(.04);
      gPad->SetBottomMargin(.14); gPad->SetTopMargin(topMargin);
      double maximum = 1.;
      for (const auto* pull : variable.pulls) {
        maximum = std::max(maximum, pull->GetMaximum());
        if (const auto* gaussian = pull->GetFunction(TString(pull->GetName()) + "_Gaussian"))
          maximum = std::max(maximum, gaussian->GetParameter(0));
      }
      auto* fits = new TLegend(.16, 1. - topMargin + .015, .98, .985);
      fits->SetBorderSize(0); fits->SetFillStyle(0); fits->SetTextFont(42);
      fits->SetTextSize(std::min(.030, .42 / (steps + 1)));
      fits->SetHeader("Gaussian fits: mean, width");
      std::vector<TH1D*> displayed;
      for (int step = 0; step < steps; ++step) {
        auto* pull = variable.pulls[step];
        const int color = colors[step % colorCount];
        // DrawCopy gives the pad ownership of the histogram and its fit, so
        // older iterations remain available after the tuning function returns.
        auto* copy = static_cast<TH1D*>(pull->DrawCopy(step == 0 ? "hist" : "hist same"));
        copy->SetName(key + "_" + variable.name + TString::Format("_P_Pull_iter%d", step));
        copy->SetStats(false); copy->SetMinimum(0.); copy->SetMaximum(1.18 * maximum);
        copy->SetLineColor(color); copy->SetLineWidth(2);
        copy->GetXaxis()->SetTitleSize(.045); copy->GetYaxis()->SetTitleSize(.045);
        copy->GetXaxis()->SetLabelSize(.04); copy->GetYaxis()->SetLabelSize(.04);
        copy->GetYaxis()->SetTitleOffset(1.35);
        displayed.push_back(copy);
        const auto* gaussian = copy->GetFunction(TString(pull->GetName()) + "_Gaussian");
        const TString label = gaussian
          ? TString::Format("iter%d: %.3f #pm %.3f, %.3f #pm %.3f", step,
              gaussian->GetParameter(1), gaussian->GetParError(1),
              gaussian->GetParameter(2), gaussian->GetParError(2))
          : TString::Format("iter%d: Gaussian fit unavailable", step);
        fits->AddEntry(copy, label, "l");
      }
      // Draw fits last so no later histogram hides an earlier fitted curve.
      for (int step = 0; step < steps; ++step)
        if (auto* gaussian = displayed[step]->GetFunction(TString(variable.pulls[step]->GetName()) + "_Gaussian")) {
          gaussian->SetLineColor(colors[step % colorCount]);
          gaussian->SetLineStyle(2); gaussian->SetLineWidth(2);
          gaussian->ResetBit(TF1::kNotDraw); gaussian->Draw("same");
        }
      fits->Draw();
    }
    canvas->cd(2);
    gPad->SetLeftMargin(.20); gPad->SetRightMargin(.04);
    gPad->SetBottomMargin(.14); gPad->SetTopMargin(.22);
    gPad->SetGridx(); gPad->SetGridy();
    double maximum = spec.trueResolution;
    for (const auto& step : result.iterations) maximum = std::max(maximum, step.estimates[variable.group]);
    TH1D frame(key + "_ResolutionFrame", ";Iteration (0 = initial guess);ResEstm", steps, xmin, xmax);
    frame.SetDirectory(nullptr); frame.SetStats(false);
    frame.SetMinimum(0.); frame.SetMaximum(1.3 * maximum);
    frame.GetXaxis()->SetTitleSize(.045); frame.GetYaxis()->SetTitleSize(.045);
    frame.GetXaxis()->SetLabelSize(.04); frame.GetYaxis()->SetLabelSize(.04);
    frame.GetYaxis()->SetTitleOffset(1.8);
    for (int step = 0; step < steps; ++step) frame.GetXaxis()->SetBinLabel(step + 1, TString::Format("%d", step));
    frame.DrawCopy("axis");
    auto* estimates = new TGraph(); estimates->SetName(key + "_" + spec.name);
    estimates->SetLineColor(kBlue + 1); estimates->SetMarkerColor(kBlue + 1);
    estimates->SetLineWidth(2); estimates->SetMarkerStyle(20); estimates->SetMarkerSize(1.15);
    for (int step = 0; step < steps; ++step)
      estimates->SetPoint(step, step, result.iterations[step].estimates[variable.group]);
    estimates->Draw("LP same");
    auto* truth = new TLine(xmin, spec.trueResolution, xmax, spec.trueResolution);
    truth->SetLineColor(kRed + 1); truth->SetLineStyle(2); truth->SetLineWidth(2); truth->Draw();
    auto* legend = new TLegend(.17, .80, .97, .985);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextFont(42); legend->SetTextSize(.034);
    TString label = variable.name; label.ReplaceAll("_", " ");
    legend->AddEntry(estimates, spec.name + " (" + label + ")", "lp");
    legend->AddEntry(truth, TString::Format("True = %.5g", spec.trueResolution), "l");
    legend->Draw();
    TLatex values; values.SetTextFont(42); values.SetTextSize(.032);
    values.SetTextColor(kBlue + 1); values.SetTextAlign(21);
    for (int step = 0; step < steps; ++step) {
      const double estimate = result.iterations[step].estimates[variable.group];
      double labelY = estimate + .06 * maximum;
      if (labelY <= spec.trueResolution && spec.trueResolution - labelY < .05 * maximum)
        labelY = spec.trueResolution + .06 * maximum;
      values.DrawLatex(step, labelY, TString::Format("%.5g", estimate));
    }
    const TString snapshot = directory + TString::Format("/CombinedMonitor_Iteration%d_Page%zu", iteration, page + 1);
    canvas->SaveAs(snapshot + ".pdf");
    canvas->SaveAs(directory + "/CombinedMonitor_" + variable.name + ".pdf");
    canvases.push_back(canvas);
  }
  const TString latest = directory + "/CombinedMonitor.pdf";
  if (canvases.size() == 1) canvases.front()->SaveAs(latest);
  else
    for (std::size_t page = 0; page < canvases.size(); ++page)
      canvases[page]->Print(latest + (page == 0 ? "(" : (page + 1 == canvases.size() ? ")" : "")));
  return canvases;
}

typedef std::function<std::vector<double>(const KFResolutionEvent&, std::vector<double>&)> KFResolutionFit;

KFResolutionTuningResult KFRunResolutionTuning(
    const TString& name, const std::vector<KFResolutionEvent>& events,
    const std::vector<KFResolutionGroup>& groups, const KFResolutionFit& fit,
    int maxIterations, double initialScale, const TString& directory,
    const std::vector<KFTestVariable>& measuredVariables = {},
    const std::vector<KFTestVariable>& unknownVariables = {},
    const std::vector<KFTestDiagnostics::CanvasPage>& measuredPages = {}) {
  KFCheckResolutionTuning(static_cast<int>(events.size()), maxIterations, initialScale);
  if (groups.empty()) throw std::invalid_argument("No momentum resolutions to tune");
  gSystem->mkdir(directory, true);
  TFile output(directory + "/tuning.root", "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create resolution-tuning output");
  std::ofstream table((directory + "/ResolutionIterations.tsv").Data());
  if (!table) throw std::runtime_error("Cannot create resolution-tuning table");
  table << "iteration\tparameter\tResEstm\tResTrue\tpull_mean\tpull_width\twidth_error\tvalid\tstable\tfinite_events\tfailed_events\n";
  table << std::setprecision(12);
  KFResolutionTuningResult result;
  result.groups = groups;
  for (const auto& group : groups) result.estimates.push_back(initialScale * group.trueResolution);
  static unsigned long serial = 0;
  const TString prefix = name + TString::Format("_%lu", ++serial);
  const double missing = std::numeric_limits<double>::quiet_NaN();
  const std::size_t measuredCount = events.front().sample.variances.size();
  std::vector<KFTestVariable> measuredSpecs = measuredVariables;
  if (measuredSpecs.empty())
    for (std::size_t row = 0; row < measuredCount; ++row)
      measuredSpecs.emplace_back(TString::Format("row%zu", row), TString::Format("measurement %zu", row),
                                KFTestRange(100, -1., 1.));
  if (measuredSpecs.size() != measuredCount)
    throw std::invalid_argument("Measured pull labels do not match the cached measurements");
  std::vector<std::vector<std::unique_ptr<TH1D>>> momentumPullHistory(measuredCount);
  std::cout << "\n" << name << ": " << events.size() << " fixed events, initial scale = " << initialScale
            << ", maximum fit iterations = " << maxIterations << std::endl;
  for (int iteration = 1; iteration <= maxIterations; ++iteration) {
    KFResolutionIteration step;
    step.step = iteration; step.estimates = result.estimates;
    step.means.assign(groups.size(), missing); step.widths.assign(groups.size(), missing);
    step.widthErrors.assign(groups.size(), missing);
    step.valid.assign(groups.size(), false); step.stable.assign(groups.size(), false);
    std::vector<std::unique_ptr<TH1D>> pulls;
    std::vector<std::unique_ptr<TH1D>> measuredPulls, unknownPulls;
    auto bookDiagnosticPulls = [&](const std::vector<KFTestVariable>& specs,
                                   const TString& kind, std::vector<std::unique_ptr<TH1D>>& histograms) {
      for (const auto& spec : specs) {
        // Pulls are dimensionless. Keep the particle and coordinate explicit
        // in the PDF labels, including theta/phi when symbol fonts differ.
        TString coordinate = spec.name;
        coordinate.ReplaceAll("_", " ");
        histograms.emplace_back(new TH1D(prefix + TString::Format("_Iteration%d_", iteration) + kind + "_" + spec.name,
          ";" + (kind == "Measured" ? TString("Correction pull: ") : TString("Truth pull: ")) + coordinate + ";Entries",
          160, -8., 8.));
        histograms.back()->SetDirectory(nullptr); histograms.back()->Sumw2();
      }
    };
    bookDiagnosticPulls(measuredSpecs, "Measured", measuredPulls);
    bookDiagnosticPulls(unknownVariables, "Unknown", unknownPulls);
    for (const auto& group : groups) {
      pulls.emplace_back(new TH1D(prefix + TString::Format("_Iteration%d_", iteration) + group.name,
        ";Correction pull: " + group.title + ";Entries", 160, -8., 8.));
      pulls.back()->SetDirectory(nullptr); pulls.back()->Sumw2();
    }
    for (const auto& event : events) {
      std::vector<double> covariance = event.sample.variances;
      for (std::size_t group = 0; group < groups.size(); ++group)
        for (int row : groups[group].rows) {
          const double sigma = result.estimates[group] * std::abs(event.truth(row, 0));
          covariance.at(row) = sigma * sigma;
        }
      std::vector<double> correction;
      try { correction = fit(event, covariance); }
      catch (const std::exception&) { ++step.failedEvents; continue; }
      if (correction.size() != covariance.size() + unknownVariables.size()) { ++step.failedEvents; continue; }
      ++step.finiteEvents;
      for (std::size_t group = 0; group < groups.size(); ++group)
        for (int row : groups[group].rows)
          if (std::isfinite(correction[row])) pulls[group]->Fill(correction[row]);
      for (std::size_t row = 0; row < measuredPulls.size(); ++row)
        if (std::isfinite(correction[row])) measuredPulls[row]->Fill(correction[row]);
      for (std::size_t row = 0; row < unknownPulls.size(); ++row)
        if (std::isfinite(correction[measuredCount + row])) unknownPulls[row]->Fill(correction[measuredCount + row]);
    }
    auto* canvas = new TCanvas(prefix + TString::Format("_Pulls_Iteration%d", iteration),
                               "Momentum correction pulls", 480 * groups.size(), 380);
    canvas->Divide(groups.size(), 1);
    bool allStable = true;
    std::cout << "  Fit iteration " << iteration << ": " << step.finiteEvents
              << " finite, " << step.failedEvents << " failed fits\n";
    for (std::size_t group = 0; group < groups.size(); ++group) {
      canvas->cd(group + 1);
      KFTestDiagnostics::DrawCorrectionPull(pulls[group].get());
      const auto* gaussian = pulls[group]->GetFunction(TString(pulls[group]->GetName()) + "_Gaussian");
      if (gaussian) {
        step.valid[group] = true;
        step.means[group] = gaussian->GetParameter(1);
        step.widths[group] = gaussian->GetParameter(2);
        step.widthErrors[group] = gaussian->GetParError(2);
        step.stable[group] = step.widths[group] > 1 - margin && step.widths[group] < 1.+margin;
      }
      allStable = allStable && step.stable[group];
      std::cout << "    " << groups[group].name << " = " << step.estimates[group]
                << ", true = " << groups[group].trueResolution << ", pull width = "
                << step.widths[group] << (step.stable[group] ? " (stable)" : " (not stable)") << '\n';
      table << iteration << '\t' << groups[group].name << '\t' << step.estimates[group]
            << '\t' << groups[group].trueResolution << '\t' << step.means[group]
            << '\t' << step.widths[group] << '\t' << step.widthErrors[group]
            << '\t' << step.valid[group] << '\t' << step.stable[group]
            << '\t' << step.finiteEvents << '\t' << step.failedEvents << '\n';
      output.cd(); pulls[group]->Write();
    }
    canvas->SaveAs(directory + TString::Format("/MomentumPulls_Iteration%d.pdf", iteration));
    output.cd(); canvas->Write();
    // Fit every displayed coordinate, but keep these results separate from
    // step.widths/step.stable: only the original momentum groups may tune V.
    std::ofstream diagnosticTable((directory + TString::Format("/PullFits_Iteration%d.tsv", iteration)).Data());
    if (!diagnosticTable) throw std::runtime_error("Cannot create diagnostic pull-fit table");
    diagnosticTable << "kind\tvariable\tentries\tvalid\tmean\tmean_error\twidth\twidth_error\n";
    diagnosticTable << std::setprecision(12);
    auto drawDiagnosticPages = [&](std::vector<std::unique_ptr<TH1D>>& histograms,
        const std::vector<KFTestVariable>& specs, const TString& kind,
        std::vector<KFTestDiagnostics::CanvasPage> pages) {
      if (pages.empty())
        for (std::size_t begin = 0; begin < histograms.size(); begin += 9) {
          const int count = std::min<std::size_t>(9, histograms.size() - begin);
          const int columns = std::min(3, count);
          KFTestDiagnostics::CanvasPage page{columns, (count + columns - 1) / columns, {}};
          for (int item = 0; item < count; ++item) page.variables.push_back(begin + item);
          pages.push_back(page);
        }
      for (std::size_t page = 0; page < pages.size(); ++page) {
        const auto& layout = pages[page];
        auto* diagnostic = new TCanvas(prefix + TString::Format("_%sPulls_Iteration%d_Page%zu", kind.Data(), iteration, page + 1),
                                       kind + " diagnostic pulls", 400 * layout.columns, 340 * layout.rows);
        diagnostic->Divide(layout.columns, layout.rows);
        for (std::size_t pad = 0; pad < layout.variables.size(); ++pad) {
          const std::size_t row = layout.variables[pad];
          auto* histogram = histograms.at(row).get();
          diagnostic->cd(pad + 1); KFTestDiagnostics::DrawCorrectionPull(histogram);
          const auto* gaussian = histogram->GetFunction(TString(histogram->GetName()) + "_Gaussian");
          diagnosticTable << kind << '\t' << specs[row].name << '\t' << histogram->GetEntries() << '\t' << bool(gaussian)
            << '\t' << (gaussian ? gaussian->GetParameter(1) : missing)
            << '\t' << (gaussian ? gaussian->GetParError(1) : missing)
            << '\t' << (gaussian ? gaussian->GetParameter(2) : missing)
            << '\t' << (gaussian ? gaussian->GetParError(2) : missing) << '\n';
          output.cd(); histogram->Write();
        }
        diagnostic->SaveAs(directory + TString::Format("/%sPulls_Iteration%d_Page%zu.pdf", kind.Data(), iteration, page + 1));
        output.cd(); diagnostic->Write();
      }
    };
    drawDiagnosticPages(measuredPulls, measuredSpecs, "Measured", measuredPages);
    drawDiagnosticPages(unknownPulls, unknownVariables, "Unknown", {});
    // Finish the numerical exports before rendering the monitor.
    diagnosticTable.close();
    table.flush();
    if (!diagnosticTable || !table) throw std::runtime_error("Cannot write resolution-tuning tables");
    std::cout.flush();
    result.iterations.push_back(step);
    std::vector<KFResolutionMonitorVariable> monitorVariables;
    // Keep every fitted momentum histogram, not just the current round.
    // PhiPhi's four kaons come first; each retains its own pull fit and page.
    for (bool shared : {true, false})
      for (std::size_t group = 0; group < groups.size(); ++group)
        if ((groups[group].rows.size() > 1) == shared)
          for (int row : groups[group].rows) {
            auto& history = momentumPullHistory[row];
            history.push_back(std::move(measuredPulls[row]));
            output.cd();
            history.back()->Write(prefix + "_" + measuredSpecs[row].name + TString::Format("_P_Pull_iter%d", iteration - 1));
            KFResolutionMonitorVariable variable{measuredSpecs[row].name, group, {}};
            for (const auto& pull : history) variable.pulls.push_back(pull.get());
            monitorVariables.push_back(std::move(variable));
          }
    for (auto* monitor : KFDrawCombinedResolutionMonitor(prefix, result, monitorVariables, directory)) {
      output.cd(); monitor->Write();
    }
    if (allStable) { result.converged = true; break; }
    // Update sigma, not variance. A width of 2 doubles ResEstm. Keep a parameter
    // inside the stable band unchanged, and recheck it after the other updates.
    // The last recorded estimate is always one that was actually fitted.
    if (iteration < maxIterations)
      for (std::size_t group = 0; group < groups.size(); ++group)
        if (step.valid[group] && !step.stable[group]) {
          const double next = result.estimates[group] * step.widths[group];
          if (std::isfinite(next) && next > 0.) result.estimates[group] = next;
        }
  }
  auto* history = KFDrawResolutionHistory(prefix, result, directory);
  output.cd(); history->Write();
  for (auto* object : *history->GetListOfPrimitives()) {
    auto* pad = dynamic_cast<TPad*>(object);
    if (!pad) continue;
    for (auto* primitive : *pad->GetListOfPrimitives())
      if (auto* graph = dynamic_cast<TGraph*>(primitive)) graph->Write();
  }
  std::cout << name << (result.converged ? ": converged" : ": iteration limit reached; not converged")
            << " after " << result.iterations.size() << " fit iterations.\n";
  return result;
}

KFResolutionTuningResult TuneXiResolutions(int nEvents = 10000, int maxIterations = 4,
                                          double initialScale = .5, UInt_t seed = 12345,
                                          const TString& directory = "figs/TuneXiResolutions") {
  const int count = KFCheckResolutionTuning(nEvents, maxIterations, initialScale);
  SetStyle();
  const std::vector<KFResolutionGroup> groups = {
    {"ResPEstm", "p", ResP, {0}}, {"ResPi1Estm", "#pi_{#Lambda}", ResPi1, {3}},
    {"ResPi2Estm", "#pi_{#Xi}", ResPi2, {6}}};
  auto variables = KFCascadeMeasurements();
  const auto unknownVariables = KFCascadeUnknowns();
  KFInitialTuningVariables(variables, groups, initialScale);
  TRandom3 random(seed); KFResolutionRandomScope randomScope(random);
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18)), target(0., 0., 0., mP);
  TLorentzVector system = beam + target;
  double productionMasses[2] = {mXi, mK}, xiMasses[2] = {mL, mPi}, lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace production, xiDecay, lambdaDecay;
  if (!production.SetDecay(system, 2, productionMasses))
    throw std::runtime_error("Kinematically forbidden Xi production");
  std::vector<KFResolutionEvent> events; events.reserve(count);
  for (int event = 0; event < count; ++event) {
    production.Generate(); TLorentzVector xi = *production.GetDecay(0);
    if (!xiDecay.SetDecay(xi, 2, xiMasses)) throw std::runtime_error("Kinematically forbidden Xi decay");
    xiDecay.Generate(); TLorentzVector lambda = *xiDecay.GetDecay(0);
    if (!lambdaDecay.SetDecay(lambda, 2, lambdaMasses))
      throw std::runtime_error("Kinematically forbidden Lambda decay");
    lambdaDecay.Generate();
    TMatrixD truth = KFSphericalParameters({*lambdaDecay.GetDecay(0), *lambdaDecay.GetDecay(1), *xiDecay.GetDecay(1)});
    events.push_back({truth, KFSmear(truth, variables, random)});
  }
  auto result = KFRunResolutionTuning("TuneXiResolutions", events, groups,
    [unknownVariables](const KFResolutionEvent& event, std::vector<double>& covariance) {
      const TMatrixD& measured = event.sample.measurements;
      CascadeFitter fitter(KFFourVector(measured, 0, mP), KFFourVector(measured, 3, mPi),
                            KFFourVector(measured, 6, mPi));
      fitter.SetInvMass(mL, mXi);
      const TLorentzVector trueXi = KFFourVector(event.truth, 0, mP) + KFFourVector(event.truth, 3, mPi) +
                                     KFFourVector(event.truth, 6, mPi);
      return KFFitResolutionSample(fitter, covariance, KFSphericalParameters({trueXi}), unknownVariables);
    }, maxIterations, initialScale, directory, variables, unknownVariables,
       {{3, 3, {0, 3, 6, 1, 4, 7, 2, 5, 8}}});
  ResPEstm = result.estimates[0]; ResPi1Estm = result.estimates[1]; ResPi2Estm = result.estimates[2];
  return result;
}

KFResolutionTuningResult TuneMissingEtaResolutions(int nEvents = 10000, int maxIterations = 4,
                                                  double initialScale = .5, UInt_t seed = 12345,
                                                  const TString& directory = "figs/TuneMissingEtaResolutions") {
  const int count = KFCheckResolutionTuning(nEvents, maxIterations, initialScale);
  SetStyle();
  const std::vector<KFResolutionGroup> groups = {
    {"ResPEstm", "p", ResP, {0}}, {"ResPi1Estm", "#pi^{-}", ResPi1, {3}},
    {"ResKBeamEstm", "K^{-} beam", ResKBeam, {6}}};
  const auto unknownVariables = KFParticleVariables("lambda", "#Lambda", .5);
  auto variables = KFJoinVariables({
    KFParticleVariables("proton", "p", .2, ResP, ResThP, ResPhP, ScalePP),
    KFParticleVariables("pion", "#pi", .1, ResPi1, ResTh, ResPh, ScalePPi1),
    KFParticleVariables("beam", "K^{-}", .05, ResKBeam, .01, .01)});
  KFInitialTuningVariables(variables, groups, initialScale);
  TRandom3 random(seed); KFResolutionRandomScope randomScope(random);
  TVector3 direction(.1, 0., 1.); direction.SetMag(pE72);
  const TLorentzVector beam(direction, std::hypot(mK, pE72)), target(0., 0., 0., mP);
  TLorentzVector system = beam + target;
  double productionMasses[2] = {mL, mEta}, lambdaMasses[2] = {mP, mPi};
  TGenPhaseSpace production, lambdaDecay;
  if (!production.SetDecay(system, 2, productionMasses))
    throw std::runtime_error("Kinematically forbidden Lambda eta production");
  std::vector<KFResolutionEvent> events; events.reserve(count);
  for (int event = 0; event < count; ++event) {
    production.Generate(); TLorentzVector lambda = *production.GetDecay(0);
    if (!lambdaDecay.SetDecay(lambda, 2, lambdaMasses))
      throw std::runtime_error("Kinematically forbidden Lambda decay");
    lambdaDecay.Generate();
    TMatrixD truth = KFSphericalParameters({*lambdaDecay.GetDecay(0), *lambdaDecay.GetDecay(1), beam});
    events.push_back({truth, KFSmear(truth, variables, random)});
  }
  auto result = KFRunResolutionTuning("TuneMissingEtaResolutions", events, groups,
    [target, unknownVariables](const KFResolutionEvent& event, std::vector<double>& covariance) {
      const TMatrixD& measured = event.sample.measurements;
      const TLorentzVector p = KFFourVector(measured, 0, mP), pion = KFFourVector(measured, 3, mPi);
      MissingMassFitter fitter(KFFourVector(measured, 6, mK), target, mEta, p, pion, p + pion);
      fitter.SetInvMass(mL);
      const TLorentzVector trueLambda = KFFourVector(event.truth, 0, mP) + KFFourVector(event.truth, 3, mPi);
      return KFFitResolutionSample(fitter, covariance, KFSphericalParameters({trueLambda}), unknownVariables);
    }, maxIterations, initialScale, directory, variables, unknownVariables,
       {{3, 3, {0, 3, 6, 1, 4, 7, 2, 5, 8}}});
  ResPEstm = result.estimates[0]; ResPi1Estm = result.estimates[1]; ResKBeamEstm = result.estimates[2];
  return result;
}

KFResolutionTuningResult TunePhiPhiResolutions(int nEvents = 10000, int maxIterations = 4,
                                              double initialScale = .5, UInt_t seed = 12345,
                                              double beamMomentum = -1.,
                                              const TString& directory = "figs/TunePhiPhiResolutions") {
  const int count = KFCheckResolutionTuning(nEvents, maxIterations, initialScale);
  SetStyle();
  const double momentum = beamMomentum < 0. ? pPbar : beamMomentum;
  if (!std::isfinite(momentum) || momentum <= 0.) throw std::invalid_argument("Invalid antiproton beam momentum");
  TVector3 direction(.1, .05, 1.); direction.SetMag(momentum);
  const TLorentzVector beam(direction, std::hypot(momentum, mP)), target(0., 0., 0., mP);
  const double sqrtS = (beam + target).M(), minimum = 2. * mK;
  if (!std::isfinite(mPhi) || !std::isfinite(mK) || mK <= 0. ||
      !std::isfinite(WidthPhi) || WidthPhi < 0. || mPhi <= minimum ||
      sqrtS <= 4. * mK || (WidthPhi == 0. && sqrtS <= 2. * mPhi))
    throw std::runtime_error("Phi phi production or decay is kinematically forbidden");
  const std::vector<KFResolutionGroup> groups = {
    {"ResPbarBeamEstm", "#bar{p} beam", ResPbarBeam, {0}},
    {"ResPEstm", "four K tracks (pooled)", ResP, {3, 6, 9, 12}}};
  auto variables = KFPbarBalanceMeasurements(true);
  KFInitialTuningVariables(variables, groups, initialScale);
  TRandom3 random(seed); KFResolutionRandomScope randomScope(random);
  TGenPhaseSpace production, decay1, decay2;
  std::vector<KFResolutionEvent> events; events.reserve(count);
  for (int event = 0; event < count; ++event) {
    double phiMass1 = mPhi, phiMass2 = mPhi;
    if (WidthPhi > 0.) {
      const double s = sqrtS * sqrtS, maxMomentum = .5 * std::sqrt(s - 4. * minimum * minimum);
      bool accepted = false;
      for (int attempt = 0; attempt < 1000000; ++attempt) {
        phiMass1 = random.BreitWigner(mPhi, WidthPhi); phiMass2 = random.BreitWigner(mPhi, WidthPhi);
        if (!std::isfinite(phiMass1) || !std::isfinite(phiMass2) || phiMass1 <= minimum ||
            phiMass2 <= minimum || phiMass1 + phiMass2 >= sqrtS) continue;
        const double sum = phiMass1 + phiMass2, difference = phiMass1 - phiMass2;
        const double pStar = std::sqrt((s - sum * sum) * (s - difference * difference)) / (2. * sqrtS);
        if (random.Rndm() * maxMomentum < pStar) { accepted = true; break; }
      }
      if (!accepted) throw std::runtime_error("Unable to generate an allowed phi mass pair");
    }
    TLorentzVector system = beam + target;
    double phiMasses[2] = {phiMass1, phiMass2}, kaonMasses[2] = {mK, mK};
    if (!production.SetDecay(system, 2, phiMasses))
      throw std::runtime_error("Kinematically forbidden phi phi production");
    production.Generate(); TLorentzVector phi1 = *production.GetDecay(0), phi2 = *production.GetDecay(1);
    if (!decay1.SetDecay(phi1, 2, kaonMasses) || !decay2.SetDecay(phi2, 2, kaonMasses))
      throw std::runtime_error("Kinematically forbidden phi decay");
    decay1.Generate(); decay2.Generate();
    TMatrixD truth = KFSphericalParameters({beam, *decay1.GetDecay(0), *decay1.GetDecay(1),
                                                 *decay2.GetDecay(0), *decay2.GetDecay(1)});
    events.push_back({truth, KFSmear(truth, variables, random)});
  }
  auto result = KFRunResolutionTuning("TunePhiPhiResolutions", events, groups,
    [target](const KFResolutionEvent& event, std::vector<double>& covariance) {
      const TMatrixD& measured = event.sample.measurements;
      BalanceFitter fitter(KFFourVector(measured, 0, mP), target,
        {KFFourVector(measured, 3, mK), KFFourVector(measured, 6, mK),
         KFFourVector(measured, 9, mK), KFFourVector(measured, 12, mK)});
      return KFFitResolutionSample(fitter, covariance);
    }, maxIterations, initialScale, directory, variables, {},
       {{3, 1, {0, 1, 2}}, {4, 3, {3, 6, 9, 12, 4, 7, 10, 13, 5, 8, 11, 14}}});
  ResPbarBeamEstm = result.estimates[0]; ResPEstm = result.estimates[1];
  return result;
}

void TestKinfit() {
  SetStyle();
  std::cout << "TestXiFit([nEvents])\n"
               "TestXi1530Fit([nEvents], [signalFraction=1./3.], [updateVariance=false])\n"
               "TestFourVectorFit([nEvents])\n"
               "TestMissingEtaFit([nEvents])\n"
               "TestPbarElasticFit([nEvents], [beamMomentum])\n"
               "TestPbarPhiPhiFit([nEvents], [beamMomentum])\n"
               "TestPbar4PiPi0Fit([nEvents], [beamMomentum])\n"
               "TestPbarElasticPi0Fit([nEvents], [beamMomentum])\n"
               "TestPbarBGSeparation([nEvents], [beamMomentum], [signalFraction=.5], [massWindow=.100], [useKaonMassWindow=true])\n"
               "TestBalanceFit([nEvents], [beamMomentum]) runs both pbar tests (default 1 GeV/c)\n"
               "TunePhiPhiResolutions([nEvents=10000], [maxIterations=4], [initialScale=.5], [seed=12345])\n"
               "TuneXiResolutions([nEvents=10000], [maxIterations=4], [initialScale=.5], [seed=12345])\n"
               "TuneMissingEtaResolutions([nEvents=10000], [maxIterations=4], [initialScale=.5], [seed=12345])\n"
               "Variable setup: KFTestVariable(name, axisTitle, {bins,min,max}, sigma, period)\n"
               "Residual convention: initial/fitted minus truth; matrix row order must match.\n";
}
