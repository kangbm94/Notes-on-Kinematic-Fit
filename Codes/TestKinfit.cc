#include "src/FourVectorFitter.cc"
#include "src/MissingMassFitter.cc"
#include "src/CascadeFitter.cc"
#include "TestKinfit.hh"

#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>

// Load with: .L TestKinfit.cc+
// Existing fitter sources are used unchanged. Besides their existing API, they
// must provide GetBestMeasurements(), GetInitialMeasurements(), GetBestUnknowns()
// and GetInitialUnknowns(), all returning physical, unscaled TMatrixD columns.
// Edit the variable lists below, not MakeHist/Smear/FillHist/DrawHist.

int nev = 100000;
double mXi = 1.321, mXi0 = 1.315, mXi1530 = 1.535;
double mL = 1.115, mP = .938, mK = .493, mPi = .139, mPi0 = .135, mEta = .5478;
double pK18 = 1.8;
double ResP = .05, ResPi1 = .05, ResPi2 = .05;
double ResThP = .005, ResTh = .005, ResThPi2 = .005;
double ResPhP = .01, ResPh = .02, ResPhPi2 = .02;
double ScalePP = 1., ScalePPi1 = 1., ScalePPi2 = 1.;

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
      {name + "_phi", "#phi_{" + title + "} [rad]", {100, -.5, .5},
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
  // Assumed CascadeFitter unknown row order: Lambda first, Xi second.
  // If your implementation orders them differently, change this list AND the
  // trueUnknown column in KFFitCascade together. Dimension mismatches throw.
  return KFJoinVariables({KFParticleVariables("lambda", "#Lambda", .5),
                           KFParticleVariables("xi", "#Xi", .5)});
}

void SetStyle() {
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gStyle->SetTitleSize(.05, "XY"); gStyle->SetLabelSize(.045, "XY");
  gStyle->SetPadLeftMargin(.15); gStyle->SetPadBottomMargin(.15);
  gROOT->ForceStyle();
}

// Every decay used below is two-body, so its Generate() weight is constant.
void KFSetDecay(TGenPhaseSpace& generator, const TLorentzVector& parent,
                 double firstMass, double secondMass) {
  double masses[2] = {firstMass, secondMass};
  TLorentzVector copy = parent;
  if (!generator.SetDecay(copy, 2, masses))
    throw std::runtime_error("Kinematically forbidden test decay");
}

struct KFCascadeEvent {
  TLorentzVector proton, pi1, pi2;
};

KFCascadeEvent KFGenerateCascade(const TLorentzVector& xi, double pionMass = -1.) {
  TGenPhaseSpace evXi, evLambda;
  KFSetDecay(evXi, xi, mL, pionMass < 0. ? mPi : pionMass);
  evXi.Generate();
  const TLorentzVector lambda = *evXi.GetDecay(0);
  KFSetDecay(evLambda, lambda, mP, mPi);
  evLambda.Generate();
  return {*evLambda.GetDecay(0), *evLambda.GetDecay(1), *evXi.GetDecay(1)};
}

void KFAddCascadeObservables(KFTestDiagnostics& test) {
  test.AddObservable({"LambdaMass", "M(p#pi_{1}) [GeV/c^{2}]", {500, 1.04, 1.2}});
  test.AddObservable({"XiMass", "M(p#pi_{1}#pi_{2}) [GeV/c^{2}]", {500, mXi - .075, mXi + .075}});
}

void KFFitCascade(KFTestDiagnostics& test, const KFCascadeEvent& event,
                   bool updateVariance = false) {
  const auto& p = event.proton;
  const auto& pi1 = event.pi1;
  const auto& pi2 = event.pi2;
  const TMatrixD trueMeasured = KFSphericalParameters({p, pi1, pi2});
  const TMatrixD trueUnknown = KFSphericalParameters({p + pi1, p + pi1 + pi2});
  const KFTestSample sample = test.Smear(trueMeasured);
  const TLorentzVector measuredP = KFFourVector(sample.measurements, 0, mP);
  const TLorentzVector measuredPi1 = KFFourVector(sample.measurements, 3, mPi);
  const TLorentzVector measuredPi2 = KFFourVector(sample.measurements, 6, mPi);

  CascadeFitter fitter(measuredP, measuredPi1, measuredPi2);
  fitter.SetInvMass(mL, mXi);
  fitter.SetMaximumStep(5); // Same test setting as the original macro.
  fitter.UpdateVariance(updateVariance);
  sample.SetVariance(fitter);
  const double chi2 = fitter.DoKinematicFit();
  if (!test.FillHist(fitter, chi2, trueMeasured, trueUnknown)) return;

  const TMatrixD fitted = fitter.GetBestMeasurements();
  const TLorentzVector fittedP = KFFourVector(fitted, 0, mP);
  const TLorentzVector fittedPi1 = KFFourVector(fitted, 3, mPi);
  const TLorentzVector fittedPi2 = KFFourVector(fitted, 6, mPi);
  test.FillObservable("LambdaMass", (measuredP + measuredPi1).M(), (fittedP + fittedPi1).M());
  test.FillObservable("XiMass", (measuredP + measuredPi1 + measuredPi2).M(),
                       (fittedP + fittedPi1 + fittedPi2).M());
}


// Generic signal/background comparison, retaining the original cut studies.
// Tail integrals include overflow, and empty selections do not divide by zero.
void KFCompareSamples(const KFTestDiagnostics& signal, const KFTestDiagnostics& background,
                       const TString& directory) {
  gSystem->mkdir(directory, true);
  auto* summary = new TCanvas("Xi1530Comparison", "Xi1530 comparison", 1200, 800);
  summary->Divide(3, 1);
  const TH1D* sig[] = {signal.Chi2Histogram(), signal.PValueHistogram(),
                       signal.ObservableHistogram("XiMass")};
  const TH1D* bkg[] = {background.Chi2Histogram(), background.PValueHistogram(),
                       background.ObservableHistogram("XiMass")};
  for (int i = 0; i < 3; ++i) {
    summary->cd(i + 1);
    std::unique_ptr<TH1D> total(static_cast<TH1D*>(sig[i]->Clone("total")));
    total->SetDirectory(nullptr); total->Add(bkg[i]); total->SetLineColor(kBlack);
    total->SetMaximum(1.2 * std::max(1., total->GetMaximum()));
    TH1* all = total->DrawCopy("hist");
    TH1* s = sig[i]->DrawCopy("hist same"); s->SetLineColor(kBlue);
    TH1* b = bkg[i]->DrawCopy("hist same"); b->SetLineColor(kRed);
    auto* legend = new TLegend(.55, .68, .90, .90);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(all, "All events", "l"); legend->AddEntry(s, "Signal", "l");
    legend->AddEntry(b, "Combinatorial background", "l"); legend->Draw();
  }
  summary->SaveAs(directory + "/comparison.pdf");

  std::vector<double> pcuts, peff, ppurity, mcuts, meff, mpurity;
  const int np = sig[1]->GetNbinsX();
  const double ns = sig[1]->Integral(0, np + 1);
  for (int bin = 1; bin <= np; ++bin) {
    const double s = sig[1]->Integral(bin, np + 1);
    const double b = bkg[1]->Integral(bin, np + 1);
    if (ns <= 0. || s + b <= 0.) continue;
    pcuts.push_back(sig[1]->GetXaxis()->GetBinLowEdge(bin));
    peff.push_back(s / ns); ppurity.push_back(s / (s + b));
  }
  const int nm = sig[2]->GetNbinsX();
  const double nsm = sig[2]->Integral(0, nm + 1);
  const double maxWindow = std::min(mXi - sig[2]->GetXaxis()->GetXmin(),
                                     sig[2]->GetXaxis()->GetXmax() - mXi);
  for (int k = 1; k <= nm / 2; ++k) {
    const double window = maxWindow * k / (nm / 2);
    double s = 0., b = 0.;
    for (int bin = 1; bin <= nm; ++bin)
      if (std::abs(sig[2]->GetBinCenter(bin) - mXi) <= window) {
        s += sig[2]->GetBinContent(bin); b += bkg[2]->GetBinContent(bin);
      }
    if (nsm <= 0. || s + b <= 0.) continue;
    mcuts.push_back(window); meff.push_back(s / nsm); mpurity.push_back(s / (s + b));
  }
  auto* cuts = new TCanvas("Xi1530Cuts", "Selection studies", 1200, 800);
  cuts->Divide(3, 1);
  auto drawPair = [&](int pad, const std::vector<double>& x,
                       const std::vector<double>& first, const std::vector<double>& second,
                       const char* title) {
    if (x.empty()) return;
    cuts->cd(pad);
    auto* g1 = new TGraph(static_cast<int>(x.size()), x.data(), first.data());
    auto* g2 = new TGraph(static_cast<int>(x.size()), x.data(), second.data());
    g1->SetTitle(title); g1->SetLineColor(kBlue); g2->SetLineColor(kRed);
    g1->SetMinimum(0.); g1->SetMaximum(1.05); g1->Draw("AL"); g2->Draw("L same");
    auto* legend = new TLegend(.58, .72, .90, .90);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(g1, "Efficiency", "l");
    legend->AddEntry(g2, "Purity", "l");
    legend->Draw();
  };
  drawPair(1, pcuts, peff, ppurity, ";P-value cut;Efficiency / purity");
  drawPair(2, mcuts, meff, mpurity, ";Mass half-window [GeV/c^{2}];Efficiency / purity");
  // The efficiency coordinates of the two selections differ; draw separately.
  cuts->cd(3);
  TGraph* pcurve = nullptr;
  if (!peff.empty()) {
    pcurve = new TGraph(static_cast<int>(peff.size()), peff.data(), ppurity.data());
    pcurve->SetTitle(";Efficiency;Purity"); pcurve->SetLineColor(kBlue);
    pcurve->SetMinimum(0.); pcurve->SetMaximum(1.05); pcurve->Draw("AL");
    pcurve->GetXaxis()->SetLimits(0., 1.);
  }
  if (!meff.empty()) {
    auto* mcurve = new TGraph(static_cast<int>(meff.size()), meff.data(), mpurity.data());
    mcurve->SetTitle(";Efficiency;Purity"); mcurve->SetLineColor(kRed);
    mcurve->SetMinimum(0.); mcurve->SetMaximum(1.05); mcurve->Draw(pcurve ? "L same" : "AL");
    auto* legend = new TLegend(.56, .28, .90, .46);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    if (pcurve) legend->AddEntry(pcurve, "P-value cut", "l");
    legend->AddEntry(mcurve, "Mass-window cut", "l"); legend->Draw();
  }
  cuts->SaveAs(directory + "/selection.pdf");
}

void TestXiFit(int nEvents = -1) {
  KFTestDiagnostics test("XiCascade", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFAddCascadeObservables(test);
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TGenPhaseSpace evVert;
  KFSetDecay(evVert, beam + target, mXi, mK);
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    evVert.Generate();
    KFFitCascade(test, KFGenerateCascade(*evVert.GetDecay(0)), true);
  }
  test.DrawHist("figs/XiFit");
  test.WriteHist("figs/XiFit/histograms.root");
}
void TestXi1530Fit(int nEvents = -1) {
  KFTestDiagnostics signal("Xi1530Signal", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFTestDiagnostics background("Xi1530Background", KFCascadeMeasurements(), KFCascadeUnknowns());
  KFAddCascadeObservables(signal); KFAddCascadeObservables(background);
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TGenPhaseSpace evVert, evXi1530;
  KFSetDecay(evVert, beam + target, mXi1530, mK);
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    int attempts = 0;
    do {
      evVert.Generate();
      if (++attempts > 1000000) throw std::runtime_error("K+ selection has no generated acceptance");
    } while (evVert.GetDecay(1)->CosTheta() < .97);
    const TLorentzVector xi1530 = *evVert.GetDecay(0);
    KFSetDecay(evXi1530, xi1530, mXi, mPi0); evXi1530.Generate();
    KFFitCascade(signal, KFGenerateCascade(*evXi1530.GetDecay(0)), false);

    // Wrong hypothesis: Xi1530 -> Xi0 pi-, Xi0 -> Lambda pi0. The charged pion
    // from Xi1530 is paired with Lambda; no generated Xi- is present here.
    KFSetDecay(evXi1530, xi1530, mXi0, mPi); evXi1530.Generate();
    KFCascadeEvent wrong = KFGenerateCascade(*evXi1530.GetDecay(0), mPi0);
    wrong.pi2 = *evXi1530.GetDecay(1);
    KFFitCascade(background, wrong, true);
  }
  signal.DrawHist("figs/Xi1530Fit"); background.DrawHist("figs/Xi1530Fit");
  signal.WriteHist("figs/Xi1530Fit/signal.root");
  background.WriteHist("figs/Xi1530Fit/background.root");
  KFCompareSamples(signal, background, "figs/Xi1530Fit");
}

void TestFourVectorFit(int nEvents = -1) {
  const auto allVariables = KFCascadeMeasurements();
  const std::vector<KFTestVariable> lambdaMeasured(allVariables.begin(), allVariables.begin() + 6);
  KFTestDiagnostics lambdaTest("LambdaSequential", lambdaMeasured,
                                 KFParticleVariables("lambda", "#Lambda", .5));
  const std::vector<KFTestVariable> pion2Variables(allVariables.begin() + 6, allVariables.end());
  auto xiMeasured = KFJoinVariables({KFParticleVariables("lambda", "#Lambda", .5),
                                      pion2Variables});
  KFTestDiagnostics xiTest("XiSequential", xiMeasured, KFParticleVariables("xi", "#Xi", .5));
  lambdaTest.AddObservable({"Mass", "M(p#pi_{1}) [GeV/c^{2}]", {500, 1.04, 1.2}});
  xiTest.AddObservable({"Mass", "M(#Lambda#pi_{2}) [GeV/c^{2}]", {500, mXi - .075, mXi + .075}});
  const TLorentzVector beam(0., 0., pK18, std::hypot(mK, pK18));
  const TLorentzVector target(0., 0., 0., mP);
  TGenPhaseSpace evVert;
  KFSetDecay(evVert, beam + target, mXi, mK);
  unsigned long skippedSecondFits = 0;
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    evVert.Generate();
    const auto event = KFGenerateCascade(*evVert.GetDecay(0));
    const TLorentzVector trueLambda = event.proton + event.pi1;
    const TLorentzVector trueXi = trueLambda + event.pi2;
    const KFTestSample sample = KFSmear(KFSphericalParameters({event.proton, event.pi1, event.pi2}),
                                        allVariables, *gRandom);
    const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pi1 = KFFourVector(sample.measurements, 3, mPi);
    const TLorentzVector pi2 = KFFourVector(sample.measurements, 6, mPi);
    FourVectorFitter lambdaFit(p, pi1, p + pi1);
    lambdaFit.SetInvMass(mL);
    KFTestSample lambdaSample;
    lambdaSample.variances.assign(sample.variances.begin(), sample.variances.begin() + 6);
    lambdaSample.SetVariance(lambdaFit);
    const double lambdaChi2 = lambdaFit.DoKinematicFit();
    if (!lambdaTest.FillHist(lambdaFit, lambdaChi2,
                             KFSphericalParameters({event.proton, event.pi1}),
                             KFSphericalParameters({trueLambda}))) {
      ++skippedSecondFits; continue;
    }
    const TMatrixD lambdaBest = lambdaFit.GetBestMeasurements();
    const TLorentzVector pFit = KFFourVector(lambdaBest, 0, mP);
    const TLorentzVector pi1Fit = KFFourVector(lambdaBest, 3, mPi);
    const TLorentzVector lambda = KFFourVector(lambdaFit.GetBestUnknowns(), 0, mL);
    lambdaTest.FillObservable("Mass", (p + pi1).M(), (pFit + pi1Fit).M());
    const TMatrixD lambdaCovariance = lambdaFit.GetUnmeasuredCovariance();
    TMatrixD xiCovariance(6, 6);
    xiCovariance.Zero();
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) xiCovariance(j, k) = lambdaCovariance(j, k);
      xiCovariance(j + 3, j + 3) = sample.variances[j + 6];
    }
    FourVectorFitter xiFit(lambda, pi2, lambda + pi2);
    xiFit.SetInvMass(mXi);
    try {
      KFSetCovariance(xiFit, xiCovariance);
    } catch (const std::invalid_argument&) {
      ++skippedSecondFits; continue;
    }
    const double xiChi2 = xiFit.DoKinematicFit();
    if (xiTest.FillHist(xiFit, xiChi2, KFSphericalParameters({trueLambda, event.pi2}),
                        KFSphericalParameters({trueXi}))) {
      const TMatrixD xiBest = xiFit.GetBestMeasurements();
      xiTest.FillObservable("Mass", (lambda + pi2).M(),
                             (KFFourVector(xiBest, 0, mL) + KFFourVector(xiBest, 3, mPi)).M());
    }
  }
  lambdaTest.DrawHist("figs/FourVectorFit"); xiTest.DrawHist("figs/FourVectorFit");
  lambdaTest.WriteHist("figs/FourVectorFit/lambda.root");
  xiTest.WriteHist("figs/FourVectorFit/xi.root");
  std::cout << "Second-stage fits skipped after invalid Lambda result/covariance: "
            << skippedSecondFits << std::endl;
}

// Example for the current MissingMassFitter interface, without introducing a
// mandatory dependency on its source. Load your corrected implementation first,
// then call TestMissingEtaFit<MissingMassFitter>(10000).
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
  TGenPhaseSpace evVert, evLambda;
  KFSetDecay(evVert, beam + target, mL, mEta);
  for (int i = 0; i < (nEvents < 0 ? nev : nEvents); ++i) {
    evVert.Generate();
    const TLorentzVector lambda = *evVert.GetDecay(0);
    KFSetDecay(evLambda, lambda, mP, mPi); evLambda.Generate();
    const TMatrixD truth = KFSphericalParameters({*evLambda.GetDecay(0), *evLambda.GetDecay(1), beam});
    const KFTestSample sample = test.Smear(truth);
    const TLorentzVector p = KFFourVector(sample.measurements, 0, mP);
    const TLorentzVector pion = KFFourVector(sample.measurements, 3, mPi);
    const TLorentzVector measuredBeam = KFFourVector(sample.measurements, 6, mK);
    MissingMassFitter fitter(measuredBeam, target, mEta, p, pion, p + pion);
    fitter.SetInvMass(mL); sample.SetVariance(fitter);
    const double chi2 = fitter.DoKinematicFit();
    if (!test.FillHist(fitter, chi2, truth, KFSphericalParameters({lambda}))) continue;
    const TMatrixD best = fitter.GetBestMeasurements();
    const TLorentzVector fittedLambda = KFFourVector(fitter.GetBestUnknowns(), 0, mL);
    test.FillObservable("LambdaMass", (p + pion).M(),
                         (KFFourVector(best, 0, mP) + KFFourVector(best, 3, mPi)).M());
    test.FillObservable("MissingMass", (measuredBeam + target - (p + pion)).M(),
                         (KFFourVector(best, 6, mK) + target - fittedLambda).M());
  }
  test.DrawHist("figs/MissingEtaFit");
  test.WriteHist("figs/MissingEtaFit/histograms.root");
}

void TestKinfit() {
  SetStyle();
  std::cout << "TestXiFit([nEvents])\n"
               "TestXi1530Fit([nEvents])\n"
               "TestFourVectorFit([nEvents])\n"
               "TestMissingEtaFit([nEvents]) after loading that fitter\n"
               "Variable setup: KFTestVariable(name, axisTitle, {bins,min,max}, sigma, period)\n"
               "Residual convention: initial/fitted minus truth; matrix row order must match.\n";
}
