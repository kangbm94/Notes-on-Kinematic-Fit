#include "src/FourVectorFitter.cc"
#include "src/CascadeFitter.cc"
#include "TestKinfit.hh"
// Code to evaluate the Kineamtic Fit performance for the decay of the Xi to Lambda and pion
// Interaction between 1.8 GeV/c Kaon- and proton target is assumed.
// Author: Kang Byungmin, kangbmw2@naver.com
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

int nev = 1000;
TGenPhaseSpace EvVert;
TGenPhaseSpace EvXi;
TGenPhaseSpace EvLd;
double pK18 = 1.8;
TLorentzVector KM(0, 0, pK18, hypot(mK, pK18));
TLorentzVector PT(0, 0, 0, mP);
TLorentzVector Vertex = KM + PT;
double VertMass[2] = {mXi, mK};
double XiDecayMass[2] = {mL, mPi};
double LdDecayMass[2] = {mP, mPi};
void TestKinfit()
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
//	MakeHists();
	cout << "TestCascadeFit()" << endl;
	cout << "TestFourVectorFit()" << endl;
}
void TestCascadeFit()
{
	MakeHistCascade();
	for (int i = 0; i < nev; ++i)
	{
		EvVert.SetDecay(Vertex, 2, VertMass);
		EvVert.Generate();
		iev = i;
		if (iev % 1000 == 0)
		{
			cout << Form("Event %d", iev) << endl;
		}
		auto Xi = *EvVert.GetDecay(0);
		auto KP = *EvVert.GetDecay(1);
		EvXi.SetDecay(Xi, 2, XiDecayMass);
		EvXi.Generate();
		auto Ld = *EvXi.GetDecay(0);
		auto Pi2 = *EvXi.GetDecay(1);
		EvLd.SetDecay(Ld, 2, LdDecayMass);
		EvLd.Generate();
		auto P = *EvLd.GetDecay(0);
		auto Pi1 = *EvLd.GetDecay(1);

		auto TVLd = Ld.Vect();
		PLd = TVLd.Mag();
		ThLd = Ld.Theta();
		PhLd = Ld.Phi();

		auto TVP = P.Vect();
		auto TVPi1 = Pi1.Vect();
		auto TVPi2 = Pi2.Vect();
		PP = TVP.Mag();
		ThP = TVP.Theta();
		PhP = TVP.Phi();
		PPi1 = TVPi1.Mag();
		ThPi1 = TVPi1.Theta();
		PhPi1 = TVPi1.Phi();
		PPi2 = TVPi2.Mag();
		ThPi2 = TVPi2.Theta();
		PhPi2 = TVPi2.Phi();

		TVector3 V1 = TVLd;
		V1.SetTheta(ThLdMeas);
		V1.SetPhi(PhLdMeas);
		TVector3 V2(0, 0, 0);
		PPMeas = gRandom->Gaus(PP * ScalePP, PP * ResP);
		ThPMeas = gRandom->Gaus(ThP, ResThP);
		PhPMeas = gRandom->Gaus(PhP, ResPhP);
		PPi1Meas = gRandom->Gaus(PPi1 * ScalePPi1, PPi1 * ResPi1);
		ThPi1Meas = gRandom->Gaus(ThPi1, ResTh);
		PhPi1Meas = gRandom->Gaus(PhPi1, ResPh);
		PPi2Meas = gRandom->Gaus(PPi2 * ScalePPi2, PPi2 * ResPi2);
		ThPi2Meas = gRandom->Gaus(ThPi2, ResThPi2);
		PhPi2Meas = gRandom->Gaus(PhPi2, ResPhPi2);

		TVector3 TVPMeas(0, 0, PPMeas);
		TVPMeas.SetTheta(ThPMeas);
		TVPMeas.SetPhi(PhPMeas);
		TVector3 TVPi1Meas(0, 0, PPi1Meas);
		TVPi1Meas.SetTheta(ThPi1Meas);
		TVPi1Meas.SetPhi(PhPi1Meas);
		TVector3 TVPi2Meas(0, 0, PPi2Meas);
		TVPi2Meas.SetTheta(ThPi2Meas);
		TVPi2Meas.SetPhi(PhPi2Meas);
		TLorentzVector PMeas(TVPMeas, hypot(PPMeas, mP));
		TLorentzVector Pi1Meas(TVPi1Meas, hypot(PPi1Meas, mPi));
		TLorentzVector Pi2Meas(TVPi2Meas, hypot(PPi2Meas, mPi));
		double rp = PPMeas * ResP;
		double rpi1 = PPi1Meas * ResPi1;
		double rpi2 = PPi2Meas * ResPi2;

		auto LdRecon = PMeas + Pi1Meas;
		auto LdReconMeas = LdRecon;
		InvMLd = LdRecon.Mag();
		auto TVLdRecon = LdRecon.Vect();
		PLdMeas = TVLdRecon.Mag();
		ThLdMeas = TVLdRecon.Theta();
		PhLdMeas = TVLdRecon.Phi();
		LdRecon = TLorentzVector(TVLdRecon,hypot(mL,TVLdRecon.Mag()));
		auto XiRecon = LdRecon + Pi2Meas;
		InvMXi = (LdRecon + Pi2).Mag();
		auto TVXiMeas = XiRecon.Vect();
		auto TVXi = Xi.Vect();
		PXiMeas = TVXiMeas.Mag();
		ThXiMeas = TVXiMeas.Theta();
		PhXiMeas = TVXiMeas.Phi();
		PXi = TVXi.Mag();
		ThXi = TVXi.Theta();
		PhXi = TVXi.Phi();

		TLorentzVector MM = Xi - P - Pi1 - Pi2;
		auto TVXiRecon = TVP + TVPi1 + TVPi2;
		double PXiRecon = TVXiRecon.Mag();
		double EXi = hypot(PXiRecon,mXi);
		double ELd = hypot(PLd,mL);
		double EP = hypot(PP,mP);
		double EPi1 = hypot(PPi1,mPi);
		double EPi2 = hypot(PPi2,mPi); 
//		cout<<Form("dE = %f",E_Xi-E_P-E_Pi1-E_Pi2)<<endl;
		CascadeFitter KFCascade(PMeas, Pi1Meas, Pi2Meas);
//		CascadeFitter KFCascade(P, Pi1, Pi2);
		KFCascade.SetMaximumStep(5);
		KFCascade.UpdateVariance(true);
//		KFCascade.ScaleParameters(false);
		double Variance[9] = {
			rp * rp, ResThP * ResThP, ResPhP * ResPhP,
			rpi1 * rpi1, ResTh * ResTh, ResPh * ResPh,
			rpi2 * rpi2, ResThPi2 * ResThPi2, ResPhPi2 * ResPhPi2
		};
		KFCascade.SetInvMass(mL,mXi);
		KFCascade.SetVariance(Variance);
		Chi2 = KFCascade.DoKinematicFit();
		double Pval = KFCascade.GetPValue();
		auto Pull = KFCascade.GetPull();
		pulls = Pull;
		auto LVCont = KFCascade.GetFittedLV();
		auto LVPCor = LVCont.at(0);
		auto LVPi1Cor = LVCont.at(1);
		auto LVPi2Cor = LVCont.at(2);
		auto LVLdCor = LVCont.at(3);
		auto LVXiCor = LVCont.at(4);

		auto TVPCor = LVPCor.Vect();
		auto TVPi1Cor = LVPi1Cor.Vect();
		auto TVPi2Cor = LVPi2Cor.Vect();
		auto TVLdCor = LVLdCor.Vect();
		auto TVXiCor = LVXiCor.Vect();

		InvMLdCor = (LVPCor+LVPi1Cor).Mag();
		auto LdCorFix = TLorentzVector((LVPCor+LVPi1Cor).Vect(),hypot(mL,(LVPCor+LVPi1Cor).Vect().Mag()));
		InvMXiCor = (LdCorFix+LVPi2Cor).Mag();

		PXiCor = TVXiCor.Mag();
		ThXiCor = TVXiCor.Theta();
		PhXiCor = TVXiCor.Phi();
		PLdCor = TVLdCor.Mag();
		ThLdCor = TVLdCor.Theta();
		PhLdCor = TVLdCor.Phi();
		PPCor = TVPCor.Mag();
		ThPCor = TVPCor.Theta();
		PhPCor = TVPCor.Phi();
		PPi1Cor = TVPi1Cor.Mag();
		ThPi1Cor = TVPi1Cor.Theta();
		PhPi1Cor = TVPi1Cor.Phi();
		PPi2Cor = TVPi2Cor.Mag();
		ThPi2Cor = TVPi2Cor.Theta();
		PhPi2Cor = TVPi2Cor.Phi();
		constsAfter = KFCascade.GetKinematicConstraints();
		constsIni = KFCascade.GetInitialConstraints();
//		if(Pval<0.01)continue;
		FillHistCascade();

	}
	TCanvas* cKF = new TCanvas("cKF","DaughterResidual",1200,800);
	cKF->Divide(3,3);
	for(int i=0;i<9;++i){
		cKF->cd(i+1);
		hCorResi[i]->Draw("same");
		hMeasResi[i]->Draw("same");
	}
	TCanvas* cKF2 = new TCanvas("cKF2","MotherResidual",1200,800);
	cKF2->Divide(3,2);
	for(int i=0;i<6;++i){
		cKF2->cd(i+1);
		hCorResi[i+9]->Draw("same");
		hMeasResi[i+9]->Draw("same");
	}
	TCanvas* cKF3 = new TCanvas("cKF3","InvMass",1200,800);
	cKF3->Divide(2,1);
	cKF3->cd(1);
	hInvMLdCor->Draw("same");
	hInvMLd->Draw("same");
	cKF3->cd(2);
	hInvMXiCor->Draw("same");
	hInvMXi->Draw("same");
	TCanvas* cKF4 = new TCanvas("cKF4","Pull",1200,800);
	cKF4->Divide(3,3);
	for(int i=0;i<9;++i){
		cKF4->cd(i+1);
		hPull[i]->Draw();
	}	
	TCanvas* cKF5 = new TCanvas("cKF5","Constraints",1200,800);
	cKF5->Divide(4,2);
	for(int i=0;i<8;++i){
		cKF5->cd(i+1);
		hConstIni[i]->Draw();
		hConst[i]->Draw("same");
	}	

}
void TestFourVectorFit()
{	
	MakeHists();
	TH1D* hViolation[8];
	TH1D* hCorrection[8];
	for(int i=0;i<8;++i){
		hViolation[i] = new TH1D(Form("hViolation%d",i),Form("hViolation%d",i),1000,-1,1);
		hCorrection[i] = new TH1D(Form("hCorrection%d",i),Form("hCorrection%d",i),1000,-1,1);
		hViolation[i]->SetLineColor(kRed);
		hCorrection[i]->SetLineColor(kBlue);
	}
	for (int i = 0; i < nev; ++i)
	{
		EvVert.SetDecay(Vertex, 2, VertMass);
		EvVert.Generate();
		iev = i;
		if (iev % 1000 == 0)
		{
			cout << Form("Event %d", iev) << endl;
		}
		auto Xi = *EvVert.GetDecay(0);
		auto KP = *EvVert.GetDecay(1);
		EvXi.SetDecay(Xi, 2, XiDecayMass);
		EvXi.Generate();
		auto Ld = *EvXi.GetDecay(0);
		auto Pi2 = *EvXi.GetDecay(1);
		EvLd.SetDecay(Ld, 2, LdDecayMass);
		EvLd.Generate();
		auto P = *EvLd.GetDecay(0);
		auto Pi1 = *EvLd.GetDecay(1);

		auto TVLd = Ld.Vect();
		PLd = TVLd.Mag();
		ThLd = Ld.Theta();
		PhLd = Ld.Phi();

		auto TVP = P.Vect();
		auto TVPi1 = Pi1.Vect();
		auto TVPi2 = Pi2.Vect();
		PP = TVP.Mag();
		ThP = TVP.Theta();
		PhP = TVP.Phi();
		PPi1 = TVPi1.Mag();
		ThPi1 = TVPi1.Theta();
		PhPi1 = TVPi1.Phi();
		PPi2 = TVPi2.Mag();
		ThPi2 = TVPi2.Theta();
		PhPi2 = TVPi2.Phi();

		ThLdMeas = gRandom->Gaus(ThLd, ResThV);
		PhLdMeas = gRandom->Gaus(PhLd, ResPhV);
		TVector3 V2(0, 0, 0);
		PPMeas = gRandom->Gaus(PP * ScalePP, PP * ResP);
		ThPMeas = gRandom->Gaus(ThP, ResThP);
		PhPMeas = gRandom->Gaus(PhP, ResPhP);
		PPi1Meas = gRandom->Gaus(PPi1 * ScalePPi1, PPi1 * ResPi1);
		ThPi1Meas = gRandom->Gaus(ThPi1, ResTh);
		PhPi1Meas = gRandom->Gaus(PhPi1, ResPh);
		PPi2Meas = gRandom->Gaus(PPi2 * ScalePPi2, PPi2 * ResPi2);
		ThPi2Meas = gRandom->Gaus(ThPi2, ResTh);
		PhPi2Meas = gRandom->Gaus(PhPi2, ResPh);

		TVector3 TVPMeas(0, 0, PPMeas);
		TVPMeas.SetTheta(ThPMeas);
		TVPMeas.SetPhi(PhPMeas);
		TVector3 TVPi1Meas(0, 0, PPi1Meas);
		TVPi1Meas.SetTheta(ThPi1Meas);
		TVPi1Meas.SetPhi(PhPi1Meas);
		TVector3 TVPi2Meas(0, 0, PPi2Meas);
		TVPi2Meas.SetTheta(ThPi2Meas);
		TVPi2Meas.SetPhi(PhPi2Meas);
		TLorentzVector PMeas(TVPMeas, hypot(PPMeas, mP));
		TLorentzVector Pi1Meas(TVPi1Meas, hypot(PPi1Meas, mPi));
		TLorentzVector Pi2Meas(TVPi2Meas, hypot(PPi2Meas, mPi));
		double rp = PPMeas * ResP;
		double rpi1 = PPi1Meas * ResPi1;
		double rpi2 = PPi2Meas * ResPi2;

		auto LdRecon = PMeas + Pi1Meas;
		auto LdReconMeas = LdRecon;
//		LdReconMeas.SetTheta(ThLdMeas);
//		LdReconMeas.SetPhi(PhLdMeas);

		auto XiRecon = LdRecon + Pi2Meas;
		InvMLd = LdRecon.Mag();
		InvMXi = (LdRecon + Pi2).Mag();
		auto TVXiMeas = XiRecon.Vect();
		auto TVXi = Xi.Vect();
		PXiMeas = TVXiMeas.Mag();
		ThXiMeas = TVXiMeas.Theta();
		PhXiMeas = TVXiMeas.Phi();
		PXi = TVXi.Mag();
		ThXi = TVXi.Theta();
		PhXi = TVXi.Phi();

		//		FourVectorFitter KFLd(PMeas,Pi1Meas,LdReconMeas);
		//		KFLd.UseVertex(true,V1,V2);
		//		double Variance[8] = {ResThV*ResThV,ResPhV*ResPhV,rp*rp,ResThP*ResThP,ResPhP*ResPhP,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};

		FourVectorFitter KFLd(PMeas, Pi1Meas, LdRecon);
//		KFLd.UpdateVariance(true);
		double Variance[6] = {rp * rp, ResThP * ResThP, ResPhP * ResPhP, rpi1 * rpi1, ResTh * ResTh, ResPh * ResPh};
		KFLd.SetInvMass(mL);
		KFLd.SetVariance(Variance);
		Chi2 = KFLd.DoKinematicFit();
		pvalLd = KFLd.GetPValue();
		auto PullLd = KFLd.GetPull();
		auto cont = KFLd.GetFittedLV();
		auto PCor = cont.at(0);
		auto Pi1Cor = cont.at(1);
		auto LdCor = cont.at(2);
		auto TVLdCor = LdCor.Vect();
		PLdCor = TVLdCor.Mag();
		ThLdCor = TVLdCor.Theta();
		PhLdCor = TVLdCor.Phi();
		auto LdFit = PCor + Pi1Cor;
		auto TVLdFit = LdFit.Vect();
		auto TVLdMeas = TVPMeas + TVPi1Meas;
		PLdMeas = TVLdMeas.Mag();
		ThLdMeas = TVLdMeas.Theta();
		PhLdMeas = TVLdMeas.Phi();
		InvMLdCor = LdFit.Mag();
		PLdFit = TVLdFit.Mag();
		ThLdFit = TVLdFit.Theta();
		PhLdFit = TVLdFit.Phi();
		auto TVPCor = PCor.Vect();
		PPCor = TVPCor.Mag();
		ThPCor = TVPCor.Theta();
		PhPCor = TVPCor.Phi();
		auto TVPi1Cor = Pi1Cor.Vect();
		PPi1Cor = TVPi1Cor.Mag();
		ThPi1Cor = TVPi1Cor.Theta();
		PhPi1Cor = TVPi1Cor.Phi();
		auto VLd = KFLd.GetUnmeasuredCovariance();

		double VPL = VLd(0,0);
		double VThL = VLd(1,1);
		double VPhL = VLd(2,2);


		auto XiFit = LdCor + Pi2Meas;
		FourVectorFitter KFXi(LdCor, Pi2Meas, XiFit);
//		KFXi.UpdateVariance(true);
		double rL = ResPLd * PLdCor;
		double VarianceXi[6] = {VPL, VThL, VPhL, rpi2 * rpi2, ResTh * ResTh, ResPh * ResPh};
		KFXi.SetVariance(VarianceXi);
		TMatrixD CovXi(6,6);
		for(int j=0;j<6;++j){
			for(int k=0;k<6;++k){
				CovXi[j][k] = 0;
			}
		}
		for(int j=0;j<3;++j){
			for(int k=0;k<3;++k){
				if(j==k)continue;
				CovXi[j][k] = VLd[j][k];
			}
		}
		KFXi.AddOffdiagonals(CovXi);
		KFXi.SetInvMass(mXi);
		Chi2Xi = KFXi.DoKinematicFit();
		pvalXi = KFXi.GetPValue();
		auto PullXi = KFXi.GetPull();
		auto contXi = KFXi.GetFittedLV();
		auto LdCorCor = contXi.at(0);
		auto Pi2Cor = contXi.at(1);
		auto XiCor = LdCorCor + Pi2Cor;

		auto TVXiCor = XiCor.Vect();
		PXiCor = TVXiCor.Mag();
		ThXiCor = TVXiCor.Theta();
		PhXiCor = TVXiCor.Phi();
		auto TVPi2Cor = Pi2Cor.Vect();
		PPi2Cor = TVPi2Cor.Mag();
		ThPi2Cor = TVPi2Cor.Theta();
		PhPi2Cor = TVPi2Cor.Phi();
		auto TVLdCorCor = LdCorCor.Vect();
		PLdCorCor = TVLdCorCor.Mag();
		ThLdCorCor = TVLdCorCor.Theta();
		PhLdCorCor = TVLdCorCor.Phi();
		FillHists();
		InvMXiCor = (LdCor + Pi2).Mag();

		auto FIni = KFLd.GetInitialConstraints();
		auto FAfter = KFLd.GetKinematicConstraints();
		double Ipx = FIni.at(0);
		double Ipy = FIni.at(1);
		double Ipz = FIni.at(2);
		double IE = FIni.at(3);
		double Apx = FAfter.at(0);
		double Apy = FAfter.at(1);
		double Apz = FAfter.at(2);
		double AE = FAfter.at(3);
		auto FIniXi = KFXi.GetInitialConstraints();
		auto FAfterXi = KFXi.GetKinematicConstraints();
		for(int j=0;j<4;++j){
			hViolation[j]->Fill(FIni.at(j));
			hViolation[j+4]->Fill(FIniXi.at(j));
			hCorrection[j]->Fill(FAfter.at(j));
			hCorrection[j+4]->Fill(FAfterXi.at(j));
		}
		for (int j = 0; j < 6; ++j)
		{
			hPullLd[j]->Fill(PullLd.at(j));
			hPullXi[j]->Fill(PullXi.at(j));
		}
	}
	TCanvas *cLd = new TCanvas("cLd", "cLd", 1200, 800);
	cLd->Divide(3, 3);
	for (int i = 0; i < 9; ++i)
	{
		cLd->cd(i + 1);
		hLdFitResi[i]->Draw("same");
		hLdMeasResi[i]->Draw("same");
		TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
		leg->AddEntry(hLdMeasResi[i], "Measured", "l");
		leg->AddEntry(hLdFitResi[i], "Fitted", "l");
		leg->Draw();
	}
	TCanvas *cXi = new TCanvas("cXi", "cXi", 1200, 800);
	cXi->Divide(3, 3);
	for (int i = 0; i < 9; ++i)
	{
		cXi->cd(i + 1);
		hXiFitResi[i]->Draw("same");
		hXiMeasResi[i]->Draw("same");
		TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
		leg->AddEntry(hXiMeasResi[i], "Measured", "l");
		leg->AddEntry(hXiFitResi[i], "Fitted", "l");
		leg->Draw();
	}
	TCanvas *cPval = new TCanvas("cPval", "cPval", 1200, 800);
	cPval->Divide(2, 1);
	cPval->cd(1);
	hPvalLd->Draw();
	cPval->cd(2);
	hPvalXi->Draw();
	TCanvas *cPullLd = new TCanvas("cPullLd", "cPullLd", 1200, 800);
	cPullLd->Divide(3, 2);
	for (int i = 0; i < 6; ++i)
	{
		cPullLd->cd(i + 1);
		hPullLd[i]->Draw();
	}
	TCanvas *cPullXi = new TCanvas("cPullXi", "cPullXi", 1200, 800);
	cPullXi->Divide(3, 2);
	for (int i = 0; i < 6; ++i)
	{
		cPullXi->cd(i + 1);
		hPullXi[i]->Draw();
	}
	TCanvas *cViolation = new TCanvas("cViolation", "cViolation", 1200, 800);
	cViolation->Divide(4, 2);
	for (int i = 0; i < 8; ++i)
	{
		cViolation->cd(i + 1);
		hCorrection[i]->Draw("same");
		hViolation[i]->Draw("same");
	}
}
