#include "src/FourVectorFitter.cc"
#include "TestKinfit.hh"
// Code to evaluate the Kineamtic Fit performance for the decay of the Xi to Lambda and pion
// Interaction between 1.8 GeV/c Kaon- and proton target is assumed.
// Author: Kang Byungmin, kangbmw2@naver.com
// https://github.com/kangbm94/Notes-on-Kinematic-Fit


int nev = 10000;
void TestKinfit(){
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	double pK18 = 1.8;
	TLorentzVector KM(0,0,pK18,hypot(mK,pK18));
	TLorentzVector PT(0,0,0,mP);
	TLorentzVector Vertex = KM + PT;
	double VertMass[2] = {mXi,mK};
	double XiDecayMass[2] = {mL,mPi};
	double LdDecayMass[2] = {mP,mPi};
	TGenPhaseSpace EvVert;
	TGenPhaseSpace EvXi;
	TGenPhaseSpace EvLd;
	MakeHists();
	for(int i=0;i<nev;++i){
		EvVert.SetDecay(Vertex,2,VertMass);
		EvVert.Generate();
		iev = i;
		if(iev%1000==0){
			cout<<Form("Event %d",iev)<<endl;
		}
		auto Xi = *EvVert.GetDecay(0);
		auto KP = *EvVert.GetDecay(1);
		EvXi.SetDecay(Xi,2,XiDecayMass);
		EvXi.Generate();
		auto Ld = *EvXi.GetDecay(0);
		auto Pi2 = *EvXi.GetDecay(1);
		EvLd.SetDecay(Ld,2,LdDecayMass);
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

		ThLdMeas = gRandom->Gaus(ThLd,ResThV);
		PhLdMeas = gRandom->Gaus(PhLd,ResPhV);
		TVector3 V1 = TVLd;
		V1.SetTheta(ThLdMeas);
		V1.SetPhi(PhLdMeas);
		TVector3 V2(0,0,0);
		PPMeas = gRandom->Gaus(PP,PP*ResP);
		ThPMeas = gRandom->Gaus(ThP,ResThP);
		PhPMeas = gRandom->Gaus(PhP,ResPhP);
		PPi1Meas = gRandom->Gaus(PPi1,PPi1*ResPi1);
		ThPi1Meas = gRandom->Gaus(ThPi1,ResTh);
		PhPi1Meas = gRandom->Gaus(PhPi1,ResPh);
		PPi2Meas = gRandom->Gaus(PPi2,PPi2*ResPi2);
		ThPi2Meas = gRandom->Gaus(ThPi2,ResTh);
		PhPi2Meas = gRandom->Gaus(PhPi2,ResPh);
		
		TVector3 TVPMeas(0,0,PPMeas);
		TVPMeas.SetTheta(ThPMeas);
		TVPMeas.SetPhi(PhPMeas);
		TVector3 TVPi1Meas(0,0,PPi1Meas);
		TVPi1Meas.SetTheta(ThPi1Meas);
		TVPi1Meas.SetPhi(PhPi1Meas);
		TVector3 TVPi2Meas(0,0,PPi2Meas);
		TVPi2Meas.SetTheta(ThPi2Meas);
		TVPi2Meas.SetPhi(PhPi2Meas);
		TLorentzVector PMeas(TVPMeas,hypot(PPMeas,mP));
		TLorentzVector Pi1Meas(TVPi1Meas,hypot(PPi1Meas,mPi));
		TLorentzVector Pi2Meas(TVPi2Meas,hypot(PPi2Meas,mPi));
		double rp = PPMeas*ResP;
		double rpi1 = PPi1Meas*ResPi1;
		double rpi2 = PPi2Meas*ResPi2;

		auto LdRecon = PMeas + Pi1Meas;
		auto LdReconMeas = LdRecon;
		LdReconMeas.SetTheta(ThLdMeas);
		LdReconMeas.SetPhi(PhLdMeas);
		


		auto XiRecon = LdRecon + Pi2Meas;
		InvMLd = LdRecon.Mag();
		InvMXi = (LdRecon+Pi2).Mag();
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


		FourVectorFitter KFLd(PMeas,Pi1Meas,LdRecon);
		double Variance[6] = {rp*rp,ResTh*ResTh,ResPh*ResPh,rpi1*rpi1,ResTh*ResTh,ResPh*ResPh};
		KFLd.SetInvMass(mL);
		KFLd.SetVariance(Variance);
		Chi2=KFLd.DoKinematicFit();
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
		auto LdFit = PCor+Pi1Cor;
		auto TVLdFit = LdFit.Vect();
		auto TVLdMeas = TVPMeas + TVPi1Meas;
		PLdMeas = TVLdMeas.Mag();
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
		
		auto XiFit = LdCor + Pi2Meas;
		FourVectorFitter KFXi(LdCor,Pi2Meas,XiFit);
		double rL = ResPLd*PLdCor;
		double VarianceXi[6] = {rp*rp,ResThLd*ResThLd,ResPhLd*ResPhLd,rpi2*rpi2,ResTh*ResTh,ResPh*ResPh};
		KFXi.SetVariance(VarianceXi);
		KFXi.SetInvMass(mXi);	
		Chi2Xi=KFXi.DoKinematicFit();
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
		double Ipx=FIni.at(0);
		double Ipy=FIni.at(1);
		double Ipz=FIni.at(2);
		double IE=FIni.at(3);
		double Apx=FAfter.at(0);
		double Apy=FAfter.at(1);
		double Apz=FAfter.at(2);
		double AE=FAfter.at(3);
	}
	TCanvas* cLd = new TCanvas("cLd","cLd",1200,800);
	cLd->Divide(3,3);
	for(int i=0;i<9;++i){
		cLd->cd(i+1);
		hLdFitResi[i]->Draw("same");
		hLdMeasResi[i]->Draw("same");
		TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
		leg->AddEntry(hLdMeasResi[i],"Measured","l");
		leg->AddEntry(hLdFitResi[i],"Fitted","l");
		leg->Draw();
	}
	TCanvas* cXi = new TCanvas("cXi","cXi",1200,800);
	cXi->Divide(3,3);
	for(int i=0;i<9;++i){
		cXi->cd(i+1);
		hXiFitResi[i]->Draw("same");
		hXiMeasResi[i]->Draw("same");
		TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
		leg->AddEntry(hXiMeasResi[i],"Measured","l");
		leg->AddEntry(hXiFitResi[i],"Fitted","l");
		leg->Draw();
	}
	TCanvas* cPval = new TCanvas("cPval","cPval",1200,800);
	cPval->Divide(2,1);
	cPval->cd(1);
	hPvalLd->Draw();
	cPval->cd(2);
	hPvalXi->Draw();
}

