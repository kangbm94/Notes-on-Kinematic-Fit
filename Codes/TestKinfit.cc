#include "src/FourVectorFitter.cc"
#include "src/CascadeFitter.cc"
#include "TestKinfit.hh"
// Code to evaluate the Kineamtic Fit performance for the decay of the Xi to Lambda and pion
// Interaction between 1.8 GeV/c Kaon- and proton target is assumed.
// Author: Kang Byungmin, kangbmw2@naver.com
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

int nev = 100000;
TGenPhaseSpace EvVert;
TGenPhaseSpace EvXi1530;
TGenPhaseSpace EvXi0;
TGenPhaseSpace EvXi;
TGenPhaseSpace EvLd;
TGenPhaseSpace EvLd_Xi0;
double pK18 = 1.8;
TLorentzVector KM(0, 0, pK18, hypot(mK, pK18));
TLorentzVector PT(0, 0, 0, mP);
TLorentzVector Vertex = KM + PT;
double VertMass[2] = {mXi, mK};
double XiDecayMass[2] = {mL, mPi};
double Xi0DecayMass[2] = {mL, mPi0};
double LdDecayMass[2] = {mP, mPi};
double VertMassXi1530[2] = {mXi1530,mK};
double Xi1530_to_XiPi0[2] = {mXi,mPi0};
double Xi1530_to_Xi0Pi[2] = {mXi0,mPi};
void TestXi1530Fit(){
	MakeHistXi();
	MakeHistXi("Xi0");
	TString figdir = "figs/Xi1530Fit/";
	gSystem->mkdir(figdir,1);
	EvVert.SetDecay(Vertex, 2, VertMassXi1530);
	for (int i = 0; i < nev; ++i){
		EvVert.Generate();
		iev = i;
		if (iev % 1000 == 0){
			cout << Form("Event %d", iev) << endl;
		}
		auto KP = *EvVert.GetDecay(1);
		while(KP.CosTheta()<0.97){
			EvVert.Generate();
			KP = *EvVert.GetDecay(1);
		}
		auto Xi1530 = *EvVert.GetDecay(0);
		EvXi1530.SetDecay(Xi1530, 2, Xi1530_to_XiPi0);
		EvXi1530.Generate();
		auto Xi = *EvXi1530.GetDecay(0);
		EvXi.SetDecay(Xi,2,XiDecayMass);
		EvXi.Generate();
		auto Ld = *EvXi.GetDecay(0);
		auto Pi2 = *EvXi.GetDecay(1);
		EvLd.SetDecay(Ld, 2, LdDecayMass);
		EvLd.Generate();
		auto P = *EvLd.GetDecay(0);
		auto Pi1 = *EvLd.GetDecay(1);
		
		auto TVXi = Xi.Vect();
		PXi = TVXi.Mag();
		ThXi = TVXi.Theta();
		PhXi = TVXi.Phi();
		auto TVLd = Ld.Vect();
		PLd = TVLd.Mag();
		ThLd = Ld.Theta();
		PhLd = Ld.Phi();
		auto TVP = P.Vect();
		PP = TVP.Mag();
		ThP = TVP.Theta();
		PhP = TVP.Phi();
		auto TVPi1 = Pi1.Vect();
		PPi1 = TVPi1.Mag();
		ThPi1 = TVPi1.Theta();
		PhPi1 = TVPi1.Phi();
		auto TVPi2 = Pi2.Vect();
		PPi2 = TVPi2.Mag();
		ThPi2 = TVPi2.Theta();
		PhPi2 = TVPi2.Phi();
		ShakeParamsXi();

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
		double rp = PP * ResP;
		double rpi1 = PPi1 * ResPi1;
		double rpi2 = PPi2 * ResPi2;

		auto LdMeas = PMeas + Pi1Meas;
		InvMLd = LdMeas.Mag();
		auto TVLdMeas = LdMeas.Vect();
		PLdMeas = TVLdMeas.Mag();
		ThLdMeas = TVLdMeas.Theta();
		PhLdMeas = TVLdMeas.Phi();
		LdMeas = TLorentzVector(TVLdMeas,hypot(mL,TVLdMeas.Mag()));//Invariant Mass is fixed to physical value
		auto XiMeas = LdMeas + Pi2Meas;
		InvMXi = XiMeas.Mag();
		auto TVXiMeas = XiMeas.Vect();
		PXiMeas = TVXiMeas.Mag();
		ThXiMeas = TVXiMeas.Theta();
		PhXiMeas = TVXiMeas.Phi();

		CascadeFitter KFXi(PMeas, Pi1Meas, Pi2Meas);
		KFXi.SetMaximumStep(5);
//		KFXi.UpdateVariance(true);
//		KFXi.ScaleParameters(false);
		double Variance[9] = {
			rp * rp, ResThP * ResThP, ResPhP * ResPhP,
			rpi1 * rpi1, ResTh * ResTh, ResPh * ResPh,
			rpi2 * rpi2, ResThPi2 * ResThPi2, ResPhPi2 * ResPhPi2
		};//Diagonal components for 9 X 9 matrix
		KFXi.SetInvMass(mL,mXi);
		KFXi.SetVariance(Variance);
		Chi2 = KFXi.DoKinematicFit();
		Pval = KFXi.GetPValue();
		auto Pull = KFXi.GetPull();
		pulls = Pull;
		auto LVCont = KFXi.GetFittedLV();
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

		InvMLdCor = LVLdCor.Mag();
		InvMXiCor = LVXiCor.Mag();
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
		constsAfter = KFXi.GetKinematicConstraints();
		constsIni = KFXi.GetInitialConstraints();
//		if(Pval<0.01)continue;
		AssignXiParam();
		FillHistsXi();
	
		//Mis-identified combinatorial background 
		EvXi1530.SetDecay(Xi1530, 2, Xi1530_to_Xi0Pi);
		EvXi1530.Generate();
		auto Xi0 = *EvXi1530.GetDecay(0);
		Pi2 = *EvXi1530.GetDecay(1); // pi- from Xi1530 is misidentified as pi- from Xi-> L pi- decay.
		EvXi.SetDecay(Xi0,2,Xi0DecayMass);
		EvXi.Generate();
		Ld = *EvXi.GetDecay(0);
		EvLd.SetDecay(Ld, 2, LdDecayMass);
		EvLd.Generate();
		P = *EvLd.GetDecay(0);
		Pi1 = *EvLd.GetDecay(1);

		TVLd = Ld.Vect();
		PLd = TVLd.Mag();
		ThLd = Ld.Theta();
		PhLd = Ld.Phi();

		TVP = P.Vect();
		TVPi1 = Pi1.Vect();
		TVPi2 = Pi2.Vect();
		PP = TVP.Mag();
		ThP = TVP.Theta();
		PhP = TVP.Phi();
		PPi1 = TVPi1.Mag();
		ThPi1 = TVPi1.Theta();
		PhPi1 = TVPi1.Phi();
		PPi2 = TVPi2.Mag();
		ThPi2 = TVPi2.Theta();
		PhPi2 = TVPi2.Phi();

		ShakeParamsXi();
		
		TVPMeas = TVector3(0, 0, PPMeas);
		TVPMeas.SetTheta(ThPMeas);
		TVPMeas.SetPhi(PhPMeas);
		TVPi1Meas = TVector3(0, 0, PPi1Meas);
		TVPi1Meas.SetTheta(ThPi1Meas);
		TVPi1Meas.SetPhi(PhPi1Meas);
		TVPi2Meas = TVector3(0, 0, PPi2Meas);
		TVPi2Meas.SetTheta(ThPi2Meas);
		TVPi2Meas.SetPhi(PhPi2Meas);
		PMeas = TLorentzVector(TVPMeas, hypot(PPMeas, mP));
		Pi1Meas = TLorentzVector(TVPi1Meas, hypot(PPi1Meas, mPi));
		Pi2Meas = TLorentzVector(TVPi2Meas, hypot(PPi2Meas, mPi));
		rp = PP * ResP;
		rpi1 = PPi1 * ResPi1;
		rpi2 = PPi2 * ResPi2;

		LdMeas = PMeas + Pi1Meas;
		InvMLd = LdMeas.Mag();
		TVLdMeas = LdMeas.Vect();
		PLdMeas = TVLdMeas.Mag();
		ThLdMeas = TVLdMeas.Theta();
		PhLdMeas = TVLdMeas.Phi();
		LdMeas = TLorentzVector(TVLdMeas,hypot(mL,TVLdMeas.Mag()));
		XiMeas = LdMeas + Pi2Meas;
		InvMXi = XiMeas.Mag();
		TVXiMeas = XiMeas.Vect();
		TVXi = Xi.Vect();
		PXiMeas = TVXiMeas.Mag();
		ThXiMeas = TVXiMeas.Theta();
		PhXiMeas = TVXiMeas.Phi();
		PXi = TVXi.Mag();
		ThXi = TVXi.Theta();
		PhXi = TVXi.Phi();
		CascadeFitter KFXi0(PMeas, Pi1Meas, Pi2Meas);
		KFXi0.SetMaximumStep(5);
		KFXi0.UpdateVariance(true);
		double VarianceXi0[9] = {
			rp * rp, ResThP * ResThP, ResPhP * ResPhP,
			rpi1 * rpi1, ResTh * ResTh, ResPh * ResPh,
			rpi2 * rpi2, ResThPi2 * ResThPi2, ResPhPi2 * ResPhPi2
		};
		KFXi0.SetInvMass(mL,mXi);
		KFXi0.SetVariance(VarianceXi0);
		Chi2 = KFXi0.DoKinematicFit();
		Pval = KFXi0.GetPValue();
		Pull = KFXi0.GetPull();
		pulls = Pull;
		LVCont = KFXi.GetFittedLV();
		LVPCor = LVCont.at(0);
		LVPi1Cor = LVCont.at(1);
		LVPi2Cor = LVCont.at(2);
		LVLdCor = LVCont.at(3);
		LVXiCor = LVCont.at(4);

		TVPCor = LVPCor.Vect();
		TVPi1Cor = LVPi1Cor.Vect();
		TVPi2Cor = LVPi2Cor.Vect();
		TVLdCor = LVLdCor.Vect();
		TVXiCor = LVXiCor.Vect();

		InvMLdCor = LVLdCor.Mag();
		InvMXiCor = LVXiCor.Mag();

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
		constsAfter = KFXi.GetKinematicConstraints();
		constsIni = KFXi.GetInitialConstraints();
		AssignXiParam("Xi0");
		FillHistsXi("Xi0");
	}
	TString key;
	vector<TString> mothers = {"Xi_","Xi0"};	
	for(auto mother:mothers){
		TCanvas* cKF0 = new TCanvas(mother +"cKF0",mother +"Chi2_Pull",1200,800);
		cKF0->Divide(2,1);
		cKF0->cd(1);
		hMap[mother + "Chi2"]->Draw();
		cKF0->cd(2);
//		hMap[mother + "Pval"]->Draw();
//		hMap[mother + "Pval"]->GetYaxis()->SetRangeUser(0,1.2*hMap[mother + "Pval"]->GetMaximum());
		auto hpval = (TH1D*)hMap[mother + "Pval"]->Clone();
		hpval->Draw();
		hpval->GetYaxis()->SetRangeUser(0,1.2*hMap[mother + "Pval"]->GetMaximum());
		cKF0->SaveAs(figdir + cKF0->GetName() + ".pdf");
		TCanvas* cKF1 = new TCanvas(mother +"cKF1",mother +"DaughterResidual",1200,800);
		cKF1->Divide(3,3);
		TLegend* leg_kf = new TLegend(0.6,0.7,0.9,0.9);
		leg_kf->SetBorderSize(0);
		leg_kf->SetFillStyle(0);
		for(int i=0;i<9;++i){
			cKF1->cd(i+1);
			key = mother + Form("%s",Title[i].Data()) + "Cor";
			auto h1 = hMap[key];
			key = mother + Form("%s",Title[i].Data());
			auto h2 = hMap[key];
			if(i == 0){
				leg_kf->AddEntry(h2,"Before KF","l");
				leg_kf->AddEntry(h1,"After KF","l");
			}
			double maxi;
			if(h1->GetMaximum()>h2->GetMaximum()){
				h1->Draw();
				h2->Draw("same");
			}
			else{
				h2->Draw();
				h1->Draw("same");
			}
		}
		cKF1->cd(1);
		leg_kf->Draw();
		cout<<"DrawingLeg"<<endl;
		cKF1->SaveAs(figdir + cKF1->GetName() + ".pdf");
		TCanvas* cKF2 = new TCanvas(mother +"cKF2",mother +"MotherResidual",1200,800);
		cKF2->Divide(3,2);
		for(int i=0;i<6;++i){
			cKF2->cd(i+1);
			key = mother + Form("%s",Title[i+9].Data()) + "Cor";
			auto h1 = hMap[key];
			key = mother + Form("%s",Title[i+9].Data());
			auto h2 = hMap[key];
//			(i == 0) ? leg_kf->Draw() : void();
			if(h1->GetMaximum()>h2->GetMaximum()){
				h1->Draw();
				h2->Draw("same");
			}
			else{
				h2->Draw();
				h1->Draw("same");
			}
		}
		cKF2->cd(1);
		leg_kf->Draw();
		cout<<"DrawingLeg"<<endl;
		cKF2->SaveAs(figdir + cKF2->GetName() + ".pdf");
		TCanvas* cKF3 = new TCanvas(mother +"cKF3",mother +"InvMass",1200,800);
		cKF3->Divide(2,1);
		cKF3->cd(1);
		hMap[mother + "InvMLdCor"]->Draw();
		hMap[mother + "InvMLd"]->Draw("same");
		cKF3->cd(2);
		hMap[mother + "InvMXiCor"]->Draw();
		hMap[mother + "InvMXi"]->Draw("same");
		cKF3->SaveAs(figdir + cKF3->GetName() + ".pdf");
		TCanvas* cKF4 = new TCanvas(mother +"cKF4",mother +"Pull",1200,800);
		cKF4->Divide(3,3);
		for(int i=0;i<9;++i){
			cKF4->cd(i+1);
			key = mother + "Pull_"+ Form("%s",Title[i].Data());
			hMap[key]->Draw();
		}	
		cKF4->SaveAs(figdir + cKF4->GetName() + ".pdf");
		TCanvas* cKF5 = new TCanvas(mother +"cKF5",mother +"Constraints",1200,800);
		cKF5->Divide(4,2);
		for(int i=0;i<8;++i){
			cKF5->cd(i+1);
			key = mother + "Consts_"+ Form("%s",Title[i].Data());
			hMap[key]->Draw();
			key = mother + "IniConsts_"+ Form("%s",Title[i].Data());
			hMap[key]->Draw("same");
		}	
		cKF5->SaveAs(figdir + cKF5->GetName() + ".pdf");
	}
	{
		cout<<"Summary"<<endl;
		TCanvas* c_sum = new TCanvas("c_sum","Events_sum",1200,800);
		c_sum->Divide(2,1);
		c_sum->cd(1);
		TH1D* hSumChi2 = (TH1D*)hMap["Xi_Chi2"]->Clone();
		hSumChi2->Add(hMap["Xi0Chi2"]);
		hSumChi2->Draw("hist");
		hSumChi2->SetLineColor(kBlack);
		TH1D* hXiChi2 = (TH1D*)hMap["Xi_Chi2"]->Clone();
		hXiChi2->SetLineColor(kAzure);
		hXiChi2->Draw("same");
		TH1D* hXi0Chi2 = (TH1D*)hMap["Xi0Chi2"]->Clone();
		hXi0Chi2->SetLineColor(kRed);
		hXi0Chi2->Draw("same");
		c_sum->cd(2);
		gPad->SetLogy();
		TH1D* hSumPval = (TH1D*)hMap["Xi_Pval"]->Clone();
		hSumPval->Add(hMap["Xi0Pval"]);
		hSumPval->Draw("hist");
		hSumPval->SetLineColor(kBlack);
		double maxi_pval = hSumPval->GetMaximum();
		hSumPval->GetYaxis()->SetRangeUser(1,1.5*maxi_pval);
		TH1D* hXiPval = (TH1D*)hMap["Xi_Pval"]->Clone();
		hXiPval->SetLineColor(kAzure);
		hXiPval->Draw("same");
		TH1D* hXi0Pval = (TH1D*)hMap["Xi0Pval"]->Clone();
		hXi0Pval->SetLineColor(kRed);
		hXi0Pval->Draw("same");
		TLegend* legs = new TLegend(0.6,0.7,0.9,0.9);
		legs->AddEntry(hSumPval,"All Events","l");
		legs->AddEntry(hXiPval,"#Xi Events","l");
		legs->AddEntry(hXi0Pval,"#Xi^{0} Events","l");
		legs->SetBorderSize(0);
		legs->SetFillStyle(0);
		legs->Draw();
		c_sum->SaveAs(figdir + c_sum->GetName() + ".pdf");
		TCanvas* c_imxi = new TCanvas("c_imxi","XiIM",1200,800);
		TH1D* hSumInvMXi = (TH1D*)hMap["Xi_InvMXi"]->Clone();
		hSumInvMXi->Add(hMap["Xi0InvMXi"]);
		hSumInvMXi->Draw("hist");
		hSumInvMXi->SetLineColor(kBlack);
		TH1D* hXiInvMXi = (TH1D*)hMap["Xi_InvMXi"]->Clone();
		hXiInvMXi->SetLineColor(kAzure);
		hXiInvMXi->Draw("same");
		TH1D* hXi0InvMXi = (TH1D*)hMap["Xi0InvMXi"]->Clone();
		hXi0InvMXi->SetLineColor(kRed);
		hXi0InvMXi->Draw("same");
		legs->Draw();
		c_imxi->SaveAs(figdir + c_imxi->GetName() + ".pdf");
		double p_cuts[1000];
		double p_purities[1000];
		double p_efficiencies[1000];
		for(int ib=1;ib<= hSumPval->GetNbinsX();++ib){
			double total = 0;
			double nxi = hXiPval->GetEntries();
			double cnt_xi=0,cnt_xi0=0;
			double p_cut = ib*1./(hSumPval->GetNbinsX());
			for(int i=ib;i <= hSumPval->GetNbinsX();++i){
				total += hSumPval->GetBinContent(i);
				cnt_xi += hXiPval->GetBinContent(i);
			}
			double purity = cnt_xi / total;
			double efficiency = cnt_xi / nxi;
			p_cuts[ib-1] = p_cut;
			p_purities[ib-1] = purity;
			p_efficiencies[ib-1] = efficiency;
		}
		double m_cuts[500];
		double m_purities[1000];
		double m_efficiencies[1000];
		double dm = 0.00015;
		for(int ib=1;ib<= (hSumInvMXi->GetNbinsX())/2;++ib){
			double total = 0;
			double nxi = hXiInvMXi->GetEntries();
			double cnt_xi=0,cnt_xi0=0;
			double bw = (ib+0.1) * 2 * dm;
			for(int i = 1; i <= hSumInvMXi->GetNbinsX();++i){
				if( abs(hXiInvMXi -> GetBinCenter(i)-mXi)< bw){
					cnt_xi+= hXiInvMXi -> GetBinContent(i);
					total += hSumInvMXi -> GetBinContent(i);
				}
			}
			double purity = cnt_xi / total;
			double efficiency = cnt_xi / nxi;
			m_cuts[ib-1] = bw;
			m_purities[ib-1] = purity;
			m_efficiencies[ib-1] = efficiency;
		}
		TCanvas* c_pval = new TCanvas("c_pval","c_pval",1200,800);
		TGraph* gPurityPval = new TGraph(hSumPval->GetNbinsX(),p_cuts,p_purities);
		gPurityPval->SetLineWidth(4);
		gPurityPval->SetLineColor(kRed);
		gPurityPval->GetYaxis()->SetLimits(0,1);
		gPurityPval->GetYaxis()->SetRangeUser(0,1);
		gPurityPval->GetXaxis()->SetTitle("P-value cut");
		gPurityPval->Draw("AL");	
		TGraph* gEfficiencyPval = new TGraph(hSumPval->GetNbinsX(),p_cuts,p_efficiencies);
		gEfficiencyPval->SetLineWidth(4);
		gEfficiencyPval->SetLineColor(kAzure);
		gEfficiencyPval->Draw("Lsame");	
		TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.8);
		leg1->AddEntry(gPurityPval,"Purity","l");
		leg1->AddEntry(gEfficiencyPval,"Efficiency","l");
		leg1->SetBorderSize(0);
		leg1->SetFillStyle(0);
		leg1->Draw();
		c_pval->SaveAs(figdir + c_pval->GetName() + ".pdf");
		TCanvas* c_mass = new TCanvas("c_mass","c_mass",1200,800);
		TGraph* gPurityMass = new TGraph(hSumInvMXi->GetNbinsX()/2,m_cuts,m_purities);
		gPurityMass->GetYaxis()->SetLimits(0,1);
		gPurityMass->GetYaxis()->SetRangeUser(0,1);
		gPurityMass->SetLineWidth(4);
		gPurityMass->SetLineColor(kRed);
		gPurityMass->GetYaxis()->SetRange(0,1);
		gPurityMass->GetXaxis()->SetTitle("Mass Window [GeV/c^{2}]");
		gPurityMass->Draw("ALsame");	
		TGraph* gEfficiencyMass = new TGraph(hSumInvMXi->GetNbinsX()/2,m_cuts,m_efficiencies);
		gEfficiencyMass->SetLineWidth(4);
		gEfficiencyMass->SetLineColor(kAzure);
		gEfficiencyMass->Draw("Lsame");	
		leg1->Draw();
		c_mass->SaveAs(figdir + c_mass->GetName() + ".pdf");

		TCanvas* c_eff_pur = new TCanvas("c_eff_pur","c_eff_pur",1200,800);
		TGraph* gEffPurPval = new TGraph(hSumPval->GetNbinsX(),p_efficiencies,p_purities);
		gEffPurPval->SetLineWidth(4);
		gEffPurPval->SetLineColor(kAzure);
		gEffPurPval->GetXaxis()->SetLimits(0,1);
		gEffPurPval->GetXaxis()->SetRangeUser(0,1);
		gEffPurPval->GetXaxis()->SetTitle("Efficiency");
		gEffPurPval->GetYaxis()->SetLimits(0,1);
		gEffPurPval->GetYaxis()->SetRangeUser(0,1);
		gEffPurPval->GetYaxis()->SetTitle("Purity");
		gEffPurPval->Draw("AL");
		TGraph* gEffPurMass = new TGraph(hSumInvMXi->GetNbinsX()/2,m_efficiencies,m_purities);
		gEffPurMass->SetLineWidth(4);
		gEffPurMass->SetLineColor(kRed);
		gEffPurMass->Draw("Lsame");
		TLegend* leg2 = new TLegend(0.6,0.4,0.9,0.6);
		leg2->AddEntry(gEffPurPval,"P-value cut","l");
		leg2->AddEntry(gEffPurMass,"Mass cut","l");
		leg2->SetBorderSize(0);
		leg2->SetFillStyle(0);
		leg2->Draw();
		c_eff_pur->SaveAs(figdir + c_eff_pur->GetName() + ".pdf");
		// Save the canvases
	}
}
void TestXiFit()
{
	MakeHistXi();
	EvVert.SetDecay(Vertex, 2, VertMass);
	for (int i = 0; i < nev; ++i){
		EvVert.Generate();
		iev = i;
		if (iev % 1000 == 0){
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
		CascadeFitter KFXi(PMeas, Pi1Meas, Pi2Meas);
//		CascadeFstter KFXi(P, Pi1, Pi2);
		KFXi.SetMaximumStep(5);
		KFXi.UpdateVariance(true);
//		KFXi.ScaleParameters(false);
		double Variance[9] = {
			rp * rp, ResThP * ResThP, ResPhP * ResPhP,
			rpi1 * rpi1, ResTh * ResTh, ResPh * ResPh,
			rpi2 * rpi2, ResThPi2 * ResThPi2, ResPhPi2 * ResPhPi2
		};
		KFXi.SetInvMass(mL,mXi);
		KFXi.SetVariance(Variance);
		Chi2 = KFXi.DoKinematicFit();
		double Pval = KFXi.GetPValue();
		auto Pull = KFXi.GetPull();
		pulls = Pull;
		auto LVCont = KFXi.GetFittedLV();
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
		constsAfter = KFXi.GetKinematicConstraints();
		constsIni = KFXi.GetInitialConstraints();
//		if(Pval<0.01)continue;
		FillHistXi();






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
	EvVert.SetDecay(Vertex, 2, VertMass);
	for (int i = 0; i < nev; ++i)
	{
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
void TestKinfit()
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	SetStyle();
//	MakeHists();
	cout << "TestXi1530Fit()" << endl;
	cout << "TestXiFit()" << endl;
	cout << "TestFourVectorFit()" << endl;
	TestXi1530Fit();
}
