
double mXi = 1.321;
double mL = 1.115;
double mP = 0.938;
double mK = 0.493;
double mPi = 0.139;
double ResP = 0.1;
double ResPi1 = 0.1;
double ResPi2 = 0.1;
double ResTh = 0.01; 
double ResPh = 0.02; 
double ResThP = 0.02; 
double ResPhP = 0.02; 
double ResThV = 0.05; 
double ResPhV = 0.05; 

double ResPLd = 0.03;
double ResThLd = 0.01;
double ResPhLd = 0.01;

double ScalePP = 1.05;
double ScalePPi1 = 1.05;
double ScalePPi2 = 1.05;
double PLd,ThLd,PhLd; 
double PP,ThP,PhP; 
double PPi1,ThPi1,PhPi1; 
double PPi2,ThPi2,PhPi2; 
double PXi,ThXi,PhXi; 
double PLdMeas,ThLdMeas,PhLdMeas; 
double PLdCorMeas,ThLdCorMeas,PhLdCorMeas; 
double PXiMeas,ThXiMeas,PhXiMeas; 
double PPMeas,ThPMeas,PhPMeas; 
double PPi1Meas,ThPi1Meas,PhPi1Meas; 
double PPi2Meas,ThPi2Meas,PhPi2Meas; 
double PLdCor,ThLdCor,PhLdCor; 
double PLdFit,ThLdFit,PhLdFit; 
double PLdCorCor,ThLdCorCor,PhLdCorCor; 
double PXiCor,ThXiCor,PhXiCor; 
double PPCor,ThPCor,PhPCor; 
double PPi1Cor,ThPi1Cor,PhPi1Cor; 
double PPi2Cor,ThPi2Cor,PhPi2Cor; 
double InvMLd,InvMXi,InvMLdCor,InvMXiCor;
double Chi2,Chi2Xi;
double PullPLd,PullThLd,PullPhLd,PullPP,PullThP,PullPhP,PullPPi1,PullThPi1,PullPhPi1;
double PullPLdCor,PullThLdCor,PullPhLdCor,PullPPi2,PullThPi2,PullPhPi2;
double pvalLd,pvalXi;
int NStep,BestStep;
int Step[200];
double StepChi2[200];
double StepMassDiff[200];
int iev;
TString ElemUnit[3] = {"[GeV/c]","[rad]","[rad]"};
TString ElemLd[9] = {"P_{#Lambda}","#Theta_{#Lambda}","#Phi_{#Lambda}",
"P_{p}","#Theta_{p}","#Phi_{p}",
"P_{#pi1}","#Theta_{#pi1}","#Phi_{#pi1}"
};
TString ElemXi[9] = {"P_{#Xi}","#Theta_{#Xi}","#Phi_{#Xi}",
"P_{#Lambda}","#Theta_{#Lambda}","#Phi_{#Lambda}",
"P_{#pi2}","#Theta_{#pi2}","#Phi_{#pi2}"
};
TH1D* hLdMeasResi[9];
TH1D* hLdFitResi[9];
TH1D* hPullLd[6];
TH1D* hPvalLd;
TH1D* hXiMeasResi[9];
TH1D* hXiFitResi[9];
TH1D* hPullXi[6];
TH1D* hPvalXi;
double range[9] = {0.5,0.3,0.3,0.5,0.3,0.3,0.5,0.3,0.3};
void MakeHists(){
	for(int i=0;i<9;++i){
		hLdMeasResi[i] = new TH1D(ElemLd[i]+"MeasResi",ElemLd[i]+"MeasResi;#Delta"+ElemLd[i]+ElemUnit[i%3],100,-range[i],range[i]);
		hLdFitResi[i] = new TH1D(ElemLd[i]+"FitResi",ElemLd[i]+"FitResi;#Delta"+ElemLd[i]+ElemUnit[i%3],100,-range[i],range[i]);
		hXiMeasResi[i] = new TH1D(ElemXi[i]+"MeasResi",ElemXi[i]+"MeasResi;#Delta"+ElemXi[i]+ElemUnit[i%3],100,-range[i],range[i]);
		hXiFitResi[i] = new TH1D(ElemXi[i]+"FitResi",ElemXi[i]+"FitResi;#Delta"+ElemXi[i]+ElemUnit[i%3],100,-range[i],range[i]);
		hLdMeasResi[i]->SetLineColor(kRed);
		hLdFitResi[i]->SetLineColor(kBlue);
		hXiMeasResi[i]->SetLineColor(kRed);
		hXiFitResi[i]->SetLineColor(kBlue);
	}
	for(int i=0;i<6;++i){
		hPullLd[i] = new TH1D("PullLd"+ElemLd[i],"PullLd"+ElemLd[i]+";Pull",100,-5,5);
		hPullXi[i] = new TH1D("PullXi"+ElemXi[i],"PullXi"+ElemXi[i]+";Pull",100,-5,5);
	}
	hPvalLd = new TH1D("PvalLd","PvalLd;Pval#Lambda",100,0,1);
	hPvalXi = new TH1D("PvalXi","PvalXi;Pval#Xi",100,0,1);
}
void FillHists(){
	hLdMeasResi[0]->Fill((PLdMeas-PLd));
	hLdMeasResi[1]->Fill((ThLdMeas-ThLd));
	hLdMeasResi[2]->Fill((PhLdMeas-PhLd));
	hLdMeasResi[3]->Fill((PPMeas-PP));
	hLdMeasResi[4]->Fill((ThPMeas-ThP));
	hLdMeasResi[5]->Fill((PhPMeas-PhP));
	hLdMeasResi[6]->Fill((PPi1Meas-PPi1));
	hLdMeasResi[7]->Fill((ThPi1Meas-ThPi1));
	hLdMeasResi[8]->Fill((PhPi1Meas-PhPi1));
	hXiMeasResi[0]->Fill((PXiMeas-PXi));
	hXiMeasResi[1]->Fill((ThXiMeas-ThXi));
	hXiMeasResi[2]->Fill((PhXiMeas-PhXi));
	hXiMeasResi[3]->Fill((PLdMeas-PLd));
	hXiMeasResi[4]->Fill((ThLdMeas-ThLd));
	hXiMeasResi[5]->Fill((PhLdMeas-PhLd));
	hXiMeasResi[6]->Fill((PPi2Meas-PPi2));
	hXiMeasResi[7]->Fill((ThPi2Meas-ThPi2));
	hXiMeasResi[8]->Fill((PhPi2Meas-PhPi2));
	hLdFitResi[0]->Fill((PLdFit-PLd));
	hLdFitResi[1]->Fill((ThLdFit-ThLd));
	hLdFitResi[2]->Fill((PhLdFit-PhLd));
	hLdFitResi[3]->Fill((PPCor-PP));
	hLdFitResi[4]->Fill((ThPCor-ThP));
	hLdFitResi[5]->Fill((PhPCor-PhP));
	hLdFitResi[6]->Fill((PPi1Cor-PPi1));
	hLdFitResi[7]->Fill((ThPi1Cor-ThPi1));
	hLdFitResi[8]->Fill((PhPi1Cor-PhPi1));
	hXiFitResi[0]->Fill((PXiCor-PXi));
	hXiFitResi[1]->Fill((ThXiCor-ThXi));
	hXiFitResi[2]->Fill((PhXiCor-PhXi));
	hXiFitResi[3]->Fill((PLdCorCor-PLd));
	hXiFitResi[4]->Fill((ThLdCorCor-ThLd));
	hXiFitResi[5]->Fill((PhLdCorCor-PhLd));
	hXiFitResi[6]->Fill((PPi2Cor-PPi2));
	hXiFitResi[7]->Fill((ThPi2Cor-ThPi2));
	hXiFitResi[8]->Fill((PhPi2Cor-PhPi2));
	hPvalLd->Fill(pvalLd);
	hPvalXi->Fill(pvalXi);
}
