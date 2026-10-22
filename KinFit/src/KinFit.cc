#include "../include/KinFit.hh"
#include "TMatrixDEigen.h"
#include "TString.h"
#include <iostream>
#ifndef KinFit_cc
#define KinFit_cc
#define Debug 0
#define ShowChi2 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
using namespace std;
void
KinematicFitter::SetVariance(double* var){
#if Debug
	cout<<"KinematicFitter::SetVariance"<<endl;
#endif
	double variance[400];
	double varianceInv[400];
	TMatrixD ScaleUp(nMeas,nMeas);
	TMatrixD ScaleDn(nMeas,nMeas);
	for(int i = 0;i<nMeas*nMeas;++i){
		variance[i]=0;
		varianceInv[i]=0;
	}
	for(int i = 0;i<nMeas;++i){
#if Debug
		std::cout<<"Variance["<<i<<"] = "<<var[i]<<std::endl;
#endif
		variance[i+i*nMeas] = var[i];
		varianceInv[i+i*nMeas] = 1./var[i];
		if(ScaleParams){
			ScaleUp(i,i) = sqrt(1./var[i]);
			ScaleDn(i,i) = sqrt(var[i]);
			variance[i+i*nMeas] = 1.;
			varianceInv[i+i*nMeas] = 1.;
		}
		else{
			ScaleUp(i,i) = 1.; 
			ScaleDn(i,i) = 1.; 
		}
	}

	TMatrixD Variance(nMeas,nMeas,variance);
	TMatrixD VarianceInv(nMeas,nMeas,varianceInv);
	TMatrixD ZeroMat(nMeas,nMeas);
	Variancies.push_back(Variance);
	VarianceInvs.push_back(VarianceInv);
	dVMats.push_back(ZeroMat);
	ScalingMats.push_back(ScaleUp);
	ScalingMats.push_back(ScaleDn);

	if(nUnkn > 0){
		TMatrixD ZeroMatU(nUnkn,nUnkn);
		VarianciesU.push_back(ZeroMatU);
	}
};
void
KinematicFitter::AddOffdiagonals(TMatrixD Cov){
	if(ScaleParams){
		auto ScaleUp = ScalingMats.at(0);
		Cov = Cov * ScaleUp;
		Cov = ScaleUp * Cov;
	}
	for(int ir = 0; ir < Cov.GetNrows();++ir){
		for(int ic = 0; ic < Cov.GetNcols();++ic){
			if(ir == ic) continue;
			Variancies.at(0)(ir,ic) += Cov(ir,ic);
		}
	}
#if Debug
	if(ScaleParams){
		cout<<"Covariance Matrix After scaling";
		Variancies.at(0).Print();
	}
	else{
		cout<<"Covariance Matrix";
		Variancies.at(0).Print();
	}
#endif
	bool CheckPositiveDefinite = false;
	int nitr = 0;
	while(!CheckPositiveDefinite and nitr < 5){
		bool PositiveDefinate = true;
		TMatrixDEigen VEig(Variancies.at(0));
		TVectorD eigenvalue = VEig.GetEigenValuesRe();
		for(int ie=0;ie<eigenvalue.GetNrows();++ie){
			if(eigenvalue(ie)<=0){
				PositiveDefinate = false;
				break;
			}
		}
		CheckPositiveDefinite = PositiveDefinate;
		if(!PositiveDefinate){
			cout<<"KinematicFit:: Variance not positive-definite"<<endl;
			cout<<"Reducing Offdiagonals..."<<endl;
			for(int irow=0;irow < Variancies.at(0).GetNrows();++irow){
				for(int icol=0;icol < Variancies.at(0).GetNcols();++icol){
					if(irow == icol) continue;
					Variancies.at(0)(irow,icol) *=0.9;
				}
			}
		}
		nitr++;
	}
#if Debug
	cout<<"Covariance, det ="<<Variancies.at(0).Determinant();
	Variancies.at(0).Print();
#endif
}
void KinematicFitter::ProcessStep(){
	if(nUnkn > 0){
		ProcessStepWUnkn();
	}
	else{
		ProcessStepWOUnkn();
	}
}
void KinematicFitter::ProcessStepWUnkn(){
	#if Debug
	std::cout<<"KinematicFitter::ProcessStep"<<std::endl;
	std::cout<<"Loading Varaibles"<<std::endl;
	std::cout<<"nMeas = "<<Measurements.size()<<std::endl;
	std::cout<<"nUnkn = "<<Unknowns.size()<<std::endl;
	std::cout<<"ScaleUp = "<<ScalingMats.size()<<std::endl;
	#endif
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step); 
	auto Meas0 = Measurements.at(0);
	auto Unkn0 = Unknowns.at(0);
	auto ScaleUp = ScalingMats.at(0);
	auto ScaleDn = ScalingMats.at(1);
	auto MS0 = ScaleUp * Meas0;// MS = M/sigma 
	auto MS = ScaleUp * Meas;
	/*
	if(UpdateVariancies and step ){
		Meas0 = Measurements.at(step-1);	
		Unkn0 = Unknowns.at(step-1);	
	}
	*/
	auto VMat = Variancies.at(step);
	auto VInv = VMat;
#if Debug
	std::cout<<"Inverting VMat"<<std::endl;
#endif
	VInv.SetTol(1e-26);
	if(VInv.Determinant() == 0){
		cout<<"Warning: Variance Matrix is singular!"<<endl;
	}
	VInv.Invert();
#if Debug
	std::cout<<"Setting Constraints"<<std::endl;
#endif

	SetConstraints();
	TMatrixD FMat = FMats.at(step);
	TMatrixD dFdM = dFdMs.at(step);
	TMatrixD dFdMS = dFdM*ScaleDn;
#if Debug
	if(step == 0){
		cout<<"Covariance, det ="<<VMat.Determinant();
		VMat.Print();
	}
#endif
	
#if Debug>1
	cout<<"dFdM";
	dFdM.Print();
	cout<<"dFdMScaled";
	dFdMS.Print();
#endif
	TMatrixD dFdU = dFdUs.at(step);
//	vector<TMatrixD> d2Fd2U = d2Fd2Us.at(step);
	double det_dFdM=0,det_dFdU=0;
	auto dFdMT = TransposeMatrix(dFdMS);
	auto dFdUT = TransposeMatrix(dFdU);
	auto rMat = FMat + dFdMS*(MS0-MS);
	auto sMat =dFdMS*VMat*dFdMT;
	auto sInv = sMat;
	if(sInv.Determinant() == 0){
		cout<<"Warning: S Matrix is singular!"<<endl;
	}
	sInv.Invert();
	auto FuSIFu =	dFdUT*sInv*dFdU;
	if(FuSIFu.Determinant() == 0){
		cout<<"Warning: FuSIFu Matrix is singular!"<<endl;
	}
	FuSIFu.Invert();
	auto dU = (FuSIFu) * (dFdUT* (sInv) * rMat) ;
	dU = dU -dU - dU;
	auto Unkn_next = Unkn + dU;

	auto Lambda = (sInv* (rMat+ dFdU*dU));
	auto LambdaT = TransposeMatrix(Lambda);
	auto dM =  VMat*dFdMT*Lambda;
	auto Meas_next = MS0 - dM; 
	
	auto dMT = TransposeMatrix(dM);
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(ScaleDn*Meas_next);	
	rMats.push_back(rMat);
	sMats.push_back(sMat);
	auto GMat = dFdMT*sInv*dFdMS;
	auto HMat = dFdMT*sInv*dFdU;
	auto UMat = dFdUT*sInv*dFdU;
	if(UMat.Determinant() == 0){
		cout<<"Warning: UMat is singular!"<<endl;
	}
	UMat.Invert();

	auto HMatT = TransposeMatrix(HMat);
	auto HUH = (HMat*UMat*HMatT);
	auto dV = VMat*(GMat -(HUH) )*VMat*(GMat - HUH)*VMat;
	dVMats.push_back(dV);
	auto VMat_next = VMat- VMat * (GMat - HUH)*VMat - VMat * (GMat - HUH)*VMat + dV;
	auto VInv_next = VMat_next;
	VInv_next.SetTol(1e-26);
	//VInv_next.Invert();
	
	double Chi2 = (dMT* (VInv)*dM)(0,0);// + 2 * (LambdaT * FMat )(0,0);
	Chi2s.push_back(Chi2);
	TMatrixD CHI2_U = TMatrixD(nUnkn,nUnkn);
	for(int iu=0;iu<nUnkn;++iu){
		CHI2_U(iu,iu)=Chi2;
	}
#if Debug
	TString StepIndi = Form("[Step::%d]",step);
	cout<<StepIndi<<"##########Constraint Matrix##########"<<endl;
	cout<<"F Mat";
	FMat.Print();
	cout<<"dFdM Mat";
	dFdMS.Print();
	cout<<"dFdU Mat";
	dFdU.Print();

	std::cout<<"s Mat, Det = "<<sMat.Determinant();
	sMat.Print();
	std::cout<<"S-1 Mat :";
	sInv.Print();
	std::cout<<"FuS-1Fu-1 :";
	FuSIFu.Print();
	std::cout<<"FuS-1Fu-1 dFdU S-1 Mat :";
	(FuSIFu * dFdUT * sInv).Print();
	std::cout<<"r Mat :";
	rMat.Print();
	std::cout<<"Unkn Mat :";
	Unkn.Print();
	std::cout<<"dU Mat :";
	dU.Print();
	std::cout<<"Unkn_next :";
	Unkn_next.Print();

	std::cout<<"dFdMT Mat :";
	dFdMT.Print();
	std::cout<<"Lambda Mat :";
	Lambda.Print();
	std::cout<<"VMat :";
	VMat.Print();
	std::cout<<"Meas Mat :";
	Meas.Print();
	std::cout<<"dM Mat :";
	dM.Print();
	std::cout<<"Meas_next :";
	Meas_next.Print();

#endif 
#if Debug>1
	cout<<StepIndi<<"##########Variance Matrix##########"<<endl;
	cout<<"V Mat : Determinant = "<<VMat.Determinant();
	VMat.Print();
	cout<<"V Mat_next : Determinant = "<<VMat_next.Determinant();
	VMat_next.Print();
#if UseHessian
	cout<<"U Cov";
	UVMat.Print();
#endif
	cout<<"V*VInv : Determinant = "<<(VMat*VInv).Determinant();
	(VMat*VInv).Print();
	cout<<"V*VInv_next : Determinant = "<<(VMat_next*VInv_next).Determinant();
	(VMat_next*VInv_next).Print();
	cout<<"dV : Determinant "<<dV.Determinant();
	dV.Print();
	TCanvas*canv = new TCanvas("canv","canv",600,600);
	gStyle->SetOptStat(0);
	TH2D* HdV = new TH2D("HDV","HDV",8,0,8,8,0,8);
	for(int ir=0;ir<nMeas;++ir){
		for(int il=0;il<nMeas;++il){
			HdV->SetBinContent(ir+1,8-il,dV(ir,il));
		}
	}
	cout<<"Chi2 = "<<Chi2<<endl;
	for(auto pul:Pull){
		cout<<Form("Pull %d = %g ,",ip,pul);
		ip++;
	}
	cout<<endl;
		

	
	HdV->Draw("colz");
	gPad->SetMargin(0.1,0.2,0.1,0.1);
	canv->Modified();
	canv->Update();
	gSystem->ProcessEvents();
	cin.ignore();
#endif
	if(UpdateVariancies and step < 1){
	}else{
		VMat_next = VMat;
		VInv_next = VInv;
	}
	if(VMat_next.Determinant() == 0){
		VMat_next = VMat;
		VInv_next = VInv;
	}
//	VMat_next = VMat - dV;
	Variancies.push_back(VMat_next);
	VarianceInvs.push_back(VInv_next);
	SampleStepPoint(step);
	step++;
//	CalcVariance(step);	
	int ip = 0;
	vector<double>Pull;
	for(int i = 0; i<nMeas;++i){
		double dm = dM(i,0); 	
		double dv = dV(i,i);
		dv = sqrt(dv);
		Pull.push_back( dm / dv);
	}
	vector<double>UPull;
//	VarianciesU.push_back(UMat);
	auto JacobiandMdU =  (FuSIFu) * (dFdUT* (sInv)*dFdMS);
	auto JacobiandMdUT = TransposeMatrix(JacobiandMdU);	
	auto VarianceU = JacobiandMdU * VMat * JacobiandMdUT;	
	VarianciesU.push_back(VarianceU);
	for(int i = 0; i<nUnkn;++i){
		double du = (Unkn-Unkn0)(i,0); 	
		double dv = (VarianceU(i,i));
//		double dv = (UMat(i,i));
		dv = sqrt(dv);
		UPull.push_back( du / dv);
	}

	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
}
void KinematicFitter::ProcessStepWOUnkn(){
	#if Debug
	std::cout<<"KinematicFitter::ProcessStep"<<std::endl;
	std::cout<<"Loading Varaibles"<<std::endl;
	std::cout<<"nMeas = "<<Measurements.size()<<std::endl;
	std::cout<<"nUnkn = "<<Unknowns.size()<<std::endl;
	std::cout<<"ScaleUp = "<<ScalingMats.size()<<std::endl;
	#endif
	auto Meas = Measurements.at(step); 
	auto Meas0 = Measurements.at(0);
	auto ScaleUp = ScalingMats.at(0);
	auto ScaleDn = ScalingMats.at(1);
	auto MS0 = ScaleUp * Meas0;// MS = M/sigma 
	auto MS = ScaleUp * Meas;
	/*
	if(UpdateVariancies and step ){
		Meas0 = Measurements.at(step-1);	
	}
	*/
	auto VMat = Variancies.at(step);
	auto VInv = VMat;
#if Debug
	std::cout<<"Inverting VMat"<<std::endl;
#endif
	VInv.SetTol(1e-26);
	if(VInv.Determinant() == 0){
		cout<<"Warning: Variance Matrix is singular!"<<endl;
	}
	VInv.Invert();
#if Debug
	std::cout<<"Setting Constraints"<<std::endl;
#endif

	SetConstraints();
	TMatrixD FMat = FMats.at(step);
	TMatrixD dFdM = dFdMs.at(step);
	TMatrixD dFdMS = dFdM*ScaleDn;
#if Debug
	if(step == 0){
		cout<<"Covariance, det ="<<VMat.Determinant();
		VMat.Print();
	}
#endif
	
#if Debug>1
	cout<<"dFdM";
	dFdM.Print();
	cout<<"dFdMScaled";
	dFdMS.Print();
#endif
//	vector<TMatrixD> d2Fd2U = d2Fd2Us.at(step);
	double det_dFdM=0;
	auto dFdMT = TransposeMatrix(dFdMS);
	auto rMat = FMat + dFdMS*(MS0-MS);
	auto sMat =dFdMS*VMat*dFdMT;
	auto sInv = sMat;
	if(sInv.Determinant() == 0){
		cout<<"Warning: S Matrix is singular!"<<endl;
	}
	sInv.Invert();

	auto Lambda = (sInv* rMat);
	auto LambdaT = TransposeMatrix(Lambda);
	auto dM =  VMat*dFdMT*Lambda;
	auto Meas_next = MS0 - dM; 
	
	auto dMT = TransposeMatrix(dM);
	TMatrixD Unkn_next(0,0);
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(ScaleDn*Meas_next);	
	rMats.push_back(rMat);
	sMats.push_back(sMat);
	auto GMat = dFdMT*sInv*dFdMS;
	auto dV = VMat*(GMat )*VMat*(GMat )*VMat;
	dVMats.push_back(dV);
	auto VMat_next = VMat- VMat * GMat *VMat - VMat * GMat *VMat + dV;
	auto VInv_next = VMat_next;
	VInv_next.SetTol(1e-26);
	//VInv_next.Invert();
	
	double Chi2 = (dMT* VInv*dM)(0,0);// + 2 * (LambdaT * FMat )(0,0);
	Chi2s.push_back(Chi2);
#if Debug
	TString StepIndi = Form("[Step::%d]",step);
	cout<<StepIndi<<"##########Constraint Matrix##########"<<endl;
	cout<<"F Mat";
	FMat.Print();
	cout<<"dFdM Mat";
	dFdMS.Print();
	cout<<"dFdU Mat";
	dFdU.Print();

	std::cout<<"s Mat, Det = "<<sMat.Determinant();
	sMat.Print();
	std::cout<<"S-1 Mat :";
	sInv.Print();
	std::cout<<"FuS-1Fu-1 :";
	FuSIFu.Print();
	std::cout<<"FuS-1Fu-1 dFdU S-1 Mat :";
	(FuSIFu * dFdUT * sInv).Print();
	std::cout<<"r Mat :";
	rMat.Print();
	std::cout<<"Unkn Mat :";
	Unkn.Print();
	std::cout<<"dU Mat :";
	dU.Print();
	std::cout<<"Unkn_next :";
	Unkn_next.Print();

	std::cout<<"dFdMT Mat :";
	dFdMT.Print();
	std::cout<<"Lambda Mat :";
	Lambda.Print();
	std::cout<<"VMat :";
	VMat.Print();
	std::cout<<"Meas Mat :";
	Meas.Print();
	std::cout<<"dM Mat :";
	dM.Print();
	std::cout<<"Meas_next :";
	Meas_next.Print();

#endif 
#if Debug>1
	cout<<StepIndi<<"##########Variance Matrix##########"<<endl;
	cout<<"V Mat : Determinant = "<<VMat.Determinant();
	VMat.Print();
	cout<<"V Mat_next : Determinant = "<<VMat_next.Determinant();
	VMat_next.Print();
#if UseHessian
	cout<<"U Cov";
	UVMat.Print();
#endif
	cout<<"V*VInv : Determinant = "<<(VMat*VInv).Determinant();
	(VMat*VInv).Print();
	cout<<"V*VInv_next : Determinant = "<<(VMat_next*VInv_next).Determinant();
	(VMat_next*VInv_next).Print();
	cout<<"dV : Determinant "<<dV.Determinant();
	dV.Print();
	TCanvas*canv = new TCanvas("canv","canv",600,600);
	gStyle->SetOptStat(0);
	TH2D* HdV = new TH2D("HDV","HDV",8,0,8,8,0,8);
	for(int ir=0;ir<nMeas;++ir){
		for(int il=0;il<nMeas;++il){
			HdV->SetBinContent(ir+1,8-il,dV(ir,il));
		}
	}
	cout<<"Chi2 = "<<Chi2<<endl;
	for(auto pul:Pull){
		cout<<Form("Pull %d = %g ,",ip,pul);
		ip++;
	}
	cout<<endl;
		

	
	HdV->Draw("colz");
	gPad->SetMargin(0.1,0.2,0.1,0.1);
	canv->Modified();
	canv->Update();
	gSystem->ProcessEvents();
	cin.ignore();
#endif
	if(UpdateVariancies and step < 1){
	}else{
		VMat_next = VMat;
		VInv_next = VInv;
	}
	if(VMat_next.Determinant() == 0){
		VMat_next = VMat;
		VInv_next = VInv;
	}
//	VMat_next = VMat - dV;
	Variancies.push_back(VMat_next);
	VarianceInvs.push_back(VInv_next);
	SampleStepPoint(step);
	step++;
//	CalcVariance(step);	
	int ip = 0;
	vector<double>Pull;
	for(int i = 0; i<nMeas;++i){
		double dm = dM(i,0); 	
		double dv = dV(i,i);
		dv = sqrt(dv);
		Pull.push_back( dm / dv);
	}
	Pulls.push_back(Pull);
}
void KinematicFitter::Finalize(){
#if Debug
	cout<<"Finalizing..."<<endl;
#endif
	for(int is=0;is<step;++is){
		double chi2 = Chi2s.at(is);
#if ShowChi2
		cout<<Form("Step %d, chisqr %g",is,chi2)<<endl;
#endif
		if(chi2<Best_Chi2 and chi2>0){
			Best_Chi2 = chi2;
			best_step = is;
			best_pull = Pulls.at(is);
			if(nUnkn>0)best_Upull = UPulls.at(is);
			auto BestF = FMats.at(is);
			vector<double>consts;
			for(int i=0;i<nConst;++i){
				consts.push_back(BestF(i,0));
			}
			best_constraints = consts;
		}
	}
	if(best_step == 0){
		best_step = step-1;
		Best_Chi2 = Chi2s.at(step-1);
		best_pull = Pulls.at(step-1); 
		if(nUnkn>0)best_Upull = UPulls.at(step-1);
		auto BestF = FMats.at(step-1);
		vector<double>consts;
		for(int i=0;i<nConst;++i){
			consts.push_back(BestF(i,0));
		}
		best_constraints = consts;
	}
	auto InitialF = FMats.at(0);
	for(int i=0;i<nConst;++i){
		initial_constraints.push_back(InitialF(i,0));
	}
	SampleStepPoint(best_step);
}


double KinematicFitter::DoKinematicFit(bool Do ){
	int Cnt = 0;
	int niti=0;
	if(!Do){
		Finalize();
		return -1;
	}
	while(1){
#if Debug
		cout<<"Processed step "<<step<<endl;
#endif
		ProcessStep();
		double Chi2 = Chi2s.at(step);
		double Chi2_prev = Chi2s.at(step-1);
		if(Cnt > 3){
			break;
		}
		if(step and Chi2<0) break;
		if(abs(Chi2_prev - Chi2 )< Chi2_cut and step > 1){
			Cnt++;
		}
		else{
			Cnt= 0;
		}
		if(step > MaxStep){
			break;
		}
	}
	Finalize();
	return Best_Chi2;
}
void
KinematicFitter::Clear(){
	step = 0;
	best_step = 0;
	Best_Chi2 = 1e18;
	ScalingMats.clear();
	Variancies.clear();
	VarianceInvs.clear();
	Measurements.clear();
	Unknowns.clear();
	Pulls.clear();
	UPulls.clear();
	Chi2s.clear();
	dVMats.clear();
	FMats.clear();
	dFdMs.clear();
	dFdUs.clear();
	rMats.clear();
	sMats.clear();
	best_constraints.clear();
	initial_constraints.clear();
}

TMatrixD
KinematicFitter::TransposeMatrix(TMatrixD M){
	int row = M.GetNrows();
	int col = M.GetNcols();
	double elem[500];
	for(int r=0;r<row;++r){
	for(int c=0;c<col;++c){
		elem[c*row + r]= M(r,c); 
	}	
	}
	return TMatrixD(col,row,elem);
}
void
KinematicFitter::RotateVariance(TMatrixD J){
	auto VMat = Variancies.at(0);
	Variancies.clear();
	VMat = TransposeMatrix(J)*VMat*J;
	Variancies.push_back(VMat);
}
#endif
