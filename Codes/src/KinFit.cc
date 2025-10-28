#include "../include/KinFit.hh"
#include "TMatrixDEigen.h"
#include <iostream>
#ifndef KinFit_cc
#define KinFit_cc
#define Debug 0
#define ShowChi2 0
#define DebugHessian 0
#define Hessian 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
using namespace std;
void
KinematicFitter::SetVariance(double* var){
	double variance[200];
	double varianceInv[200];
	TMatrixD ScaleUp(nMeas,nMeas);
	TMatrixD ScaleDn(nMeas,nMeas);
	for(int i = 0;i<nMeas*nMeas;++i){
		variance[i]=0;
		varianceInv[i]=0;
	}
	for(int i = 0;i<nMeas;++i){
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
	TMatrixD ZeroMatU(nUnkn,nUnkn);
	Variancies.push_back(Variance);
	VarianceInvs.push_back(VarianceInv);
	VarianciesU.push_back(ZeroMatU);
//	TMatrixD UHessian(nUnkn,nUnkn);
//	UHessians.push_back(UHessian);
	dVMats.push_back(ZeroMat);
//	CalcVariance(0);
	ScalingMats.push_back(ScaleUp);
	ScalingMats.push_back(ScaleDn);
#if Debug
	cout<<"Setting Variance... det ="<<Variance.Determinant();
	Variance.Print();
#endif
};
void
KinematicFitter::AddOffdiagonals(TMatrixD Cov){
#if Debug
	cout<<"Covariance";
	Cov.Print();
#endif
	if(ScaleParams){
		auto ScaleUp = ScalingMats.at(0);
		Cov = Cov * ScaleUp;
		Cov = ScaleUp * Cov;
	}
#if Debug
	if(ScaleParams){
		cout<<"Covariance After scaling";
		Cov.Print();
	}
#endif
	Variancies.at(0)+= Cov;
	bool CheckPositiveDefinite = false;
	int nitr = 0;
	while(!CheckPositiveDefinite and nitr < 5){
		bool PositiveDefinate = true;
		TMatrixDEigen VEig(Variancies.at(0));
		TVectorD eigenvalue = VEig.GetEigenValuesRe();
		for(int ie=0;ie<eigenvalue.GetNrows();++ie){
			if(eigenvalue(ie)<0){
				PositiveDefinate = false;
				break;
			}
			CheckPositiveDefinite = true;
		}
		if(!PositiveDefinate){
			cout<<"KinematicFit:: Variance not positive-definite"<<endl;
			cout<<"Reducing Offdiagonals..."<<endl;
			#if Debug
			Variancies.at(0).Print();	
			#endif
			for(int irow=0;irow < Variancies.at(0).GetNrows();++irow){
				for(int icol=0;icol < Variancies.at(0).GetNcols();++icol){
					if(irow == icol) continue;
					Variancies.at(0)(irow,icol) *=0.9;
				}
			}
		}
		nitr++;
	}
//	cout<<"CovarianceMat : ";
#if Debug
	cout<<"Variance, det ="<<Variancies.at(0).Determinant();
	Variancies.at(0).Print();
#endif
}
void KinematicFitter::ProcessStep(){
#if Debug
	cout<<"Processing step "<<step<<endl;
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
	VInv.SetTol(1e-26);
	VInv.Invert();
	SetConstraints();

	TMatrixD FMat = FMats.at(step);
	TMatrixD dFdM = dFdMs.at(step);
	TMatrixD dFdMS = dFdM*ScaleDn;
	FMat.SetTol(1e-26);
	dFdM.SetTol(1e-26);
	dFdMS.SetTol(1e-26);
#if Debug
	cout<<"Variance, det ="<<VMat.Determinant();
	VMat.Print();
#endif
	
#if Debug 
	cout<<"dFdM";
	dFdM.Print();
	cout<<"dFdMScaled";
	dFdMS.Print();
#endif
	TMatrixD dFdU = dFdUs.at(step);
	dFdU.SetTol(1e-26);
//	vector<TMatrixD> d2Fd2U = d2Fd2Us.at(step);
	double det_dFdM=0,det_dFdU=0;
	auto dFdMT = TransposeMatrix(dFdMS);
	auto dFdUT = TransposeMatrix(dFdU);
	auto rMat = FMat + dFdMS*(MS0-MS);
	dFdMT.SetTol(1e-26);
	dFdUT.SetTol(1e-26);
	rMat.SetTol(1e-26);
	auto sMat =dFdMS*VMat*dFdMT;
	sMat.SetTol(1e-26);
	auto sInv = sMat;
#if Debug
	cout<<"sMat, det ="<<sMat.Determinant();
	sMat.Print();
#endif
	sInv.Invert();
	auto FuSIFu =	dFdUT*sInv*dFdU;
	FuSIFu.SetTol(1e-26);
#if Debug
	cout<<"FuSIFu, det ="<<FuSIFu.Determinant();
	FuSIFu.Print();
#endif
	FuSIFu.Invert();
	auto dU = (FuSIFu) * (dFdUT* (sInv) * rMat) ;
	dU = dU -dU - dU;
	dU =  dU;
	auto Unkn_next = Unkn +  dU;

	auto Lambda = (sInv* (rMat+ dFdU*dU));
	auto LambdaT = TransposeMatrix(Lambda);
		
	auto dM =  VMat*dFdMT*Lambda;
	auto Meas_next = MS0 -  dM; 
	
	auto dMT = TransposeMatrix(dM);
//	auto UHessian = TMatrixD(nUnkn,nUnkn);
	for(int ic =0; ic < nConst;++ic){
//		UHessian+= 2*Lambda(ic,0)*d2Fd2U.at(ic);
	}
#if DebugHessian
	cout<<Form("Step %d Hessian",step);
	UHessian.Print();
	cout<<Form("Step %d UVariance",step);
	UVMat.Print();
	cin.ignore();
#endif
	Unknowns.push_back(Unkn_next);	
	Measurements.push_back(ScaleDn*Meas_next);	
	rMats.push_back(rMat);
	sMats.push_back(sMat);
	auto GMat = dFdMT*sInv*dFdMS;
	auto HMat = dFdMT*sInv*dFdU;
	auto UMat = dFdUT*sInv*dFdU;
	UMat.SetTol(1e-26);	
#if Debug
	cout<<"UMat, det ="<<UMat.Determinant();
	UMat.Print();
#endif

	UMat.Invert();
//	UHessians.push_back(UHessian);
	auto HMatT = TransposeMatrix(HMat);
	auto HUH = (HMat*UMat*HMatT);
	auto dV = VMat*(GMat -(HUH) )*VMat*(GMat - HUH)*VMat;
	dVMats.push_back(dV);
	auto VMat_next = VMat- VMat * (GMat - HUH)*VMat - VMat * (GMat - HUH)*VMat + dV;
	auto VInv_next = VMat_next;
//	auto UVMat = UHessian;
//	UVMat.Invert();
	VInv_next.SetTol(1e-26);
#if Debug
	cout<<"Variance, det ="<<VMat.Determinant();
	VMat.Print();
#endif
	VInv_next.Invert();
	
	double Chi2 = (dMT* (VInv)*dM)(0,0) + 2 * (LambdaT * FMat )(0,0);
	Chi2s.push_back(Chi2);
	TMatrixD CHI2_U = TMatrixD(nUnkn,nUnkn);
	for(int iu=0;iu<nUnkn;++iu){
		CHI2_U(iu,iu)=Chi2;
	}
//	UVMat=UVMat*CHI2_U;


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
#if Debug>1
	TString StepIndi = Form("[Step::%d]",step);
	cout<<StepIndi<<"##########Variance Matrix##########"<<endl;
	cout<<"V Mat : Determinant = "<<VMat.Determinant();
	VMat.Print();
	cout<<"V Mat_next : Determinant = "<<VMat_next.Determinant();
	VMat_next.Print();
//	cout<<"U Cov";
//	UVMat.Print();
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
		

	cout<<StepIndi<<"##########Constraint Matrix##########"<<endl;
	cout<<"F Mat";
	FMat.Print();
	cout<<"dFdM Mat";
	dFdMS.Print();
	cout<<"dFdU Mat";
	dFdU.Print();
	cout<<"r Mat";
	rMat.Print();
	cout<<"s Mat, Det = "<<sMat.Determinant();
	sMat.Print();
	cout<<"s Inv";
	sInv.Print();
	

	cout<<StepIndi<<"##########Updated Matrix##########"<<endl;
	cout<<"dU Mat";
	dU.Print();
	cout<<"UnknNext Mat";
	Unkn_next.Print();
	cout<<"dM Mat";
	dM.Print();
	cout<<"MeasNext Mat";
	Meas_next.Print();
	cout<<"Lambda Mat";
	Lambda.Print();
	cout<<"GMat ";
	GMat.Print();
	cout<<"HMat";
	HMat.Print();
	cout<<"UMat";
	UMat.Print();
	cout<<"HUHTMat";
	(HMat*UMat*HMatT).Print();
	
	HdV->Draw("colz");
	gPad->SetMargin(0.1,0.2,0.1,0.1);
	canv->Modified();
	canv->Update();
	gSystem->ProcessEvents();
	cin.ignore();
#endif
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
			best_Upull = UPulls.at(is);
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
		best_Upull = UPulls.at(step-1);
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
		ProcessStep();
#if Debug
		cout<<"Processed step "<<step<<endl;
#endif
		double Chi2 = Chi2s.at(step);
		double Chi2_prev = Chi2s.at(step-1);
		if(Cnt > 3 and step>10){
			break;
		}
		if(step and Chi2<0) break;
		if(abs(Chi2_prev - Chi2 )< Chi2_cut and step > 0){
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
	UHessians.clear();
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