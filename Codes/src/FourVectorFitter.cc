//#include "KinFit.cc"
#include "../include/FourVectorFitter.hh"
#include "KinFit.cc"
#ifndef FourVectorFitter_cc
#define FourVectorFitter_cc
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void FourVectorFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	MeasDir = status;
	Clear();
	if(MeasDir){//Production and Decay vertex of R is known, hence the direction is known. However, the momentum magnitude is not directly measured. I used the sum of P and Q as an initial value.
		auto dir = Vert1 - Vert2;	
		auto th = dir.Theta();
		auto ph = dir.Phi();
		R.SetTheta(th);
		R.SetPhi(ph);
	}
	Initialize();
}
void FourVectorFitter::Initialize(){
	Clear();
	if(MeasDir){
		nMeas = 8;nUnkn = 1; nConst = 4; 
	}
	else{
		nMeas = 6;nUnkn = 3; nConst = 4; //For 1-C fit, when Vertex information is useless.
	}
	mP = P.Mag();
	TVector3 TV_P = P.Vect();
	double p_P = TV_P.Mag();	
	double th_P = TV_P.Theta();
	double ph_P = TV_P.Phi();
	
	mQ = Q.Mag();
	TVector3 TV_Q = Q.Vect(); 
	double p_Q = TV_Q.Mag();	
	double th_Q = TV_Q.Theta();
	double ph_Q = TV_Q.Phi();

	TVector3 TV_R = R.Vect(); 
	double p_R = TV_R.Mag();	
	double th_R = TV_R.Theta();
	double ph_R = TV_R.Phi();

	double meas[8];
	double unkn[3];
	if(MeasDir){
		double temp[] = {th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_R,0,0};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	else{
		double temp[] = {p_P,th_P,ph_P,p_Q,th_Q,ph_Q,0,0};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_R,th_R,ph_R};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
	vector<double>Pull;
	Pull.resize(nMeas);
	vector<double>UPull;
	UPull.resize(nUnkn);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
	Chi2s.push_back(-1);
	MassDiffs.push_back(1e9);
}
void FourVectorFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Q= Meas(5,0); 
		th_Q= Meas(6,0); 
		ph_Q= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 
	}

	double f1 = 
		-p_R*sin(th_R)*cos(ph_R) 
		+p_P*sin(th_P)*cos(ph_P) 
		+p_Q*sin(th_Q)*cos(ph_Q) ;//Constraint on x momentum
	double f2 = 
		-p_R*sin(th_R)*sin(ph_R) 
		+p_P*sin(th_P)*sin(ph_P) 
		+p_Q*sin(th_Q)*sin(ph_Q) ;//Constraint on y momentum
	double f3 =  
		-p_R*cos(th_R) 
		+p_P*cos(th_P)
		+p_Q*cos(th_Q);//Constraint on z momentum 
	double f4	=
		- sqrt(p_R*p_R+mR*mR)
		+ sqrt(p_P*p_P+mP*mP)
		+ sqrt(p_Q*p_Q+mQ*mQ);//Constraint on Energy

	double df1du1 = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
	double df1du2 = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1 in case of 3-C fit.However, I didnt want to change the token... Mathematically it should be df1 / d (Th_R)
	double df1du3 = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
	double df1dm1 = sin(th_P)*cos(ph_P);//...
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df1dm4 = sin(th_Q)*cos(ph_Q);
	double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
	double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	
	double df2du1 = -sin(th_R)*sin(ph_R);
	double df2du2 = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
	double df2du3 = -p_R*sin(th_R)*cos(ph_R);
	double df2dm1 = sin(th_P)*sin(ph_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	double df2dm4 = sin(th_Q)*sin(ph_Q);
	double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
	double df2dm6 = p_Q*sin(th_Q)*cos(ph_Q);

	double df3du1 = -cos(th_R);
	double df3du2 = p_R*sin(th_R);
	double df3du3 = 0;
	double df3dm1 = cos(th_P);
	double df3dm2 = -p_Q*sin(th_P);
	double df3dm3 = 0;
	double df3dm4 = cos(th_Q);
	double df3dm5 = -p_Q*sin(th_Q);
	double df3dm6 = 0;
	
	double ER = sqrt(p_R*p_R+mR*mR);
	double df4du1 = -p_R/ER;
	double df4du2 = 0;
	double df4du3 = 0;
	double df4dm1 = p_P/sqrt(p_P*p_P+mP*mP);
	double df4dm2 = 0;
	double df4dm3 = 0;
	double df4dm4 = p_Q/sqrt(p_Q*p_Q+mQ*mQ);
	double df4dm5 = 0;
	double df4dm6 = 0;


	double df1du1du1 = 0;
	double df1du1du2 = -cos(th_R)*cos(ph_R);
	double df1du1du3 = sin(th_R)*sin(ph_R);
	
	double df1du2du1 = df1du1du2;
	double df1du2du2 = p_R*sin(th_R)*cos(ph_R);
	double df1du2du3 = p_R*cos(th_R)*sin(ph_R);
	
	double df1du3du1 = df1du1du3;
	double df1du3du2 = df1du2du3;
	double df1du3du3 = p_R*sin(th_R)*cos(ph_R);

	double df2du1du1 = 0;
	double df2du1du2 = -cos(th_R)*sin(ph_R);
	double df2du1du3 = -sin(th_R)*cos(ph_R);

	double df2du2du1 = df2du1du2;
	double df2du2du2 = p_R*sin(th_R)*sin(ph_R);
	double df2du2du3 = -p_R*cos(th_R)*cos(ph_R);
	
	double df2du3du1 = df2du1du3;
	double df2du3du2 = df2du2du3;
	double df2du3du3 = p_R*sin(th_R)*sin(ph_R);
	double df3du1du1 = 0;
	double df3du1du2 = sin(th_R);
	double df3du1du3 = 0;
	
	double df3du2du1 = df3du1du2;
	double df3du2du2 = p_R*cos(th_R);
	double df3du2du3 = 0;

	double df3du3du1 = 0;
	double df3du3du2 = 0;
	double df3du3du3 = 0;
	
	double df4du1du1 = -mR*mR/ER/ER/ER;
	double df4du1du2 = 0;
	double df4du1du3 = 0;
	
	double df4du2du1 = 0;
	double df4du2du2 = 0;
	double df4du2du3 = 0;

	double df4du3du1 = 0;
	double df4du3du2 = 0;
	double df4du3du3 = 0;



	double fs[]={f1,f2,f3,f4};
	double dfdms[200] ;
	double dfdus[200] ;
	if(MeasDir){
		double temp[] = {
			df1du2,df1du3,df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,
			df2du2,df2du3,df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,
			df3du2,df3du3,df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,
			df4du2,df4du3,df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,
			df2du1,
			df3du1,
			df4du1
		};	
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[]={
			df1du1du1
		};
		double temp2[]={
			df2du1du1
		};
		double temp3[]={
			df3du1du1
		};
		double temp4[]={
			df4du1du1
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);
	}
	else{
		double temp[] = {
			df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,
			df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,
			df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,
			df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
		double temp1[] = {
			df1du1du1,df1du1du2,df1du1du3,
			df1du2du1,df1du2du2,df1du2du3,
			df1du3du1,df1du3du2,df1du3du3
		};
		double temp2[] = {
			df2du1du1,df2du1du2,df2du1du3,
			df2du2du1,df2du2du2,df2du2du3,
			df2du3du1,df2du3du2,df2du3du3
		};
		double temp3[] = {
			df3du1du1,df3du1du2,df3du1du3,
			df3du2du1,df3du2du2,df3du2du3,
			df3du3du1,df3du3du2,df3du3du3
		};
		double temp4[] = {
			df4du1du1,df4du1du2,df4du1du3,
			df4du2du1,df4du2du2,df4du2du3,
			df4du3du1,df4du3du2,df4du3du3
		};
		TMatrixD d2F1dU(nUnkn,nUnkn,temp1);
		TMatrixD d2F2dU(nUnkn,nUnkn,temp2);
		TMatrixD d2F3dU(nUnkn,nUnkn,temp3);
		TMatrixD d2F4dU(nUnkn,nUnkn,temp4);
		vector<TMatrixD> d2FdU = {
			d2F1dU,
			d2F2dU,
			d2F3dU,
			d2F4dU
		};
		d2Fd2Us.push_back(d2FdU);


	}

	



	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void FourVectorFitter::SampleStepPoint(int steps){
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_R,th_R,ph_R;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Q= Meas(5,0); 
		th_Q= Meas(6,0); 
		ph_Q= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 
	}
	double Ppx = p_P*sin(th_P)*cos(ph_P);
	double Ppy = p_P*sin(th_P)*sin(ph_P);
	double Ppz = p_P*cos(th_P);
	double Qpx = p_Q*sin(th_Q)*cos(ph_Q);
	double Qpy = p_Q*sin(th_Q)*sin(ph_Q);
	double Qpz = p_Q*cos(th_Q);
	double Rpx = p_R*sin(th_R)*cos(ph_R);
	double Rpy = p_R*sin(th_R)*sin(ph_R);
	double Rpz = p_R*cos(th_R);
	TLorentzVector PP(Ppx,Ppy,Ppz,hypot(mP,p_P));
	TLorentzVector QQ(Qpx,Qpy,Qpz,hypot(mQ,p_Q));
	TLorentzVector RR(Rpx,Rpy,Rpz,hypot(mR,p_R));
	auto V = PP + QQ;
	double MassDiff = V.Mag()-mR;
	PCor = PP;
	QCor = QQ;
	RCor = RR;
	MassDiffs.push_back(MassDiff);
}
FourVectorFitter::FourVectorFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	//Does Kinematic fitting for R -> P+Q decay.
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
TMatrixD
FourVectorFitter::JacobianSphToCart(double p, double th, double ph){
	// x = p sin(th) cos(ph)
	// y = p sin(th) sin(ph)
	// z = p cos(th)
	//V_c = J^T V J|->
	//		dxdp, dxdth,dxdph
	//J	=	dydp, dydth,dydph
	//		dzdp, dzdth,dzdph

	double dxdp = sin(th)*cos(ph);
	double dydp = sin(th)*sin(ph);
	double dzdp = cos(th);

	double dxdth = p*cos(th)*cos(ph);
	double dydth = p*cos(th)*sin(ph);
	double dzdth = -p*sin(th);

	double dxdph = -p*sin(th)*sin(ph);
	double dydph = p*sin(th)*cos(ph);
	double dzdph = 0;
	double mat[9] = 
	{ dxdp, dxdth, dxdph,
		dydp, dydth, dydph,
		dzdp, dzdth, dzdph
	};
	/*
	double mat[9] = 
	{ dxdp, dydp, dzdp,
		dxdth, dydth, dzdth,
		dxdph, dydph, dzdph
	};
	*/
	return TMatrixD(3,3,mat);

}
void
FourVectorFitter::CalcVariance(int istep){
	auto Meas = Measurements.at(istep); 
	auto Unkn = Unknowns.at(istep);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q;
	if(MeasDir){
		p_R= Unkn(0,0); 
		th_R= Meas(0,0); 
		ph_R= Meas(1,0); 
		p_P= Meas(2,0); 
		th_P= Meas(3,0); 
		ph_P= Meas(4,0);
		p_Q= Meas(5,0); 
		th_Q= Meas(6,0); 
		ph_Q= Meas(7,0); 
	}
	else{
		p_R= Unkn(0,0); 
		th_R= Unkn(1,0); 
		ph_R= Unkn(2,0); 
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);
		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 
	}
	TMatrixD Jsc_P = JacobianSphToCart(p_P,th_P,ph_P);
	TMatrixD Jsc_Q = JacobianSphToCart(p_Q,th_Q,ph_Q);
	
	double El_Jsc_PQ[6*6]= {0};
	for(int ic =0;ic<3;++ic){
	for(int ir =0;ir<3;++ir){
		int col_P = ic, row_P = ir;
		int col_Q = ic+3, row_Q = ir+3;
		El_Jsc_PQ[row_P+6*col_P] = Jsc_P(ic,ir);
		El_Jsc_PQ[row_Q+6*col_Q] = Jsc_Q(ic,ir);
	}
	}
	TMatrixD Jsc_PQ = TMatrixD(6,6,El_Jsc_PQ);
	/*
	Jsc_P.Print();
	Jsc_Q.Print();
	Jsc_PQ.Print();
	cin.ignore();
	*/
	TMatrixD Jsc_PQ_T = TransposeMatrix(Jsc_PQ);
	TMatrixD VMat = Variancies.at(istep);
	TMatrixD dV = dVMats.at(istep);
	TMatrixD VMat_C = Jsc_PQ_T*(VMat-dV)*Jsc_PQ;
//	TMatrixD VMat_C = Jsc_PQ_T*(VMat)*Jsc_PQ;
	double El_contract[18]={
		1,0,0,1,0,0,
		0,1,0,0,1,0,
		0,0,1,0,0,1
	};
	TMatrixD ContT(3,6,El_contract);
	TMatrixD Cont = TransposeMatrix(ContT);
	TMatrixD UVMat_C = ContT*VMat_C*Cont;
	TMatrixD Jcs_R = JacobianSphToCart(p_R,th_R,ph_R);
	Jcs_R.Invert();
	auto Jcs_RT = TransposeMatrix(Jcs_R);
	TMatrixD UVMat = Jcs_RT*UVMat_C*Jcs_R;
//	VarianciesU.push_back(UVMat);
}















void
FourVectorFitter::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
FourVectorFitter::ToDecayPlane(){
	auto Zaxis =(P + Q).Vect();
	auto vP = P.Vect();
	auto vQ = Q.Vect();
	auto Yaxis = vP.Cross(vQ);
//	double YNorm = 1./(Yaxis.Mag());
//	Yaxis = YNorm * Yaxis;
	double Th_F = Zaxis.Theta();
	double Ph_F = Zaxis.Phi();
	double RotZ[9] ={
		cos(-Ph_F),	-sin(-Ph_F),	0,	
		sin(-Ph_F),	cos(-Ph_F),		0,
		0,					0,						1
	};
	double RotY[9] ={
		cos(Th_F),	0,					-sin(Th_F),
		0,					-1,					0,
		sin(Th_F),	0,					cos(Th_F)
	};
	TMatrixD RZ(3,3,RotZ);
	TMatrixD RY(3,3,RotY);
	TMatrixD R_F = RY * RZ;
	Yaxis = R_F * Yaxis;
	double Th_Y = Yaxis.Theta();
	double Ph_Y = Yaxis.Phi();
	double RotX[9] ={
		1,				0,					0,
		0,				cos(-Ph_Y),	-sin(-Ph_Y),
		0,				sin(-Ph_Y),	cos(-Ph_Y)
	};
	TMatrixD RX(3,3,RotX);
	Yaxis = RX * Yaxis;

}
#endif
