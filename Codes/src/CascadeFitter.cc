//#include "KinFit.cc"
#include "../include/CascadeFitter.hh"
#include "KinFit.cc"
#ifndef CascadeFitter_cc
#define CascadeFitter_cc
#define DebugCKF 0
#define LambdaFit 0
#define InvMassFit 0
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void CascadeFitter::UseVertex(bool status,TVector3 Vert1,TVector3 Vert2){
	UseVertexFlag = status;
	Clear();
	if(UseVertexFlag){//Production and Decay vertex of R is known, hence the direction is known. However, the momentum magnitude is not directly measured. I used the sum of P and Q as an initial value.
	}
	Initialize();
}
void CascadeFitter::Initialize(){
#if DebugCKF
	cout<<"Initializing..."<<endl;
#endif
	Clear();
	if(UseVertexFlag){
		nMeas = 11;nUnkn = 2; nConst = 9; 
	}
	else{
		nMeas = 9;nUnkn = 3; nConst = 5;
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

	mR = R.Mag();
	TVector3 TV_R = R.Vect(); 
	double p_R = TV_R.Mag();	
	double th_R = TV_R.Theta();
	double ph_R = TV_R.Phi();

	TVector3 TV_L = (P+Q).Vect();
	double p_L = TV_L.Mag();
	double th_L = TV_L.Theta();
	double ph_L = TV_L.Phi();

	TVector3 TV_X = (P+Q+R).Vect();
	double p_X = TV_X.Mag();
	double th_X = TV_X.Theta();
	double ph_X = TV_X.Phi();
	double meas[20];
	double unkn[20];
	if(UseVertexFlag){
	}
	else{
		double temp[] = {p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_R,th_R,ph_R};
		for(int i=0;i<nMeas;++i)meas[i]=temp[i];
		double temp2[] = {p_X,th_X,ph_X};
		for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
	}
	TMatrixD Meas0(nMeas,1,meas);	
	TMatrixD Unkn0(nUnkn,1,unkn);
#if DebugCKF
	cout<<"Meas0 : ";
	Meas0.Print();
	cout<<"Unkn0 : ";
	Unkn0.Print();
#endif
	vector<double>Pull;
	Pull.resize(nMeas);
	vector<double>UPull;
	UPull.resize(nUnkn);
	Measurements.push_back(Meas0);
	Unknowns.push_back(Unkn0);
	Pulls.push_back(Pull);
	UPulls.push_back(UPull);
	Chi2s.push_back(-1);
	MassDiffsL.push_back(1e9);
	MassDiffsX.push_back(1e9);
}
void CascadeFitter::SetConstraints(){
	auto Meas = Measurements.at(step); 
	auto Unkn = Unknowns.at(step);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_X,th_X,ph_X;	
	if(UseVertexFlag){
	}
	else{
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);

		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 

		p_R= Meas(6,0);
		th_R= Meas(7,0);
		ph_R= Meas(8,0);

		p_X= Unkn(0,0);
		th_X= Unkn(1,0);
		ph_X= Unkn(2,0);
	}
	double E_P = hypot(p_P,mP);
	double dEdp_P = p_P/E_P;
	double p_Px = p_P*sin(th_P)*cos(ph_P);
	double p_Py = p_P*sin(th_P)*sin(ph_P);
	double p_Pz = p_P*cos(th_P);
	
	double E_Q = hypot(p_Q,mQ);
	double dEdp_Q = p_Q/E_Q;
	double p_Qx = p_Q*sin(th_Q)*cos(ph_Q);
	double p_Qy = p_Q*sin(th_Q)*sin(ph_Q);
	double p_Qz = p_Q*cos(th_Q);


	double E_R = hypot(p_R,mR);
	double dEdp_R = p_R/E_R;
	double p_Rx = p_R*sin(th_R)*cos(ph_R);
	double p_Ry = p_R*sin(th_R)*sin(ph_R);
	double p_Rz = p_R*cos(th_R);

	double p_Lx = p_Px + p_Qx;
	double p_Ly = p_Py + p_Qy;
	double p_Lz = p_Pz + p_Qz;
	double p_L = hypot(p_Lx,hypot(p_Ly,p_Lz));
	double E_L = hypot(p_L,mL);
	double dEdp_L = p_L/E_L;
	//For L derivatives...
	double cosPQ = (p_Px*p_Qx + p_Py*p_Qy + p_Pz*p_Qz)/p_P/p_Q;

	double dp_Ldp_P = (p_P+p_Q*cosPQ)/p_L;
	double dp_Ldth_P = p_P * (p_Qx *cos(th_P)*cos(ph_P)+p_Qy*cos(th_P)*sin(ph_P) - p_Qz *sin(th_P) )/p_L;
	double dp_Ldph_P = (-p_Qx *p_Py + p_Qy*p_Px)/p_L;

	double dp_Ldp_Q = (p_Q+p_P*cosPQ)/p_L;
	double dp_Ldth_Q = p_Q * (p_Px *cos(th_Q)*cos(ph_Q)+p_Py*cos(th_Q)*sin(ph_Q) - p_Pz *sin(th_Q) )/p_L;
	double dp_Ldph_Q = -dp_Ldph_P;//If P and Q are rotated the same ammount on phi, P+Q should have same magnitude. Hence, dP_L/dph_Q + dP_L/dph_P = 0
	//
		


	double E_X = hypot(p_X,mX);
	double dEdp_X = p_X/E_X;
	double p_Xx = p_X*sin(th_X)*cos(ph_X);
	double p_Xy = p_X*sin(th_X)*sin(ph_X);
	double p_Xz = p_X*cos(th_X);

	double f1 =
		-p_X*sin(th_X)*cos(ph_X) 
		+p_P*sin(th_P)*cos(ph_P)
		+p_Q*sin(th_Q)*cos(ph_Q)
		+p_R*sin(th_R)*cos(ph_R);
	double f2 =
		-p_X*sin(th_X)*sin(ph_X)
		+p_P*sin(th_P)*sin(ph_P)
		+p_Q*sin(th_Q)*sin(ph_Q)
		+p_R*sin(th_R)*sin(ph_R);
	double f3 =
		-p_X*cos(th_X)
		+p_P*cos(th_P)
		+p_Q*cos(th_Q)
		+p_R*cos(th_R);
	double f4	=
		-E_L + E_P + E_Q;//Constraint on Lambda Energy
	double f5 = 
		-E_X + E_P +E_Q + E_R;//Constraint on Xi Energy
	
	double df1du1 = -sin(th_X)*cos(ph_X);
	double df1du2 = -p_X*cos(th_X)*cos(ph_X);
	double df1du3 = p_X*sin(th_X)*sin(ph_X);
	
	double df1dm1 = sin(th_P)*cos(ph_P);
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df1dm4 = sin(th_Q)*cos(ph_Q);
	double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);
	double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	double df1dm7 = sin(th_R)*cos(ph_R);
	double df1dm8 = p_R*cos(th_R)*cos(ph_R);
	double df1dm9 = -p_R*sin(th_R)*sin(ph_R);


	double df2du1 = -sin(th_X)*sin(ph_X);
	double df2du2 = -p_X*cos(th_X)*sin(ph_X);
	double df2du3 = -p_X*sin(th_X)*cos(ph_X);

	double df2dm1 = sin(th_P)*sin(ph_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	double df2dm4 = sin(th_Q)*sin(ph_Q);
	double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);
	double df2dm6 = p_Q*sin(th_Q)*cos(ph_Q);
	double df2dm7 = sin(th_R)*sin(ph_R);
	double df2dm8 = p_R*cos(th_R)*sin(ph_R);
	double df2dm9 = p_R*sin(th_R)*cos(ph_R);


	double df3du1 = -cos(th_X);
	double df3du2 = p_X*sin(th_X);
	double df3du3 = 0;

	double df3dm1 = cos(th_P);
	double df3dm2 = -p_P*sin(th_P);
	double df3dm3 = 0;
	double df3dm4 = cos(th_Q);
	double df3dm5 = -p_Q*sin(th_Q);
	double df3dm6 = 0;
	double df3dm7 = cos(th_R);
	double df3dm8 = -p_R*sin(th_R);
	double df3dm9 = 0;

	double df4du1 = 0;
	double df4du2 = 0;
	double df4du3 = 0;

	double df4dm1 = 
	-dEdp_L * dp_Ldp_P
	+ dEdp_P;
	double df4dm2 = -dEdp_L * dp_Ldth_P;
	double df4dm3 = -dEdp_L * dp_Ldph_P;
	double df4dm4 =
	-dEdp_L * dp_Ldp_Q
	+ dEdp_Q;
	double df4dm5 = -dEdp_L * dp_Ldth_Q;
	double df4dm6 = -dEdp_L * dp_Ldph_Q;
	double df4dm7 = 0;
	double df4dm8 = 0;
	double df4dm9 = 0;

	double df5du1 = -dEdp_X;
	double df5du2 = 0;
	double df5du3 = 0;

	double df5dm1 = dEdp_P;//p_X is an Unknown, so p_X = p_P + p_x + p_R is a constraint, not a definition. In contrast, p_L = p_P + p_Q is a definition.
	double df5dm2 = 0;
	double df5dm3 = 0;
	double df5dm4 = dEdp_Q;
	double df5dm5 = 0;
	double df5dm6 = 0;
	double df5dm7 = dEdp_R;
	double df5dm8 = 0;
	double df5dm9 = 0;

	double fs[]={f1,f2,f3,f4,f5};
	double dfdms[200] ;
	double dfdus[200] ;
	if(UseVertexFlag){
	}
	else{
		double temp[] = {
			df1dm1,df1dm2,df1dm3,df1dm4	,df1dm5,df1dm6,df1dm7,df1dm8,df1dm9,
			df2dm1,df2dm2,df2dm3,df2dm4	,df2dm5,df2dm6,df2dm7,df2dm8,df2dm9,
			df3dm1,df3dm2,df3dm3,df3dm4	,df3dm5,df3dm6,df3dm7,df3dm8,df3dm9,
			df4dm1,df4dm2,df4dm3,df4dm4	,df4dm5,df4dm6,df4dm7,df4dm8,df4dm9,
			df5dm1,df5dm2,df5dm3,df5dm4	,df5dm5,df5dm6,df5dm7,df5dm8,df5dm9,
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,df1du2,df1du3,
			df2du1,df2du2,df2du3,
			df3du1,df3du2,df3du3,
			df4du1,df4du2,df4du3,
			df5du1,df5du2,df5du3
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
	}

	TMatrixD FMat(nConst,1,fs);
	TMatrixD dFdM(nConst,nMeas,dfdms);
	TMatrixD dFdU(nConst,nUnkn,dfdus);
	FMats.push_back(FMat);//Constraint Matrices for Each step
	dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
	dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void CascadeFitter::SampleStepPoint(int steps){
	auto Meas = Measurements.at(steps); 
	auto Unkn = Unknowns.at(steps); 
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_X,th_X,ph_X;	
	if(UseVertexFlag){
	}
	else{
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);

		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 

		p_R= Meas(6,0);
		th_R= Meas(7,0);
		ph_R= Meas(8,0);

		p_X= Unkn(0,0);
		th_X= Unkn(1,0);
		ph_X= Unkn(2,0);
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
	double Lpx = Ppx + Qpx;
	double Lpy = Ppy + Qpy;
	double Lpz = Ppz + Qpz;
	double Xpx = p_X*sin(th_X)*cos(ph_X);
	double Xpy = p_X*sin(th_X)*sin(ph_X);
	double Xpz = p_X*cos(th_X);

	double p_L = hypot(Lpx,hypot(Lpy,Lpz));
	TLorentzVector PP(Ppx,Ppy,Ppz,hypot(mP,p_P));
	TLorentzVector QQ(Qpx,Qpy,Qpz,hypot(mQ,p_Q));
	TLorentzVector RR(Rpx,Rpy,Rpz,hypot(mR,p_R));
	TLorentzVector LL(Lpx,Lpy,Lpz,hypot(mL,p_L));
	TLorentzVector XX(Xpx,Xpy,Xpz,hypot(mX,p_X));
	auto L_PQ = PP+QQ;
	auto X_PQR = PP+QQ+RR;

	double MassDiffL = L_PQ.Mag()-mL;
	double MassDiffX = X_PQR.Mag()-mX;
	PCor = PP;
	QCor = QQ;
	RCor = RR;
	LCor = LL;
	XCor = XX;
	MassDiffsL.push_back(MassDiffL);
	MassDiffsX.push_back(MassDiffX);
}
CascadeFitter::CascadeFitter(TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
	P=P_;
	Q=Q_;
	R=R_;
	Initialize();
};
TMatrixD
CascadeFitter::JacobianSphToCart(double p, double th, double ph){
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
CascadeFitter::CalcVariance(int istep){//
	auto Meas = Measurements.at(istep); 
	auto Unkn = Unknowns.at(istep);
	double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_L,th_L,ph_L,p_X,th_X,ph_X;	
	if(UseVertexFlag){
	}
	else{
		p_P= Meas(0,0); 
		th_P= Meas(1,0); 
		ph_P= Meas(2,0);

		p_Q= Meas(3,0); 
		th_Q= Meas(4,0); 
		ph_Q= Meas(5,0); 

		p_R= Meas(6,0);
		th_R= Meas(7,0);
		ph_R= Meas(8,0);

		p_L= Unkn(0,0);
		th_L= Unkn(1,0);
		ph_L= Unkn(2,0);

		p_X= Unkn(3,0);
		th_X= Unkn(4,0);
		ph_X= Unkn(5,0);
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
CascadeFitter::Rotate(){
	auto VMat = Variancies.at(0);
	Initialize();
	Variancies.push_back(VMat);
	TMatrixD J;
	RotateVariance(J);
}
void
CascadeFitter::ToDecayPlane(){
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
