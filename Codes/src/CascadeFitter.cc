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
		nMeas = 9;nUnkn = 6; nConst = 8; //For 1-C fit, when Vertex information is useless.
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
		double temp2[] = {p_L,th_L,ph_L,p_X,th_X,ph_X};
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
	
	double E_L = hypot(p_L,mL);
	double dEdp_L = p_L/E_L;
	double p_Lx = p_L*sin(th_L)*cos(ph_L);
	double p_Ly = p_L*sin(th_L)*sin(ph_L);
	double p_Lz = p_L*cos(th_L);
	
	double E_X = hypot(p_X,mX);
	double dEdp_X = p_X/E_X;
	double p_Xx = p_X*sin(th_X)*cos(ph_X);
	double p_Xy = p_X*sin(th_X)*sin(ph_X);
	double p_Xz = p_X*cos(th_X);



	double f1 = 
		-p_L*sin(th_L)*cos(ph_L)
		+p_P*sin(th_P)*cos(ph_P)
		+p_Q*sin(th_Q)*cos(ph_Q);//Constraint on Lambda x momentum
	double f2 =
		-p_L*sin(th_L)*sin(ph_L)
		+p_P*sin(th_P)*sin(ph_P)
		+p_Q*sin(th_Q)*sin(ph_Q);//Constraint on Lambda y momentum	
	double f3 =
		-p_L*cos(th_L)
		+p_P*cos(th_P)
		+p_Q*cos(th_Q);//Constraint on Lambda z momentum
	double f4	=
#if InvMassFit
		mL*mL - (pow(E_P + E_Q,2)- pow(hypot(p_Px+p_Qx,hypot(p_Py+p_Qy,p_Pz+p_Qz)),2));//Constraint on Lambda Energy
#else
		-E_L + E_P + E_Q;//Constraint on Lambda Energy
#endif
//		- E_L*E_L + pow(E_P + E_Q,2);//Constraint on Lambda Energy
#if LambdaFit
	double f5 = 
		-p_X*sin(th_X)*cos(ph_X)
		+p_L*sin(th_L)*cos(ph_L) 
		+p_R*sin(th_R)*cos(ph_R); //Constraint on Xi x momentum
	double f6 = 
		-p_X*sin(th_X)*sin(ph_X)
		+p_L*sin(th_L)*sin(ph_L) 
		+p_R*sin(th_R)*sin(ph_R);//Constraint on Xi y momentum
	double f7 =  
		-p_X*cos(th_X)
		+p_L*cos(th_L)
		+p_R*cos(th_R);//Constraint on Xi z momentum 
	double f8	=
#if InvMassFit
		mX*mX - (pow(E_L + E_R,2)- pow(hypot(p_Lx+p_Rx,hypot(p_Ly+p_Ry,p_Lz+p_Rz)),2));//Constraint on Xi Energy
#else
		-E_X + E_L + E_R;//Constraint on Xi Energy
#endif
//		-E_X+E_L+E_R;//Constraint on Xi-> Lambda pi Energy
//		-E_X*E_X+pow(E_L+E_R,2);//Constraint on Xi-> Lambda pi Energy
	// u1 = p_L, u2 = th_L, u3 = ph_L, u4 = p_X, u5 = th_X, u6 = ph_X
	// m1 = p_P, m2 = th_P, m3 = ph_P, m4 = p_Q, m5 = th_Q, m6 = ph_Q, m7 = p_R, m8 = th_R, m9 = ph_R
#else
	double f5 = 
		-p_X*sin(th_X)*cos(ph_X)
		+p_P*sin(th_P)*cos(ph_P) 
		+p_Q*sin(th_Q)*cos(ph_Q) 
		+p_R*sin(th_R)*cos(ph_R); //Constraint on Xi x momentum
	double f6 = 
		-p_X*sin(th_X)*sin(ph_X)
		+p_P*sin(th_P)*sin(ph_P) 
		+p_Q*sin(th_Q)*sin(ph_Q) 
		+p_R*sin(th_R)*sin(ph_R);//Constraint on Xi y momentum
	double f7 =  
		-p_X*cos(th_X)
		+p_P*cos(th_P)
		+p_Q*cos(th_Q)
		+p_R*cos(th_R);//Constraint on Xi z momentum 
	double f8	=
#if InvMassFit
		mX*mX - (pow(E_P + E_Q + E_R,2)- pow(hypot(p_Px+p_Qx+p_Rx,hypot(p_Py+p_Qy+p_Ry,p_Pz+p_Qz+p_Rz)),2));//Constraint on Xi Energy
#else
		-E_X + E_P +E_Q + E_R;//Constraint on Xi Energy
#endif
//		-E_X*E_X+pow(E_P+E_Q+E_R,2);//Constraint on Xi Energy
#endif
	double df1du1 = -sin(th_L)*cos(ph_L);
	double df1du2 = -p_L*cos(th_L)*cos(ph_L);
	double df1du3 = p_L*sin(th_L)*sin(ph_L);
	double df1du4 = 0;
	double df1du5 = 0;
	double df1du6 = 0;
	
	double df1dm1 = sin(th_P)*cos(ph_P);
	double df1dm2 = p_P*cos(th_P)*cos(ph_P);
	double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df1dm4 = sin(th_Q)*cos(ph_Q);
	double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);
	double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	double df1dm7 = 0;
	double df1dm8 = 0;
	double df1dm9 = 0;


	double df2du1 = -sin(th_L)*sin(ph_L);
	double df2du2 = -p_L*cos(th_L)*sin(ph_L);
	double df2du3 = -p_L*sin(th_L)*cos(ph_L);
	double df2du4 = 0;
	double df2du5 = 0;
	double df2du6 = 0;

	double df2dm1 = sin(th_P)*sin(ph_P);
	double df2dm2 = p_P*cos(th_P)*sin(ph_P);
	double df2dm3 = p_P*sin(th_P)*cos(ph_P);
	double df2dm4 = sin(th_Q)*sin(ph_Q);
	double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);
	double df2dm6 = p_Q*sin(th_Q)*cos(ph_Q);
	double df2dm7 = 0;
	double df2dm8 = 0;
	double df2dm9 = 0;


	double df3du1 = -cos(th_L);
	double df3du2 = p_L*sin(th_L);
	double df3du3 = 0;
	double df3du4 = 0;
	double df3du5 = 0;
	double df3du6 = 0;

	double df3dm1 = cos(th_P);
	double df3dm2 = -p_P*sin(th_P);
	double df3dm3 = 0;
	double df3dm4 = cos(th_Q);
	double df3dm5 = -p_Q*sin(th_Q);
	double df3dm6 = 0;
	double df3dm7 = 0;
	double df3dm8 = 0;
	double df3dm9 = 0;
#if InvMassFit
	double df4du1 = 0;
	double df4du2 = 0;
	double df4du3 = 0;
	double df4du4 = 0;
	double df4du5 = 0;
	double df4du6 = 0;
	
	double df4dm1 = -2*( (E_P+E_Q)*dEdp_P
		- (p_Px+p_Qx)*cos(th_P)*cos(ph_P)
		- (p_Py+p_Qy)*cos(th_P)*sin(ph_P)
		- (p_Pz+p_Qz)*cos(th_P) 
	);  
	double df4dm2 = -2*(
		-(p_Px+p_Qx)*p_P*cos(th_P)*cos(ph_P)
		-(p_Py+p_Qy)*p_P*cos(th_P)*sin(ph_P)
		+(p_Pz+p_Qz)*p_P*sin(th_P)
	);
	double df4dm3 = -2*(
		-(p_Px+p_Qx)*p_P*sin(th_P)*sin(ph_P)
		+(p_Py+p_Qy)*p_P*sin(th_P)*cos(ph_P)
	);
	
	
	double df4dm4 = -2*( (E_P+E_Q)*dEdp_Q
		- (p_Px+p_Qx)*sin(th_Q)*cos(ph_Q)
		- (p_Py+p_Qy)*sin(th_Q)*sin(ph_Q)
		- (p_Pz+p_Qz)*cos(th_Q) 
	);
	double df4dm5 = -2*(
		-(p_Px+p_Qx)*p_Q*cos(th_Q)*cos(ph_Q)
		-(p_Py+p_Qy)*p_Q*cos(th_Q)*sin(ph_Q)
		+(p_Pz+p_Qz)*p_Q*sin(th_Q)
	);
	double df4dm6 = -2*(
		-(p_Px+p_Qx)*p_Q*sin(th_Q)*sin(ph_Q)
		+(p_Py+p_Qy)*p_Q*sin(th_Q)*cos(ph_Q)
	);
	double df4dm7 = 0;
	double df4dm8 = 0;
	double df4dm9 = 0;
#else
	double df4du1 = -dEdp_L;
	double df4du2 = 0;
	double df4du3 = 0;
	double df4du4 = 0;
	double df4du5 = 0;
	double df4du6 = 0;

	double df4dm1 = dEdp_P;
	double df4dm2 = 0;
	double df4dm3 = 0;
	double df4dm4 = dEdp_Q;
	double df4dm5 = 0;
	double df4dm6 = 0;
	double df4dm7 = 0;
	double df4dm8 = 0;
	double df4dm9 = 0;
#endif
//	double df4du1 = -2*dEdp_L*E_L	;
/*
*/

#if LambdaFit
	double df5du1 = sin(th_L)*cos(ph_L);
	double df5du2 = p_L*cos(th_L)*cos(ph_L);
	double df5du3 = -p_L*sin(th_L)*sin(ph_L);
	double df5du4 = -sin(th_X)*cos(ph_X);
	double df5du5 = -p_X*cos(th_X)*cos(ph_X);
	double df5du6 = p_X*sin(th_X)*sin(ph_X);

	double df5dm1 = 0;
	double df5dm2 = 0;// d/ dth_P
	double df5dm3 = 0;
	double df5dm4 = 0;
	double df5dm5 = 0;// d/ dth_Q
	double df5dm6 = 0;
	double df5dm7 = sin(th_R)*cos(ph_R);//df5 / d(P_R)
	double df5dm8 = p_R*cos(th_R)*cos(ph_R);// d/ dth_R
	double df5dm9 = -p_R*sin(th_R)*sin(ph_R);//df5 / d(Ph_R)

	double df6du1 = sin(th_L)*sin(ph_L);
	double df6du2 = p_L*cos(th_L)*sin(ph_L);
	double df6du3 = p_L*sin(th_L)*cos(ph_L);
	double df6du4 = -sin(th_X)*sin(ph_X);
	double df6du5 = -p_X*cos(th_X)*sin(ph_X);
	double df6du6 = -p_X*sin(th_X)*cos(ph_X);

	double df6dm1 = 0;
	double df6dm2 = 0;// d/ dth_P
	double df6dm3 = 0;
	double df6dm4 = 0;
	double df6dm5 = 0;// d/ dth_Q
	double df6dm6 = 0;
	double df6dm7 = sin(th_R)*sin(ph_R);
	double df6dm8 = p_R*cos(th_R)*sin(ph_R);// d/ dth_R
	double df6dm9 = p_R*sin(th_R)*cos(ph_R);


	double df7du1 = cos(th_L);
	double df7du2 = -p_L*sin(th_L);
	double df7du3 = 0;
	double df7du4 = -cos(th_X);
	double df7du5 = p_X*sin(th_X);
	double df7du6 = 0;

	double df7dm1 = 0; 
	double df7dm2 = 0; 
	double df7dm3 = 0;
	double df7dm4 = 0; 
	double df7dm5 = 0; 
	double df7dm6 = 0;
	double df7dm7 = cos(th_R);
	double df7dm8 = -p_R*sin(th_R);
	double df7dm9 = 0;
#if InvMassFit	
	double df8du1 = -2*( (E_L+E_X)*dEdp_L
		- (p_Lx+p_Xx)*cos(th_L)*cos(ph_L)
		- (p_Ly+p_Xy)*cos(th_L)*sin(ph_L)
		- (p_Lz+p_Xz)*cos(th_L) 
	);  
	double df8du2 = -2*(
		-(p_Lx+p_Rx)*p_L*cos(th_L)*cos(ph_L)
		-(p_Ly+p_Ry)*p_L*cos(th_L)*sin(ph_L)
		+(p_Lz+p_Rz)*p_L*sin(th_L)
	);
	double df8du3 = -2*(
		-(p_Lx+p_Rx)*p_L*sin(th_L)*sin(ph_L)
		+(p_Ly+p_Ry)*p_L*sin(th_L)*cos(ph_L)
	);
	double df8du4 = 0;
	double df8du5 = 0;
	double df8du6 = 0;
	
	double df8dm1 = 0;
	double df8dm2 = 0;
	double df8dm3 = 0;
	double df8dm4 = 0;
	double df8dm5 = 0;
	double df8dm6 = 0;
	double df8dm7 = -2*( (E_L+E_R)*dEdp_R
		- (p_Lx+p_Rx)*sin(th_R)*cos(ph_R)
		- (p_Ly+p_Ry)*sin(th_R)*sin(ph_R)
		- (p_Lz+p_Rz)*cos(th_R) 
	);
	double df8dm8 = -2*(
		-(p_Lx+p_Rx)*p_R*cos(th_R)*cos(ph_R)
		-(p_Ly+p_Ry)*p_R*cos(th_R)*sin(ph_R)
		+(p_Lz+p_Rz)*p_R*sin(th_R)
	);
	double df8dm9 = -2*(
		-(p_Lx+p_Rx)*p_R*sin(th_R)*sin(ph_R)
		+(p_Ly+p_Ry)*p_R*sin(th_R)*cos(ph_R)
	);
#else
	double df8du1 = dEdp_L;
	double df8du2 = 0;
	double df8du3 = 0;
	double df8du4 = -dEdp_X;
	double df8du5 = 0;
	double df8du6 = 0;

	double df8dm1 = 0;
	double df8dm2 = 0;
	double df8dm3 = 0;
	double df8dm4 = 0;
	double df8dm5 = 0;
	double df8dm6 = 0;
	double df8dm7 = dEdp_R;
	double df8dm8 = 0;
	double df8dm9 = 0;	
#endif
#else
	double df5du1 = 0; 
	double df5du2 = 0; 
	double df5du3 = 0; 
	double df5du4 = -sin(th_X)*cos(ph_X);
	double df5du5 = -p_X*cos(th_X)*cos(ph_X);
	double df5du6 = p_X*sin(th_X)*sin(ph_X);

	double df5dm1 = sin(th_P)*cos(ph_P);
	double df5dm2 = p_P*cos(th_P)*cos(ph_P);
	double df5dm3 = -p_P*sin(th_P)*sin(ph_P);
	double df5dm4 = sin(th_Q)*cos(ph_Q);
	double df5dm5 = p_Q*cos(th_Q)*cos(ph_Q);
	double df5dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
	double df5dm7 = sin(th_R)*cos(ph_R);//df5 / d(P_R)
	double df5dm8 = p_R*cos(th_R)*cos(ph_R);// d/ dth_R
	double df5dm9 = -p_R*sin(th_R)*sin(ph_R);//df5 / d(Ph_R)

	double df6du1 = 0; 
	double df6du2 = 0; 
	double df6du3 = 0; 
	double df6du4 = -sin(th_X)*sin(ph_X);
	double df6du5 = -p_X*cos(th_X)*sin(ph_X);
	double df6du6 = -p_X*sin(th_X)*cos(ph_X);

	double df6dm1 = sin(th_P)*sin(ph_P);
	double df6dm2 = p_P*cos(th_P)*sin(ph_P);
	double df6dm3 = p_P*sin(th_P)*cos(ph_P);
	double df6dm4 = sin(th_Q)*sin(ph_Q);
	double df6dm5 = p_Q*cos(th_Q)*sin(ph_Q);
	double df6dm6 = p_Q*sin(th_Q)*cos(ph_Q);
	double df6dm7 = sin(th_R)*sin(ph_R);
	double df6dm8 = p_R*cos(th_R)*sin(ph_R);// d/ dth_R
	double df6dm9 = p_R*sin(th_R)*cos(ph_R);


	double df7du1 = 0; 
	double df7du2 = 0; 
	double df7du3 = 0; 
	double df7du4 = -cos(th_X);
	double df7du5 = p_X*sin(th_X);
	double df7du6 = 0;

	double df7dm1 = cos(th_P);
	double df7dm2 = -p_P*sin(th_P);
	double df7dm3 = 0;
	double df7dm4 = cos(th_Q); 
	double df7dm5 = -p_Q*sin(th_Q);
	double df7dm6 = 0;
	double df7dm7 = cos(th_R);
	double df7dm8 = -p_R*sin(th_R);
	double df7dm9 = 0;
	
	double df8du1 = 0;
	double df8du2 = 0;
	double df8du3 = 0;
	double df8du4 = -dEdp_X;
	double df8du5 = 0;
	double df8du6 = 0;

	double df8dm1 = dEdp_P;
	double df8dm2 = 0;
	double df8dm3 = 0;
	double df8dm4 = dEdp_Q;
	double df8dm5 = 0;
	double df8dm6 = 0;
	double df8dm7 = dEdp_R;
	double df8dm8 = 0;
	double df8dm9 = 0;	
#endif
	double fs[]={f1,f2,f3,f4,f5,f6,f7,f8};
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
			df6dm1,df6dm2,df6dm3,df6dm4	,df6dm5,df6dm6,df6dm7,df6dm8,df6dm9,
			df7dm1,df7dm2,df7dm3,df7dm4	,df7dm5,df7dm6,df7dm7,df7dm8,df7dm9,
			df8dm1,df8dm2,df8dm3,df8dm4	,df8dm5,df8dm6,df8dm7,df8dm8,df8dm9
		};
		for(int i=0;i<nMeas*nConst;++i){
			dfdms[i]=temp[i];
		};
		double tempu[] = {
			df1du1,df1du2,df1du3,df1du4,df1du5,df1du6,
			df2du1,df2du2,df2du3,df2du4,df2du5,df2du6,
			df3du1,df3du2,df3du3,df3du4,df3du5,df3du6,
			df4du1,df4du2,df4du3,df4du4,df4du5,df4du6,
			df5du1,df5du2,df5du3,df5du4,df5du5,df5du6,
			df6du1,df6du2,df6du3,df6du4,df6du5,df6du6,
			df7du1,df7du2,df7du3,df7du4,df7du5,df7du6,
			df8du1,df8du2,df8du3,df8du4,df8du5,df8du6
		};
		for(int i=0;i<nUnkn*nConst;++i){
			dfdus[i]=tempu[i];
		};
#if Hessian
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
#endif
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
	double Ppx = p_P*sin(th_P)*cos(ph_P);
	double Ppy = p_P*sin(th_P)*sin(ph_P);
	double Ppz = p_P*cos(th_P);
	double Qpx = p_Q*sin(th_Q)*cos(ph_Q);
	double Qpy = p_Q*sin(th_Q)*sin(ph_Q);
	double Qpz = p_Q*cos(th_Q);
	double Rpx = p_R*sin(th_R)*cos(ph_R);
	double Rpy = p_R*sin(th_R)*sin(ph_R);
	double Rpz = p_R*cos(th_R);
	double Lpx = p_L*sin(th_L)*cos(ph_L);
	double Lpy = p_L*sin(th_L)*sin(ph_L);
	double Lpz = p_L*cos(th_L);
	double Xpx = p_X*sin(th_X)*cos(ph_X);
	double Xpy = p_X*sin(th_X)*sin(ph_X);
	double Xpz = p_X*cos(th_X);
	
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
	//Does Kinematic fitting for R -> P+Q decay.
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
