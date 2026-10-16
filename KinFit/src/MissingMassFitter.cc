//#include "KinFit.cc"
#include "../include/MissingMassFitter.hh"
#include "KinFit.cc"
#ifndef MissingMassFitter_cc
#define MissingMassFitter_cc
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void MissingMassFitter::Initialize(){
  Clear();
  nMeas = 9;
  nUnkn = 3;
  nConst = 5; //For 2-C fit, when Vertex information is useless.
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

  mB = B.Mag();
  TVector3 TV_B = B.Vect();
  double p_B = TV_B.Mag();
  double th_B = TV_B.Theta();
  double ph_B = TV_B.Phi();

  Missing = B + Target - R;

  double meas[9];
  double unkn[3];
  double temp[] = {p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_B,th_B,ph_B};
  for(int i=0;i<nMeas;++i)meas[i]=temp[i];
  double temp2[] = {p_R,th_R,ph_R};
  for(int i=0;i<nUnkn;++i)unkn[i]=temp2[i];
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
  Initialized = true;
}
void MissingMassFitter::SetConstraints(){
  auto Meas = Measurements.at(step); 
  auto Unkn = Unknowns.at(step);
  double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_B,th_B,ph_B;
  p_R= Unkn(0,0); 
  th_R= Unkn(1,0); 
  ph_R= Unkn(2,0); 
  p_P= Meas(0,0); 
  th_P= Meas(1,0); 
  ph_P= Meas(2,0);
  p_Q= Meas(3,0); 
  th_Q= Meas(4,0); 
  ph_Q= Meas(5,0); 
  p_B= Meas(6,0); 
  th_B= Meas(7,0); 
  ph_B= Meas(8,0); 
  
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
  double f4  =
    - sqrt(p_R*p_R+mR*mR)
    + sqrt(p_P*p_P+mP*mP)
    + sqrt(p_Q*p_Q+mQ*mQ);//Constraint on Energy

  double pX_x = p_B*sin(th_B)*cos(ph_B) -p_R*sin(th_R)*cos(ph_R);
  double pX_y = p_B*sin(th_B)*sin(ph_B) -p_R*sin(th_R)*sin(ph_R);
  double pX_z = p_B*cos(th_B) -p_R*cos(th_R);
  double pX = sqrt(pX_x*pX_x+pX_y*pX_y+pX_z*pX_z);
  double eB = sqrt(p_B*p_B+mB*mB);
  double eX = eB + Target.Mag() - sqrt(p_R*p_R+mR*mR);


  double f5 = -eX + sqrt(pX*pX+MM*MM);//Constraint on Missing Mass 
  
  double df1du1 = -sin(th_R)*cos(ph_R);//df1 / d(P_R)
  double df1du2 = -p_R*cos(th_R)*cos(ph_R);//It could be df1/dm1 in case of 3-C fit.However, I didnt want to change the token... Mathematically it should be df1 / d (Th_R)
  double df1du3 = p_R*sin(th_R)*sin(ph_R);//df1 / d(Ph_R)
  double df1dm1 = sin(th_P)*cos(ph_P);//...
  double df1dm2 = p_P*cos(th_P)*cos(ph_P);// d/ dth_P
  double df1dm3 = -p_P*sin(th_P)*sin(ph_P);
  double df1dm4 = sin(th_Q)*cos(ph_Q);
  double df1dm5 = p_Q*cos(th_Q)*cos(ph_Q);// d/ dth_Q
  double df1dm6 = -p_Q*sin(th_Q)*sin(ph_Q);
  double df1dm7 = 0;
  double df1dm8 = 0;
  double df1dm9 = 0;
  
  double df2du1 = -sin(th_R)*sin(ph_R);
  double df2du2 = -p_R*cos(th_R)*sin(ph_R);// d/ dth_R
  double df2du3 = -p_R*sin(th_R)*cos(ph_R);
  double df2dm1 = sin(th_P)*sin(ph_P);
  double df2dm2 = p_P*cos(th_P)*sin(ph_P);// d/ dth_P
  double df2dm3 = p_P*sin(th_P)*cos(ph_P);
  double df2dm4 = sin(th_Q)*sin(ph_Q);
  double df2dm5 = p_Q*cos(th_Q)*sin(ph_Q);// d/ dth_Q
  double df2dm6 = p_Q*sin(th_Q)*cos(ph_Q);
  double df2dm7 = 0;
  double df2dm8 = 0;
  double df2dm9 = 0;

  double df3du1 = -cos(th_R);
  double df3du2 = p_R*sin(th_R);
  double df3du3 = 0;
  double df3dm1 = cos(th_P);
  double df3dm2 = -p_P*sin(th_P);
  double df3dm3 = 0;
  double df3dm4 = cos(th_Q);
  double df3dm5 = -p_Q*sin(th_Q);
  double df3dm6 = 0;
  double df3dm7 = 0;
  double df3dm8 = 0;
  double df3dm9 = 0;
  
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
  double df4dm7 = 0;
  double df4dm8 = 0;
  double df4dm9 = 0;

  double dpX_xdpR = -sin(th_R)*cos(ph_R);
  double dpX_ydpR = -sin(th_R)*sin(ph_R);
  double dpX_zdpR = -cos(th_R);
  double dpX_xdthR = -p_R*cos(th_R)*cos(ph_R);
  double dpX_ydthR = -p_R*cos(th_R)*sin(ph_R);
  double dpX_zdthR = p_R*sin(th_R);
  double dpX_xdphR = p_R*sin(th_R)*sin(ph_R);
  double dpX_ydphR = -p_R*sin(th_R)*cos(ph_R);
  double dpX_zdphR = 0;
  double deXdpR = -p_R/ER;
  double deXdthR = 0;
  double deXdphR = 0;
 
  double dpX_xdpB = sin(th_B)*cos(ph_B);
  double dpX_ydpB = sin(th_B)*sin(ph_B);
  double dpX_zdpB = cos(th_B);
  double dpX_xdthB = p_B*cos(th_B)*cos(ph_B);
  double dpX_ydthB = p_B*cos(th_B)*sin(ph_B);
  double dpX_zdthB =- p_B*sin(th_B);
  double dpX_xdphB = -p_B*sin(th_B)*sin(ph_B);
  double dpX_ydphB = p_B*sin(th_B)*cos(ph_B);
  double dpX_zdphB = 0;
  double deXdpB = p_B/eB;
  double deXdthB = 0;
  double deXdphB = 0;

  double df5du1 = p_R/sqrt(p_R*p_R+mR*mR) // deX/ dpR
  + (pX_x*dpX_xdpR + pX_y*dpX_ydpR + pX_z*dpX_zdpR)/sqrt(pX*pX+MM*MM);
  double df5du2 = (pX_x*dpX_xdthR + pX_y*dpX_ydthR + pX_z*dpX_zdthR)/sqrt(pX*pX+MM*MM);
  double df5du3 = (pX_x*dpX_xdphR + pX_y*dpX_ydphR + pX_z*dpX_zdphR)/sqrt(pX*pX+MM*MM);
  double df5dm1 = 0;
  double df5dm2 = 0;
  double df5dm3 = 0;
  double df5dm4 = 0;
  double df5dm5 = 0;
  double df5dm6 = 0;
  double df5dm7 = -p_B/eB + (pX_x*dpX_xdpB + pX_y*dpX_ydpB + pX_z*dpX_zdpB)/sqrt(pX*pX+MM*MM);
  double df5dm8 = (pX_x*dpX_xdthB + pX_y*dpX_ydthB + pX_z*dpX_zdthB)/sqrt(pX*pX+MM*MM);
  double df5dm9 = (pX_x*dpX_xdphB + pX_y*dpX_ydphB + pX_z*dpX_zdphB)/sqrt(pX*pX+MM*MM);

  double fs[]={f1,f2,f3,f4,f5};
  double dfdms[200] ;
  double dfdus[200] ;
  double temp[] = {
    df1dm1,df1dm2,df1dm3,df1dm4,df1dm5,df1dm6,df1dm7,df1dm8,df1dm9,
    df2dm1,df2dm2,df2dm3,df2dm4,df2dm5,df2dm6,df2dm7,df2dm8,df2dm9,
    df3dm1,df3dm2,df3dm3,df3dm4,df3dm5,df3dm6,df3dm7,df3dm8,df3dm9,
    df4dm1,df4dm2,df4dm3,df4dm4,df4dm5,df4dm6,df4dm7,df4dm8,df4dm9,
    df5dm1,df5dm2,df5dm3,df5dm4,df5dm5,df5dm6,df5dm7,df5dm8,df5dm9
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


  TMatrixD FMat(nConst,1,fs);
  TMatrixD dFdM(nConst,nMeas,dfdms);
  TMatrixD dFdU(nConst,nUnkn,dfdus);
  FMats.push_back(FMat);//Constraint Matrices for Each step
  dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
  dFdUs.push_back(dFdU);// same, but for unmeasured params.
}
void MissingMassFitter::SampleStepPoint(int steps){
  auto Meas = Measurements.at(steps); 
  auto Unkn = Unknowns.at(steps); 
  double p_P,th_P,ph_P,p_Q,th_Q,ph_Q,p_R,th_R,ph_R,p_B,th_B,ph_B;
  p_R= Unkn(0,0); 
  th_R= Unkn(1,0); 
  ph_R= Unkn(2,0); 
  p_P= Meas(0,0); 
  th_P= Meas(1,0); 
  ph_P= Meas(2,0);
  p_Q= Meas(3,0); 
  th_Q= Meas(4,0); 
  ph_Q= Meas(5,0); 
  p_B = Meas(6,0);
  th_B = Meas(7,0);
  ph_B = Meas(8,0);
  TVector3 PVector = p_P* TVector3(sin(th_P)*cos(ph_P),sin(th_P)*sin(ph_P),cos(th_P));
  TVector3 QVector = p_Q* TVector3(sin(th_Q)*cos(ph_Q),sin(th_Q)*sin(ph_Q),cos(th_Q));
  TVector3 RVector = p_R* TVector3(sin(th_R)*cos(ph_R),sin(th_R)*sin(ph_R),cos(th_R));
  TVector3 BVector = p_B* TVector3(sin(th_B)*cos(ph_B),sin(th_B)*sin(ph_B),cos(th_B));
  TLorentzVector PP(PVector,hypot(mP,p_P));
  TLorentzVector QQ(QVector,hypot(mQ,p_Q));
  TLorentzVector RR(RVector,hypot(mR,p_R));
  TLorentzVector BB(BVector,hypot(mB,p_B));
  auto V = PP + QQ;
  double MassDiff = V.Mag()-mR;
  PCor = PP;
  QCor = QQ;
  RCor = RR;
  BCor = BB;
  MissingCor = BB + Target - RR;
  MassDiffs.push_back(MassDiff);
}
MissingMassFitter::MissingMassFitter(
    TLorentzVector B_,TLorentzVector Target_, double MissingMass,
    TLorentzVector P_,TLorentzVector Q_, TLorentzVector R_){
  //Does Kinematic fitting for R -> P+Q decay.
  P=P_;
  Q=Q_;
  R=R_;
  B=B_;
  Target=Target_;
  MM = MissingMass;
  Initialize();
};
TMatrixD
MissingMassFitter::JacobianSphToCart(double p, double th, double ph){
  // x = p sin(th) cos(ph)
  // y = p sin(th) sin(ph)
  // z = p cos(th)
  //V_c = J^T V J|->
  //    dxdp, dxdth,dxdph
  //J  =  dydp, dydth,dydph
  //    dzdp, dzdth,dzdph

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
MissingMassFitter::CalcVariance(int istep){//legacy
    return;
  auto Meas = Measurements.at(istep); 
  auto Unkn = Unknowns.at(istep);
  double p_R,th_R,ph_R,p_P,th_P,ph_P,p_Q,th_Q,ph_Q;
    p_R= Unkn(0,0); 
    th_R= Unkn(1,0); 
    ph_R= Unkn(2,0); 
    p_P= Meas(0,0); 
    th_P= Meas(1,0); 
    ph_P= Meas(2,0);
    p_Q= Meas(3,0); 
    th_Q= Meas(4,0); 
    ph_Q= Meas(5,0); 
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
//  TMatrixD VMat_C = Jsc_PQ_T*(VMat)*Jsc_PQ;
  double El_contract[18]={//reduce matrix dimension
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
//  VarianciesU.push_back(UVMat);
}

void
MissingMassFitter::Rotate(){//legacy
  return;
  auto VMat = Variancies.at(0);
  Initialize();
  Variancies.push_back(VMat);
  TMatrixD J;
  RotateVariance(J);
}
#endif
