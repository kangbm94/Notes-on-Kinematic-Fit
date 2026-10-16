//#include "KinFit.cc"
#include "../include/BalanceFitter.hh"
#include "KinFit.cc"
#ifndef BalanceFitter_cc
#define BalanceFitter_cc
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit

void BalanceFitter::Initialize(){
  Clear();
  nMeas = 3;nUnkn = 0; nConst = 4;
  double meas[100];
  meas[0] = B.P();
  meas[1] = B.Theta();
  meas[2] = B.Phi();
  for(int i=0;i<nPs;++i){
    nMeas += 3;
    meas[3*(i+1)+0] = Ps.at(i).P();
    meas[3*(i+1)+1] = Ps.at(i).Theta();
    meas[3*(i+1)+2] = Ps.at(i).Phi();
  }

  
  TMatrixD Meas0(nMeas,1,meas);  
  vector<double>Pull;
  Pull.resize(nMeas);
  Measurements.push_back(Meas0);
  Pulls.push_back(Pull);
  Chi2s.push_back(-1);
  Initialized = true;
}
void BalanceFitter::SetConstraints(){
  auto Meas = Measurements.at(step); 

  double p_B = Meas(0,0);
  double th_B = Meas(1,0);
  double ph_B = Meas(2,0);
  vector<double> p_Ps;
  vector<double> th_Ps;
  vector<double> ph_Ps;
  for(int i=0;i<nPs;++i){
    p_Ps.push_back(Meas(3*i+3,0));
    th_Ps.push_back(Meas(3*i+4,0));
    ph_Ps.push_back(Meas(3*i+5,0));
  }
  
  
  double f1 = -p_B*sin(th_B)*cos(ph_B);
  for(int i=0;i<nPs;++i){
    f1 += p_Ps.at(i)*sin(th_Ps.at(i))*cos(ph_Ps.at(i));
  }//Constraint on x momentum  
  double f2 = -p_B*sin(th_B)*sin(ph_B); 
  for(int i=0;i<nPs;++i){
    f2 += p_Ps.at(i)*sin(th_Ps.at(i))*sin(ph_Ps.at(i));
  }//Constraint on y momentum
  double f3 = -p_B*cos(th_B);
  for(int i=0;i<nPs;++i){
    f3 += p_Ps.at(i)*cos(th_Ps.at(i));
  }//Constraint on z momentum
  double f4  = -sqrt(p_B*p_B+mB*mB) - mTarget;
  for(int i=0;i<nPs;++i){
    f4 += sqrt(p_Ps.at(i)*p_Ps.at(i)+mPs.at(i)*mPs.at(i));
  }//Constraint on energy


  vector<double> df1dms;
  vector<double> df2dms;
  vector<double> df3dms;
  vector<double> df4dms;

  df1dms.push_back(-sin(th_B)*cos(ph_B));// d/ dp_B
  df1dms.push_back(-p_B*cos(th_B)*cos(ph_B));// d/ dth_B
  df1dms.push_back(p_B*sin(th_B)*sin(ph_B));// d/ dph_B
  for(int i=0;i<nPs;++i){
    df1dms.push_back(sin(th_Ps.at(i))*cos(ph_Ps.at(i)));// d/ dp_P
    df1dms.push_back(p_Ps.at(i)*cos(th_Ps.at(i))*cos(ph_Ps.at(i)));// d/ dth_P
    df1dms.push_back(-p_Ps.at(i)*sin(th_Ps.at(i))*sin(ph_Ps.at(i)));// d/ dph_P
  }

  df2dms.push_back(-sin(th_B)*sin(ph_B));// d/ dp_B
  df2dms.push_back(-p_B*cos(th_B)*sin(ph_B));// d/ dth_B
  df2dms.push_back(-p_B*sin(th_B)*cos(ph_B));// d/ dph_B
  for(int i=0;i<nPs;++i){
    df2dms.push_back(sin(th_Ps.at(i))*sin(ph_Ps.at(i)));// d/ dp_P
    df2dms.push_back(p_Ps.at(i)*cos(th_Ps.at(i))*sin(ph_Ps.at(i)));// d/ dth_P
    df2dms.push_back(p_Ps.at(i)*sin(th_Ps.at(i))*cos(ph_Ps.at(i)));// d/ dph_P
  }

  df3dms.push_back(-cos(th_B));// d/ dp_B
  df3dms.push_back(p_B*sin(th_B));// d/ dth_B
  df3dms.push_back(0);// d/ dph_B
  for(int i=0;i<nPs;++i){
    df3dms.push_back(cos(th_Ps.at(i)));// d/ dp_P
    df3dms.push_back(-p_Ps.at(i)*sin(th_Ps.at(i)));// d/ dth_P
    df3dms.push_back(0);// d/ dph_P
  }

  df4dms.push_back(-p_B/sqrt(p_B*p_B+mB*mB));// d/ dp_B
  df4dms.push_back(0);// d/ dth_B
  df4dms.push_back(0);// d/ dph_B
  for(int i=0;i<nPs;++i){
    double p_P = p_Ps.at(i);
    double m_P = mPs.at(i);
    df4dms.push_back(p_P/sqrt(p_P*p_P+m_P*m_P));// d/ dp_P
    df4dms.push_back(0);// d/ dth_P
    df4dms.push_back(0);// d/ dph_P
  }

  double fs[]={f1,f2,f3,f4};
  double dfdms[200] ;
  double dfdus[200] ;

  vector<vector<double>> dfdm_vec = {df1dms,df2dms,df3dms,df4dms};

  TMatrixD FMat(nConst,1,fs);
  TMatrixD dFdM(nConst,nMeas,dfdms);
  for(int i=0;i<nConst;++i){
    for(int j=0;j<nMeas;++j){
      dFdM(i,j) = dfdm_vec.at(i).at(j);
    }
  }
  FMats.push_back(FMat);//Constraint Matrices for Each step
  dFdMs.push_back(dFdM);//Constraint matrix differentiated by measurement params.
}
void BalanceFitter::SampleStepPoint(int steps){
  auto Meas = Measurements.at(steps); 
  double p_B = Meas(0,0);
  double th_B = Meas(1,0);
  double ph_B = Meas(2,0);
  TVector3 VB = p_B*TVector3(sin(th_B)*cos(ph_B),sin(th_B)*sin(ph_B),cos(th_B));
  TLorentzVector BB(VB,hypot(p_B,mB));
  BCor = BB;
  for(int i=0;i<nPs;++i){
    double p_P = Meas(3*i+3,0);
    double th_P = Meas(3*i+4,0);
    double ph_P = Meas(3*i+5,0);
    TVector3 VP = p_P*TVector3(sin(th_P)*cos(ph_P),sin(th_P)*sin(ph_P),cos(th_P));
    TLorentzVector PP(VP,hypot(p_P,mPs.at(i)));
    PCor.at(i) = PP;
  }
}
BalanceFitter::BalanceFitter(TLorentzVector B_,TLorentzVector Target_,
  vector<TLorentzVector> Ps_){
  //Does Kinematic fitting for R -> P+Q decay.
  B = B_;
  Target = Target_;
  Ps = Ps_;
  nPs = Ps.size();
  mB = B.M();
  mTarget = Target.M();
  for(int i=0;i<nPs;++i){
    mPs.push_back(Ps.at(i).M());
    PCor.push_back(Ps.at(i));
  }
  Initialize();
};
TMatrixD
BalanceFitter::JacobianSphToCart(double p, double th, double ph){
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
BalanceFitter::CalcVariance(int istep){//legacy
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
BalanceFitter::Rotate(){
  auto VMat = Variancies.at(0);
  Initialize();
  Variancies.push_back(VMat);
  TMatrixD J;
  RotateVariance(J);
}
#endif
