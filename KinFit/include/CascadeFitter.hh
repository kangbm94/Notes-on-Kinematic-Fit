#ifndef CascadeFitter_h
#define CascadeFitter_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class CascadeFitter: virtual public KinematicFitter{
	// P: proton, Q: pion, R: pion, L: Lambda, X: Xi
	//L -> P + Q;
	//X -> L + R;
	protected:	
		TLorentzVector P;
		TVector3 Pres;
		TLorentzVector PCor;
		double mP;

		TLorentzVector Q;
		TVector3 Qres;
		TLorentzVector QCor;
		double mQ;
		
		TLorentzVector R;
		TVector3 Rres;
		TLorentzVector RCor;
		double mR;

		TLorentzVector L;
		TVector3 Lres;
		TLorentzVector LCor;
		vector<double> MassDiffsL;
		double mL;
		
		
		TLorentzVector X;
		TVector3 Xres;
		TLorentzVector XCor;
		vector<double> MassDiffsX;
		double mX;

		bool UseVertexFlag = 0;
	public:
		CascadeFitter(){}
		CascadeFitter(TLorentzVector P_,TLorentzVector Q_,TLorentzVector R_);
		void SetInvMass(double ML,double MX){
			mL = ML;
			mX = MX;
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {PCor,QCor,RCor,LCor,XCor};	
			return ret;
		}
		void UseVertex(bool status,TVector3 Vert1,TVector3 Vert2);
		void ToDecayPlane();
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif
