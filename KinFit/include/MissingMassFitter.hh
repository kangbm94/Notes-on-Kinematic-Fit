#ifndef MissingMassFitter_h
#define MissingMassFitter_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class MissingMassFitter: virtual public KinematicFitter{
	//R -> P + Q;
    //B + Target -> R + X  B : beam
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
		vector<double> MassDiffs;
		double mR;
		
        TLorentzVector B;
        TLorentzVector BCor;
		double mB;
        TLorentzVector Target;
        
        double MM;
        TLorentzVector Missing;
        TLorentzVector MissingCor;
	public:
		MissingMassFitter(){}
		MissingMassFitter(
            TLorentzVector B_,TLorentzVector Target_, double MissingMass,
            TLorentzVector P_,TLorentzVector Q_,TLorentzVector R_);
		void SetInvMass(double IM){
			mR = IM;
		}
		void SetMissingMass(double MissingMass){
			MM = MissingMass;
		}
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {BCor,MissingCor,PCor,QCor,RCor};
			return ret;
		}
	protected:
		virtual void Initialize();
		virtual void SampleStepPoint(int steps);
		virtual void SetConstraints();
		virtual void Rotate();
		TMatrixD JacobianSphToCart(double p, double th, double ph);		
		virtual void CalcVariance(int istep);
};
#endif

