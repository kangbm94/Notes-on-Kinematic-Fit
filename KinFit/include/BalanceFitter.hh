#ifndef BalanceFitter_h
#define BalanceFitter_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include "KinFit.hh"
#include <TVector3.h>
#include <TLorentzVector.h>
class BalanceFitter: virtual public KinematicFitter{
	//B + Target -> P1 + P2 + ... + Pn
	protected:	

		TLorentzVector B;
		TLorentzVector BCor;
		double mB;

		TLorentzVector Target;//Stationary target
		double mTarget;

		vector<TLorentzVector> Ps;
		vector<TLorentzVector> PCor;
		vector<double> mPs;

		int nPs;

	public:
		BalanceFitter(){}
		BalanceFitter(TLorentzVector B_,TLorentzVector Target_,
			vector<TLorentzVector> Ps_);
		vector<TLorentzVector> GetFittedLV(){
			vector<TLorentzVector> ret = {BCor};
			for(auto& p : PCor) ret.push_back(p);
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
