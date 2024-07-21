#ifndef KinFit_h
#define KinFit_h
// Author: Kang Byungmin, kangbmw2@naver.com
// For the mathematics of the fitting, please refer to:
// https://github.com/kangbm94/Notes-on-Kinematic-Fit
#include <vector>
#include <TMatrixD.h>
#include <Math/ProbFunc.h>
using namespace std;
class KinematicFitter{
	protected:
		// R -> P + Q;
		int step = 0;
		int step1st = 0;
		int best_step = 0;
		int nConst;
		int nMeas;
		int nUnkn;
		double Chi2_cut = 0.1;
		double Best_Chi2 = 1e18;
		vector<double> best_pull ;
		vector<double> best_Upull ;
		int MaxStep = 100;
		double Best_Chi22nd = 1e18;
		double damping = 1.00;
		vector<TMatrixD>ScalingMats;
		
		vector<vector<double>> Pulls;
		vector<vector<double>> UPulls;
		vector<double> Lambdas;
		vector<double> Chi2s;
		vector<TMatrixD> Measurements;
		vector<TMatrixD> Unknowns;
		vector<TMatrixD> Variancies;
		vector<TMatrixD> VarianceInvs;
		vector<TMatrixD> dVMats;
		vector<TMatrixD> UHessians;//Maybe not used

		vector<TMatrixD> FMats;
		vector<TMatrixD> dFdMs;//\pdf{Constraint}{Measurement}
		vector<TMatrixD> dFdUs;//\pdf{Constraint }{Unknown }
		vector<TMatrixD> VarianciesU;
		vector<vector<TMatrixD>> d2Fd2Us;//\pdf^2 Constraints / dU_idU_j. This is a 3D Tensor object that is contracted with the lambda vector, resulting in a Hesian matrix.
		vector<TMatrixD> rMats;
		vector<TMatrixD> sMats;
		vector<double>	best_constraints;
		vector<double>	initial_constraints;
		bool UpdateVariancies = false;
		bool ScaleParams = true;

	public:
		KinematicFitter(){};
		//Setters
		void SetVariance(double* var);
		TMatrixD GetVariance(int i){
			return Variancies.at(i);
		}
		void AddOffdiagonals(TMatrixD Cov);//AddsCovaraince.
		void SetChi2DifCut(double cut){
			Chi2_cut = cut; }
		void SetMaximumStep(int max){
			MaxStep = max;
		}
		void UpdateVariance(bool status = true){
			UpdateVariancies = status;
		}
		void ScaleParameters(bool status = true){
			ScaleParams = status;
		}
		//Getters
		double GetLambda(int ent = -1){
			if(ent == -1) ent = step;
			return Lambdas.at(ent);
		}
		void Clear();
		virtual double DoKinematicFit(bool Do = true);
		double GetPValue(){
			return 1-ROOT::Math::chisquared_cdf(Best_Chi2,GetNDF());	
		};
		int GetNStep(){
			return step;
		}
		int GetBestStep(){
			return best_step;
		}
		vector<double>GetPull(){
			return best_pull;		
		}
		vector<double>GetUPull(){
			return best_Upull;		
		}
		TMatrixD GetUnmeasuredCovariance(){
			auto UCov = VarianciesU.at(best_step);
			return UCov;
			//return Best_Chi2*(UCov.Invert());
		}
		vector<double>GetStepChi2(){
			return Chi2s;
		}
		vector<vector<double>>GetStepPull(){
			return Pulls;		
		}
		int GetNDF(){
			return nConst - nUnkn;	
		}
		vector<double> GetKinematicConstraints(){
			return best_constraints;
		}
		vector<double> GetInitialConstraints(){
			return initial_constraints;
		}
	protected:
		//User Part: Set variables and constraints as you want.
		virtual void Initialize(){};
		virtual void SampleStepPoint(int steps){};
		virtual void SetConstraints(){};
		virtual void CalcVariance(int istep){};
		virtual void Rotate(){};//General methods for parameter transformation	
		//Core Part. Do not modify unless necessary
		void Finalize();	
		TMatrixD TransposeMatrix(TMatrixD M);		
		void ProcessStep();
		void RotateVariance(TMatrixD Jacobian);
};
#endif