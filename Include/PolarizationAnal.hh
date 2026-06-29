#include "Constants.hh"
#ifndef PolarizationAnal_hh
#define PolarizationAnal_hh
#define SecondDerivative 1
const int nAngle = 3;
const double CMDomain[nAngle + 1] = {1, 0.97, 0.95, 0.9};
double SysW_PhA = 0, SysW_Th = 0, SysW_PhC = 0, SysW_Pol = 0;
int nbinCM = 10;
double lbCM = 0.7, hbCM = 1;
//double XiParam[2] = {-0.401,-2.1};
double XiParam[2] = {-0.390,-1.2};//https://pdg.lbl.gov/2024/tables/contents_tables_baryons.html
//double XiParam[2] = {-0.390,-90};
double AXi = XiParam[0];
double BXi = sqrt(1 - AXi*AXi) * sin(XiParam[1]*acos(-1)/180);
double CXi = sqrt(1 - AXi*AXi) * cos(XiParam[1]*acos(-1)/180);

double LdParam[2] = {0.747,-6.5};
double ALd = LdParam[0];
double BLd = sqrt(1 - ALd*ALd) * sin(LdParam[1]*acos(-1)/180);
double CLd = sqrt(1 - ALd*ALd) * cos(LdParam[1]*acos(-1)/180);
		
static double KA = AXi*ALd,KB = acos(-1)/4*BXi*ALd,KC = acos(-1)/4*CXi*ALd;
namespace PolaSystematics{
    static int nbin = 10;
    double pld, pxi;
    double wl, wx;
    static void
    Chi2Lambda(int &npar, double *gin, double &f, double *par, int flag){
        double chi2 = 0;
        double P = par[0];
        double x0 = -1., x1 = 1.;
        double bw = (x1 - x0) / nbin;
        for (int i = 0; i < nbin; ++i){
            double x_i = x0 + (i + 0.5) * bw;
            double y_i = (1 + pld * x_i);
            double e_i = y_i * wl;
            chi2 += (y_i - (1 + P * x_i)) * (y_i - (1 + P * x_i)) / (e_i * e_i);
        }
        f = chi2;
    }
    static void
    Chi2Xi(int &npar, double *gin, double &f, double *par, int flag){
        double chi2 = 0;
        double P = par[0];
        double x0 = -1., x1 = 1.;
        double bw = (x1 - x0) / nbin;
        for (int i = 0; i < nbin; ++i){
            double x_i = x0 + (i + 0.5) * bw;
            double y_cth_i = (1 + pxi * AXi * x_i);
            double e_cth_i = y_cth_i * wx;
            double y_cphc_i = (1 + pxi * KC * x_i);
            double e_cphc_i = y_cphc_i * wx;
            chi2 += (y_cth_i - (1 + P * AXi * x_i)) * (y_cth_i - (1 + P * AXi * x_i)) / (e_cth_i * e_cth_i);
            chi2 += (y_cphc_i - (1 + P * KC * x_i)) * (y_cphc_i - (1 + P * KC * x_i)) / (e_cphc_i * e_cphc_i);
        }
        f = chi2;
    }

    void CalcPolaSys(double pl_, double px_, double& ple, double& pxe){
        wl = SysW_PhA;
        wx = SysW_Pol;
        pld = pl_;
        pxi = px_;
        TMinuit* Minuit = new TMinuit(1);
        Minuit->SetPrintLevel(-1);
        Minuit->SetFCN(Chi2Lambda);
        Minuit->DefineParameter(0, "P", pl_, 0.01, 0, 0);
        Minuit->Migrad();
        Minuit->Command("MINOS");
        double par, parerr;
        Minuit->GetParameter(0, par, parerr);
        ple = parerr;

        Minuit->SetFCN(Chi2Xi);
        Minuit->DefineParameter(0, "P", pxi, 0.01, 0, 0);
        Minuit->Migrad();
        Minuit->Command("MINOS");
        Minuit->GetParameter(0, par, parerr);
        pxe = parerr;
    }

}


namespace{
	double weight_cap = 100; 
	static vector<double> FitTh;
	static vector<double> FitPhA;//Z alpha
	static vector<double> FitPhB;//X beta
	static vector<double> FitPhC;//Y gamma
	static vector<double> FitThWeight;
	static vector<double> FitPhAWeight;
	static vector<double> FitPhBWeight;
	static vector<double> FitPhCWeight;
	static bool FitPhAFlag ;

	static vector<double> FitXTh;
	static vector<double> FitXPhA;
	static vector<double> FitXPhB;
	static vector<double> FitXPhC;
	static vector<double> FitYTh;
	static vector<double> FitYPhA;
	static vector<double> FitYPhB;
	static vector<double> FitYPhC;
	static vector<double> FitETh;
	static vector<double> FitEPhA;
	static vector<double> FitEPhB;
	static vector<double> FitEPhC;
	static double Fitdensity = 0;
	static double w_sum = 0; 
	static double w2_sum = 0; 

	static double Score2Sum_Lambda = 0;
	static double Hessian_Lambda = 0;
	static double Score2Sum_Xi = 0;
	static double Hessian_Xi = 0;

	static double Neff_Lambda = 0;

	static vector<double> FitThX;
	static vector<double> FitThXWeight;
	static vector<double> FitThZ;
	static vector<double> FitThZWeight;
	static double Score2Sum_ThX = 0;
	static double Hessian_ThX = 0;
	static double Score2Sum_ThZ = 0;
	static double Hessian_ThZ = 0;
	
	static void
	Chi2Function(	int& npar, double* gin, double& f, double* par, int iflag){
		double chi2 = 0;
		double P_xi = par[0];
		int nb = FitXTh.size();
		double bw = 2. / nb;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitYTh.size();++i){
			if(FitETh[i]==0) FitETh[i]=1.;
			if(FitEPhB[i]==0) FitEPhB[i]=1.;
			if(FitEPhC[i]==0) FitEPhC[i]=1.;
			double Th0 = FitXTh[i] - bw/2;
			double Th1 = FitXTh[i] + bw/2;
			double CorXTh = FitXTh[i];
			double CorXPhC = FitXPhC[i];
#if SlopeCorrection
			CorXTh = (FitXTh[i] + 1./3 *P_xi * AXi* (Th0*Th0 + Th0*Th1 + Th1*Th1)) / (1 + 1./2* P_xi * AXi * (FitXTh[i]));
			CorXPhC = (FitXPhC[i] + 1./3 *P_xi * KC* (Th0*Th0 + Th0*Th1 + Th1*Th1)) / (1 + 1./2* P_xi * KC * (FitXPhC[i]));
#endif
			chi2 += pow((Fitdensity*(1 +P_xi* AXi*CorXTh)-FitYTh[i])/FitETh[i],2);
#if IncludeBeta
			chi2 += pow((Fitdensity*(1 +P_xi* KB*FitXPhB[i])-FitYPhB[i])/FitEPhB[i],2);
#endif
			chi2 += pow((Fitdensity *(1 +P_xi* KC*FitXPhC[i])-FitYPhC[i])/FitEPhC[i],2);
			w_sum += FitThWeight[i];
			w2_sum += FitThWeight[i] * FitThWeight[i];
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f = chi2 * n_eff / nb;
	}
	static void
	Chi2FunctionLambda(	int& npar, double* gin, double& f, double* par, int iflag){
		double chi2 = 0;
		double P_l = par[0];
		int nb = FitXPhA.size();
		double bw = 2. / nb;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitYPhA.size();++i){
			if(FitEPhA[i]==0) FitEPhA[i]=1.;
			double CorXPhA = FitXPhA[i];
			double PhA0 = FitXPhA[i] - bw/2;
			double PhA1 = FitXPhA[i] + bw/2;
#if SlopeCorrection
			CorXPhA = (FitXPhA[i] + 1./3 *P_l * (PhA0*PhA0 + PhA0*PhA1 + PhA1*PhA1)) / (1 + 1./2* P_l * (FitXPhA[i]));
#endif
			chi2 += pow((Fitdensity*(1 + P_l*CorXPhA)-FitYPhA[i])/FitEPhA[i],2);
			w_sum += FitPhAWeight[i];
			w2_sum += FitPhAWeight[i] * FitPhAWeight[i];
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f = chi2* n_eff / nb;
	}
	
/*
	Error Estimations in Unbinned Log-Likelihood Fits:
	Define the (log) Likelihood l as
	L = sum_i w_i * log(f(x_i; theta)) = sum_i w_i l;
	where x_i are the datapoints, and theta is the parameter vector.
	We define the score vector as
	s_i = d l / d theta; 

	and Jacobian matrix J is the derivative of S with respect to theta:
	J = sum_i w_i s_i * w_i s_i^T;

	The Hessian matrix H is the second derivative of l with respect to theta:
	H = sum_i w_i * d^2l/dtheta^2.

	Note that H and S are the sum of weighted variables.
	We define the covariance as
	Cov = H^-1 * J * H^-1;

	Then, the error of the parameter can be estimated as sqrt(Cov).
	In Our case, since our parameter is a single variable(pXi or pLd), the covariance is reduced to a single value, and the error can be estimated as sqrt(Score2Sum)/Hessian.
*/
	static void
	LikelihoodFunction(	int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_xi = par[0];
		double norm = 0;
		double wnorm = 0;
		int np = FitTh.size();
		Score2Sum_Xi = 0;
		Hessian_Xi = 0;
		Neff_Lambda = 0;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitTh.size();++i){
			double Th_i = FitTh[i];
			double PhA_i = FitPhA[i];
			double PhB_i = FitPhB[i];
			double PhC_i = FitPhC[i];
			double wTh_i = FitThWeight[i];
			double wPhA_i = FitPhAWeight[i];
			double wPhB_i = FitPhBWeight[i];
			double wPhC_i = FitPhCWeight[i];
			double probTh = (1 + P_xi*AXi*cos( Th_i))/2;
			double probPhA =(1 + KA*cos(PhA_i))/2;
			double probPhB =(1 + P_xi*KB*cos(PhB_i))/2;
			double probPhC =(1 + P_xi*KC*cos(PhC_i))/2;
			if(wTh_i > weight_cap) continue;

			likelihood -= wTh_i*log(probTh);
			likelihood -= wPhC_i*log(probPhC);
		
			//derivative of log (F) is 1/F dFdx.
			double dProbTh = AXi*cos(Th_i)/2;
			double dProbPhC = KC*cos(PhC_i)/2;
		
			double ddProbTh = 0;
			double ddProbPhC = 0;

			double dLikelihood_dTh = dProbTh / probTh;
			double dLikelihood_dPhC = dProbPhC / probPhC;
			double d2Likelihood_dTh2 = (ddProbTh / probTh - dProbTh * dProbTh/(probTh * probTh)) ;
			double d2Likelihood_dPhC2 = (ddProbPhC / probPhC - dProbPhC * dProbPhC/(probPhC * probPhC)) ;
			double d2Likelihood_dThdPhC = 0; //Assume that Th and PhC are independent.
			
			double Score_Xi = wTh_i * (
				-dLikelihood_dTh + -dLikelihood_dPhC
			); //Assume that weight for Th, PhB, PhC are the same.

			Hessian_Xi += wTh_i * (
				-d2Likelihood_dPhC2 - d2Likelihood_dTh2 - 2 * dLikelihood_dTh * dLikelihood_dPhC
			);
			#if SecondDerivative
			Score2Sum_Xi += wTh_i*wTh_i * (
				-d2Likelihood_dPhC2 - d2Likelihood_dTh2 - 2 * dLikelihood_dTh * dLikelihood_dPhC
			);
			#else
			Score2Sum_Xi += Score_Xi * Score_Xi;
			#endif
#if IncludeBeta
			likelihood -= wPhB_i*log(probPhB);
			Score_Xi += wPhB_i * ( (KB * cos(PhB_i)) / (probPhB) );
			Hessian_Xi += wPhB_i * ( (KB * KB * cos(PhB_i) * cos(PhB_i)) / (probPhB * probPhB) );
#endif
			w_sum += wTh_i;
			w2_sum += wTh_i * wTh_i;
		//	if(PhAFlag)likelihood -= wPh3_i*log(probPh3);
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f=likelihood* n_eff / np;
	}
	static void
	LikelihoodFunctionLambda(int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_ld = par[0];
		int np = FitPhA.size();
		Score2Sum_Lambda = 0;
		Hessian_Lambda = 0;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitPhA.size();++i){
			double PhA_i = FitPhA[i];
			double wPhA_i = FitPhAWeight[i];
			if(wPhA_i > weight_cap) continue;
			double probPhA = (1 + P_ld*cos(PhA_i))/2;
			likelihood -= wPhA_i*log(probPhA);
			w_sum += wPhA_i;
			w2_sum += wPhA_i * wPhA_i;

			double dProbPhA = cos(PhA_i)/2;
			double ddProbPhA = 0;
			
			double dLikelihood_dPhA = dProbPhA / probPhA ;
			double d2Likelihood_dPhA2 = (ddProbPhA / probPhA - dProbPhA * dProbPhA/(probPhA * probPhA)) ;
	
			double Score_Lambda = wPhA_i * ( -dLikelihood_dPhA ); //Assume that weight for Th, PhB, PhC are the same.;
			Hessian_Lambda += wPhA_i * ( -d2Likelihood_dPhA2 );
			#if SecondDerivative
			Score2Sum_Lambda += wPhA_i * wPhA_i * ( -d2Likelihood_dPhA2 );
			#else
			Score2Sum_Lambda += Score_Lambda * Score_Lambda;
			#endif
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f=likelihood* n_eff / np;
	}
	static void
	LikelihoodFunctionPolaX(int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_X = par[0];
		int np = FitPhA.size();
		Score2Sum_ThX = 0;
		Hessian_ThX = 0;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitThX.size();++i){
			double ThX_i = FitThX[i];
			double wThX_i = FitThXWeight[i];
			double probThX = (1 + P_X*AXi*cos(ThX_i))/2;
			if(wThX_i > weight_cap) continue;
			likelihood -= wThX_i*log(probThX);
			w_sum += wThX_i;
			w2_sum += wThX_i * wThX_i;
			double dProbThX = AXi*cos(ThX_i)/2;
			double ddProbThX = 0;
			double dLikelihood_dThX = dProbThX / probThX ;
			double d2Likelihood_dThX2 = (ddProbThX / probThX - dProbThX * dProbThX/(probThX * probThX)) ;
			double Score_ThX = wThX_i * ( dLikelihood_dThX);
			Hessian_ThX += wThX_i * ( -d2Likelihood_dThX2 );
			#if SecondDerivative
			Score2Sum_ThX += wThX_i * wThX_i * ( -d2Likelihood_dThX2 );
			#else
			Score2Sum_ThX += Score_ThX * Score_ThX;
			#endif
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f=likelihood* n_eff / np;
	}
	static void
	LikelihoodFunctionPolaZ(int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_Z = par[0];
		int np = FitPhA.size();
		Score2Sum_ThZ = 0;
		Hessian_ThZ = 0;
		w_sum = 0;
		w2_sum = 0;
		for(int i=0;i<FitThZ.size();++i){
			double ThZ_i = FitThZ[i];
			double wThZ_i = FitThZWeight[i];
			double probThZ = (1 + P_Z*AXi*cos(ThZ_i))/2;
			if(wThZ_i > weight_cap) continue;
			likelihood -= wThZ_i*log(probThZ);
			w_sum += wThZ_i;
			w2_sum += wThZ_i * wThZ_i;
			double dProbThZ = AXi*cos(ThZ_i)/2;
			double ddProbThZ = 0;
			double dLikelihood_dThZ = dProbThZ / probThZ ;
			double d2Likelihood_dThZ2 = (ddProbThZ / probThZ - dProbThZ * dProbThZ/(probThZ * probThZ)) ;
			double Score_ThZ = wThZ_i * ( dLikelihood_dThZ);
			Hessian_ThZ += wThZ_i * ( -d2Likelihood_dThZ2 );
			#if SecondDerivative
			Score2Sum_ThZ += wThZ_i * wThZ_i * ( -d2Likelihood_dThZ2 );
			#else
			Score2Sum_ThZ += Score_ThZ * Score_ThZ;
			#endif
		}
		double n_eff = w_sum * w_sum / w2_sum;
		f=likelihood* n_eff / np;
	}
}
class PolaConverter{
	private:
		double mXi = 1.32171;
		double mp = 938.272/1000;
		double mL = 1115.683/1000;
		TVector3 Km;
		TVector3 Kp;
		TVector3 Xi;
		TVector3 Ld;
		TVector3 P;
		TVector3 Pi1;
		TVector3 Pi2;
		TVector3 PolXi;// Y axis of the Xi polarization frame
		TVector3 PolLd;
		double cTh,Th,Ph,PhA,PhB,PhC;
		double ThX, ThZ;
		TVector3 XLd;
		TVector3 YLd;
		TVector3 ZLd;
		
		TVector3 PolXi_X;
		TVector3 PolXi_Z;

	public:
		PolaConverter(TVector3 PKm,TVector3 PKp,TVector3 PXi,TVector3 PLd,TVector3 PP, TVector3 PPi1, TVector3 PPi2,bool PiMode = 0){
			Km=PKm,Kp=PKp,Xi=PXi,Ld=PLd,P=PP,Pi1=PPi1,Pi2=PPi2;	
			double du = 0.045,dv = 0.096;
			/* 
			Km.RotateY(du);
			Km.RotateX(dv);
			Kp.RotateY(du);
			Kp.RotateX(dv);
			*/
			PolXi= NormalCross(Km,Kp);
			TLorentzVector LVXi(Xi,hypot(Xi.Mag(),mXi));
			TLorentzVector LVLd(Ld,hypot(Ld.Mag(),mL));
			TLorentzVector LVP(P,hypot(P.Mag(),mP));
			TLorentzVector LVPi1(Pi1,hypot(Pi1.Mag(),mPi));
			TLorentzVector LVPi2(Pi2,hypot(Pi2.Mag(),mPi));
			auto BoostXi = LVXi.BoostVector();
			auto LVLdCM = LVLd;
			auto LVPi2CM = LVPi2;
			LVLdCM.Boost(-BoostXi);
			auto LdCM = LVLdCM.Vect();
			LVPi2CM.Boost(-BoostXi);
			auto Pi2CM = LVPi2CM.Vect();
			if(PiMode)LdCM = -Pi2CM;
			cTh = cos(PolXi.Angle(LdCM));
			Th=acos(cTh);
			ZLd = LdCM*(1./LdCM.Mag());
//			auto LdV = LVLd.Vect();	
//			ZLd = LdV*(1./LdV.Mag());
			XLd = NormalCross(PolXi,ZLd);
			YLd = NormalCross(ZLd,XLd);
			double PLx = BXi*sin(Th)/(1+AXi*cTh);
			double PLy = CXi*sin(Th)/(1+AXi*cTh);
			double PLz = (AXi+cTh)/(1+AXi*cTh);
	
			PolLd = PLx*XLd+PLy*YLd+PLz*ZLd;
			auto BoostLd = LVLd.BoostVector();
			auto LVPCM = LVP;
			LVPCM.Boost(-BoostLd);
			auto PCM = LVPCM.Vect();//proton CM momentum

			auto LVPi1CM = LVPi1;
			LVPi1CM.Boost(-BoostLd);
			auto Pi1CM = LVPi1CM.Vect();
			if(PiMode)PCM = -Pi1CM;
			Ph = acos(NormalDot(PolLd,PCM));
			PhA = acos(NormalDot(ZLd,PCM));//alpha angle
			PhB = acos(NormalDot(XLd,PCM));//beta angle
			PhC = acos(NormalDot(YLd,PCM));//gamma angle


			PolXi_Z = Km.Unit();
			PolXi_X = NormalCross(PolXi,PolXi_Z);// x = y cross z
			ThX = PolXi_X.Angle(LdCM);
			ThZ = PolXi_Z.Angle(LdCM);

			if(isnan(cTh)
			 or isnan(PhB)
			 or isnan(PhC)
			 or isnan(PhA)
			){
				cout<<"Something Wrong!"<<endl;
				cout<<Form("cTh,PhA,PhB,PhC= (%.2g,%.2g,%.2g,%.2g)",cTh,PhA,PhB,PhC)<<endl; 
				cout<<Form("PKm = %f %f %f",PKm.X(),PKm.Y(),PKm.Z())<<endl;	
				cout<<Form("PKp = %f %f %f",PKp.X(),PKp.Y(),PKp.Z())<<endl;
				cout<<Form("PolXi = %f %f %f",PolXi.X(),PolXi.Y(),PolXi.Z())<<endl;
				cout<<Form("LdCM = %f %f %f, mag%g",LdCM.X(),LdCM.Y(),LdCM.Z(),LdCM.Mag())<<endl;
				
				cout<<Form("PXi = %f %f %f",PXi.X(),PXi.Y(),PXi.Z())<<endl;
				cout<<Form("PLd = %f %f %f",PLd.X(),PLd.Y(),PLd.Z())<<endl;
				cout<<Form("PP = %f %f %f",PP.X(),PP.Y(),PP.Z())<<endl;
			}
		}
		double GetTheta(){
			return Th;
		}
		double GetPhi(){
			return Ph;
		}
		double GetPhiA(){
			return PhA;
		}
		double GetPhiB(){
			return PhB;
		}
		double GetPhiC(){
			return PhC;
		}
		double GetThetaX(){
			return ThX;
		}
		double GetThetaZ(){
			return ThZ;
		}
		TVector3 GetPolAxis(){
			return PolXi;
		}
		TVector3 GetPolLambda(){
			return PolLd;
		}
		TVector3 GetPolAxisX(){
			return PolXi_X;
		}
		TVector3 GetPolAxisZ(){
			return PolXi_Z;
		}
};
class PolaFitter{
	private:
		TMinuit* Minuit = new TMinuit(1);
		TMinuit* MinuitChi2 = new TMinuit(1);
		TMinuit* MinuitLambda = new TMinuit(1);
		TMinuit* MinuitLambdaChi2 = new TMinuit(1);
		TMinuit* MinuitThX = new TMinuit(1);
		TMinuit* MinuitThZ = new TMinuit(1);
		bool FitXZ = 0;
	public:
		PolaFitter(){}
		PolaFitter(vector<double> Th,vector<double> PhA,vector<double> PhB,vector<double> PhC){
			FitTh = Th;
			FitPhA = PhA;
			FitPhB = PhB;
			FitPhC = PhC;
			FitPhAFlag = false;
		}
		void AddThXZ(vector<double> ThX,vector<double> ThZ){
			FitThX = ThX;
			FitThZ = ThZ;
			FitXZ = 1;
		}
		void SetWeight(vector<double> ThWeight,vector<double> PhAWeight,vector<double> PhBWeight,vector<double> PhCWeight){
			FitThWeight = ThWeight;
			FitPhAWeight = PhAWeight;
			FitPhBWeight = PhBWeight;
			FitPhCWeight = PhCWeight;
			FitThXWeight = ThWeight;
			FitThZWeight = ThWeight;
		}
		void SetWeightCap(double cap){
			weight_cap = cap;
		}
		void ClearGraphData(){
			FitXTh.clear();
			FitXPhA.clear();
			FitXPhB.clear();
			FitXPhC.clear();
			
			FitYTh.clear();
			FitYPhA.clear();
			FitYPhB.clear();
			FitYPhC.clear();
		
			FitETh.clear();
			FitEPhA.clear();
			FitEPhB.clear();
			FitEPhC.clear();
		} 
		void Clear(){
			FitTh.clear();
			FitPhB.clear();
			FitPhC.clear();
			FitPhA.clear();
			FitThX.clear();
			FitThZ.clear();
			FitThWeight.clear();
			FitPhBWeight.clear();
			FitPhCWeight.clear();
			FitPhAWeight.clear();
			FitThXWeight.clear();
			FitThZWeight.clear();
			ClearGraphData();
		}
		void SetHistos(TH1* hTh, TH1* hPhA, TH1* hPhB, TH1* hPhC){
			ClearGraphData();
			cout<<"Setting Histos"<<endl;
			int np = hTh->GetNbinsX();
			Fitdensity = hTh->Integral()/np;
			cout<<"Density = "<<Fitdensity<<endl;
			for(int ip=1;ip<np+1;++ip){
				double x,y,ex,ey;
				x = hTh->GetBinCenter(ip);
				y = hTh->GetBinContent(ip);
				ex = hTh->GetBinWidth(ip);
				ey = hTh->GetBinError(ip);
				FitXTh.push_back(x);
				FitYTh.push_back(y);
				FitETh.push_back(ey);

				x = hPhA->GetBinCenter(ip);
				y = hPhA->GetBinContent(ip);
				ex = hPhA->GetBinWidth(ip);
				ey = hPhA->GetBinError(ip);
				FitXPhA.push_back(x);
				FitYPhA.push_back(y);
				FitEPhA.push_back(ey);

				x = hPhB->GetBinCenter(ip);
				y = hPhB->GetBinContent(ip);
				ex = hPhB->GetBinWidth(ip);
				ey = hPhB->GetBinError(ip);
				FitXPhB.push_back(x);
				FitYPhB.push_back(y);
				FitEPhB.push_back(ey);

				x = hPhC->GetBinCenter(ip);
				y = hPhC->GetBinContent(ip);
				ex = hPhC->GetBinWidth(ip);
				ey = hPhC->GetBinError(ip);
				FitXPhC.push_back(x);
				FitYPhC.push_back(y);
				FitEPhC.push_back(ey);
			}
			cout<<"Setting Done"<<endl;
		}
		void SetGraphs(TGraphErrors* gTh,TGraphErrors* gPhA, TGraphErrors* gPhB, TGraphErrors* gPhC){
			cout<<"Setting Graphs"<<endl;
			ClearGraphData();
			int np = gTh->GetN();
			for(int ip=0;ip<np;++ip){
				double x,y;
				gTh->GetPoint(ip,x,y);
				double ex = gTh->GetErrorX(ip);
				double ey = gTh->GetErrorY(ip);
				FitXTh.push_back(x);
				FitYTh.push_back(y);
				FitETh.push_back(ey);
				
				gPhA->GetPoint(ip,x,y);
				ex = gPhA->GetErrorX(ip);
				ey = gPhA->GetErrorY(ip);
				FitXPhA.push_back(x);
				FitYPhA.push_back(y);
				FitEPhA.push_back(ey);
				
				gPhB->GetPoint(ip,x,y);
				ex = gPhB->GetErrorX(ip);
				ey = gPhB->GetErrorY(ip);
				FitXPhB.push_back(x);
				FitYPhB.push_back(y);
				FitEPhB.push_back(ey);
				
				gPhC->GetPoint(ip,x,y);
				ex = gPhC->GetErrorX(ip);
				ey = gPhC->GetErrorY(ip);
				FitXPhC.push_back(x);
				FitYPhC.push_back(y);
				FitEPhC.push_back(ey);
			}
			Fitdensity = 0;
			double YSum = 0;
			for(auto y:FitYTh){
				YSum+=y;
			}
			Fitdensity = YSum / np;
			cout<<"Den = "<<Fitdensity<<endl;
			cout<<"Setting Done"<<endl;
		}
		void IncludePhA(bool flag = true){
			FitPhAFlag = flag;
		}
		void DoFit(){
			if(FitThWeight.size()!= FitTh.size()){
				FitThWeight = vector<double>(FitTh.size(),1);
				FitPhBWeight = vector<double>(FitPhB.size(),1);
				FitPhCWeight = vector<double>(FitPhC.size(),1);
				FitPhAWeight = vector<double>(FitPhA.size(),1);
			}
			Minuit->SetPrintLevel(-1);
			Minuit->SetFCN(LikelihoodFunction);
			Minuit->DefineParameter(0,"P_xi",0.1,0.01,0,0);
			Minuit->Migrad();
			Minuit->Command("MINOS");
		}
		void DoChi2Fit(){
			MinuitChi2->SetPrintLevel(-1);
			MinuitChi2->SetFCN(Chi2Function);
			MinuitChi2->DefineParameter(0,"P_xi",0.0,0.01,0,0);
			MinuitChi2->Migrad();
		}
		void DoLambdaFit(){
			MinuitLambda->SetPrintLevel(-1);
			MinuitLambda->SetFCN(LikelihoodFunctionLambda);
			MinuitLambda->DefineParameter(0,"P_ld",-0.,0.001,0,0);
			MinuitLambda->Migrad();
			MinuitLambda->Command("MINOS");
		}
		void DoChi2FitLambda(){
			MinuitLambdaChi2->SetPrintLevel(-1);
			MinuitLambdaChi2->SetFCN(Chi2FunctionLambda);
			MinuitLambdaChi2->DefineParameter(0,"P_ld",-0.,0.001,0,0);
			MinuitLambdaChi2->Migrad();
		}
		void DoFitXZ(){
			MinuitThX->SetPrintLevel(-1);
			MinuitThX->SetFCN(LikelihoodFunctionPolaX);
			MinuitThX->DefineParameter(0,"P_x",0.1,0.01,0,0);
			MinuitThX->Migrad();
			MinuitThX->Command("MINOS");

			MinuitThZ->SetPrintLevel(-1);
			MinuitThZ->SetFCN(LikelihoodFunctionPolaZ);
			MinuitThZ->DefineParameter(0,"P_z",0.1,0.01,0,0);
			MinuitThZ->Migrad();
			MinuitThZ->Command("MINOS");
		}


		void GetLambdaResults(double& P_ld,double& P_ld_err){
			double err_MINOS;
			MinuitLambda->GetParameter(0,P_ld,err_MINOS);
			P_ld_err = sqrt(Score2Sum_Lambda)/Hessian_Lambda;
			double n_eff = w_sum * w_sum / w2_sum;
			double rough_err = sqrt(3./n_eff);
			cout<<"P_ld = "<<P_ld<<" +- "<<P_ld_err<<" (MINOS error = "<<err_MINOS<<", rough error = "<<rough_err<<")"<<endl;
		}
		void GetResults(double& P_xi,double& P_xi_err){
			double err_MINOS;
			Minuit->GetParameter(0,P_xi,err_MINOS);
			P_xi_err = sqrt(Score2Sum_Xi)/Hessian_Xi;
			cout<<"P_xi = "<<P_xi<<" +- "<<P_xi_err<<" (MINOS error = "<<err_MINOS<<")"<<endl;
		}
		void GetResultsThX(double& P_x,double& P_x_err){
			double err_MINOS;
			MinuitThX->GetParameter(0,P_x,err_MINOS);
			P_x_err = sqrt(Score2Sum_ThX)/Hessian_ThX;
			cout<<"P_x = "<<P_x<<" +- "<<P_x_err<<" (MINOS error = "<<err_MINOS<<")"<<endl;
		}
		void GetResultsThZ(double& P_z,double& P_z_err){
			double err_MINOS;
			MinuitThZ->GetParameter(0,P_z,err_MINOS);
			P_z_err = sqrt(Score2Sum_ThZ)/Hessian_ThZ;
			cout<<"P_z = "<<P_z<<" +- "<<P_z_err<<" (MINOS error = "<<err_MINOS<<")"<<endl;
		}
		void GetResultsChi2(double& P_xi,double& P_xi_err){
			MinuitChi2->GetParameter(0,P_xi,P_xi_err);
		}
		void GetResultsChi2Lambda(double& P_ld,double& P_ld_err){
			MinuitLambdaChi2->GetParameter(0,P_ld,P_ld_err);
		}
		double GetDensity(){
			return Fitdensity;
		}
};
class PolaAnalyzer{
	private:
		bool debug = 0;
		bool CheckXZ = 0;
		double weight_cap = 20;
		vector<double> Th;
		vector<double> PhA;//Z alpha
		vector<double> PhB;//X beta
		vector<double> PhC;//Y gamma
		vector<double> ThX;
		vector<double> ThZ;
		vector<double> Weight;
		vector<double> Scale;
		map<TString, vector<double>> par_vects;
		vector<TString> pars = {
			"MM", "CThCM", "CThCMKpXi","SqrtS", "PT", "ET", "PXi", "EXi",
			"Coplanarity", "Colinearity", "MMCKm_KpXi", "MMpKm_KpXi","CosKpXi","dPpKmKpXi",
			"MMCor","BECor","MMCKm_KpXiCor","MMpKm_KpXiCor","SqrtSCor","dPpKmKpXiCor",
			"WeightKurama","EffAna","TargetScale","CrossSection",
			"ThX","ThZ","FiducialFlag","KFXiPval","ResidualTracks"
		};
		TString title_base = "";

	public:
		PolaAnalyzer(){
			Clear();
			for(auto p:pars){
				par_vects[p] = vector<double>();
			}
		}
		void Clear();
		~PolaAnalyzer(){
			Clear();
			par_vects.clear();
		}
		void FillEvent(double Th_, double PhA_, double PhB_, double PhC_,
			vector<double> par_values,double weight_);
		PolaAnalyzer* GetSelectedEvents(TString par_name, double p_min, double p_max);
		void SetTitle(TString t){title_base=t;};
		TString GetTitle(){return title_base;};
		void XZPolaFlag(bool flag = 1){ CheckXZ = flag;}
		void MakeFile(TString fname);
		void LoadFile(TString fname);
		void SetWeightCap(double cap){weight_cap = cap;}
		bool CheckParameter(TString par_name);
		bool CheckIndex(int index);
		int GetEntries(){return Th.size();}
		double GetWeight(int index){
			if(!CheckIndex(index)) return -9999;
			else return Weight[index];
		}
		double GetTheta(int index){
			if(!CheckIndex(index)) return -9999;
			else return Th[index];
		}
		double GetPhiA(int index){
			if(!CheckIndex(index))return -9999;
			else return PhA[index];
		}
		double GetPhiB(int index){
			if(!CheckIndex(index))return -9999;
			else return PhB[index];
		}
		double GetPhiC(int index){
			if(!CheckIndex(index))return -9999;
			else return PhC[index];
		}
		double GetValue(TString par_name,int index){
			if(!CheckParameter(par_name) or !CheckIndex(index)) return -9999;
			else return par_vects[par_name][index];
		}
		double GetCumulatedWeight();
		double GetCumulatedCrossSection();
		double GetCumulatedTargetScale();
	private:

		double p_xi = -9999, p_xi_err = -9999, p_xi_sys;
		double p_ld = -9999, p_ld_err = -9999, p_ld_sys;
		double p_xi_x = -9999, p_xi_x_err = -9999;
		double p_xi_z = -9999, p_xi_z_err = -9999;
		TH1D* hTh = nullptr;
		TH1D* hPhA = nullptr;
		TH1D* hPhB = nullptr;
		TH1D* hPhC = nullptr;

		TH1D* hThEff = nullptr;
		TH1D* hPhAEff = nullptr;
		TH1D* hPhBEff = nullptr;
		TH1D* hPhCEff = nullptr;
		
		TGraphErrors* gTh = nullptr;
		TGraphErrors* gPhA = nullptr;
		TGraphErrors* gPhB = nullptr;
		TGraphErrors* gPhC = nullptr;
		
		TH1D* hThX = nullptr;
		TH1D* hThZ = nullptr;
		TH1D* hThXEff = nullptr;
		TH1D* hThZEff = nullptr;



	public:
		void Analyze(int nbin = 10){
			PolaFitter PolF(Th, PhA, PhB, PhC);
			vector<double> Weight_vector;
			for(int i=0;i<Th.size();++i){
				double weight = Weight[i];
				double eff_ana = par_vects["EffAna"][i];
				weight *= 1./eff_ana;
#if KuramaEff
				double kurama_w = par_vects["WeightKurama"][i];
				weight *= kurama_w;
#endif
				Weight_vector.push_back(weight);
			}
			PolF.SetWeight(Weight_vector, Weight_vector, Weight_vector, Weight_vector);
			PolF.SetWeightCap(weight_cap);
			if(CheckXZ){
				PolF.AddThXZ(ThX, ThZ);
				PolF.DoFitXZ();
				PolF.GetResultsThX(p_xi_x,p_xi_x_err);
				PolF.GetResultsThZ(p_xi_z,p_xi_z_err);
				//cout<<Form("P_X = %g +- %g",p_xi_x,p_xi_x_err)<<endl;
				//cout<<Form("P_Z = %g +- %g",p_xi_z,p_xi_z_err)<<endl;
			}
			PolF.DoFit();
			PolF.GetResults(p_xi, p_xi_err);
			PolF.DoLambdaFit();
			PolF.GetLambdaResults(p_ld, p_ld_err);
			PolF.Clear();
			//cout<<Form("P_Xi = %g +- %g",p_xi,p_xi_err)<<endl;
			//cout<<Form("P_Ld = %g +- %g",p_ld,p_ld_err)<<endl;
#if CalcSystematics
			PolaSystematics::CalcPolaSys(p_ld, p_xi, p_ld_sys, p_xi_sys);
#endif
			hTh = new TH1D("hTh" + title_base,";cos#theta;Corrected Yield",nbin,-1,1);
			hPhA = new TH1D("hPhA" + title_base,";cos#phi_{#alpha};Corrected Yield",nbin,-1,1);
			hPhB = new TH1D("hPhB" + title_base,";cos#phi_{#beta};Corrected Yield",nbin,-1,1);
			hPhC = new TH1D("hPhC" + title_base,";cos#phi_{#gamma};Corrected Yield",nbin,-1,1);

			hThEff = new TH1D("hThEff" + title_base,";cos#theta;Efficiency",nbin,-1,1);
			hPhAEff = new TH1D("hPhAEff" + title_base,";cos#phi_{#alpha};Efficiency",nbin,-1,1);
			hPhBEff = new TH1D("hPhBEff" + title_base,";cos#phi_{#beta};Efficiency",nbin,-1,1);
			hPhCEff = new TH1D("hPhCEff" + title_base,";cos#phi_{#gamma};Efficiency",nbin,-1,1);

			hThX = new TH1D("hThX" + title_base,";cos#theta_{X};Corrected Yield",nbin,-1,1);
			hThZ = new TH1D("hThZ" + title_base,";cos#theta_{Z};Corrected Yield",nbin,-1,1);
			hThXEff = new TH1D("hThXEff" + title_base,";cos#theta_{X};Efficiency",nbin,-1,1);
			hThZEff = new TH1D("hThZEff" + title_base,";cos#theta_{Z};Efficiency",nbin,-1,1);

			for(int i=0;i<Th.size();++i){
				if(Weight[i] > weight_cap){
					continue;
				}
				hTh->Fill(cos(Th[i]),Weight_vector[i]);
				hPhA->Fill(cos(PhA[i]),Weight_vector[i]);
				hPhB->Fill(cos(PhB[i]),Weight_vector[i]);
				hPhC->Fill(cos(PhC[i]),Weight_vector[i]);
				hThEff->Fill(cos(Th[i]));
				hPhAEff->Fill(cos(PhA[i]));
				hPhBEff->Fill(cos(PhB[i]));
				hPhCEff->Fill(cos(PhC[i]));
				if(CheckXZ){
					hThX->Fill(cos(ThX[i]),Weight_vector[i]);
					hThZ->Fill(cos(ThZ[i]),Weight_vector[i]);
					hThXEff->Fill(cos(ThX[i]));
					hThZEff->Fill(cos(ThZ[i]));
				}
			}
			for(int i=1;i<=nbin;++i){
				double contTh = hTh->GetBinContent(i);
				double contPhA = hPhA->GetBinContent(i);
				double contPhB = hPhB->GetBinContent(i);
				double contPhC = hPhC->GetBinContent(i);
				double effTh = hThEff->GetBinContent(i);
				double effPhA = hPhAEff->GetBinContent(i);
				double effPhB = hPhBEff->GetBinContent(i);
				double effPhC = hPhCEff->GetBinContent(i);
				if(effTh>0)hThEff->SetBinContent(i,effTh/contTh);
				if(effPhA>0)hPhAEff->SetBinContent(i,effPhA/contPhA);
				if(effPhB>0)hPhBEff->SetBinContent(i,effPhB/contPhB);
				if(effPhC>0)hPhCEff->SetBinContent(i,effPhC/contPhC);
				if(CheckXZ){
					double contThX = hThX->GetBinContent(i);
					double contThZ = hThZ->GetBinContent(i);
					double effThX = hThXEff->GetBinContent(i);
					double effThZ = hThZEff->GetBinContent(i);
					if(effThX>0)hThXEff->SetBinContent(i,effThX/contThX);
					if(effThZ>0)hThZEff->SetBinContent(i,effThZ/contThZ);
				}
			}
		}
		TH1D* GetThHist(){return hTh;};
		TH1D* GetPhAHist(){return hPhA;};
		TH1D* GetPhBHist(){return hPhB;};
		TH1D* GetPhCHist(){return hPhC;};
		double GetPXi(){
			if(!isnan(p_xi_err) and p_xi_err >0){
				return p_xi;
			}
			else return NAN;
		};
		double GetPXiErr(){return p_xi_err;};
		double GetPLd(){
			if(!isnan(p_ld_err) and p_ld_err >0){
				return p_ld;
			}
			else return NAN;
		};
		double GetPLdErr(){return p_ld_err;};
		double GetPXiSys(){return p_xi_sys;};
		double GetPLdSys(){return p_ld_sys;};
		double GetPXi_X(){
			if(!isnan(p_xi_x_err) and p_xi_x_err >0){
				return p_xi_x;
			}
			else return NAN;
		}
		double GetPXi_XErr(){return p_xi_x_err;};
		double GetPXi_Z(){
			if(!isnan(p_xi_z_err) and p_xi_z_err >0){
				return p_xi_z;
			}
			else return NAN;
		}
		double GetPXi_ZErr(){return p_xi_z_err;};
		TH1D* GetHist(int i){
			switch (i){
			case 0:
				return hTh;
			case 1:
				return hPhA;
			case 2:
				return hPhB;
			case 3:
				return hPhC;
			case 4:
				return hThX;
			case 5:
				return hThZ;
			default:
				return hTh;
			}
		}
		TH1D* GetEffHist(int i){
			switch (i){
			case 0:
				return hThEff;
			case 1:
				return hPhAEff;
			case 2:
				return hPhBEff;
			case 3:
				return hPhCEff;
			default:
				return hThEff;
			}
		}
		double GetSlope(int i){
			switch (i){
			case 0:
				return p_xi*AXi;
			case 1:
				return p_ld;
			case 2:
				return p_xi*KB;
			case 3:
				return p_xi*KC;
			case 4:
				return p_xi_x*AXi;
			case 5:
				return p_xi_z*AXi;
			default:
				return -9999;
			}
		}
		double GetDensity(){
			double conts = hTh->Integral();
			return conts / hTh->GetNbinsX();
		}

};
void PolaAnalyzer::Clear(){
	Th.clear();
	PhA.clear();
	PhB.clear();
	PhC.clear();
	ThX.clear();
	ThZ.clear();
	Weight.clear();
	
	par_vects.clear();
	for(auto p:pars){
		par_vects[p] = vector<double>();
	}
	CheckXZ = 0;
};
void PolaAnalyzer::FillEvent(double Th_, double PhA_, double PhB_, double PhC_,
	vector<double> par_values,double weight_){
	Th.push_back(Th_);
	PhA.push_back(PhA_);
	PhB.push_back(PhB_);
	PhC.push_back(PhC_);
	if(weight_ > weight_cap){
		weight_ = 0;
	}
	Weight.push_back(weight_);
	for(int i=0;i<pars.size();++i){
		//cout<<Form("Filling %s = %g",pars[i].Data(),par_values[i])<<endl;
		par_vects[pars[i]].push_back(par_values[i]);
		if(pars[i] == "ThX"){
			ThX.push_back(par_values[i]);
		}
		if(pars[i] == "ThZ"){
			ThZ.push_back(par_values[i]);	
		}
	}
};
PolaAnalyzer* PolaAnalyzer::GetSelectedEvents(TString par_name, double p_min, double p_max){
	PolaAnalyzer* Cont = new PolaAnalyzer();
	Cont->SetWeightCap(weight_cap);
	if(par_name == ""){
		Cont = this;
		return Cont;
	}
	if(par_vects.find(par_name)==par_vects.end()){
		cout<<"Parameter "<<par_name<<" not found!"<<endl;
		cout<<"ParLists :"<<endl;
		for(auto p:par_vects){
			cout<<p.first<<endl;
		}
		return nullptr;
	}
	int idata = 0;
	for(int i=0;i<Th.size();++i){
		double p_val = par_vects[par_name][i];
		vector<double> par_values;
		for(auto p:pars){
			par_values.push_back(par_vects[p][i]);
		}
		if(p_val >= p_min and p_val < p_max){
			Cont->FillEvent(Th[i],PhA[i],PhB[i],PhC[i],
				par_values,
				Weight[i]);
			++idata;
		}
	}
	Cont->XZPolaFlag(CheckXZ);
	Cont->SetTitle(title_base+Form("_%s_%.3g_%.3g",par_name.Data(),p_min,p_max));
	//cout<<Form("Selected %d events for %s in [%.2g,%.2g)",idata,par_name.Data(),p_min,p_max)<<endl;
	return Cont;
}
void PolaAnalyzer::MakeFile(TString fname){
	TFile* file = new TFile(fname,"RECREATE");
	TNamed* titlen = new TNamed("Title",title_base.Data());
	titlen->Write();
	TTree* tree = new TTree("tree","tree");
	double t_Th, t_PhA, t_PhB, t_PhC, t_Weight;
	double params[100] = {0};
	tree->Branch("Th",&t_Th);
	tree->Branch("PhA",&t_PhA);
	tree->Branch("PhB",&t_PhB);
	tree->Branch("PhC",&t_PhC);
	tree->Branch("Weight",&t_Weight);
	for(int i=0;i<pars.size();++i){
		tree->Branch(pars[i].Data(),&params[i]);
	}
	for(int i=0;i<Th.size();++i){
		t_Th = Th[i];
		t_PhA = PhA[i];
		t_PhB = PhB[i];
		t_PhC = PhC[i];
		t_Weight = Weight[i];
		for(int j=0;j<pars.size();++j){
			params[j] = par_vects[pars[j]][i];
		}
		tree->Fill();
	}
	tree->Write();
	file->Close();
	cout<<Form("File %s created with %d entries.",fname.Data(),(int)Th.size())<<endl;
}
void PolaAnalyzer::LoadFile(TString fname){
	Clear();
	TFile* file = new TFile(fname,"READ");
	TTree* tree = (TTree*)file->Get("tree");
	TNamed* titlen = (TNamed*)file->Get("Title");
	title_base = titlen->GetTitle();
	double t_Th, t_PhA, t_PhB, t_PhC, t_Weight;
	double params[100] = {0};
	tree->SetBranchAddress("Th",&t_Th);
	tree->SetBranchAddress("PhA",&t_PhA);
	tree->SetBranchAddress("PhB",&t_PhB);
	tree->SetBranchAddress("PhC",&t_PhC);
	tree->SetBranchAddress("Weight",&t_Weight);
	for(int i=0;i<pars.size();++i){
		tree->SetBranchAddress(pars[i].Data(),&params[i]);
	}
	int nentries = tree->GetEntries();
	for(int i=0;i<nentries;++i){
		tree->GetEntry(i);
		if(i%10000==0){
			cout<<"Loading entry "<<i<<"/"<<nentries<<endl;
		}
		vector<double> par_values;
		for(int j=0;j<pars.size();++j){
			par_values.push_back(params[j]);
		}
		FillEvent(t_Th,t_PhA,t_PhB,t_PhC,
			par_values,t_Weight);
	}
	file->Close();
	cout<<Form("File %s loaded with %d entries.",fname.Data(),(int)Th.size())<<endl;
}
bool PolaAnalyzer::CheckParameter(TString par_name){
	if(par_vects.find(par_name)==par_vects.end()){
		cout<<"Parameter "<<par_name<<" not found!"<<endl;
	}
	return par_vects.find(par_name)!=par_vects.end();
}
bool PolaAnalyzer::CheckIndex(int index){
	if(index < 0 or index >= Th.size()){
		cout<<"Index "<<index<<" out of range!"<<endl;
	}
	return index >= 0 and index < Th.size();
}
double PolaAnalyzer::GetCumulatedWeight(){
	double sum = 0;
	for(auto w:Weight){
		sum+=w;
	}
	return sum;
}
double PolaAnalyzer::GetCumulatedCrossSection(){
	double sum = 0;
	for(int i=0;i<Weight.size();++i){
		sum+=Weight[i]*par_vects["CrossSection"][i];//Weight is Xi recon weight, which is not included in KK cross section.
	}
	return sum;
}
double PolaAnalyzer::GetCumulatedTargetScale(){
	double sum = 0;
	for(int i=0;i<Weight.size();++i){
		sum+=Weight[i]*par_vects["TargetScale"][i];
	}
	return sum;
}
class PolaGen{
	private:
		double P_Xi = 0.;//Xi Polarization
		double CMAngle = 0;
		TF1 fXiPol = TF1("fXiPol","1+[0]*x",-1,1);
		TF1 fLdPol = TF1("fLdPol","1+[0]*x",-1,1);
		TF1 fLdPolX = TF1("fLdPolX","1+[0]*x",-1,1);
		TVector3 TVKm;
		TVector3 TVKp;
		TVector3 TVXi;
		TVector3 TVLd;
		TVector3 TVP;
		TVector3 TVPi1;
		TVector3 TVPi2;
		TVector3 TVPolXi;
		TVector3 TVPolLd;
		
		TLorentzVector LVKm;
		TLorentzVector LVKp;
		TLorentzVector LVXi;
		TLorentzVector LVLd;
		TLorentzVector LVP;
		TLorentzVector LVPi1;
		TLorentzVector LVPi2;
		TVector3 Plane;
		
		double cTh,Th,Ph,Ph1,Ph2,Ph3;
		int reset_count;
	public:
		PolaGen(){
			Initialize(0,25,1.82);
		}
		PolaGen(double PolXi,double Angle,double PKm){
			Initialize(PolXi,Angle,PKm);
		}
		void Initialize(double PolXi,double Angle,double PKm){
			P_Xi = PolXi;
			fXiPol.SetParameter(0,P_Xi*AXi);
			
			TVKm = TVector3(0,0,PKm);
			CMAngle = Angle;
			TVector3 TVTarget(0,0,0);
			LVKm = TLorentzVector(TVKm,hypot(mk,TVKm.Mag()));
			TLorentzVector LVTarget(TVTarget,hypot(mp,TVTarget.Mag()));
			auto LVCM = LVKm + LVTarget;
			auto BetaCM = LVCM.BoostVector();
			auto MCM = LVCM.M();
			auto pKpCM = Pmx(MCM,mk,mXi);
			TVector3 TVKpCM = pKpCM*TVector3(sin(CMAngle*TMath::DegToRad()),0,cos(CMAngle*TMath::DegToRad()));	
			TVector3 TVXiCM = -TVKpCM;
			auto LVKpCM = TLorentzVector(TVKpCM,hypot(mk,TVKpCM.Mag()));
			auto LVXiCM = TLorentzVector(TVXiCM,hypot(mXi,TVXiCM.Mag()));
			LVKp = LVKpCM;
			LVXi = LVXiCM;
			LVKp.Boost(BetaCM);
			LVXi.Boost(BetaCM);
			TVXi = LVXi.Vect();
			Plane = LVKm.Vect().Cross(LVXi.Vect());
			Plane = Plane*(1./Plane.Mag());
		}
		void Generate(){
			double pLdCM = Pmx(mXi,mLd,mpi);
			double ThetaLdCM = acos(fXiPol.GetRandom(-1,1));
			double PhiLdCM = 2*TMath::Pi()*gRandom->Rndm();
			double pLdCM_x = pLdCM*sin(ThetaLdCM)*cos(PhiLdCM);
			double pLdCM_y = pLdCM*sin(ThetaLdCM)*sin(PhiLdCM);
			double pLdCM_z = pLdCM*cos(ThetaLdCM);
			TVector3 ZXiCM = Plane;//Z axis of XiCM frame
			TVector3 XXiCM = TVXi.Cross(ZXiCM);
			XXiCM = XXiCM*(1./XXiCM.Mag());
			TVector3 YXiCM = ZXiCM.Cross(XXiCM);
			TVector3 TVLdCM = pLdCM_x*XXiCM + pLdCM_y*YXiCM + pLdCM_z*ZXiCM;
			LVLd = TLorentzVector(TVLdCM,hypot(mLd,TVLdCM.Mag()));
			TVector3 TVPi2CM = -TVLd;
			LVPi2 = TLorentzVector(TVPi2CM,hypot(mpi,TVPi2CM.Mag()));
			auto BetaXi = LVXi.BoostVector();
			LVLd.Boost(BetaXi);
			LVPi2.Boost(BetaXi);
			TVLd = LVLd.Vect();
			TVPi2 = LVPi2.Vect();
			TVector3 TVLdCMDir = TVLdCM*(1./TVLdCM.Mag());
			TVector3 PolarityLdZ = TVLdCMDir;
			TVector3 PolarityLdX =
			Plane.Cross(TVLdCMDir);// Magnitute is sin(ThetaLdCM)
			TVector3 PolarityLdY = TVLdCMDir.Cross(PolarityLdX);
			TVector3 PolarityLd =
			+P_Xi * CXi *PolarityLdX
			+P_Xi * BXi *PolarityLdY
			+(AXi + P_Xi*cos(ThetaLdCM))*PolarityLdZ;
			PolarityLd = PolarityLd * (1./(1 +P_Xi* AXi*cos(ThetaLdCM)));

			double P_Ld = PolarityLd.Mag();
			double pPCM = Pmx(mLd,mp,mpi);
			
			auto PLd_dir = PolarityLd.Unit();
			double rcth = gRandom->Uniform(-1,1);
			double rsth = sqrt(1-rcth*rcth);
			double rph = gRandom->Uniform(-acos(-1),acos(-1));
			TVector3 r_dir(rsth*cos(rph),rsth*sin(rph),rcth);
#if 1
//			fLdPol.SetParameter(0,(PolarityLd*TVLdCMDir) *ALd);
			fLdPol.SetParameter(0,P_Ld *ALd);//Proton momentum spectrum allignment with Lambda polarity should be P_Ld*ALd, Just like L momentum alligns with Xi polarity is P_Xi*AXi. 
			double ThetaPCM = acos(fLdPol.GetRandom(-1,1));
			double PhiPCM = 2*TMath::Pi()*gRandom->Rndm();//
			bool retry = 1;
			double pPCM_x = pPCM*sin(ThetaPCM)*cos(PhiPCM);
			double pPCM_y = pPCM*sin(ThetaPCM)*sin(PhiPCM);
			double pPCM_z = pPCM*cos(ThetaPCM);
#else
			double pPCM_x,pPCM_y,pPCM_z;
			double r = gRandom->Uniform(0,1);
			if(r <1+ P_Ld*ALd){
				pPCM_x = PLd_dir.x();
				pPCM_y = PLd_dir.y();
				pPCM_z = PLd_dir.z();
			}
			else{
				pPCM_x = r_dir.x();
				pPCM_y = r_dir.y();
				pPCM_z = r_dir.z();
			}
#endif
//			fLdPolX.SetParameter(0,P_Xi *CXi);
#if 1
			TVector3 ZLdSpin = PLd_dir;
			TVector3 XLdSpin = Plane.Cross(ZLdSpin);//if P_Xi == 0, TVLdCM is alligned with ZLdCM
#else
			TVector3 ZLdCM = PolarityLd;
			TVector3 XLdCM = ZLdCM.Cross(TVLdCMDir);//if P_Xi == 0, TVLdCM is alligned with ZLdCM
#endif
			int cnt = 999;
			XLdSpin = XLdSpin.Unit();
			TVector3 YLdSpin = ZLdSpin.Cross(XLdSpin);
			TVector3 TVPCM = pPCM_x*XLdSpin + pPCM_y*YLdSpin + pPCM_z*ZLdSpin;
			LVP = TLorentzVector(TVPCM,hypot(mP,TVPCM.Mag()));
			TVector3 TVPi1CM = -TVPCM;
			LVPi1 = TLorentzVector(TVPi1CM,hypot(mPi,TVPi1CM.Mag()));
			auto BetaLd = LVLd.BoostVector();
			LVP.Boost(BetaLd);
			LVPi1.Boost(BetaLd);
			TVP = LVP.Vect();
			TVPi1 = LVPi1.Vect();
			TVector3 LdZCM = TVLdCMDir;
			TVector3 LdXCM = (TVLdCMDir.Cross(Plane)).Unit();
			TVector3 LdYCM = LdZCM.Cross(LdXCM);
			Ph3 = acos(NormalDot(LdZCM,TVPCM));
			Ph1 = acos(NormalDot(LdXCM,TVPCM));
			Ph2 = acos(NormalDot(LdYCM,TVPCM));
			if(isnan(TVP.Mag())){
				cout<<Form("Plane: %f %f %f",Plane.X(),Plane.Y(),Plane.Z())<<endl;
				cout<<Form("TVLd: %f %f %f",TVLd.X(),TVLd.Y(),TVLd.Z())<<endl;
//				cout<<Form("PolarityLdZ: %g,%g,%g ",PolarityLdZ.X(),PolarityLdZ.Y(),PolarityLdZ.Z())<<endl;
//				cout<<Form("PolarityLdX: %g,%g,%g ",PolarityLdX.X(),PolarityLdX.Y(),PolarityLdX.Z())<<endl;
//				cout<<Form("PolarityLdY: %g,%g,%g ",PolarityLdY.X(),PolarityLdY.Y(),PolarityLdY.Z())<<endl;
				cout<<Form("PolarityLd: %f %f %f",PolarityLd.X(),PolarityLd.Y(),PolarityLd.Z())<<endl;
				cout<<Form("TVLd*PolarityLd: %f",TVLdCM*PolarityLd*(1./TVLdCM.Mag()))<<endl;
				cout<<Form("TVLdCM: %f %f %f",TVLdCM.X(),TVLdCM.Y(),TVLdCM.Z())<<endl;
				cout<<Form("TVPCM: %f %f %f",TVPCM.X(),TVPCM.Y(),TVPCM.Z())<<endl;	
			}
		}
		TLorentzVector GetKm(){
			return LVKm;
		}
		TLorentzVector GetKp(){
			return LVKp;
		}
		TLorentzVector GetXi(){
			return LVXi;
		}
		TLorentzVector GetLd(){
			return LVLd;
		}
		TLorentzVector GetP(){
			return LVP;
		}
		TLorentzVector GetPi1(){
			return LVPi1;
		}
		TLorentzVector GetPi2(){
			return LVPi2;
		}
		TVector3 GetPolXi(){
			return TVPolXi;
		}
		TVector3 GetPolLd(){
			return TVPolLd;
		}
		int GetResetCount(){
			return reset_count;
		}
};
#endif
