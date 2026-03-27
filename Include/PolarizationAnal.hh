#include "Math.hh"
#include "Constants.hh"
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

class PolaAnal{
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
		TVector3 PolXi;
		TVector3 PolLd;
		double cTh,Th,Ph,PhA,PhB,PhC;
		TVector3 XLd;
		TVector3 YLd;
		TVector3 ZLd;

	public:
		PolaAnal(TVector3 PKm,TVector3 PKp,TVector3 PXi,TVector3 PLd,TVector3 PP, TVector3 PPi1, TVector3 PPi2,bool PiMode = 0, bool debug = 0){
			Km=PKm,Kp=PKp,Xi=PXi,Ld=PLd,P=PP,Pi1=PPi1,Pi2=PPi2;	
			PolXi= NormalCross(Km,Kp);
//			TVector3 Dir0 = NormalCross(Xi,Km);
//			PolXi= NormalCross(Xi,Kp);
//			if(Dir0*PolXi<0)PolXi = -PolXi;
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
			if((isnan(cTh)
			 or isnan(PhB)
			 or isnan(PhC)
			 or isnan(PhA)) and debug
			){
				cout<<"Something Wrong! "<<endl;
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
		TVector3 GetPolAxis(){
			return PolXi;
		}
		TVector3 GetPolLambda(){
			return PolLd;
		}

};
#if 0
TF1* fitTheta = new TF1("fitTheta","[0] + [1]*x",-1,1);
TF1* fitPhi1 = new TF1("fitPhi1","[0] + [1]*x",-1,1);
TF1* fitPhi2 = new TF1("fitPhi2","[0] + [1]*x",-1,1);
TF1* UbfitTheta = new TF1("UbfitTheta","[0] + [1]*x",-1,1);
TF1* UbfitPhi1 = new TF1("UbfitPhi1","[0] + [1]*x",-1,1);
TF1* UbfitPhi2 = new TF1("UbfitPhi2","[0] + [1]*x",-1,1);
#endif
namespace
{
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
	static void
	Chi2Function(	int& npar, double* gin, double& f, double* par, int iflag){
		double chi2 = 0;
		double P_xi = par[0];
		int nb = FitXTh.size();
		double bw = 2. / nb;
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
			chi2 += pow((Fitdensity +P_xi*Fitdensity* AXi*CorXTh-FitYTh[i])/FitETh[i],2);
#if IncludeBeta
			chi2 += pow((Fitdensity +P_xi*Fitdensity* KB*FitXPhB[i]-FitYPhB[i])/FitEPhB[i],2);
#endif
			chi2 += pow((Fitdensity +P_xi*Fitdensity* KC*FitXPhC[i]-FitYPhC[i])/FitEPhC[i],2);
		}
		f = chi2;
	}
	static void
	Chi2FunctionLambda(	int& npar, double* gin, double& f, double* par, int iflag){
		double chi2 = 0;
		double P_l = par[0];
		int nb = FitXPhA.size();
		double bw = 2. / nb;
		for(int i=0;i<FitYPhA.size();++i){
			if(FitEPhA[i]==0) FitEPhA[i]=1.;
			double CorXPhA = FitXPhA[i];
			double PhA0 = FitXPhA[i] - bw/2;
			double PhA1 = FitXPhA[i] + bw/2;
#if SlopeCorrection
			CorXPhA = (FitXPhA[i] + 1./3 *P_l * (PhA0*PhA0 + PhA0*PhA1 + PhA1*PhA1)) / (1 + 1./2* P_l * (FitXPhA[i]));
#endif
			chi2 += pow((Fitdensity + P_l*Fitdensity*CorXPhA-FitYPhA[i])/FitEPhA[i],2);
		}
		f = chi2;
	}
	static void
	LikelihoodFunction(	int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_xi = par[0];
		double norm = 0;
		double wnorm = 0;
		int np = FitTh.size();
		for(int i=0;i<FitTh.size();++i){
			double Th_i = FitTh[i];
			double PhA_i = FitPhA[i];
			double PhB_i = FitPhB[i];
			double PhC_i = FitPhC[i];
			double wTh_i = FitThWeight[i];
			double wPhA_i = FitPhAWeight[i];
			double wPhB_i = FitPhBWeight[i];
			double wPhC_i = FitPhCWeight[i];
			double probTh = 0.5*(1 + P_xi*AXi*cos( Th_i));
			double probPhA = 0.5*(1 + KA*cos(PhA_i));
			double probPhB = 0.5*(1 + P_xi*KB*cos(PhB_i));
			double probPhC = 0.5*(1 + P_xi*KC*cos(PhC_i));

			likelihood -= wTh_i*log(probTh);
#if IncludeBeta
			likelihood -= wPhB_i*log(probPhB);
#endif
			likelihood -= wPhC_i*log(probPhC);

			norm += wTh_i;
			wnorm += wTh_i*wTh_i;
	/*
			probTh*=1./wTh_i;
			probPh1*=1./wPh1_i;
			probPh2*=1./wPh2_i;

			likelihood -= log(probTh);
			likelihood -= log(probPh1);
			likelihood -= log(probPh2);
			*/
		//	if(PhAFlag)likelihood -= wPh3_i*log(probPh3);
		}
		double n_eff = sqrt(np /wnorm);
		likelihood *= n_eff;
		f=likelihood;
	}
	static void
	LikelihoodFunctionLambda(	int& npar, double* gin, double& f, double* par, int iflag){
		double likelihood = 0;
		double P_ld = par[0];
		int np = FitPhA.size();
		double norm = 0, wnorm = 0;
		for(int i=0;i<FitPhA.size();++i){
			double PhA_i = FitPhA[i];
			double wPhA_i = FitPhAWeight[i];
			double probPhA = 0.5*(1 + P_ld*cos(PhA_i));
			likelihood -= wPhA_i*log(probPhA);
			norm += wPhA_i;
			wnorm += wPhA_i*wPhA_i;
		}
//		double n_eff = np * 1./norm;
		double n_eff = sqrt(np /wnorm);
		likelihood *= n_eff;
		f=likelihood;
	}
}
class PolaFitter{
	private:
		TMinuit* Minuit = new TMinuit(1);
		TMinuit* MinuitChi2 = new TMinuit(1);
		TMinuit* MinuitLambda = new TMinuit(1);
		TMinuit* MinuitLambdaChi2 = new TMinuit(1);
	public:
		PolaFitter(){}
		PolaFitter(vector<double> Th,vector<double> PhA,vector<double> PhB,vector<double> PhC){
			FitTh = Th;
			FitPhA = PhA;
			FitPhB = PhB;
			FitPhC = PhC;
			FitPhAFlag = false;
		}
		void SetWeight(vector<double> ThWeight,vector<double> PhAWeight,vector<double> PhBWeight,vector<double> PhCWeight){
			FitThWeight = ThWeight;
			FitPhAWeight = PhAWeight;
			FitPhBWeight = PhBWeight;
			FitPhCWeight = PhCWeight;
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
			FitThWeight.clear();
			FitPhBWeight.clear();
			FitPhCWeight.clear();
			FitPhAWeight.clear();
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
			cout<<"Chi2FitDone"<<endl;
		}
		void DoLambdaFit(){
			MinuitLambda->SetPrintLevel(-1);
			MinuitLambda->SetFCN(LikelihoodFunctionLambda);
			MinuitLambda->DefineParameter(0,"P_ld",-0.,0.001,0,0);
			MinuitLambda->Migrad();
			MinuitLambda->Command("MINOS");
		}
		void DoChi2FitLambda(){
//			cout<<"DoingChi2FitLambda"<<endl;
			MinuitLambdaChi2->SetPrintLevel(-1);
			MinuitLambdaChi2->SetFCN(Chi2FunctionLambda);
			MinuitLambdaChi2->DefineParameter(0,"P_ld",-0.,0.001,0,0);
			MinuitLambdaChi2->Migrad();
			cout<<"LambdaChi2FitDone"<<endl;
		}
		void GetLambdaResults(double& P_ld,double& P_ld_err){
			MinuitLambda->GetParameter(0,P_ld,P_ld_err);
		}
		void GetResults(double& P_xi,double& P_xi_err){
			Minuit->GetParameter(0,P_xi,P_xi_err);
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



#if 0

double SimFitPola( double dens,vector<double> pols, vector<vector<double>> ents, vector<vector<double>> errs,double& error){
	int nbin = pols.size();
	for(int j=0;j<nbin;++j){
		for(int i=0;i<3;++i){
			ents[i][j]*=1./dens;	
			errs[i][j]*=1./dens;	
		}
	}
	auto ThEnt = ents[0];
	auto ThErr = errs[0];
	auto PhBEnt = ents[1];
	auto PhBErr = errs[1];
	auto PhCEnt = ents[2];
	auto PhCErr = errs[2];
	double SigKX=0;
	double SigKY=0;
	double ErrK=0;
	//f_i = 1 + P*K_j*X_j, K_j = DecayParam,X_j= DecayAngle
	//Chi2 = Sum[(y_i - f_i )^2/sig_i^2]
	//dChi2/dP = 0 -> Sum[   (y_i-f_i)*d(f_i)/d(dP)/sig_i^2 = 0;
	//df_i/dP = K_jX_j
	//dChi2/dP = Sum[(y_i-1-PK_jX_j)*K_jX_j/sig_i^2
	//P = Sum[y_i-1]*K_jX_j/sig_i^2 / Sum[ K_jX_jK_jX_j]/sig_i^2
	//

	for(int i=0;i<nbin;++i){
		SigKX += AXi*AXi*pols[i]*pols[i]/ThErr[i]/ThErr[i];
		SigKX += KB*KB*pols[i]*pols[i]/PhBErr[i]/PhBErr[i];
		SigKX += KC*KC*pols[i]*pols[i]/PhCErr[i]/PhCErr[i];
		
		SigKY += AXi*(ThEnt[i]-1)*pols[i]/ThErr[i]; 
		SigKY += KB*(PhBEnt[i]-1)*pols[i]/PhBErr[i]; 
		SigKY += KC*(PhCEnt[i]-1)*pols[i]/PhCErr[i]; 
		
		double dPdYA = AXi*pols[i]/ThErr[i]/ThErr[i];
		double dPdYB = KB*pols[i]/PhBErr[i]/PhBErr[i];
		double dPdYC = KC*pols[i]/PhCErr[i]/PhCErr[i];
		ErrK += dPdYA*dPdYA*ThErr[i]*ThErr[i];
		ErrK += dPdYB*dPdYB*PhBErr[i]*PhBErr[i];
		ErrK += dPdYC*dPdYC*PhCErr[i]*PhCErr[i];
	}
	double P = SigKY/SigKX;
	for(int i=0;i<nbin;++i){
		double f_a = 1 + AXi*P*pols[i];
		double f_b = 1 + AXi*P*pols[i];
		double f_c = 1 + AXi*P*pols[i];
		double e_a = (ThEnt[i]-f_a)/ThErr[i];
		double e_b = (PhBEnt[i]-f_b)/PhBErr[i];
		double e_c = (PhCEnt[i]-f_c)/PhCErr[i];
		error+= e_a*e_a;
		error+= e_b*e_b;
		error+= e_c*e_c;
	}

	return P;
}
#endif
