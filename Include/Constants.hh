#ifndef Constants_hh
#define Constants_hh
#include <map>
#include <vector>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <stdio.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TStyle.h>
#include <TEfficiency.h>
#ifndef SavePDF
#define SavePDF 0
#endif
#define Skip1D 1
//#include "group.hh"
int phi_conf = 0;


TString figdir = "figs/";
double PI = acos(-1);
double mLd = 1115.683/1000;
double mXi = 1.32171;
double mXi1530 = 1.535;
double mp = 0.938272;
double mn = 0.939565;
double md = 1.875612;
double mpi = 0.13957018;
double mk = 0.493677;
double mks = 0.497611;
double mphi = 1.019461;
double mPi0 = 0.1349768;
double m_eta = 0.547862;
double m_etaprime = 0.95778;
double gphi = 4.249*0.001;
double ctau_L = 78.9;//Lambda ctau[mm]
double ctau_Xi = 49.1;//Xi ctau[mm]
double ctau_Xi0 = 87.1;//Xi0 ctau[mm]
double amu = 0.931494; //MeV/c^2
double Kmass = 0.493677;
double ximass = 1.32172;
double lmass = 1.115683;
double pi0mass = 0.1349768;
TVector3 beam_dir(-0.063,-0.0016,1);
beam_dir=beam_dir.Unit();
TVector3 beam_avg = 1.82*beam_dir;

double phi_max = PI/2.;







const auto mH = 1.008 * amu;
const auto m12C = 11.174864;
const auto m11B = 10.252548;
const auto m11C = 10.252548;
const auto m11Be = 10.264568;
const auto mN = 14.0067*amu;
const auto mO = 15.9994*amu;
double N_avo = 6.022214076*TMath::Power(10, 23); //mol^-1
double ub_sr = TMath::Power(10, -30.); //ub/sr(= 10^-30 cm^-2)
double nb_sr = TMath::Power(10, -33.); //nb/sr(= 10^-33 cm^-2)
double CalcFermiMass(double p){
    TLorentzVector LV11B(0,0,-p, hypot(p, m11B));
    TLorentzVector LV12C(0,0,0, m12C);
    auto LVP = LV12C - LV11B;
    return LVP.M();
}

















double dif_cross_Burgun_th[] = {
 -0.975,-0.925,-0.875,-0.825,-0.775,
 -0.725,-0.675,-0.625,-0.575,-0.525,
 -0.475,-0.425,-0.375,-0.325,-0.275,
 -0.225,-0.175,-0.125,-0.075,-0.025
 ,0.025,0.075,0.125,0.175,0.225,
 0.275,0.325,0.375,0.425,0.475,
 0.525,0.575,0.625,0.675,0.725,
 0.775,0.825,0.875,0.925,0.975
};
double dif_cross_Burgun_thw[]={
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025,
	0.025,0.025,0.025,0.025,0.025
};
double dif_cross_Burgun_sig[] = {
36.97,47.5232,33.5323,37.0786,11.5782,
11.8711,16.3677,9.8798,8.59776,11.642,
13.0426,7.90643,4.97277,4.81638,2.9217,
3.39786,3.5702,2.87688,0,6.65527,
0.953314,0,0,8.52703,8.54695,
9.99635,13.437,9.03108,4.90802,12.0285,
13.7986,7.44123,6.94116,6.91426,3.03228,
4.94687,14.1762,5.25369,8.65653,12.5495
};
double dif_cross_Burgun_err[] = {
 7.53985,8.61369,6.86346,7.43923,3.90092,
 4.26053,5.70793,4.07823,2.87389,3.77241,
 4.48466,3.66583,2.97151,2.7673,2.4834,
 1.92655,1.77613,0,0,2.73144,
 0.953314,0,0,3.24744,3.37794,
 3.97463,4.30436,3.29426,2.60692,4.22068,
 4.23164,2.81611,2.99442,3.23449,2.71251,
 2.92469,4.46175,3.21258,3.4915,4.53646
};
double dif_cross_Dauber_th[] = {
	-0.95,-0.85,-0.75,-0.65,-0.55,
	-0.45,-0.35,-0.25,-0.15,-0.05,
	0.05,0.15,0.25,0.35,0.45,
	0.55,0.65,0.75,0.85,0.95
};
double dif_cross_Dauber_thw[]={
	0.05,0.05,0.05,0.05,0.05,
	0.05,0.05,0.05,0.05,0.05,
	0.05,0.05,0.05,0.05,0.05,
	0.05,0.05,0.05,0.05,0.05
};
double dif_cross_Dauber_sig[] = {
	65.0459,42.5087,31.5716,20.3895,19.2428,
	7.4778,12.49,2.62953,4.41017,1.54748,
	5.99044,5.87841,7.45412,10.2977,8.1518,
	4.33366,9.78125,6.6781,14.3463,5.91757
};
double dif_cross_Dauber_err[] = {
	8.04615,6.2318,5.84562,3.73617,4.0085,
	2.62133,3.68152,1.33252,2.56395,0,
	2.38634,2.49837,3.0822,2.9474,3.10861,
	1.96736,2.71879,2.44645,3.59681,2.32258
};
double dif_cross_E05_th[] = {
0.999087 , 0.991802 , 0.977336 , 0.955918 , 0.927851 , 0.893528 , 0.853443 
};
double dif_cross_E05_thw[] = {
0.001825 , 0.005448 , 0.008988 , 0.012393 , 0.015626 , 0.0186335 , 0.0213865 ,
};
double dif_cross_E05_sig[] = {
8.45833 , 10.0728 , 9.49404 , 9.75844 , 9.52568 , 8.81777 , 7.99634 
};
double dif_cross_E05_err[] = {
1.55715 , 1.61452 , 1.55357 , 1.61496 , 1.57771 , 1.56018 , 1.58422  
};
double dif_cross_E05_labth[] = {
    1,3,5,7,9,11,13
};
double dif_cross_E05_labthw[] = {
    1,1,1,1,1,1,1
};
double dif_cross_E05_labsig[] = {
    50.68,60.08,56.10,56.86,54.52,49.34,43.58
};
double dif_cross_E05_laberr[] = {
    9.33,9.63,9.19,9.41,9.03,8.73,8.60
};


//double xi_acceptance_in_ch2 = 0.992;
//double xi_acceptance_in_carbon = 0.987;

int run_mtx_mod = 5707;



double xi_acceptance_in_ch2 = 0.9958;//portion of w>50 events
double xi_acceptance_in_carbon = 0.9895;

double BrLd = 0.641;
double Kctau = 3.7114;
double L_target_ch2 = 2.35; //cm
double L_target_dia = 2.02; //cm
double L_target_wall = 0.0; //cm
//L_target_wall = 0.; //cm
double km_survival = 0.9577;
double err_km_survival = 0.0001;
double eff_bh12_ch2 = 0.9742;
double eff_bft_ch2 = 0.9783;
double eff_k18_ch2 = 0.9661;
double eff_daq_ch2 = 0.927687; //Avg.
double eff_sdcin_ch2 = 0.9987;
double eff_sdcout_ch2 = 0.9887;
double R_target_ch2 = 0.879566;

double eff_HTOFmp2 = 0.801228;
double err_eff_HTOFmp2 = 0.0316632;
double eff_TrigB_ch2 = 0.8942;
double err_eff_TrigB_ch2 = 0.010;
double eff_TrigB_prod = 0.8186;
double err_eff_TrigB_prod = 0.0039;

double ps_ch2 = 3.0;
double density_ch2 = 0.95;  //g*cm^-3
//density_ch2 = 1.13;
double density_dia = 3.223;
double density_targetwall = 1.42;
double W_ch2 = 12.0107+1.008*2;
double W_dia = 12.0107;
double W_targetwall = (22*m12C + 10 * mH + 2 * mN + 5 * mO) *1/(29.); //W_per_Nucleus other than H.

double nbeam_dia_TrigA = 18546891916; //BeamA-PS
double nbeam_dia_TrigB = 132486003858; //BeamB
double nbeam_ch2_TrigA = 1810021992; //BeamA-PS
double nbeam_ch2_TrigB = 5429939620; //BeamB

double trigged_A_dia = 623443932;
double trigged_B_dia = 187512377;
double trigged_A_ch2 = 31580956;
double trigged_B_ch2 = 4887339;


double nbeam_dia_TrigA_norB = nbeam_dia_TrigA*(1.-trigged_B_dia/trigged_A_dia);
double nbeam_ch2_TrigA_norB = nbeam_ch2_TrigA*(1.-trigged_B_ch2/trigged_A_ch2);

double trigA_factor_dia = W_dia/(N_avo*density_dia*L_target_dia*nbeam_dia_TrigA)/ub_sr;
double trigB_factor_dia = W_dia/(N_avo*density_dia*L_target_dia*nbeam_dia_TrigB)/ub_sr;
double trigA_norB_factor_dia = W_dia/(N_avo*density_dia*L_target_dia*nbeam_dia_TrigA_norB)/nb_sr;

double trigA_factor_ch2 = W_ch2/(N_avo*density_ch2*L_target_ch2*nbeam_ch2_TrigA)/ub_sr;
double trigB_factor_ch2 = W_ch2/(N_avo*density_ch2*L_target_ch2*nbeam_ch2_TrigB)/ub_sr;
double trigA_norB_factor_ch2 = W_ch2/(N_avo*density_ch2*L_target_ch2*nbeam_ch2_TrigA_norB)/nb_sr;

map <int,double> IntensityMap;
void LoadBeamIntensity(TFile* file){
	TTree* tree = (TTree*)file->Get("tree");
	int RunNumber;
	double Intensity;
	int ent = tree->GetEntries();
	tree->SetBranchAddress("RunNumber",&RunNumber);
	tree->SetBranchAddress("Intensity",&Intensity);
	for(int i=0;i<ent;++i){
		tree->GetEntry(i);
	}
}




//vector<TString> particle = {"P", "Pi1", "Pi2" , "L", "Xi"};
//vector<TString> Corvariable = {"CosTh", "Mom", "Ph", "Dist", "CosPsi", "CosOpen"};
//vector<TString> PolVariable = {"PolCTh", "PolCPhA", "PolCPhB", "PolCPhC"};
//vector<TString> CorCases = {"AcptEff", "InTargetEff", "TrackedEff", "TrigBAcptEff", "TrigBInTargetEff", "TrigBTrackedEff", "TrigBNeutralTrackedEff", "TrigBNeutralInTargetEff"};
double CosPsi(TVector3 v)
{
    return v.y() / v.Mag();
}
double DistT(TVector3 V)
{
    return hypot(V.x(), V.z() + 143);
}
double Dist(TVector3 V)
{
    TVector3 V_(V.X(),V.Y(),V.Z()+143);
    return V_.Mag();
}
double Omega(TVector3 v){
    return atan2(v.X(),v.Z());
}

double GetVariable(TVector3 P1, TVector3 P2, TVector3 V, TString v){
    double ret = 0;
    if (v == "CosTh"){
        ret = P1.CosTheta();
    }
    if (v == "Mom"){
        ret = P1.Mag();
    }
    if (v == "Ph"){
        ret = P1.Phi();
    }
    if (v == "DistT"){
        ret = DistT(V);
    }
    if (v == "Dist"){
        ret = Dist(V);
    }
    if (v == "CosPsi"){
        ret = CosPsi(P1);
    }
    if (v == "CosOpen"){
        // P1 is the mother, P2 and P1 - P2 are the daughters of two-body decay.
        ret = cos(P2.Angle(P1- P2));
    }
    if (v == "CosKink"){
        ret = cos(P1.Angle(P2+P1));
    }
    return ret;
}






//Common Methods//
double Pmx(double M,double m1, double m2){
	double ppp = (M+m1+m2)*(M+m1-m2)*(M-m1+m2)*(M-m1-m2)/(4.0*M*M);
	if (ppp>0) return std::sqrt(ppp);
	else       return -1.;
}
double WeightThPh(TEfficiency* eff, double Th, double Ph){
	int Bin = eff->FindFixBin(cos(Th),Ph);
	double effi = eff->GetEfficiency(Bin);
	return 1./effi;
}
double BlindTh(double Th, double rad){
	double cTh = cos(Th);
	if(abs(cTh)> rad ) return 1.;
	else{
		return 2.*acos(-1)/(2*acos(-1)-2* sqrt(rad*rad - cTh*cTh));
	}
}
double BlindPh(double Ph, double rad){
	double Ph_ = abs(Ph) - acos(-1)/2;
	if(abs(Ph_) > rad) return 1.;
	else{
		return 2./(2- sqrt(rad*rad - Ph_*Ph_));
	}
}
double Weight1D(TEfficiency* eff, double x){
	int Bin = eff->FindFixBin(x);
	double effi = eff->GetEfficiency(Bin);
	return 1./effi;
}
double WeightDistPsi(TEfficiency* eff, double dist, double psi){
	int Bin = eff->FindFixBin(dist,cos(psi));
	double effi = eff->GetEfficiency(Bin);
	return 1./effi;
}
double Weight2D(TEfficiency* eff, double x, double y){
	int Bin = eff->FindFixBin(x,y);
	double effi = eff->GetEfficiency(Bin);
	return 1./effi;
}

TGraphErrors* DivisionHistsSignal(TH1D* h_num, TH1D* h_den,double offset = 0, double ey_th = 0.2){
    int nbin = h_num->GetNbinsX();
    TGraphErrors* gr = new TGraphErrors();
    cout<<"Dividing hist :"<< h_num->GetName() << " / " << h_den->GetName() << endl;
    int np = 0;
    for(int ibin = 1; ibin <= nbin; ++ibin){
        double x = h_num->GetBinCenter(ibin) + offset * h_num->GetBinWidth(ibin);
        double y_num = h_num->GetBinContent(ibin);
        double y_den = h_den->GetBinContent(ibin);
        double ey_num = h_num->GetBinError(ibin);
        double ey_den = h_den->GetBinError(ibin);
        double y = 0;
        double ey = 0;
        if (y_den != 0){
            y = y_num / y_den;
            ey = sqrt( (ey_num*ey_num)/(y_den*y_den) + (y_num*y_num*ey_den*ey_den)/(y_den*y_den*y_den*y_den) );
        }
        if(ey < ey_th){
            gr->SetPoint(np, x, y);
            gr->SetPointError(np, 0, ey);
            np++;
        }
    }
    return gr;
}
TH1D* MakeDividedHistSignal(TH1D* h_num, TH1D* h_den, TString name, double ey_th = 0.2){
    TGraphErrors* gr = DivisionHistsSignal(h_num, h_den, 0, ey_th);
    TH1D* h_div = new TH1D(name, name, h_num->GetNbinsX(), h_num->GetXaxis()->GetXmin(), h_num->GetXaxis()->GetXmax());
    int npoints = gr->GetN();
    for(int ipt = 0; ipt < npoints; ++ipt){
        double x, y;
        gr->GetPoint(ipt, x, y);
        int bin = h_div->GetXaxis()->FindBin(x);
        double ey = gr->GetErrorY(ipt);
        h_div->SetBinContent(bin, y);
        h_div->SetBinError(bin, ey);
    }
    delete gr;
    return h_div;
}
TGraphErrors* DivisionHistsBGSubtracked(TH1D* h_num, TH1D* h_den,double offset = 0, double ey_th = 0.2){
    int nbin = h_num->GetNbinsX();
    TGraphErrors* gr = new TGraphErrors();
    cout<<"Dividing hist :"<< h_num->GetName() << " / " << h_den->GetName() << endl;
    int np = 0;
    for(int ibin = 1; ibin <= nbin; ++ibin){
        double x = h_num->GetBinCenter(ibin) + offset * h_num->GetBinWidth(ibin);
        double y_num = h_num->GetBinContent(ibin);
        double y_den = h_den->GetBinContent(ibin);
        double y = 0;
        double ey = 0;
        if (y_den != 0){
            y = y_num / y_den;
            double ey_num = h_num->GetBinError(ibin);
            double ey_den = h_den->GetBinError(ibin);
            double ey_bg = sqrt(ey_num * ey_num  - ey_den * ey_den);
            ey = 1./y_den * sqrt(ey_bg * ey_bg + (1-y)*(1-y)* ey_den * ey_den);
        }
        if(ey < ey_th){
            gr->SetPoint(np, x, y);
            gr->SetPointError(np, 0, ey);
            np++;
        }
    }
    return gr;
}
TH1D* MakeDividedHistBGSubtracked(TH1D* h_num, TH1D* h_den, TString name, double ey_th = 0.2){
    TGraphErrors* gr = DivisionHistsBGSubtracked(h_num, h_den, 0, ey_th);
    TH1D* h_div = new TH1D(name, name, h_num->GetNbinsX(), h_num->GetXaxis()->GetXmin(), h_num->GetXaxis()->GetXmax());
    int npoints = gr->GetN();
    for(int ipt = 0; ipt < npoints; ++ipt){
        double x, y;
        gr->GetPoint(ipt, x, y);
        int bin = h_div->GetXaxis()->FindBin(x);
        double ey = gr->GetErrorY(ipt);
        h_div->SetBinContent(bin, y);
        h_div->SetBinError(bin, ey);
    }
    delete gr;
    return h_div;
}









////Analysis functions////////







map<TString, TH1 *> hMap;
map<TString, TH2 *> hMap2D;
map<TString, TGraph*> gMap;
map<TString, TEfficiency *> effMap;
vector<TString> LCorrection;
vector<TString> XiCorrection;
vector<TString> targets = {"CH2", "Carbon"};
TH2D *hKuramaEff;
TH2D *hAcptPhi;

TString AcceptanceHistTitle1D(TString pre,TString part, TString var, TString suff){
  TString title = pre + part + var + suff;
  return title;
}
TString AcceptanceHistTitle2D(TString pre,TString part, TString var1, TString var2, TString suff){
  TString title = pre + part + var1 + var2 + suff;
  return title;
}

namespace Acpt{
    TString HistTitle1D(TString pre,TString part, TString var, TString suff){
        TString title = pre + part + var + suff;
        return title;
    }
    TString HistTitle2D(TString pre,TString part, TString var1, TString var2, TString suff){
        TString title = pre + part + var1 + var2 + suff;
        return title;
    }

};

bool CheckKey(TString key){
    bool check = true;
    if(hMap.find(key) == hMap.end() && hMap2D.find(key) == hMap2D.end() && gMap.find(key) == gMap.end() && effMap.find(key) == effMap.end()){
        check = false;
        cout<<"Key "<<key<<" not found in any map!"<<endl;
    }
    return check;
}
#if 0
void LoadAcceptanceMap(vector<TFile *> AcptFiles, TFile* KuramaAcpt){
    cout<<"Loading Acceptance Map..."<<endl;
#if KuramaPhi
	hKuramaEff = (TH2D*)KuramaAcpt->Get("BeamAcceptancePhi");
    hKuramaEff ->GetXaxis()->SetTitle("#phi_{K^{+}}");
    hKuramaEff ->GetYaxis()->SetTitle("P_{K^{+}}[GeV/c]");
    hKuramaEff->GetXaxis()->SetTitleSize(0.05);
    hKuramaEff->GetYaxis()->SetTitleSize(0.05);
#else
	hKuramaEff = (TH2D*)KuramaAcpt->Get("BeamAcceptance");
    hKuramaEff ->GetXaxis()->SetTitle("#theta_{K^{+}}[deg]");
    hKuramaEff ->GetYaxis()->SetTitle("P_{K^{+}}[GeV/c]");
    hKuramaEff->GetXaxis()->SetTitleSize(0.05);
    hKuramaEff->GetYaxis()->SetTitleSize(0.05);
#endif
#if RebinKurama > 1
    hKuramaEff->RebinX(RebinKurama);
    hKuramaEff->RebinY(RebinKurama);
    double Scale = 1./(RebinKurama * RebinKurama);
    hKuramaEff->Scale(Scale);
#endif
    vector<vector<TString>> params = {LCorrection, XiCorrection};
#if SavePDF
    TCanvas * c_Kurama = new TCanvas("Kurama","Kurama",1000,1000);
    gPad->SetMargin(0.15,0.1,0.15,0.1);
    hKuramaEff->Draw("colz");
    c_Kurama->SaveAs(figdir+"KuramaEff.png");
#endif

    for (int it = 0; it < targets.size(); it++){
        auto t = targets[it];
        cout<<"Target : "<<t<<endl;
        auto file = AcptFiles[it];
        TString key;
        vector<TString> pars = {params.at(0).at(0), params.at(1).at(0)}; 
        for (int ip= 0; ip<params.size();++ip){
            auto particle = pars[ip];
            cout<<particle<<endl;
            vector<TString> vars = {params.at(ip).at(1), params.at(ip).at(2)};  
            TString conf = params.at(ip).at(3);
            cout<<conf<<endl;
            for (int iv1 = 0; iv1 < vars.size(); ++iv1){
                auto v1 = vars[iv1];
                cout<<v1<<endl;
                key = Acpt::HistTitle1D("", particle, v1, conf);
                TString MapKey = "G4" + key;
                TString TargetKey = t + key;
#if Skip1D && 0
#else
                cout<<"Loading "<<MapKey<<endl;
                auto Eff1D = (TEfficiency *)file->Get(MapKey);
                Eff1D->SetName(TargetKey);
                effMap[TargetKey] = Eff1D;
#endif
                for (int iv2 = iv1 + 1; iv2 < vars.size(); ++iv2){
                    TString ct = "canv"+t+ Acpt::HistTitle2D("", particle, v1, vars[iv2], "");
                    auto v2 = vars[iv2];
                    cout<<v2<<endl;
                    key = Acpt::HistTitle2D("", particle, v1, v2, conf);
                    MapKey = "G4" + key;
                    TargetKey = t + key;
                    auto Eff2D = (TEfficiency *)file->Get(MapKey);
                    Eff2D->SetName(TargetKey);
                    effMap[TargetKey] = Eff2D;
                    if(!Eff2D){
                        cout<<"Error: "<<MapKey<<" not found!"<<endl;
                        continue;
                    }
                    for(int id = 0; id < ndist-1; ++id){
                        TString distKey = GetDistKey(id);
                        MapKey = "G4" + key + distKey;
                        TargetKey = t + key + distKey;
                        auto Eff2D = (TEfficiency *)file->Get(MapKey);
                        if(!Eff2D){
                            cout<<"Error: "<<MapKey<<" not found!"<<endl;
                            continue;
                        }
                        Eff2D->SetName(TargetKey);
                        effMap[TargetKey] = Eff2D;
//                        cout<<"Loaded "<<TargetKey<<endl;
                    }
#if SavePDF
                    TCanvas *c_map = new TCanvas(ct,ct,1000,1000);
                    gPad->SetMargin(0.15,0.1,0.15,0.1);
                    Eff2D->Draw("colz");
                    //Eff2D->GetPaintedGraph()->GetXaxis()->SetTitle(Form("%s_{%s}",v1.Data(),p.Data()));
                    //Eff2D->GetPaintedGraph()->GetYaxis()->SetTitle(Form("%s_{%s}",v2.Data(),p.Data()));
                    //Eff2D->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.05);
                    //Eff2D->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.05);
                    c_map->SaveAs(figdir+TargetKey+".png");
#endif
                }
            }
        }
    }
    cout<<"Done!"<<endl;
}
void LoadAcceptanceMapAll(vector<TFile *> AcptFiles, TFile* KuramaAcpt){
    cout<<"Loading Acceptance Map..."<<endl;
#if KuramaPhi
	hKuramaEff = (TH2D*)KuramaAcpt->Get("BeamAcceptancePhi");
#else
	hKuramaEff = (TH2D*)KuramaAcpt->Get("BeamAcceptance");
#endif
    vector<vector<TString>> params = {LCorrection, XiCorrection};
    for (int it = 0; it < targets.size(); it++)
    {
        auto t = targets[it];
        auto file = AcptFiles[it];
        TKey* key;
        TIter nextkey(file->GetListOfKeys());
        while ((key = (TKey*)nextkey()))
        {
            TObject* obj = key->ReadObj();
            if(obj->InheritsFrom("TEfficiency")){
                auto Eff = (TEfficiency*)obj;
                TString name = key->GetName();
                name.ReplaceAll("G4",t);
                Eff->SetName(name);
                effMap[name] = Eff;
            }
        }
        
    }
    cout<<"Done!"<<endl;
}


#endif


static TString eq_sigmaM2 = "4*TMath::Power([3], 2)*(1.+[3]/(x*x))*[0]+4*TMath::Power([3], 2)*x*x*[1]+4*x*x*(x*x+[3])*[2]";

bool KaonSelection(Double_t mass2, Double_t mom, Double_t nsigma)
{
	if(TMath::IsNaN(mass2) || mass2<0) return false;
	if(mom > 1.4) return false;
	Double_t pdgmass2 = mk*mk;

	//Measured values(sigma of M2 spectrum)
	TF1 *f_sigmaM2 = new TF1("f_sigmaM2", eq_sigmaM2.Data(), 0., 5.); 
	f_sigmaM2 -> FixParameter(0, 0.000193292);
	f_sigmaM2 -> FixParameter(1, -0.000407003);
	f_sigmaM2 -> FixParameter(2, 0.000134182);
	f_sigmaM2 -> FixParameter(3, pdgmass2);

	Double_t m2cut = nsigma*TMath::Sqrt(f_sigmaM2 -> Eval(mom)); //nsigma cut for M^2
	return (TMath::Abs(mass2 - pdgmass2) < m2cut);
}
double GetMass(int pdg){
	double mass = 0;
	switch(abs(pdg))
	{
		case 211:
			mass = mpi;
			break;
		case 321:
			mass = mk;
			break;
		case 2212:
			mass = mp;
			break;
		case 3122:
			mass = mLd;
			break;
		case 3312:
			mass = mXi;
			break;
		case 3322:
			mass = mXi1530;
			break;
		case 111:	
			mass = pi0mass;
			break;
		case 310:
			mass = mks;
			break;
		default:
			mass = 0;
			break;
	}
	return mass;

}
double solid_riemann(double *theta, double mom_kp, TH2D *hist_acceptance)
{

    double binsize = 0.1; // deg.
    double bins = 0;
    double solid = 0;
    while (true)
    {
        double mint = theta[0] + binsize * bins;
        double t = theta[0] + binsize * (0.5 + bins);
        double maxt = theta[0] + binsize * (1. + bins);
        if (t >= theta[1])
            break;
        int bin = hist_acceptance->FindBin(t, mom_kp);
        double acceptance = hist_acceptance->GetBinContent(bin);
        solid += 2. * TMath::Pi() * (TMath::Cos(mint * TMath::DegToRad()) - TMath::Cos(maxt * TMath::DegToRad())) * acceptance; // sr
        bins++;
    }
    return solid;
}
double acpt_kp(double theta, double mom_kp, TH2D *hist_acceptance){
  int bin = hist_acceptance -> FindBin(theta, mom_kp);
  double acceptance = hist_acceptance -> GetBinContent(bin);
  return acceptance;
}
TString KuramaMomBin(double mom){
    double mom_mev = mom*1000;
    int imom = mom_mev / 50;
    TString name  = Form("P_%d_%d", imom*50, (imom+1)*50);
    return name;
}   
TString KuramaName(double mom){
    double mom_mev = mom*1000;
    int imom = mom_mev / 50;
    TString name  = "BeamAcceptanceThPh_" + KuramaMomBin(mom); 
    return name;
}
double pkp_acpt_0 = 0.4;
double pkp_acpt_1 = 2;
double dpkp_acpt = 0.01;
TString KuramaMomBinNarrow(double mom){
    double mom_mev = mom*1000;
    int dp = dpkp_acpt*1000;
    int imom = mom_mev / dp;
    TString name  = Form("P_%d_%d", imom*dp, (imom+1)*dp);
    return name;
}   

TString KuramaNameNarrow(double mom){
    double mom_mev = mom*1000;
    double dp = dpkp_acpt*1000;
    int imom = mom_mev / dp;
    TString name  = "BeamAcceptanceThPh_" + KuramaMomBinNarrow(mom); 
    return name;
}
void LoadNarowAcceptance(TString filename){
    TFile* file = new TFile(filename);
    for(int imom = pkp_acpt_0*1000; imom < pkp_acpt_1*1000; imom += dpkp_acpt*1000){
        double mom = imom/1000.;
        TString name = KuramaNameNarrow(mom);
        if(hMap2D[name]){
            hMap2D[name]->Delete();
        }
        hMap2D[name] = (TH2D*)file->Get(name);
    }
}
double IsKpAccepted(TVector3 kp_mom, TH2D* hMap){
  double theta = kp_mom.Theta();
  theta = theta * TMath::RadToDeg();
  double phi_kp = kp_mom.Phi();
  double mom_kp = kp_mom.Mag(); // GeV/c
  int bin = hMap -> FindBin(theta, phi_kp);
  double acceptance = hMap -> GetBinContent(bin);
  return acceptance;
}
double CalcAcceptance(TVector3 km_mom, TVector3 kp_mom, TH2D* hMap, int nPhi = 1000){
    auto zHat = km_mom.Unit();
    auto xHat = (zHat.Cross(TVector3(0,1,0))).Unit();
    auto yHat = zHat.Cross(xHat);
    auto dir = kp_mom.Unit();
    double cth = dir.Dot(zHat);
    double sum_acceptance = 0;
    double phi_kp = acos(dir.Dot(xHat)/sin(acos(cth)));
    for(int iphi = 0; iphi < nPhi; ++iphi){
        double phi0 = -PI + 2*PI/nPhi*iphi;
        double phi1 = -PI + 2*PI/nPhi*(iphi+1);
        double phi = (phi0 + phi1)/2;
        #if PhiRegion
        phi0 = phi_min + (phi_max - phi_min)/nPhi*iphi;
        phi1 = phi_min + (phi_max - phi_min)/nPhi*(iphi+1);
        phi = (phi0 + phi1)/2;
        #endif
        TVector3 rot_dir = sin(acos(cth))*cos(phi)*xHat
                        + sin(acos(cth))*sin(phi)*yHat
                        + cth*zHat;
        double acceptance = IsKpAccepted(rot_dir*kp_mom.Mag(), hMap);
        sum_acceptance += acceptance;
    }
    return sum_acceptance * 1./nPhi;
}
double CalcAcceptance(TVector3 km_mom, TVector3 kp_mom, int nPhi = 1000){
    if(kp_mom.Mag() < pkp_acpt_0 or kp_mom.Mag() > pkp_acpt_1){
        return 0;
    }
    TH2D* heffMap = (TH2D*)hMap2D[KuramaNameNarrow(kp_mom.Mag())];
    auto zHat = km_mom.Unit();
    auto xHat = (zHat.Cross(TVector3(0,1,0))).Unit();
    auto yHat = zHat.Cross(xHat);
    auto dir = kp_mom.Unit();
    double cth = dir.Dot(zHat);
    double sum_acceptance = 0;
    double phi_kp = acos(dir.Dot(xHat)/sin(acos(cth)));
    for(int iphi = 0; iphi < nPhi; ++iphi){
        double phi0 = -PI + 2*PI/nPhi*iphi;
        double phi1 = -PI + 2*PI/nPhi*(iphi+1);
        double phi = (phi0 + phi1)/2;
        #if PhiRegion
        phi0 = phi_min + (phi_max - phi_min)/nPhi*iphi;
        phi1 = phi_min + (phi_max - phi_min)/nPhi*(iphi+1);
        phi = (phi0 + phi1)/2;
        #endif
        TVector3 rot_dir = sin(acos(cth))*cos(phi)*xHat
                        + sin(acos(cth))*sin(phi)*yHat
                        + cth*zHat;
        double acceptance = IsKpAccepted(rot_dir*kp_mom.Mag(), heffMap);
        sum_acceptance += acceptance;
    }
    return sum_acceptance * 1./nPhi;
}
TVector3 TargetMomentum(double* pars, double mom0, double charge){
    double helix_cx = pars[0];
    double helix_cy = pars[1];
    double t = atan2(helix_cy, helix_cx);
    double dz = pars[4];
    double pz = -dz*mom0*charge;
    double trans = -1./hypot(1,dz)*mom0*charge;
    double px = sin(t)*trans;
    double py = cos(t)*trans;
    TVector3 p_vec(-px, pz, py);
    return p_vec;
}
TVector3 RandomVect(){
    double cth = gRandom->Uniform(-1,1);
    double sth = sqrt(1.-cth*cth);
    double phi = gRandom->Uniform(0, 2*PI);
    double x = sth * cos(phi);
    double y = sth * sin(phi);
    double z = cth;
    return TVector3(x,y,z);
}
TVector3 pFermi(double pF_min = 0, double pF_max = 0.250){
    double pF = pow(gRandom->Uniform(pow(pF_min,3), pow(pF_max,3)), 1./3.);
    TVector3 dir = RandomVect();
    return dir * pF;
}

vector<int> CountLayers(int* nhit){
	vector<int> Multis;
	Multis.resize(10,0);
	for(int il = 0; il < 12; ++il){
		if(il == 6 or il == 7) continue;
		int nh = nhit[il];
		if(nh > 9) nh = 9;
		for(int ih = 0; ih< nh+1; ++ih){
			Multis[ih]++;
		}
	}
	return Multis;
}
int CountTracks(int* nhit, int th){
    auto layer_multi = CountLayers(nhit);
    int nt = -1;
    for(int il=0; il < layer_multi.size(); ++il){
        if(layer_multi[il] >= th){
            nt++;
        }
    }
    return nt;
}

void DisableLargestBranches(TTree* tree, int n = 20){
  struct BI { TBranch* b; Long64_t s; };
  std::vector<BI> v;

  TObjArray* bl = tree->GetListOfBranches();
  for (int i = 0; i < bl->GetEntries(); ++i) {
    auto b = (TBranch*)bl->At(i);
    v.push_back({b, b->GetTotBytes()});
  }

  std::sort(v.begin(), v.end(),
            [](auto& a, auto& b){ return a.s > b.s; });
  for (int i = 0; i < std::min(n, (int)v.size()); ++i){
		cout<<Form("Disabling %d th branch : %s",i, v[i].b->GetName())<<endl;
		tree->SetBranchStatus(v[i].b->GetName(), 0);
	}
}


#endif	
