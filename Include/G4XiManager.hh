#ifndef G4XiManager_hh
#define G4XiManager_hh
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "./Dir.hh"
#include "g4genfit.h"
#include "XiEfficiency.hh"
using namespace std;
double weight_th = 50;
class G4XiManager{
private:
    vector<TString> particle = {"P", "Pi1", "Pi2" , "L", "Xi"};
    vector<TString> Corvariable = {"CosTh", "Mom", "Ph", "Dist", "CosPsi", "CosOpen"};
    vector<CorrectionParameter> CorParams;

public:
    G4XiManager(){};
    ~G4XiManager(){};
public:
    G4XiManager(g4genfit* gf_file, TString target){
        this->gf = gf_file;
        SetTarget(target);
        InitializeEvent();
    }
    G4XiManager(g4genfit* gf_file, TString target, TString file_name){
        this->gf = gf_file;
        SetTarget(target);
        InitializeEvent();
        TFile* f_cor = new TFile(file_name, "read");
        vector<TString> pola_hists = {
            "h_PolCThGen0HistEffAcpt",
            "h_PolCPhAGen0HistEffAcpt",
            "h_PolCPhBGen0HistEffAcpt",
            "h_PolCPhCGen0HistEffAcpt"
        };
        vector<TString> pola_keys = {
            "Th", "PhA", "PhB", "PhC"
        };
        for (auto p : particle){
            for(auto v: Corvariable){
                TString key_hist = "h_" + p + v + "Gen0" + "HistEffAcpt";
                TString key_name = p + v;
                hPartMap[key_name] = (TH1D*)f_cor->Get(key_hist);
            }
        }
        for (size_t i = 0; i < pola_hists.size(); ++i) {
            TString hist_name = pola_hists[i];
            TString key_name = pola_keys[i];
            TH1D* hist = (TH1D*)f_cor->Get(hist_name);
            hPolaMap[key_name] = (TH1D*)hist->Clone();
        }

    }
    g4genfit* GetGF() { return gf; }
    void InitializeEvent();
    void ProcessEvent();
    void SetTarget(TString t) {
         target = t;
     }
     void SetCorrectionParameter(TString num, TString den, TString p_cor, vector<TString> var_cor){
         CorParams.push_back({num, den, p_cor, var_cor});
         cout<<"Added correction parameter: "<<num<<"/"<<den<<" for "<<p_cor<<" with variables ";
         for(auto v:var_cor) cout<<v<<" ";
         cout<<endl;
     }


    TString GetTarget(){ return target; }
    void Lighten();
private:
    g4genfit* gf;
    bool debug = 0;
private:
    //MetaParameters
    TString target;
    int conf;

private:
    //Phisics Variables
    TVector3 km_mom, kp_mom ;
    TLorentzVector lv_km,lv_kp,lv_cm;
    double cth_cm, sqrt_sKmp;
    TVector3 beta;
    
    TVector3 mom_pGen, mom_pi1Gen, mom_pi2Gen, mom_lGen, mom_xiGen, mom_xiprodGen;
    TVector3 mom_p, mom_pi1, mom_pi2, mom_l, mom_xi, mom_xiprod;
    TVector3 km_cm, kp_cm;
    
    TVector3 l_vertexGen, xi_vertexGen, xiprod_vertexGen;
    TVector3 l_vertex, xi_vertex, xiprod_vertex;
    map<TString, vector<TVector3>> DCGen;// DataContainer
    map<TString, vector<TVector3>> DC;

    bool AcptFlag;
    bool Xiflag;
    bool KKVtxFlag;
    double weightGen = 1;
    vector<double> weightVecGen;
    double weight = 1;
    vector<double> weightVec;


    double pola_cthweight;
    double pola_phaweight;
    double pola_phbweight;
    double pola_phcweight;
    map<TString, TH1D*> hPartMap;
    map<TString, TH1D*> hPolaMap;
    map<TString, double> weightPartMap;

public:

    bool IsAccepted(){ return AcptFlag; }
    bool IsXiFlag(){ return Xiflag; }
    bool IsKKVtxFlag(){ return KKVtxFlag; }
    double GetWeightPolaGen(TString pola, double val){
        auto bin = hPolaMap[pola]->FindBin(val);
        double eff = hPolaMap[pola]->GetBinContent(bin);
        return 1./eff;
    }
    double GetWeightParticle(TString part, TString var, double val){
        if(hPartMap[part + var] == nullptr){
            cout << "No hist for " << part + var << endl;
            return 1.;
        }
        auto bin = hPartMap[part + var]->FindBin(val);
        if(bin == -1){
            cout<<"No bin for "<< part + var << " with value "<< val << endl;
            return 1.;
        } 
        double eff = hPartMap[part + var]->GetBinContent(bin);
				if(eff < 1./weight_th) return 1;
        return 1./eff;
    }
    double GetWeight(TString part, TString num, TString den, TString v1, TString v2 = "", map<TString, vector<TVector3>> DataCont = {}){
        if(DataCont.empty()) DataCont = DC;
        TString key;
        double var1,var2;
        var1 = GetVariable(DataCont[part][0], DataCont[part][1], DataCont[part][2], v1);
        double w;
        if(v2 != ""){
            var2 = GetVariable(DataCont[part][0], DataCont[part][1], DataCont[part][2], v2);
            key = EffTitle2D(target, part, v1, v2, num, den);
#if PerDist
            if(v1 != "DistT" and v2 != "DistT") key += GetDistKey(DataCont[part][2]);
#endif
            if(effMap[key] == nullptr){
                cout<<"Efficiency histogram "<<key<<" not found!"<<endl;
                return 1.;
            }
            w = effMap[key]->GetEfficiency(effMap[key]->FindFixBin(var1, var2));
            w = 1./w;
        }
        else{
            key = EffTitle1D(target, part, v1, num, den);
            if(effMap[key] == nullptr){
                cout<<"Efficiency histogram "<<key<<" not found!"<<endl;
                return 1.;
            }
            w = effMap[key]->GetEfficiency(effMap[key]->FindFixBin(var1));
            w = 1./w;
        }
        return w;
    }
    double GetWeightGen(int d = -1){
        if(d == -1) return weightGen;
        else if(d < weightVecGen.size()) return weightVecGen[d];
        else {
            cout<<"No weight for index "<<d<<endl;
            return 1.;
        }
    }
    double GetWeight(int d = -1){
        if(d == -1) return weight;
        else if(d < weightVec.size()) return weightVec[d];
        else {
            cout<<"No weight for index "<<d<<endl;
            return 1.;
        }
    }
    map<TString, vector<TVector3>> GetDataContainer(int isgen){
        return isgen ? DCGen : DC;
    }
public:
    TVector3 GetKmMom() { return km_mom; }
    TVector3 GetKpMom() { return kp_mom; }
    double GetCosThetaCM() { return cth_cm; }


    TVector3 GetProtonMomGen() { return mom_pGen; }
    TVector3 GetPi1MomGen() { return mom_pi1Gen; }
    TVector3 GetPi2MomGen() { return mom_pi2Gen; }
    TVector3 GetLambdaMomGen() { return mom_lGen; }
    TVector3 GetXiMomGen() { return mom_xiGen; }
    TVector3 GetXiProdMomGen() { return mom_xiprodGen; }
    TVector3 GetLambdaVtxGen() { return l_vertexGen; }
    TVector3 GetXiVtxGen() { return xi_vertexGen; }
    TVector3 GetXiProdVtxGen() { return xiprod_vertexGen; }
    
    TVector3 GetProtonMomRec() { return mom_p; }
    TVector3 GetPi1MomRec() { return mom_pi1; }
    TVector3 GetPi2MomRec() { return mom_pi2; }
    TVector3 GetLambdaMomRec() { return mom_l; }
    TVector3 GetXiMomRec() { return mom_xi; }
    TVector3 GetXiProdMomRec() { return mom_xiprod; }
    TVector3 GetLambdaVtxRec() { return l_vertex; }
    TVector3 GetXiVtxRec() { return xi_vertex; }
    TVector3 GetXiProdVtxRec() { return xiprod_vertex; }
    TVector3 GetParticle(TString name){
        if(DC.count(name) == 0){
            cout<<"No particle with name "<<name<<endl;
            return TVector3(0,0,0);
        }
        return DC[name][0];
    }

};
void
G4XiManager::InitializeEvent(){
    AcptFlag = 1;
    km_mom.SetXYZ(0, 0, 0);
    kp_mom.SetXYZ(0, 0, 0);
    mom_pGen.SetXYZ(0, 0, 0);
    mom_pi1Gen.SetXYZ(0, 0, 0);
    mom_pi2Gen.SetXYZ(0, 0, 0);
    mom_lGen.SetXYZ(0, 0, 0);
    mom_xiGen.SetXYZ(0, 0, 0);
    mom_xiprodGen.SetXYZ(0, 0, 0);
   
    l_vertexGen.SetXYZ(0, 0, 0);
    xi_vertexGen.SetXYZ(0, 0, 0);
    xiprod_vertexGen.SetXYZ(0,0,0);

    
    mom_p.SetXYZ(0, 0, 0);
    mom_pi1.SetXYZ(0, 0, 0);
    mom_pi2.SetXYZ(0, 0, 0);
    mom_l.SetXYZ(0, 0, 0);
    mom_xi.SetXYZ(0, 0, 0);
    mom_xiprod.SetXYZ(0, 0, 0);
   
    l_vertex.SetXYZ(0, 0, 0);
    xi_vertex.SetXYZ(0, 0, 0);
    xiprod_vertex.SetXYZ(0,0,0);

    Xiflag = 0;
    KKVtxFlag = 0;

    weightGen = 1;
    weight = 1;
    weightVecGen.clear();
    weightVec.clear();
    DCGen.clear();
    DC.clear();
}
void
G4XiManager::ProcessEvent(){
    InitializeEvent();
    km_mom.SetXYZ(gf->KmMom_x->at(0), gf->KmMom_y->at(0), gf->KmMom_z->at(0));
    kp_mom.SetXYZ(gf->KpMom_x->at(0), gf->KpMom_y->at(0), gf->KpMom_z->at(0));
    
    lv_kp = TLorentzVector(kp_mom, hypot(kp_mom.Mag(), mk));
    lv_km = TLorentzVector(km_mom, hypot(km_mom.Mag(), mk));
    TLorentzVector lv_target = TLorentzVector(0, 0, 0, mp);
    lv_cm = lv_km + lv_target;
    sqrt_sKmp = lv_cm.M();
    beta = lv_cm.BoostVector();
    TLorentzVector lv_kmcm = lv_km;
    TLorentzVector lv_kpcm = lv_kp;
    lv_kmcm.Boost(-beta);
    lv_kpcm.Boost(-beta);
    km_cm = lv_kmcm.Vect();
    kp_cm = lv_kpcm.Vect();
    cth_cm = cos(km_cm.Angle(kp_cm));

    mom_pGen.SetXYZ(gf->G4pmom_x, gf->G4pmom_y, gf->G4pmom_z);
    mom_pi1Gen.SetXYZ(gf->G4pi1mom_x, gf->G4pi1mom_y, gf->G4pi1mom_z);
    mom_pi2Gen.SetXYZ(gf->G4pi2mom_x, gf->G4pi2mom_y, gf->G4pi2mom_z);
    mom_lGen = mom_pGen + mom_pi1Gen;
    mom_xiGen = mom_lGen + mom_pi2Gen;
    mom_xiprodGen.SetXYZ(gf->G4ximom_x, gf->G4ximom_y, gf->G4ximom_z);

    l_vertexGen.SetXYZ(gf->G4pi1vtx_x, gf->G4pi1vtx_y, gf->G4pi1vtx_z);
    xi_vertexGen.SetXYZ(gf->G4pi2vtx_x, gf->G4pi2vtx_y, gf->G4pi2vtx_z);
    xiprod_vertexGen.SetXYZ(gf->G4xivtx_x, gf->G4xivtx_y, gf->G4xivtx_z);

    DCGen = {
        {"P", {mom_pGen, mom_pi1Gen, l_vertexGen}},
        {"Pi1", {mom_pi1Gen, mom_pGen, l_vertexGen}},
        {"Pi2", {mom_pi2Gen, mom_lGen, xi_vertexGen}},
        {"L", {mom_lGen, mom_pi1Gen, l_vertexGen}},
        {"Xi", {mom_xiGen, mom_pi2Gen, xi_vertexGen}}
    };
    for(auto cp:CorParams){
        TString num = cp.num;
        TString den = cp.den;
        TString pc = cp.p_cor;
        vector<TString> var_cor = cp.var_cor;
        TString v1 = var_cor[0];
        TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
        double w =GetWeight(pc, num, den, v1, v2, DCGen);
        weightGen *= w;
        weightVecGen.push_back(w);
    }

    Xiflag = gf->Xiflag;
    if(Xiflag){
        mom_p.SetXYZ(gf->KFXiDecaysMom_x->at(0), gf->KFXiDecaysMom_y->at(0), gf->KFXiDecaysMom_z->at(0));
        mom_pi1.SetXYZ(gf->KFXiDecaysMom_x->at(1), gf->KFXiDecaysMom_y->at(1), gf->KFXiDecaysMom_z->at(1));
        mom_pi2.SetXYZ(gf->KFXiDecaysMom_x->at(2), gf->KFXiDecaysMom_y->at(2), gf->KFXiDecaysMom_z->at(2));
        mom_l.SetXYZ(gf->KFLambdaMom_x,gf->KFLambdaMom_y,gf->KFLambdaMom_z);
        mom_xi.SetXYZ(gf->KFXiMom_x,gf->KFXiMom_y,gf->KFXiMom_z);
        mom_xiprod.SetXYZ(gf->KFXiProductionVtxMom_x,gf->KFXiProductionVtxMom_y,gf->KFXiProductionVtxMom_z);
        l_vertex.SetXYZ(gf->GFLambdaDecayVtx_x,gf->GFLambdaDecayVtx_y,gf->GFLambdaDecayVtx_z);  
        xi_vertex.SetXYZ(gf->GFXiDecayVtx_x,gf->GFXiDecayVtx_y,gf->GFXiDecayVtx_z);
        xiprod_vertex.SetXYZ(gf->KFXiProductionVtx_x,gf->KFXiProductionVtx_y,gf->KFXiProductionVtx_z);
        if(mom_xiprod.Mag() != 0 and !isnan(mom_xiprod.Mag())) 
            KKVtxFlag = 1;
        else
            KKVtxFlag = 0;
        DC = {
            {"P", {mom_p, mom_pi1, l_vertex}},
            {"Pi1", {mom_pi1, mom_p, l_vertex}},
            {"Pi2", {mom_pi2, mom_l, xi_vertex}},
            {"L", {mom_l, mom_pi1, l_vertex}},
            {"Xi", {mom_xi, mom_l, xi_vertex}}
        };
        for(auto cp:CorParams){
            TString num = cp.num;
            TString den = cp.den;
            TString pc = cp.p_cor;
            vector<TString> var_cor = cp.var_cor;
            TString v1 = var_cor[0];
            TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
            double w =GetWeight(pc, num, den, v1, v2);
            weight *= w;
            weightVec.push_back(w);
        }
    }
}
void
G4XiManager::Lighten(){
    TTree* tree = gf->fChain;
    tree->SetBranchStatus("*", 0);
    vector<TString> xyz = {"_x", "_y", "_z"};
    vector<TString> vects = {
        "KmMom", "KpMom",
        "G4pmom", "G4pi1mom", "G4pi2mom","G4lmom", "G4ximom",
        "G4pi1vtx", "G4pi2vtx", "G4xivtx",
        "KFLambdaMom", "KFXiMom", "KFXiProductionVtxMom","KFXiDecaysMom",
        "GFLambdaDecayVtx", "GFXiDecayVtx", "KFXiProductionVtx"
    };
    for(auto v: vects){
        for(auto c: xyz){
            tree->SetBranchStatus((v + c).Data(), 1);
        }
    }
    tree->SetBranchStatus("Xiflag", 1);
    tree->SetBranchStatus("KFXiPval", 1);
}
#endif
