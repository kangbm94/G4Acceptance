#include "Include/g4genfitxipscat.h"
#include "Include/XiPScatEfficiency.hh"
#include "Include/PolarizationAnal.hh"
using namespace std;
#include <stdexcept>
#include <limits>
#include <TNamed.h>
//#include "G4PolaBranches.hh"
double thr_weight = 20; // cut under 5%

int nbinMom = 20;
double minMomP = 0, maxMomP = 2, minMomPi = 0, maxMomPi = 0.5;

int nbinCosTh = 20;
double minCosTh = -1, maxCosTh = 1;
int nbinAngle = 20;
double minAngle = -3.2, maxAngle = 3.2;

int nbinKpMom = 20;
double minKpMom = 0.5, maxKpMom = 1.5;
int nbinBE = 50;
double minBE = -100, maxBE = 800;

int nbinDistT = 30;
nbinDistT = 30;//After 260521
//double minDistT = 0, maxDistT = 350;
double minDistT = 0, maxDistT = 250;//After 260521
double minDistTXi = 0, maxDistTXi = 100;//After 260521
int nbinTrackLen = 30;
int nbinTrackLen2D = 30;
double minTrackLen = 0, maxTrackLen = 300;

int nbinLMom = 30;
double minMomL = 0, maxMomL = 1.5;
int nbinCosThL = 20;
double minCosThL = 0.6, maxCosThL = 1;

int nbinq = 30;
double minq = 0.2, maxq = 1.6;

double GFltarget_distcut = 25.;
double ltarget_distcut = 25.;
double GFltarget_ycut = 20.;

bool XiPOfflineFlag(g4genfitcarbon* XiEv){
  if(XiEv->XiPflag and XiEv->XipVtxDist < 5)
    return true;
  return false;
}

// Use the standard Xi definitions; proton truth flags are MC-only.
bool SuffixCheck(TString suff, g4genfitcarbon* XiEv){
  suff.ReplaceAll("TrigB", "");
  if(suff == "Gen") return true;
  if(suff == "LAcpt") return XiEv->lgood or XiEv->Xiflag;
  if(suff == "GoodL") return XiEv->lgood;
  if(suff == "XiAcpt") return XiEv->Xiflag;
  if(suff == "GoodXi") return XiEv->xigood;
  if(suff == "PTracked") return XiEv->p_tracked;
  if(suff == "Pi1Tracked") return XiEv->pi1_tracked;
  if(suff == "Pi2Tracked") return XiEv->pi2_tracked;
  if(suff == "PPi1Tracked") return XiEv->p_tracked and XiEv->pi1_tracked;
  if(suff == "AllTracked") return XiEv->p_tracked and XiEv->pi1_tracked and XiEv->pi2_tracked;
  if(suff == "GoodLAndPi2Tracked") return XiEv->lgood and XiEv->pi2_tracked;
  if(suff == "ScatPTracked") return XiEv->p_scat_tracked;
  if(suff == "ScatPGood") return XiEv->p_scat_good and XiPOfflineFlag(XiEv);
  if(suff == "GoodXiAndScatPTracked") return XiEv->xigood and XiEv->p_scat_tracked and XiPOfflineFlag(XiEv);
  if(suff == "GoodXiAndScatPGood") return XiEv->xigood and XiEv->p_scat_good and XiPOfflineFlag(XiEv);
  if(suff == "XiPflag") return XiPOfflineFlag(XiEv);
  return false;
}
bool TrigCheck(TString Trig, g4genfitcarbon* XiEv){
  bool ret = false;
  if(Trig == "") ret = 1;
  if(Trig == "TrigB") ret = (XiEv->nhHtof > 1);
  return ret;
}


TString SetAxis(TString part, TString var, int& bin, double& min, double& max){
  TString axis = "";
  if(var == "Mom"){
    axis = "Momentum [GeV/c]";
    if(part.Contains("Pi")){
      bin = nbinMom; min = minMomPi; max = maxMomPi;
    }
    else if(part == "ScatP"){
      bin = nbinMom; min = 0; max = maxMomP;
    }
    else if(part.Contains("L") or part.Contains("Xi") or part.Contains("P")){
      bin = nbinMom; min = minMomL; max = maxMomL;
    }
  }
  else if(var == "CosTh"){
    axis = "cos#theta";
    if(part.Contains("Pi") or part == "ScatP"){
      bin = nbinCosTh; min = minCosTh; max = maxCosTh;
    }
    else{
      bin = nbinCosTh; min = 0; max = maxCosTh;
    }
  }
  else if(var == "Ph"){
    axis = "#phi [rad]";
    bin = nbinAngle; min = minAngle; max = maxAngle;
  }
  else if(var == "DistT"){
    axis = "Distance from target [mm]";
    bin = nbinDistT; min = minDistT; max = maxDistT;
    if(part.Contains("Pi2") or part.Contains("Xi") or (LVtxIsXi and part.Contains("L"))){
      //bin = nbinMom; min = minMomL; max = maxMomL;After 260521
      bin = nbinDistT; min = minDistTXi; max = maxDistTXi;
    }
  }
  else if(var == "CosPsi"){
    axis = "cos#Psi";
    bin = nbinCosTh; min = minCosTh; max = maxCosTh;
  }
  else if(var == "CosOpen"){
    axis = "cos#alpha_{open}";
    bin = nbinCosTh; min = minCosTh; max = maxCosTh;
  }
  return axis;
};

void ResetXiPScatHistograms(){
  for(auto hist:hMap) delete hist.second;
  for(auto hist:hMap2D) delete hist.second;
  for(auto eff:effMap) delete eff.second;
  hMap = CheckedMap<TH1>();
  hMap2D = CheckedMap<TH1>();
  effMap = CheckedMap<TEfficiency>();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);
}
TH1* ReadXiPScatHistogram(TFile* file, TString key){
  if(file == nullptr or file->IsZombie()) throw runtime_error("Cannot read Xi-p histogram file.");
  TH1* hist = dynamic_cast<TH1*>(file->Get(key));
  if(hist == nullptr) throw runtime_error(string("Missing histogram: ") + key.Data());
  return hist;
}
void InitializeHistograms(){
  int nbin, nbinx, nbiny;
  double minx, maxx, miny, maxy;
  TString Xtitle, Ytitle;
  for(auto t:triggers){
  for(auto s:suffix){
    s = s + t;
    for(auto p:particle){
      for(int iv1=0;iv1<variable.size();++iv1){
        auto v1 = variable[iv1];
        Xtitle = SetAxis(p, v1, nbinx, minx, maxx);
        TString key = AcceptanceHistTitle1D("G4", p, v1, s);
#if Debug
        cout<<"Making "<<key<<endl;
#endif
#if LooseBin
  nbinx = nbinx/2;
#endif

	hMap[key] = new TH1D(key, ";" + Xtitle, nbinx, minx, maxx);
  for(int id = 0; id < ndist-1;++id){
    if(v1 == "DistT") continue;
    TString distkey = key + GetDistKey(id);
    hMap[distkey] = new TH1D(distkey, ";" + Xtitle, nbinx, minx, maxx);
  }
	for(int iv2=iv1+1;iv2<variable.size();++iv2){
	  auto v2 = variable[iv2];
    Ytitle = SetAxis(p, v2, nbiny, miny, maxy);

	  key = AcceptanceHistTitle2D("G4", p, v1, v2, s);
#if Debug
	  cout<<"Making "<<key<<endl;
#endif
#if LooseBin
  nbiny = nbiny/2;
#endif

	  hMap2D[key] = new TH2D(key, ";" + Xtitle + ";" + Ytitle, nbinx, minx, maxx, nbiny, miny, maxy);
    for(int id = 0; id < ndist-1;++id){
      if(v1 == "DistT" or v2 == "DistT") continue;
      TString distkey = key + GetDistKey(id);
      hMap2D[distkey] = new TH2D(distkey, ";" + Xtitle + ";" + Ytitle, nbinx, minx, maxx, nbiny, miny, maxy);
    }
        }//variable_2D
      }//variable_1D
    }//particle
  }//suffix
  }//trigger
}
void LoadOrAdd1DHistogram(TString key, TFile* file){
  TH1* hist = ReadXiPScatHistogram(file, key);
  if(hMap.find(key) == hMap.end()){
    hMap[key] = (TH1*)hist->Clone(key);
    hMap[key]->SetDirectory(nullptr);
  }
  else hMap[key]->Add(hist);
}
void LoadOrAdd2DHistogram(TString key, TFile* file){
  TH1* hist = ReadXiPScatHistogram(file, key);
  if(hMap2D.find(key) == hMap2D.end()){
    hMap2D[key] = (TH1*)hist->Clone(key);
    hMap2D[key]->SetDirectory(nullptr);
  }
  else hMap2D[key]->Add(hist);
}
void LoadHistograms(TFile* file){
  if(file == nullptr or file->IsZombie()) throw runtime_error("Cannot read map worker file.");
  TNamed* format = dynamic_cast<TNamed*>(file->Get("XiPScatMapFormat"));
  if(format == nullptr or TString(format->GetTitle()) != "XiPScat_StandardXi_20261006")
    throw runtime_error("Old Xi-p histogram worker. Regenerate maps with the revised macros.");
  cout<<"Loading Histograms from "<<file->GetName()<<endl;
  for(auto t:triggers){
  for(auto s:suffix){
    s = s + t;
    for(auto p:particle){
      for(int iv1=0;iv1<variable.size();++iv1){
        auto v1 = variable[iv1];
        TString key = AcceptanceHistTitle1D("G4", p, v1, s);
#if Debug
        cout<<"Loading "<<key<<endl;
#endif

  LoadOrAdd1DHistogram(key, file);
  for(int id = 0; id < ndist-1;++id){
    if(v1 == "DistT") continue;
    TString distkey = key + GetDistKey(id);
    LoadOrAdd1DHistogram(distkey, file);
  }
	for(int iv2=iv1+1;iv2<variable.size();++iv2){
	  auto v2 = variable[iv2];
	  key = AcceptanceHistTitle2D("G4", p, v1, v2, s);
#if Debug
	  cout<<"Loading "<<key<<endl;
#endif
	  LoadOrAdd2DHistogram(key, file);
    for(int id = 0; id < ndist-1;++id){
      if(v1 == "DistT" or v2 == "DistT") continue;
      TString distkey = key + GetDistKey(id);
      LoadOrAdd2DHistogram(distkey, file);
    }
        }//variable_2D
      }//variable_1D
    }//particle
  }//suffix
  }//trigger
}

void FillParticle(TVector3 P1, TVector3 P2, TVector3 V, TString pre, TString part, TString suff){
  double var; 
  TString key;
  for(int iv1=0;iv1<variable.size();++iv1){
    auto v1 = variable[iv1];
    key = AcceptanceHistTitle1D("G4", part, v1, suff);
    var = GetXiPVariable(P1, P2, V, part, v1);
  #if Debug
    cout<<"Filling "<<key<<endl;
  #endif

    hMap[key]->Fill(var);
    if(v1 != "DistT" ){
      TString distkey = key + GetDistKey(V);
      hMap[distkey]->Fill(var);
    }
    for(int iv2 = iv1+1;iv2<variable.size();++iv2){
      key = AcceptanceHistTitle2D("G4", part, v1, variable[iv2], suff);
      double var2 = GetXiPVariable(P1, P2, V, part, variable[iv2]);
#if Debug
    cout<<"Filling "<<key<<endl;
#endif
      hMap2D[key]->Fill(var, var2);
      if(v1 == "DistT" or variable[iv2] == "DistT") continue;
      TString distkey = key + GetDistKey(V);
#if Debug
    cout<<"Filling "<<distkey<<endl;
#endif
      hMap2D[distkey]->Fill(var, var2);
    }
  }
}


template<typename T>
void FillHist(T* XiEv){

  TVector3 G4KmMom = TVector3(XiEv->G4kmmom_x, XiEv->G4kmmom_y, XiEv->G4kmmom_z);
  TVector3 G4KpMom = TVector3(XiEv->G4kpmom_x, XiEv->G4kpmom_y, XiEv->G4kpmom_z);
  TVector3 G4KpVtx = TVector3(XiEv->G4kpvtx_x, XiEv->G4kpvtx_y, XiEv->G4kpvtx_z);

  TVector3 G4P = TVector3(XiEv->G4pmom_x, XiEv->G4pmom_y, XiEv->G4pmom_z);
  TVector3 G4Pi1 = TVector3(XiEv->G4pi1mom_x, XiEv->G4pi1mom_y, XiEv->G4pi1mom_z);
  TVector3 G4Pi2 = TVector3(XiEv->G4pi2mom_x, XiEv->G4pi2mom_y, XiEv->G4pi2mom_z);
  TVector3 G4L = TVector3(XiEv->G4lmom_x, XiEv->G4lmom_y, XiEv->G4lmom_z);
  TVector3 G4Xi = G4L + G4Pi2;
  TVector3 G4XiProd = TVector3(XiEv->G4ximom_x, XiEv->G4ximom_y, XiEv->G4ximom_z);
  TVector3 G4PScat = TVector3(XiEv->G4p_scatmom_x, XiEv->G4p_scatmom_y, XiEv->G4p_scatmom_z);
  
  TVector3 G4DecayVtxL = TVector3(XiEv->G4pvtx_x, XiEv->G4pvtx_y, XiEv->G4pvtx_z);
  TVector3 G4DecayVtxXi = TVector3(XiEv->G4pi2vtx_x, XiEv->G4pi2vtx_y, XiEv->G4pi2vtx_z);
  TVector3 G4ProdVtxXi = TVector3(XiEv->G4xivtx_x, XiEv->G4xivtx_y, XiEv->G4xivtx_z);
  if(G4Xi.Mag() < 1e-5) return;
  map<TString, vector<TVector3>> DataContainer = 
  {
    {"P",{G4P,G4Pi1,G4DecayVtxL}},
    {"Pi1",{G4Pi1,G4P,G4DecayVtxL}},
    {"Pi2",{G4Pi2,G4L,G4DecayVtxXi}},
    #if LVtxIsXi
    {"L",{G4L,G4Pi1,G4DecayVtxXi}},
    #else
    {"L",{G4L,G4Pi1,G4DecayVtxL}},
    #endif
    {"Xi",{G4Xi,G4Pi2,G4DecayVtxXi}},
    {"ScatP",{G4PScat,G4XiProd,G4ProdVtxXi}}
  };
  for(auto t:triggers){
  for(auto suff:suffix){
    bool check = SuffixCheck(suff, XiEv) and TrigCheck(t, XiEv);
    suff = suff + t;
    if(check){
      for(auto& [part, vec]:DataContainer){
        FillParticle(vec[0], vec[1], vec[2], "G4", part, suff);
      }
    }
  }
  }
}

void MakeEfficiencies(){
  for(auto t:triggers){
    for(auto s:suffix){
      if(s != "Gen") MakeEff(s+t, "Gen"+t);
    }
  }
  // The same five Xi factors, followed by one conditional proton factor.
  // XiPflag is an observed candidate sample, not an extra truth-matched factor.
  vector<pair<TString, TString>> EffConfs = {
    {"PPi1Tracked", "PTracked"}, {"PPi1Tracked", "Pi1Tracked"},
    {"AllTracked", "PTracked"}, {"AllTracked", "Pi1Tracked"},
    {"AllTracked", "PPi1Tracked"}, {"AllTracked", "Pi2Tracked"},
    {"GoodL", "PPi1Tracked"}, {"GoodLAndPi2Tracked", "GoodL"},
    {"GoodXi", "GoodLAndPi2Tracked"},
    {"GoodXiAndScatPGood", "GoodXi"},
    {"GoodXiAndScatPTracked", "GoodXi"}, // MC diagnostic only.
    {"XiAcpt", "LAcpt"}, {"GoodXi", "GoodL"}, {"GoodXi", "AllTracked"}
  };
  for(auto t:triggers){
    for(auto conf:EffConfs) MakeEff(conf.first+t, conf.second+t);
  }
}
void
SetBranches(TTree* tree){
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("nhHtof",1);
	tree->SetBranchStatus("G4kmmom_x",1);
	tree->SetBranchStatus("G4kmmom_y",1);
	tree->SetBranchStatus("G4kmmom_z",1);
	tree->SetBranchStatus("G4kpmom_x",1);
	tree->SetBranchStatus("G4kpmom_y",1);
	tree->SetBranchStatus("G4kpmom_z",1);
	tree->SetBranchStatus("G4kpvtx_x",1);
	tree->SetBranchStatus("G4kpvtx_y",1);
	tree->SetBranchStatus("G4kpvtx_z",1);
	
	tree->SetBranchStatus("G4pmom_x",1);
	tree->SetBranchStatus("G4pmom_y",1);
	tree->SetBranchStatus("G4pmom_z",1);
	tree->SetBranchStatus("G4pi1mom_x",1);
	tree->SetBranchStatus("G4pi1mom_y",1);
	tree->SetBranchStatus("G4pi1mom_z",1);
	tree->SetBranchStatus("G4pi2mom_x",1);
	tree->SetBranchStatus("G4pi2mom_y",1);
	tree->SetBranchStatus("G4pi2mom_z",1);
	tree->SetBranchStatus("G4lmom_x",1);
	tree->SetBranchStatus("G4lmom_y",1);
	tree->SetBranchStatus("G4lmom_z",1);
	tree->SetBranchStatus("G4ximom_x",1);
	tree->SetBranchStatus("G4ximom_y",1);
	tree->SetBranchStatus("G4ximom_z",1);

	tree->SetBranchStatus("G4pvtx_x",1);
	tree->SetBranchStatus("G4pvtx_y",1);
	tree->SetBranchStatus("G4pvtx_z",1);
	tree->SetBranchStatus("G4pi2vtx_x",1);
	tree->SetBranchStatus("G4pi2vtx_y",1);
	tree->SetBranchStatus("G4pi2vtx_z",1);
	tree->SetBranchStatus("G4xivtx_x",1);
	tree->SetBranchStatus("G4xivtx_y",1);
	tree->SetBranchStatus("G4xivtx_z",1);

  tree->SetBranchStatus("G4p_scatmom",1);
  tree->SetBranchStatus("G4p_scatmom_x",1);
  tree->SetBranchStatus("G4p_scatmom_y",1);
  tree->SetBranchStatus("G4p_scatmom_z",1);
  tree->SetBranchStatus("G4p_scatvtx_x",1);
  tree->SetBranchStatus("G4p_scatvtx_y",1);
  tree->SetBranchStatus("G4p_scatvtx_z",1);

	tree->SetBranchStatus("lgood",1);
	tree->SetBranchStatus("xigood",1);
	tree->SetBranchStatus("Xiflag",1);
	tree->SetBranchStatus("XiPflag",1);
  tree->SetBranchStatus("p_scat_tracked",1);
  tree->SetBranchStatus("p_scat_good",1);
	tree->SetBranchStatus("Pimflag",1);
	tree->SetBranchStatus("Emptyflag",1);

  tree->SetBranchStatus("p_tracked",1);
  tree->SetBranchStatus("pi1_tracked",1);
  tree->SetBranchStatus("pi2_tracked",1);

  tree->SetBranchStatus("KFXiDecaysMom_x",1);
  tree->SetBranchStatus("KFXiDecaysMom_y",1);
  tree->SetBranchStatus("KFXiDecaysMom_z",1);
  tree->SetBranchStatus("KFLambdaMom_x",1);
  tree->SetBranchStatus("KFLambdaMom_y",1);
  tree->SetBranchStatus("KFLambdaMom_z",1);
  tree->SetBranchStatus("KFXiMom_x",1);
  tree->SetBranchStatus("KFXiMom_y",1);
  tree->SetBranchStatus("KFXiMom_z",1);
  tree->SetBranchStatus("GFLambdaDecayVtx_x",1);
  tree->SetBranchStatus("GFLambdaDecayVtx_y",1);
  tree->SetBranchStatus("GFLambdaDecayVtx_z",1);
  tree->SetBranchStatus("GFXiDecayVtx_x",1);
  tree->SetBranchStatus("GFXiDecayVtx_y",1);
  tree->SetBranchStatus("GFXiDecayVtx_z",1);
  tree->SetBranchStatus("KFXiProductionVtx_x",1);
  tree->SetBranchStatus("KFXiProductionVtx_y",1);
  tree->SetBranchStatus("KFXiProductionVtx_z",1);
  tree->SetBranchStatus("KFXiProductionVtxMom_x",1);
  tree->SetBranchStatus("KFXiProductionVtxMom_y",1);
  tree->SetBranchStatus("KFXiProductionVtxMom_z",1);
  tree->SetBranchStatus("XiResidualsMom_x",1);
  tree->SetBranchStatus("XiResidualsMom_y",1);
  tree->SetBranchStatus("XiResidualsMom_z",1);

  tree->SetBranchStatus("XipVtxScatPMom_x",1);
  tree->SetBranchStatus("XipVtxScatPMom_y",1);
  tree->SetBranchStatus("XipVtxScatPMom_z",1);
  tree->SetBranchStatus("XipVtxDist",1);

}
bool FileCheck(TString filename){
  TFile* file = TFile::Open(filename);
  bool ok = file != nullptr and !file->IsZombie() and !file->TestBit(TFile::kRecovered);
  if(file != nullptr){ file->Close(); delete file; }
  return ok;
}
