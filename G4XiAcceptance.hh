#include "Include/g4genfitcarbon.h"
#include "Include/XiEfficiency.hh"
using namespace std;
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
double minDistT = 0, maxDistT = 350;

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


bool SuffixCheck(TString suff, g4genfitcarbon* XiEv){
  bool ret = false;
  if(suff.Contains("Gen")) ret = 1;
  if(suff.Contains("LAcpt")) ret = (XiEv->lgood or XiEv->Xiflag);
  if(suff.Contains("GoodL")) ret = XiEv->lgood;
  if(suff.Contains("XiAcpt")) ret = XiEv->Xiflag;
  if(suff.Contains("GoodXi")) ret = XiEv->xigood;
  if(suff.Contains("PTracked")) ret = XiEv->p_tracked;
  if(suff.Contains("Pi1Tracked")) ret = XiEv->pi1_tracked;
  if(suff.Contains("Pi2Tracked")) ret = XiEv->pi2_tracked;
  if(suff.Contains("PPi1Tracked")) ret = (XiEv->p_tracked and XiEv->pi1_tracked);
  if(suff.Contains("AllTracked")) ret = (XiEv->p_tracked and XiEv->pi1_tracked and XiEv->pi2_tracked);
  if(suff.Contains("GoodLAndPi2Tracked")) ret = (XiEv->lgood and XiEv->pi2_tracked);
  
  return ret;
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
    else if(part.Contains("L") or part.Contains("Xi") or part.Contains("P")){
      bin = nbinMom; min = minMomL; max = maxMomL;
    }
  }
  else if(var == "CosTh"){
    axis = "cos#theta";
    if(part.Contains("Pi")){
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

void FillParticle(TVector3 P1, TVector3 P2, TVector3 V, TString pre, TString part, TString suff){
  double var; 
  TString key;
  for(int iv1=0;iv1<variable.size();++iv1){
    auto v1 = variable[iv1];
    key = AcceptanceHistTitle1D("G4", part, v1, suff);
    var = GetVariable(P1, P2, V, v1);
  #if Debug
    cout<<"Filling "<<key<<endl;
  #endif

    hMap[key]->Fill(var);
    for(int iv2 = iv1+1;iv2<variable.size();++iv2){
      key = AcceptanceHistTitle2D("G4", part, v1, variable[iv2], suff);
      double var2 = GetVariable(P1, P2, V, variable[iv2]);
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
  TVector3 G4DecayVtxL = TVector3(XiEv->G4pvtx_x, XiEv->G4pvtx_y, XiEv->G4pvtx_z);
  TVector3 G4DecayVtxXi = TVector3(XiEv->G4pi2vtx_x, XiEv->G4pi2vtx_y, XiEv->G4pi2vtx_z);
  if(G4Xi.Mag() < 1e-5) return;
  map<TString, vector<TVector3>> DataContainer = 
  {
    {"P",{G4P,G4Pi1,G4DecayVtxL}},
    {"Pi1",{G4Pi1,G4P,G4DecayVtxL}},
    {"Pi2",{G4Pi2,G4L,G4DecayVtxXi}},
    {"L",{G4L,G4Pi1,G4DecayVtxL}},
    {"Xi",{G4Xi,G4Pi2,G4DecayVtxXi}},
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
      if(s == "Gen") continue;
      MakeEff(s, "Gen");
      MakeEff(s+t, "Gen");
      MakeEff(s+t, "Gen"+t);
    }
    } 
    //Correction Starategy:
    // 1. Gen -> P or Pi1Tracked // p or pi1 kinematics
    // 2. P or Pi1Tracked -> PPi1 Tracked // pi1 or p kinematics
    // 3-a. PPi1Tracked -> GoodL // L kinematics?
    // 4-a. GoodL -> GoodLAndPi2Tracked // pi2 kinematics
    // 5-a. GoodLAndPi2Tracked -> GoodXi // Xi kinematics?
    //
    // 3-b. PPi1Tracked -> AllTracked // pi2 kinematics
    // 4-b. AllTracked -> GoodXi // Xi kinematics? 
   
    //
    //
    //   Gen = Gen / PTracked (P) * PTracked / PPi1Tracked(Pi1) * PPi1Tracked / GoodL (L) * GoodL / GoodLAndPi2Tracked (Pi2) * GoodLAndPi2Tracked / GoodXi (Xi)
    //
    //
    vector<pair<TString, TString>> EffConfs = {// Numerator, Denominator for additional efficiency calculation.
      {"PPi1Tracked", "PTracked"},
      {"PPi1Tracked", "Pi1Tracked"},

      {"AllTracked", "PTracked"},
      {"AllTracked", "Pi1Tracked"},
      {"AllTracked", "PPi1Tracked"},
      {"AllTracked", "Pi2Tracked"},

      {"GoodL", "PPi1Tracked"},
      {"GoodLAndPi2Tracked", "GoodL"},
      {"GoodXi", "GoodLAndPi2Tracked"},

      {"XiAcpt", "LAcpt"},
      {"GoodXi", "GoodL"},
      {"GoodL", "PPi1Tracked"},
      {"GoodXi", "AllTracked"}
    };
		if(date == 260422){
			EffConfs = {
				{"GoodL","Gen"},
				{"GoodXi","GoodL"}
			};
		}
    for(auto t:triggers){
      for(auto& [num, den]:EffConfs){
        MakeEff(num+t, den+t);
      }
    }   }

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

	tree->SetBranchStatus("lgood",1);
	tree->SetBranchStatus("xigood",1);
	tree->SetBranchStatus("Xiflag",1);
	tree->SetBranchStatus("Pimflag",1);
	tree->SetBranchStatus("Emptyflag",1);

  tree->SetBranchStatus("p_tracked",1);
  tree->SetBranchStatus("pi1_tracked",1);
  tree->SetBranchStatus("pi2_tracked",1);

}
