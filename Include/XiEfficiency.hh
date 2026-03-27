#ifndef XiEfficiency_hh
#define XiEfficiency_hh
#include "Dir.hh"
vector<double> VertDistances = {0, 15, 30, 45, 60, 90, 120, 150, 250, 400};
vector<TString> particle = {"P", "Pi1", "Pi2" , "L", "Xi"};
vector<TString> variable = {"CosTh", "Mom", "Ph", "DistT", "CosPsi", "CosOpen"};
vector<TString> suffix = {"Gen",
   "LAcpt", "GoodL", "GoodXi",
    "XiAcpt", "PTracked", "Pi1Tracked", "Pi2Tracked",
    "PPi1Tracked","AllTracked", "GoodLAndPi2Tracked"
};
vector<TString> triggers = {"", "TrigB"};
struct CorrectionParameter{
    TString num, den;
    TString p_cor;
    vector<TString> var_cor;
};

vector<CorrectionParameter> CorrPars;
int ndist = VertDistances.size();
TString GetDistKey(TVector3 V){
  double dist = DistT(V);
  int id = ndist-2;
  for(int i=0;i<ndist-1;++i){
    if(dist >= VertDistances[i] && dist < VertDistances[i+1]){
      id = i; 
      break;
    }    
  }
  return Form("_d_%g_%g", VertDistances[id], VertDistances[id+1]);
} 
TString GetDistKey(int id){ 
  return Form("_d_%g_%g", VertDistances[id], VertDistances[id+1]);
}
TString EffSufix(TString Den){
  TString suf;
  suf = Den == "Gen" ? "Eff" : "Per" + Den + "Eff";
  return suf;
}
TString EffTitle1D(TString pre, TString part, TString var, TString Num, TString Den){
  TString title = AcceptanceHistTitle1D(pre, part, var, Num) + EffSufix(Den);
  return title;
}
TString EffTitle2D(TString pre, TString part, TString var1, TString var2, TString Num, TString Den){
  TString title = AcceptanceHistTitle2D(pre, part, var1, var2, Num) + EffSufix(Den);
  return title;
}
double GetWeight(TVector3 P1, TVector3 P2, TVector3 V, TString tgt, TString particle, TString conf, TString v1, TString v2 = "", TVector3 VL = TVector3(0,0,0)){
    TString MapKey;
    std::map<TString, TEfficiency*>::iterator Eff;
    if(v2 == ""){
        MapKey = AcceptanceHistTitle1D(tgt, particle, v1, conf);
        Eff = effMap.find(MapKey);
        if (Eff == effMap.end()){
            cout << "Efficiency not found for :" << endl;
            cout << MapKey << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        MapKey = AcceptanceHistTitle2D(tgt, particle, v1, v2, conf);
        #if PerDist
            MapKey = AcceptanceHistTitle2D(tgt, particle, v1, v2, conf) + GetDistKey(VL);
        #endif
        Eff = effMap.find(MapKey);
        if (Eff == effMap.end()){
            MapKey = AcceptanceHistTitle2D(tgt, particle, v2, v1, conf);
            #if PerDist
                MapKey = AcceptanceHistTitle2D(tgt, particle, v2, v1, conf) + GetDistKey(VL);
            #endif
            auto tmp = v2;
            v2 = v1;
            v1 = tmp;
            Eff = effMap.find(MapKey);
        }
        if (Eff == effMap.end()){
            cout << "Efficiency not found for :" << endl;
            cout << MapKey << endl;
            exit(EXIT_FAILURE);
        }
    } 
    double x = GetVariable(P1, P2, V, v1);
    double eff;
    if (v2 != ""){
        double y = GetVariable(P1, P2, V, v2);
        eff = Eff->second->GetEfficiency(Eff->second->FindFixBin(x, y));
    }
    else{
        eff = Eff->second->GetEfficiency(Eff->second->FindFixBin(x));
    }
    return 1. / eff;
}

void LoadEff1D(TFile* file, TString tgt, TString Num, TString Den){
    TString key;
    for(auto p:particle){
        for(int iv = 0; iv<variable.size();++iv){
            auto v = variable[iv];
            if(p == "P" or p == "Pi1" or p == "Pi2"){
                if(v == "CosOpen") continue;
            }
            key = EffTitle1D("G4", p, v, Num, Den);
            TString mapkey = EffTitle1D(tgt, p, v, Num, Den);
            effMap[mapkey] = (TEfficiency*)file->Get(key);
        }
    }
}
void LoadEff2D(TFile* file, TString tgt, TString Num, TString Den){
    TString key;
    for(auto p:particle){
        for(int iv1 = 0; iv1<variable.size();++iv1){
            auto v1 = variable[iv1];
            if(p == "P" or p == "Pi1" or p == "Pi2"){
                if(v1 == "CosOpen") continue;
            }
            for(int iv2 = iv1+1; iv2<variable.size();++iv2){
                auto v2 = variable[iv2];
                if(p == "P" or p == "Pi1" or p == "Pi2"){
                    if(v2 == "CosOpen") continue;
                }
                key = EffTitle2D("G4", p, v1, v2, Num, Den);
                TString mapkey = EffTitle2D(tgt, p, v1, v2, Num, Den);
                effMap[mapkey] = (TEfficiency*)file->Get(key);
                for(int id = 0; id < ndist-1;++id){
                  if(v1 == "DistT" or v2 == "DistT") continue;
                  key = EffTitle2D("G4", p, v1, v2, Num, Den) + GetDistKey(id);
                  mapkey = EffTitle2D(tgt, p, v1, v2, Num, Den) + GetDistKey(id);
                  effMap[mapkey] = (TEfficiency*)file->Get(key);
                }
            }
        }
    }
}
void LoadEff(TFile* file, TString tgt, TString Num, TString Den){
    LoadEff1D(file, tgt, Num, Den);
    LoadEff2D(file, tgt, Num, Den);
    for(auto eff:effMap){
        if(eff.second == nullptr){
            cout<<"Efficiency histogram "<<eff.first<<" not loaded !"<<endl;
        }
        else{
            //cout<<"Efficiency histogram "<<eff.first<<" loaded."<<endl;
        }
    }
}
void LoadEfficiencies(TFile* file, TString tgt){
    for(auto s:suffix){
      if(s == "Gen") continue;
      LoadEff(file, tgt, s, "Gen");
    }
    for(auto cp:CorrPars){
        TString num = cp.num;
        TString den = cp.den;
        LoadEff(file, tgt, num, den);
    }
}
void MakeEff1D(TString Num, TString Den){
  TString key0, key, key_eff;
  for(auto p:particle){
    for(int iv1 = 0;iv1<variable.size();++iv1){
      auto v1 = variable[iv1];
      key0 = AcceptanceHistTitle1D("G4", p, v1, Den);
      key = AcceptanceHistTitle1D("G4", p, v1, Num);
      #if Debug
        cout<<"Making Eff "<<key<<" over "<<key0<<endl;
      #endif
      auto hGen = hMap[key0];
      auto hSuff = hMap[key];
      key_eff = EffTitle1D("G4", p, v1, Num, Den);
      effMap[key_eff] = new TEfficiency(*hSuff, *hGen);
      effMap[key_eff]->SetNameTitle(key_eff, key_eff);
    }
  }
}
void MakeEff2D(TString Num, TString Den){
  TString key0, key, key_eff;
  for(auto p:particle){
    for(int iv1 = 0;iv1<variable.size();++iv1){
      auto v1 = variable[iv1];
      for(int iv2 = iv1+1;iv2<variable.size();++iv2){
        auto v2 = variable[iv2];
        key0 = AcceptanceHistTitle2D("G4", p, v1, v2, Den);
        key = AcceptanceHistTitle2D("G4", p, v1, v2, Num);
      #if Debug
        cout<<"Making Eff "<<key<<" over "<<key0<<endl;
      #endif
        auto hGen = hMap2D[key0];
        auto hSuff = hMap2D[key];
        key_eff = EffTitle2D("G4", p, v1, v2, Num, Den);
        effMap[key_eff] = new TEfficiency(*hSuff, *hGen);
        effMap[key_eff]->SetNameTitle(key_eff, key_eff);
        for(int id = 0; id < ndist-1;++id){
          if(v1 == "DistT" or v2 == "DistT") continue;
          auto hGen = hMap2D[key0 + GetDistKey(id)];
          auto hSuff = hMap2D[key + GetDistKey(id)];
      #if Debug
        cout<<"Making Eff "<<key + GetDistKey(id)<<" over "<<key0 + GetDistKey(id)<<endl;
      #endif
          key_eff = EffTitle2D("G4", p, v1, v2, Num, Den) + GetDistKey(id);
          effMap[key_eff] = new TEfficiency(*hSuff, *hGen);
          effMap[key_eff]->SetNameTitle(key_eff, key_eff);
        }
      }
    }
  }
}
void MakeEff(TString Num, TString Den){
  MakeEff1D(Num, Den);
  MakeEff2D(Num, Den);
}















#endif