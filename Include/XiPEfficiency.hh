#ifndef XiPEfficiency_hh
#define XiPEfficiency_hh
#include "XiEfficiency.hh"
int nbinMomP_XiP = 10; 
double minMomP_XiP = 0, maxMomP_XiP = 1.; 
double GeV = 1000.;
int MomBin_XiP(double p){ 
  int bin = -1; 
  for(int i=0;i<nbinMomP_XiP;++i){
    double p0 = minMomP_XiP + i*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP*GeV;
    double p1 = minMomP_XiP + (i+1)*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP*GeV;
    if(p >= p0 and p < p1){
      bin = i;
      break;
    }   
  }
  if(bin < 0) return nbinMomP_XiP -1; 
  return bin;
}
TString GetMomKey_XiP(TVector3 P){ 
  double mag = P.Mag();
  int bin = MomBin_XiP(mag);
  double p0 = minMomP_XiP + bin*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP * GeV;
  double p1 = minMomP_XiP + (bin+1)*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP * GeV;
  TString momkey = Form("_p_%d_%d", (int)p0, (int)p1);
  return momkey;
}
void LoadResidualAcceptance(vector<TFile*> AcptFiles){
    cout<<"Loading Residual Track Acceptance Map..."<<endl;
    TString t = "Carbon";
    TString key;
    int it= 0;
    auto file = AcptFiles[it];
    vector<TString> pars = {"CosTh", "Ph"}; 
    for(int ip1 = 0; ip1 < pars.size(); ++ip1){
        auto v1 = pars[ip1];
        for(int ip2 = ip1 + 1; ip2 < pars.size(); ++ip2){
            auto v2 = pars[ip2];
            key = "G4Xip" + v1 + v2 + "XipAcpt";
            TString MapKey = key + "perXi";
            cout<<"Loading "<<MapKey<<endl;
            auto Eff2D = (TEfficiency *)file->Get(MapKey);
            TString TargetKey = t + MapKey;
            Eff2D->SetName(TargetKey);
            effMap[TargetKey] = Eff2D;
            cout<<"Loaded "<<TargetKey <<endl;
            for(int imom = 0; imom < nbinMomP_XiP; ++imom){
                double p0 = minMomP_XiP + imom*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP*GeV;
                double p1 = minMomP_XiP + (imom+1)*(maxMomP_XiP - minMomP_XiP)/nbinMomP_XiP*GeV;
                TVector3 P(0,0,(p0+p1)/2.);
                TString momKey = GetMomKey_XiP(P);
                MapKey = key + momKey + "perXi";
                cout<<"Loading "<<MapKey<<endl;
                auto Eff2D = (TEfficiency *)file->Get(MapKey);
                TargetKey = t + MapKey;
                Eff2D->SetName(TargetKey);
                effMap[TargetKey] = Eff2D;
                cout<<"Loaded "<<TargetKey<<endl;
            }
        }
    }
}
void LoadResidualPurity(vector<TFile*> AcptFiles){
    TString t = "Carbon";
    TString key;
    int it= 0;
    auto file = AcptFiles[it];
    vector<TString> pars = {"Mom", "cth"}; 
    vector<TString> BH2MTs = {"","_Single_","_Multi_"};
    for(int ip1 = 0; ip1 < pars.size(); ++ip1){
        auto v1 = pars[ip1];
        for(auto mt: BH2MTs){
            key = t + "PurityGraph" + mt + v1 + "_P_Graph";
            cout<<"Loading "<<key<<endl;
            auto graph = (TGraphErrors*)file->Get(key);
            gMap[key + "Loaded"] = graph;
        }
    }
}
#endif