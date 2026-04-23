#include "G4XiAcceptance.hh"
#include "TGraphErrors.h"
#include <vector>
double weight_th = 50;
TString trig = "";
vector<TString> CheckLists = {
    "Gen"
};
//map<TString, vector<TString>> CorrectionConf {
    //{"AllTracked", {"GoodXi"}}
//};
#if CH2
#if date == 260422
		CorrPars = {
			{"GoodXi", "GoodL", "Pi2", {"CosTh","Ph"}},
			{"GoodL", "Gen", "Pi1", {"CosTh","Ph"}}
		};
#else
		CorrPars = {
		{"GoodXi", "GoodLAndPi2Tracked", "Xi", {"CosTh","CosOpen"}},
		{"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "CosPsi"}},
		{"GoodL", "PPi1Tracked", "L", {"DistT","CosOpen"}},
		{"PPi1Tracked", "PTracked", "Pi1", {"CosTh", "Ph"}},
		{"PTracked", "Gen", "P", {"CosTh", "Mom"}}
		};
#endif
#else
#if date == 260422
		CorrPars = {
			{"GoodXi", "GoodL", "Pi2", {"CosTh","Ph"}},
			{"GoodL", "Gen", "Pi1", {"CosTh","Ph"}}
		};
#else
		CorrPars = {
		{"GoodXi", "GoodLAndPi2Tracked", "Xi", {"CosTh","CosOpen"}},
		//{"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "Ph"}},
		{"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "Ph"}},
		{"GoodL", "PPi1Tracked", "L", {"CosTh","CosOpen"}},
		{"PPi1Tracked", "PTracked", "Pi1", {"CosTh", "Ph"}},
		{"PTracked", "Gen", "P", {"Mom", "DistT"}}
		};
#endif
#endif
TString Correction(TString part, TString num, TString den, TString v1, TString v2 = ""){
    TString suf;
    if(v2 == "") suf = "Cor" + part + num + den + v1;
    else suf = "Cor" + part + num + den + v1 + v2;
    return suf;
}
TString CorrectionHists(TString pre, TString part, TString var, TString suff, TString part_cor, TString num, TString den, TString v1, TString v2 = ""){
    return AcceptanceHistTitle1D(pre, part, var, suff) + Correction(part_cor,num, den, v1, v2);
}
map<TString, int> colorMap = {
    {"Gen", kBlack},
    {"PTracked", kRed},
    {"PPi1Tracked", kBlue},
    {"GoodL", kGreen+2},
    {"GoodLAndPi2Tracked", kMagenta},
    {"GoodXi", kCyan+2}
};
void InitializeTriggerCondtions(){
#if TrigB
    cout<<"Applying trigger bias correction..."<<endl;
    for(int ic=0;ic<CheckLists.size();++ic){
        CheckLists[ic] = CheckLists[ic] + "TrigB";
    }
    cout<<"Updated CheckLists: "<<endl;
    for(auto &c:CorrPars){
        cout<<"Updating CorrectionConf for "<<c.num<<" over "<<c.den<<endl;
        c.num = c.num + "TrigB"; 
        c.den = c.den + "TrigB";
    }
    for(auto c:CorrPars){
        cout<<"Updated CorrectionConf for "<<c.num<<" over "<<c.den<<endl;
    }
    cout<<"Updated CorrectionConf: "<<endl;
    map<TString, int> colorMap_temp;
    for(auto c: colorMap){
        TString key = c.first + "TrigB";
        colorMap_temp[key] = colorMap[c.first];
    }
    for(auto c: colorMap_temp){
        colorMap[c.first] = c.second;
    }
    trig = "TrigB";
    cout<<"Updated colorMap: "<<endl;
#endif
}
void InitializeCorrectionHistograms(TString tgt){
    for(auto cp:CorrPars){
        TString num = cp.num;
        TString den = cp.den;
        bool num_acpt = 0, den_acpt = 0;
        for(auto chk:CheckLists){
            if(chk == num) num_acpt = 1;
            if(chk == den) den_acpt = 1;
        }
        if(!num_acpt) CheckLists.push_back(num);
        if(!den_acpt) CheckLists.push_back(den);
    }
    for(auto chk:CheckLists){
        cout<<"CheckList: "<<chk<<endl;
    }
    int nbinx;
    double minx, maxx;
    TString Xtitle;
    for(auto chk: CheckLists){
        for(auto p:particle){
            for(int iv=0;iv<variable.size();++iv){
                auto v = variable[iv];
                TString key = AcceptanceHistTitle1D(tgt, p, v, chk);
                Xtitle = SetAxis(p, v, nbinx, minx, maxx);
                hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
                hMap[key] -> SetLineColor(colorMap[chk]);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            hMap[key] -> SetLineColor(colorMap[chk]);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1, v2);
                hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
                hMap[key] -> SetLineColor(colorMap[chk]);
            }
        }
    }
            }
        }
    }
    for(auto p:particle){
        for(int iv=0;iv<variable.size();++iv){
            auto v = variable[iv];
            TString key = AcceptanceHistTitle1D(tgt, p, v, "XiAcpt"+trig);
            Xtitle = SetAxis(p, v, nbinx, minx, maxx);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiAcptCor"+trig);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            key = AcceptanceHistTitle1D(tgt, p, v, "GoodXiCor"+trig);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        }
    }
};
class Event{
    private:
        g4genfitcarbon* gf;
        TVector3 P, Pi1, Pi2, L, Xi;
        TVector3 VL,VXi;
        map<TString, vector<TVector3>> DC;// DataContainer
        TString tgt;
    public:
        Event(g4genfitcarbon* g, TString t): gf(g), tgt(t){
            P.SetXYZ(gf->G4pmom_x, gf->G4pmom_y, gf->G4pmom_z);
            Pi1.SetXYZ(gf->G4pi1mom_x, gf->G4pi1mom_y, gf->G4pi1mom_z);
            Pi2.SetXYZ(gf->G4pi2mom_x, gf->G4pi2mom_y, gf->G4pi2mom_z);
            L.SetXYZ(gf->G4lmom_x, gf->G4lmom_y, gf->G4lmom_z);
            Xi = L + Pi2;
            VL.SetXYZ(gf->G4pvtx_x, gf->G4pvtx_y, gf->G4pvtx_z);
            VXi.SetXYZ(gf->G4pi2vtx_x, gf->G4pi2vtx_y, gf->G4pi2vtx_z);
            DC = {
                {"P",{P,Pi1,VL}},
                {"Pi1",{Pi1,P,VL}},
                {"Pi2",{Pi2,L,VXi}},
                {"L",{L,Pi1,VL}},
                {"Xi",{Xi,Pi2,VXi}},
            };
        }
        map<TString, vector<TVector3>> GetDataContainer(){
            return DC;
        }
        bool Good(){
            if(Xi.Mag() <1e-5) return false;
            #if TrigB
            if(gf->nhHtof < 2) return false;
            #endif
            return true;
        }
        double GetWeight(TString part, TString num, TString den, TString v1, TString v2 = ""){
            TString key;
            double var1,var2;
            var1 = GetVariable(DC[part][0], DC[part][1], DC[part][2], v1);
            double w;
            if(v2 != ""){
                var2 = GetVariable(DC[part][0], DC[part][1], DC[part][2], v2);
                key = EffTitle2D(tgt, part, v1, v2, num, den);
#if PerDist
                if(v1 != "DistT" and v2 != "DistT") key += GetDistKey(DC[part][2]);
#endif
                if(effMap[key] == nullptr){
                    cout<<"Efficiency histogram "<<key<<" not found!"<<endl;
                    return 1.;
                }
                w = effMap[key]->GetEfficiency(effMap[key]->FindFixBin(var1, var2));
                w = 1./w;
            }
            else{
                key = EffTitle1D(tgt, part, v1, num, den);
                if(effMap[key] == nullptr){
                    cout<<"Efficiency histogram "<<key<<" not found!"<<endl;
                    return 1.;
                }
                w = effMap[key]->GetEfficiency(effMap[key]->FindFixBin(var1));
                w = 1./w;
            }
            return w;
        }
};

void FillHistograms(g4genfitcarbon* gf, TString tgt){
    Event event(gf, tgt);
    if(!event.Good()) return;
    auto DC = event.GetDataContainer();
    TString key;
    double w;
    double weight = 1;
    for(auto cp:CorrPars){
        TString num = cp.num;
        TString den = cp.den;
        TString pc = cp.p_cor;
        vector<TString> var_cor = cp.var_cor;
        TString v1 = var_cor[0];
        TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
        weight *= event.GetWeight(pc, num, den, v1, v2);
    }
    if(weight > weight_th) return;
    for(auto p:particle){
        for(int iv=0;iv<variable.size();++iv){
            auto v = variable[iv];
            for(auto chk: CheckLists){
                if(!SuffixCheck(chk, gf)) continue;
                if(!TrigCheck(trig, gf)) continue;
                key = AcceptanceHistTitle1D(tgt, p, v, chk);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1], DC[p][2], v));
        for(auto cp:CorrPars){
            TString num = cp.num;
            if(num != chk) continue;
            TString den = cp.den;
            TString pc = cp.p_cor;
            for(int iv1 = 0; iv1 < variable.size();++iv1){
                auto v1 = variable[iv1];
                if(pc == "P" or pc == "Pi1" or pc == "Pi2"){
                    if(v1 == "CosOpen") continue;
                }
                w = event.GetWeight(pc, num, den, v1);
                key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1], DC[p][2], v), w);
                for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                    auto v2 = variable[iv2];
                    if(pc == "P" or pc == "Pi1" or pc == "Pi2"){
                        if(v2 == "CosOpen") continue;
                    }
                    w = event.GetWeight(pc, num, den, v1, v2);
                    key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1, v2);
                    hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1], DC[p][2], v), w);
                }
            }
        }
            }
            if(SuffixCheck("XiAcpt", gf)){
                key = AcceptanceHistTitle1D(tgt, p, v, "XiAcpt" + trig);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v));
                key = AcceptanceHistTitle1D(tgt, p, v, "XiAcptCor" + trig);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v), weight);
            }
            if(SuffixCheck("GoodXi", gf)){
                key = AcceptanceHistTitle1D(tgt, p, v, "GoodXiCor" + trig);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v), weight);
            }
        }
    }
}
void NormalizeHistograms(TString tgt){
    TString key;
    for(auto p:particle){
        for(int iv=0;iv<variable.size();++iv){
            double maxi = 0;
            auto v = variable[iv];
            key = AcceptanceHistTitle1D(tgt, p, v, "Gen" + trig);
            maxi = hMap[key]->GetMaximum();
            key = AcceptanceHistTitle1D(tgt, p, v, "XiAcpt"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kBlue);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiAcptCor"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kRed);
            hMap[key]->SetMarkerColor(kRed);
            hMap[key]->SetMarkerStyle(20);
            key = AcceptanceHistTitle1D(tgt, p, v, "GoodXiCor"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kGreen+2);
            hMap[key]->SetMarkerColor(kGreen+2);
            hMap[key]->SetMarkerStyle(20);
            for(auto chk: CheckLists){
                key = AcceptanceHistTitle1D(tgt, p, v, chk);
                hMap[key]->Scale(1./maxi);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1);
            hMap[key]->Scale(1./maxi);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1, v2);
                hMap[key]->Scale(1./maxi);
            }
        }
    }
            }
        }
    }
}
TGraphErrors* MakeDivisionGraph(TH1* h1, TH1* h2){
    TGraphErrors* g = new TGraphErrors();
    for(int i = 1; i <= h1->GetNbinsX();++i){
        double x = h1->GetBinCenter(i);
        double y1 = h1->GetBinContent(i);
        double y2 = h2->GetBinContent(i);
        double err1 = h1->GetBinError(i);
        double err2 = h2->GetBinError(i);
        if(y2 == 0) continue;
        double y = y1 / y2;
        double err = y * sqrt(pow(err1/y1, 2) + pow(err2/y2, 2));
        g->SetPoint(g->GetN(), x, y);
        g->SetPointError(g->GetN()-1, 0, err);
    }
    auto color = h1->GetLineColor();
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->SetLineWidth(1);
    return g;
}

double CalcChi2(TH1* Gen, TH1* Cor, double trun = 0.2){
    double chi2 = 0;
    int ent_eff = 0;
    for(int i = 1; i <= Gen->GetNbinsX();++i){
        double gen = Gen->GetBinContent(i);
        if(gen == 0) continue;
        ent_eff ++;
    }
    double den = Gen->GetEntries() / ent_eff;
    vector<double> dels;
    for(int i = 1; i <= Gen->GetNbinsX();++i){
        double gen = Gen->GetBinContent(i);
        double gen_err = Gen->GetBinError(i);
        double cor = Cor->GetBinContent(i);
        double cor_err = Cor->GetBinError(i);
        double err = hypot(gen_err, cor_err);
        dels.push_back((gen - cor));
        //dels.push_back((gen - cor)/err);
    }
    std::sort(dels.begin(), dels.end());
    int n = dels.size();
    int nc = 0;
    for(int i = 0; i < n * (1-trun);++i){
        chi2 += pow(dels[i], 2);
        nc ++;
    }
    return chi2/nc;
}

double CalcChi2(TString tgt, TString num, TString den, TString part_cor, TString v1, TString v2 = ""){
    double chi2 = 0;
    TString key_num, key_den, key_cor;
    for(auto p:particle){
        for(auto v:variable){
            if(p == "P" or p == "Pi1" or p == "Pi2"){
                if(v == "CosOpen") continue;
            }
            key_den = AcceptanceHistTitle1D(tgt, p, v, den);
            key_cor = CorrectionHists(tgt, p, v, num, part_cor, num, den, v1, v2);
            chi2 += CalcChi2(hMap[key_den], hMap[key_cor]);
        }
    }
    return chi2;
}



struct Chi2Map{
    TString num, den;
    TString part_cor;
    map<vector<TString>, double> cor_chi2;
    //{{{v1, v2}, chi2}}
};
vector<Chi2Map> chi2Maps;
vector<Chi2Map> chi2Maps_best;
// map< {Particle, Particle_c,  num, den} , map<Corrections, chi2>>
// Corrections = {var1, var2}
void MakeChi2Map(TString tgt){
    TString key_num, key_den, key_cor;
    TString num, den;
    TString part;
    TString part_cor;
    TString vc_1, vc_2;
    for(int ic = 0; ic < CorrPars.size(); ++ic){
        map<vector<TString>, double> cor_chi2;
        auto cor = CorrPars[ic];
        num = cor.num;
        den = cor.den;
        part_cor = cor.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            if(part_cor == "P" or part_cor == "Pi1" or part_cor == "Pi2"){
                if(v1 == "CosOpen") continue;
            }
            cor_chi2[{v1,""}] = CalcChi2(tgt, num, den, part_cor, v1);
            cout<<"chi2 for "<<part_cor<<" "<<num<<"/"<<den<<" with var "<<v1<<" is "<<cor_chi2[{v1,""}]<<endl;
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                if(part_cor == "P" or part_cor == "Pi1" or part_cor == "Pi2"){
                    if(v2 == "CosOpen") continue;
                }
                cor_chi2[{v1, v2}] = CalcChi2(tgt, num, den, part_cor, v1, v2);
                cout<<"chi2 for "<<part_cor<<" "<<num<<"/"<<den<<" with var "<<v1<<" and "<<v2<<" is "<<cor_chi2[{v1, v2}]<<endl;
            }
        }
        Chi2Map chi2map = {num, den, part_cor, cor_chi2};
        chi2Maps.push_back(chi2map);
        double chi2_min = -1;
        vector<TString> best_var;
        map<vector<TString>, double> cor_chi2_best;
        for(auto& [var, chi2]: cor_chi2){
            if(chi2_min < 0) chi2_min = chi2;
            if(chi2 < chi2_min){
                chi2_min = chi2;
                best_var = var;
            }
        }
        cor_chi2_best[best_var] = chi2_min;
        Chi2Map chi2map_best = {num, den, part_cor, cor_chi2_best};
        chi2Maps_best.push_back(chi2map_best);
        cout<<"Best chi2 for "<<part_cor<<" "<<num<<"/"<<den<<" is "<<chi2_min<<" with var "<<best_var[0]<<" "<<best_var[1]<<endl;
    }
}
