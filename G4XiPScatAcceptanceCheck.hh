#include "G4XiPScatAcceptance.hh"
#include "TGraphErrors.h"
#include <vector>
double weight_th = 100;
int over_cnt = 0;
TString trig = "";
vector<TString> CheckLists = {
    "Gen"
};
vector<TString> EventVars = {
    "ctau_L", "ctau_Xi", "CThKK", "SqrtS",
    "PolCTh", "PolCPhA", "PolCPhB", "PolCPhC",
    "PolCThX", "PolCThZ", "MXiPScat", "CThXiPScat", "Collinearity"
};
TString Correction(TString part, TString num, TString den, TString v1, TString v2 = ""){
    TString suf;
    if(v2 == "") suf = "Cor" + part + num + den + v1;
    else suf = "Cor" + part + num + den + v1 + v2;
    return suf;
}
TString EventTitle(TString tgt, TString var, TString suff){
    TString title = tgt + var + suff;
    return title;
}
TString CorrectedEventTitle(TString tgt, TString var, TString suff, TString part_cor, TString num, TString den, TString v1, TString v2 = ""){
    return EventTitle(tgt, var, suff) + Correction(part_cor, num, den, v1, v2); 
}
TString SetEventVarAxis(TString var, int &nbinx, double &minx, double &maxx){
    TString axis = "";
    if(var == "ctau_L"){
        axis = "c#tau of #Lambda [mm]";
        nbinx = 20; minx = 0; maxx = 320;
    }
    else if(var == "ctau_Xi"){
        axis = "c#tau of #Xi [mm]";
        nbinx = 20; minx = 0; maxx = 200;
    }
    else if(var == "CThKK"){
        axis = "cos#theta*_{K^{-}K^{+}}";
        nbinx = 20; minx = 0.8; maxx = 1;
    }
    else if(var == "SqrtS"){
        axis = "#sqrt{s} [GeV]";
        nbinx = 30; minx = 1.8; maxx = 2.4;
    }
    else if(var == "PolCTh"){
        axis = "cos#theta";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "PolCPhA"){
        axis = "cos#phi_{#alpha}";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "PolCPhB"){
        axis = "cos#phi_{#beta}";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "PolCPhC"){
        axis = "cos#phi_{#gamma}";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "PolCThX"){
        axis = "cos#theta_{X}";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "PolCThZ"){
        axis = "cos#theta_{Z}";
        nbinx = 10; minx = -1; maxx = 1;
    }
    else if(var == "MXiPScat"){
        axis = "M_{#XiP} [GeV]";
        nbinx = 20; minx = 1.3; maxx = 1.6;
    }
    else if(var == "CThXiPScat"){
        axis = "cos#theta_{#XiP}";
        nbinx = 20; minx = -1; maxx = 1;
    }
    else if(var == "Collinearity"){
        axis = "Collinearity";
        nbinx = 20; minx = 0.9; maxx = 1;
    }
    return axis;
}
//map<TString, vector<TString>> CorrectionConf {
    //{"AllTracked", {"GoodXi"}}
//};
#if CH2
    #if date == 260422
		CorrPars = {
			{"GoodXi", "GoodL", "Pi2", {"CosTh","Ph"}},
			{"GoodL", "Gen", "Pi1", {"CosTh","Ph"}}
		};
    #elif date >= 260817
		CorrPars = {
		{"GoodXi", "GoodLAndPi2Tracked", "Xi", {"Mom","CosOpen"}},
		{"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "CosPsi"}},
		{"GoodL", "PPi1Tracked", "L", {"CosTh","CosOpen"}},
		{"PPi1Tracked", "PTracked", "Pi1", {"Mom", "CosPsi"}},
		{"PTracked", "Gen", "P", {"CosTh", "Mom"}}
		};
    #else
		CorrPars = {
		{"GoodXi", "GoodLAndPi2Tracked", "Xi", {"CosOpen"}},
		{"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "CosPsi"}},
		{"GoodL", "PPi1Tracked", "L", {"CosTh"}},
		{"PPi1Tracked", "PTracked", "Pi1", {"CosTh", "CosPsi"}},
		{"PTracked", "Gen", "P", {"CosTh", "Mom"}}
		};
    #endif
    #if SecondaryCorr
        #if date < 260806
            double p0 = 0.989681, p1 = -0.0194191;
        #elif date >= 260817
            double p0 = 0.965351, p1 = -0.0226586;
        #else
            double p0 = 0.992069, p1 = -0.0135689;
        #endif
        TGraph* gr_PhA = new TGraph();
        gr_PhA->SetPoint(0, -1, p0 - p1);
        gr_PhA->SetPoint(1, 0, p0 );
        gr_PhA->SetPoint(2, 1, p0 + p1);
        SecondaryCorrection["PolCPhA"] = gr_PhA;
    #endif
#else// Carbon
    CorrPars = {
        {"GoodXiAndScatPGood", "GoodXi", "ScatP", {"CosTh", "Mom"}},
        {"GoodXi", "GoodLAndPi2Tracked", "Xi", {"CosOpen"}},
        {"GoodLAndPi2Tracked", "GoodL", "Pi2", {"CosTh","Ph"}},
        {"GoodL", "PPi1Tracked", "L", {"CosTh","CosOpen"}},
        {"PPi1Tracked", "PTracked", "Pi1", {"CosTh", "Ph"}},
        {"PTracked", "Gen", "P", {"CosTh", "Mom"}}
    };
    #if SecondaryCorr
        #if date < 260806
            double p0 = 0.989681, p1 = -0.0194191;
        #elif date >= 260817
            double p0 = 0.987175, p1 = -0.025111;
        #else
            double p0 = 0.992069, p1 = -0.0135689;
        #endif
        TGraph* gr_PhA = new TGraph();
        gr_PhA->SetPoint(0, -1, p0 - p1);
        gr_PhA->SetPoint(1, 0, p0 );
        gr_PhA->SetPoint(2, 1, p0 + p1);
        SecondaryCorrection["PolCPhA"] = gr_PhA;
    #endif
#endif
TString CorrectionHists(TString pre, TString part, TString var, TString suff, TString part_cor, TString num, TString den, TString v1, TString v2 = ""){
    return AcceptanceHistTitle1D(pre, part, var, suff) + Correction(part_cor,num, den, v1, v2);
}
map<TString, int> colorMap = {
    {"Gen", kBlack},
    {"PTracked", kRed},
    {"PPi1Tracked", kBlue},
    {"GoodL", kGreen+2},
    {"GoodLAndPi2Tracked", kMagenta},
    {"GoodXi", kCyan+2},
    {"GoodXiAndScatPTracked", kOrange+7},
    {"GoodXiAndScatPGood", kViolet+7}
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
    TH1D* h_weight = new TH1D("h_weight", "h_weight", 1000, 0, 200);
    hMap["h_weight"] = h_weight;

    for(auto cp:CorrPars){
        TString num = cp.num;
        TString den = cp.den;
        bool num_acpt = 0, den_acpt = 0;
        cout<<"Checking if "<<num<<" and "<<den<<" are in CheckLists..."<<endl;
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
            }//iv2
        }//iv1
    }//cp
            }//iv
        }//particle
        for(auto ev: EventVars){
            TString key = EventTitle(tgt, ev, chk);
            Xtitle = SetEventVarAxis(ev, nbinx, minx, maxx);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            hMap[key] -> SetLineColor(colorMap[chk]);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            hMap[key] -> SetLineColor(colorMap[chk]);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1, v2);
                hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
                hMap[key] -> SetLineColor(colorMap[chk]);
            }//iv2
        }//iv1
    }//cp
        }//ev
    }//CheckLists
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
            key = AcceptanceHistTitle1D(tgt, p, v, "XiRecAcpt"+trig);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiRecAcptCor"+trig);
            hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        }
    }
    for(auto ev: EventVars){
        TString key = EventTitle(tgt, ev, "XiAcpt"+trig);
        Xtitle = SetEventVarAxis(ev, nbinx, minx, maxx);
        hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        key = EventTitle(tgt, ev, "XiAcptCor"+trig);
        hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        key = EventTitle(tgt, ev, "GoodXiCor"+trig);
        hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        key = EventTitle(tgt, ev, "XiRecAcpt"+trig);
        hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
        key = EventTitle(tgt, ev, "XiRecAcptCor"+trig);
        hMap[key] = new TH1D(key, key + ";" + Xtitle, nbinx, minx, maxx);
    }
};
void LoadOrAddHistogram(TString key, TFile* file, TString chk = ""){
    if(hMap.find(key) == hMap.end()){
        hMap[key] = (TH1*)file->Get(key);
        if(chk == ""){
            cout<<"Warning: No CheckList specified for histogram "<<key<<"."<<endl;
        }
        else{
            hMap[key] -> SetLineColor(colorMap[chk]);
        }
    }
    else{
        hMap[key] ->Add((TH1*)file->Get(key));
    }
}
void LoadCorrectionHistograms(TFile* file,TString tgt){
    cout<<"Loading histograms from file: "<<file->GetName()<<endl;
    LoadOrAddHistogram("h_weight", file);
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
    for(auto chk: CheckLists){
        for(auto p:particle){
            for(int iv=0;iv<variable.size();++iv){
                auto v = variable[iv];
                TString key = AcceptanceHistTitle1D(tgt, p, v, chk);
                LoadOrAddHistogram(key, file, chk);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1);
            LoadOrAddHistogram(key, file, chk);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1, v2);
                LoadOrAddHistogram(key, file, chk);
            }//iv2
        }//iv1
    }//cp
            }//iv
        }//particle
        for(auto ev: EventVars){
            TString key = EventTitle(tgt, ev, chk);
            LoadOrAddHistogram(key, file, chk);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1);
            LoadOrAddHistogram(key, file, chk);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1, v2);
                LoadOrAddHistogram(key, file, chk);
            }//iv2
        }//iv1
    }//cp
        }//ev
    }//CheckLists
    vector<TString> XiAcptHists = {
        "XiAcpt"+trig, "XiAcptCor"+trig, "GoodXiCor"+trig,
        "XiRecAcpt"+trig, "XiRecAcptCor"+trig
    };
    for(auto xiconf: XiAcptHists){
        for(auto p:particle){
            for(int iv=0;iv<variable.size();++iv){
                auto v = variable[iv];
                TString key = AcceptanceHistTitle1D(tgt, p, v, xiconf);
                LoadOrAddHistogram(key, file);
            }
        }
        for(auto ev: EventVars){
            TString key = EventTitle(tgt, ev, xiconf);
            LoadOrAddHistogram(key, file);
        }
    }
}
class Event{
    private:
        g4genfitcarbon* gf;
        TVector3 Km, Kp;
        TVector3 P, Pi1, Pi2, L, Xi, XiProd, PScat;
        TVector3 VL,VXi,VXiProd;
        TVector3 PRec, Pi1Rec, Pi2Rec, LRec, XiRec, XiProdRec, PScatRec;
        TVector3 VLRec,VXiRec,VXiProdRec;
        map<TString, vector<TVector3>> DC;// DataContainer
        map<TString, vector<TVector3>> DCRec;// DataContainer
        TString tgt;
    public:
        Event(g4genfitcarbon* g, TString t): gf(g), tgt(t){
            P.SetXYZ(gf->G4pmom_x, gf->G4pmom_y, gf->G4pmom_z);
            Pi1.SetXYZ(gf->G4pi1mom_x, gf->G4pi1mom_y, gf->G4pi1mom_z);
            Pi2.SetXYZ(gf->G4pi2mom_x, gf->G4pi2mom_y, gf->G4pi2mom_z);
            L.SetXYZ(gf->G4lmom_x, gf->G4lmom_y, gf->G4lmom_z);
            Xi = L + Pi2;
            XiProd.SetXYZ(gf->G4ximom_x, gf->G4ximom_y, gf->G4ximom_z);
            PScat.SetXYZ(gf->G4p_scatmom_x, gf->G4p_scatmom_y, gf->G4p_scatmom_z); 
            VL.SetXYZ(gf->G4pvtx_x, gf->G4pvtx_y, gf->G4pvtx_z);
            VXi.SetXYZ(gf->G4pi2vtx_x, gf->G4pi2vtx_y, gf->G4pi2vtx_z);
            VXiProd.SetXYZ(gf->G4xivtx_x, gf->G4xivtx_y, gf->G4xivtx_z);
            Km.SetXYZ(gf->G4kmmom_x, gf->G4kmmom_y, gf->G4kmmom_z);
            Kp.SetXYZ(gf->G4kpmom_x, gf->G4kpmom_y, gf->G4kpmom_z);
            DC = {
                {"P",{P,Pi1,VL}},
                {"Pi1",{Pi1,P,VL}},
                {"Pi2",{Pi2,L,VXi}},
                #if LVtxIsXi
                {"L",{L,Pi1,VXi}},
                #else
                {"L",{L,Pi1,VL}},
                #endif
                {"Xi",{Xi,Pi2,VXi}},
                {"ScatP",{PScat,XiProd,VXiProd}}
            };
            if(gf->Xiflag){
                PRec.SetXYZ(gf->KFXiDecaysMom_x->at(0), gf->KFXiDecaysMom_y->at(0), gf->KFXiDecaysMom_z->at(0));
                Pi1Rec.SetXYZ(gf->KFXiDecaysMom_x->at(1), gf->KFXiDecaysMom_y->at(1), gf->KFXiDecaysMom_z->at(1));
                Pi2Rec.SetXYZ(gf->KFXiDecaysMom_x->at(2), gf->KFXiDecaysMom_y->at(2), gf->KFXiDecaysMom_z->at(2));
                LRec.SetXYZ(gf->KFLambdaMom_x, gf->KFLambdaMom_y, gf->KFLambdaMom_z);
                XiRec.SetXYZ(gf->KFXiMom_x, gf->KFXiMom_y, gf->KFXiMom_z);
                XiProdRec.SetXYZ(gf->KFXiProductionVtxMom_x, gf->KFXiProductionVtxMom_y, gf->KFXiProductionVtxMom_z);
                VLRec.SetXYZ(gf->GFLambdaDecayVtx_x, gf->GFLambdaDecayVtx_y, gf->GFLambdaDecayVtx_z);
                VXiRec.SetXYZ(gf->GFXiDecayVtx_x, gf->GFXiDecayVtx_y, gf->GFXiDecayVtx_z);
                VXiProdRec.SetXYZ(gf->KFXiProductionVtx_x, gf->KFXiProductionVtx_y, gf->KFXiProductionVtx_z);
                if(gf->XiPflag){
                    PScatRec.SetXYZ(gf->XiResidualsMom_x->at(0), gf->XiResidualsMom_y->at(0), gf->XiResidualsMom_z->at(0));
                }
            }
            else{
                PRec.SetXYZ(0,0,0);
                Pi1Rec.SetXYZ(0,0,0);
                Pi2Rec.SetXYZ(0,0,0);
                LRec.SetXYZ(0,0,0);
                XiRec.SetXYZ(0,0,0);
                XiProdRec.SetXYZ(0,0,0);
                PScatRec.SetXYZ(0,0,0);
                VLRec.SetXYZ(0,0,0);
                VXiRec.SetXYZ(0,0,0);
                VXiProdRec.SetXYZ(0,0,0);
            }
            DCRec = {
                {"P",{PRec,Pi1Rec,VLRec}},
                {"Pi1",{Pi1Rec,PRec,VLRec}},
                {"Pi2",{Pi2Rec,LRec,VXiRec}},
                #if LVtxIsXi
                {"L",{LRec,Pi1Rec,VXiRec}},
                #else
                {"L",{LRec,Pi1Rec,VLRec}},
                #endif
                {"Xi",{XiRec,Pi2Rec,VXiRec}},
                {"ScatP",{PScatRec,XiProdRec,VXiProdRec}}
            };
        }
        map<TString, vector<TVector3>> GetDataContainer(){
            return DC;
        }
        map<TString, vector<TVector3>> GetDataContainerRec(){
            return DCRec;
        }
        bool Good(){
            if(Xi.Mag() <1e-5) return false;
            #if TrigB
            if(gf->nhHtof < 2) return false;
            #endif
            return true;
        }
        double GetWeightDC(TString part, TString num, TString den, TString v1, TString v2 = "", map<TString, vector<TVector3>> DC_ = map<TString, vector<TVector3>>()){
            TString key;
            double var1,var2;
            var1 = GetVariable(DC_[part][0], DC_[part][1], DC_[part][2], v1);
            double w;
            if(v2 != ""){
                var2 = GetVariable(DC_[part][0], DC_[part][1], DC_[part][2], v2);
                key = EffTitle2D(tgt, part, v1, v2, num, den);
#if PerDist
                if(v1 != "DistT" and v2 != "DistT" and part != "Xi" and part != "Pi2") key += GetDistKey(DC_[part][2]);
#endif
#if PerXiDist
                if(v1 != "DistT" and v2 != "DistT" and (part == "Xi" or part == "Pi2")) key += GetDistKey(DC_[part][2]);
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
#if PerDist
                if(v1 != "DistT" and part != "Xi" and part != "Pi2") key += GetDistKey(DC_[part][2]);
#endif
#if PerXiDist
                if(v1 != "DistT" and (part == "Xi" or part == "Pi2")) key += GetDistKey(DC_[part][2]);
#endif
                if(effMap[key] == nullptr){
                    cout<<"Efficiency histogram "<<key<<" not found!"<<endl;
                    return 1.;
                }
                w = effMap[key]->GetEfficiency(effMap[key]->FindFixBin(var1));
                w = 1./w;
            }
            return w;
        }
        double GetWeight(TString part, TString num, TString den, TString v1, TString v2 = ""){
            return GetWeightDC(part, num, den, v1, v2, DC);
        }
        double GetWeightRec(TString part, TString num, TString den, TString v1, TString v2 = ""){
            return GetWeightDC(part, num, den, v1, v2, DCRec);
        }
        double GetEventVariable(TString var){
            double val;
            if(var == "ctau_L"){
                double gb_L = L.Mag()/mLd;
                val = (VL - VXi).Mag() / gb_L;
            }
            else if(var == "ctau_Xi"){
                double gb_Xi = Xi.Mag()/mXi;
                val = (VL - VXi).Mag() / gb_Xi;
            }
            else if(var == "CThKK"){
                TLorentzVector LVKM(Km, hypot(Km.Mag(), mk));
                TLorentzVector LVKP(Kp, hypot(Kp.Mag(), mk));
                TLorentzVector LVpTarget(0,0,0,mp);
                TLorentzVector LVpKM = LVpTarget + LVKM;
                auto Boost = LVpKM.BoostVector();
                LVKM.Boost(-Boost);
                LVKP.Boost(-Boost);
                TVector3 TVKp = LVKP.Vect();
                TVector3 TVKm = LVKM.Vect();
                val = cos(TVKp.Angle(TVKm));
            }
            else if(var == "SqrtS"){
                TLorentzVector LVKP(Kp, hypot(Kp.Mag(), mk));
                TLorentzVector LVXi(Xi, hypot(Xi.Mag(), mXi));
                auto LVKpXi = LVKP + LVXi;
                val = LVKpXi.M();
            }
            else if(var == "PolCTh"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetTheta());
            }
            else if(var == "PolCPhA"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetPhiA());
            }
            else if(var == "PolCPhB"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetPhiB());
            }
            else if(var == "PolCPhC"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetPhiC());
            }
            else if(var == "PolCThX"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetThetaX());
            }
            else if(var == "PolCThZ"){
                PolaConverter PC(Km, Kp, Xi, L, P, Pi1, Pi2);
                val = cos(PC.GetThetaZ());
            }
            else if(var == "MXiPScat"){
                TLorentzVector LVPScat(PScat, hypot(PScat.Mag(), mp));
                TLorentzVector LVXiProd(XiProd, hypot(XiProd.Mag(), mXi));
                auto LVXiPScat = LVPScat + LVXiProd;
                val = LVXiPScat.M();
            }
            else if(var == "CThXiPScat"){
                TLorentzVector LVPScat(PScat, hypot(PScat.Mag(), mp));
                TLorentzVector LVXiProd(XiProd, hypot(XiProd.Mag(), mXi));
                auto LVXiPScat = LVPScat + LVXiProd;
                auto Boost = LVXiPScat.BoostVector();
                LVPScat.Boost(-Boost);
                TVector3 TVPScat = LVPScat.Vect();
                TVector3 TVXiProd = LVXiProd.Vect();
                val = cos(TVPScat.Angle(TVXiProd));
            }
            else if(var == "Collinearity"){
                TVector3 TVL = L.Unit();
                TVector3 TVXi = Xi.Unit();
                val = TVL.Dot(TVXi);
            }
            return val;
        }
        double GetEventVariableRec(TString var){
            double val;
            if(var == "ctau_L"){
                double gb_L = LRec.Mag()/mLd;
                val = (VLRec - VXiRec).Mag() / gb_L;
            }
            else if(var == "ctau_Xi"){
                double gb_Xi = XiRec.Mag()/mXi;
                val = (VLRec - VXiRec).Mag() / gb_Xi;
            }
            else if(var == "CThKK"){
                TLorentzVector LVKM(Km, hypot(Km.Mag(), mk));
                TLorentzVector LVKP(Kp, hypot(Kp.Mag(), mk));
                TLorentzVector LVpTarget(0,0,0,mp);
                TLorentzVector LVpKM = LVpTarget + LVKM;
                auto Boost = LVpKM.BoostVector();
                LVKM.Boost(-Boost);
                LVKP.Boost(-Boost);
                TVector3 TVKp = LVKP.Vect();
                TVector3 TVKm = LVKM.Vect();
                val = cos(TVKp.Angle(TVKm));
            }
            else if(var == "SqrtS"){
                TLorentzVector LVKP(Kp, hypot(Kp.Mag(), mk));
                TLorentzVector LVXi(XiRec, hypot(XiRec.Mag(), mXi));
                auto LVKpXi = LVKP + LVXi;
                val = LVKpXi.M();
            }
            else if(var == "PolCTh"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetTheta());
            }
            else if(var == "PolCPhA"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetPhiA());
            }
            else if(var == "PolCPhB"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetPhiB());
            }
            else if(var == "PolCPhC"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetPhiC());
            }
            else if(var == "PolCThX"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetThetaX());
            }
            else if(var == "PolCThZ"){
                PolaConverter PC(Km, Kp, XiRec, LRec, PRec, Pi1Rec, Pi2Rec);
                val = cos(PC.GetThetaZ());
            }
            else if(var == "MXiPScat"){
                TLorentzVector LVPScat(PScatRec, hypot(PScatRec.Mag(), mp));
                TLorentzVector LVXiProd(XiProdRec, hypot(XiProdRec.Mag(), mXi));
                auto LVXiPScat = LVPScat + LVXiProd;
                val = LVXiPScat.M();
            }
            else if(var == "CThXiPScat"){
                TLorentzVector LVPScat(PScatRec, hypot(PScatRec.Mag(), mp));
                TLorentzVector LVXiProd(XiProdRec, hypot(XiProdRec.Mag(), mXi));
                auto LVXiPScat = LVPScat + LVXiProd;
                auto Boost = LVXiPScat.BoostVector();
                LVPScat.Boost(-Boost);
                TVector3 TVPScat = LVPScat.Vect();
                TVector3 TVXiProd = LVXiProd.Vect();
                val = cos(TVPScat.Angle(TVXiProd));
            }
            else if(var == "Collinearity"){
                TVector3 TVL = LRec.Unit();
                TVector3 TVXi = XiRec.Unit();
                val = TVL.Dot(TVXi);
            }
            return val;
        }
};
bool DebugLog = 1;
void FillHistograms(g4genfitcarbon* gf, TString tgt){
    Event event(gf, tgt);
    if(DebugLog) cout<<Form("Processing event...")<<endl;
    if(!event.Good()) return;
    if(DebugLog) cout<<Form("Loading event data...")<<endl;
    auto DC = event.GetDataContainer();
    auto DCRec = event.GetDataContainerRec();
    TString key;
    double w;
    double weight = 1;
    double weightRec = 1;
    for(auto cp:CorrPars){
        TString num = cp.num;
        TString den = cp.den;
        TString pc = cp.p_cor;
        if(DebugLog) cout<<Form("Applying correction for %s over %s...", num.Data(), den.Data())<<endl;
        vector<TString> var_cor = cp.var_cor;
        TString v1 = var_cor[0];
        TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
        if(DebugLog) cout<<"Loading Gen weight";
        weight *= event.GetWeight(pc, num, den, v1, v2);
        if(DebugLog) cout<<Form(" for %s over %s: w = %g", num.Data(), den.Data(), w)<<endl;
        if(gf->XiPflag){
            if(DebugLog) cout<<"Loading weight";
            double w = event.GetWeightRec(pc, num, den, v1, v2);
            if(DebugLog) cout<<Form(" for %s over %s: w = %g", num.Data(), den.Data(), w)<<endl;
            if(w > weight_th or isnan(w)){
                double var1 = GetVariable(DCRec[pc][0], DCRec[pc][1], DCRec[pc][2], v1);
                double var2 = v2 != "" ? GetVariable(DCRec[pc][0], DCRec[pc][1], DCRec[pc][2], v2) : 0;
                cout<<Form("Event weightRec for %s correction over threshold! w = %g, event count over threshold = %d", Correction(pc,num,den,v1,v2).Data(), w, over_cnt)<<endl;
                cout<<Form("Event variables: %s = %g, %s = %g", v1.Data(), var1, v2.Data(), var2)<<endl;
            }
            weightRec *= w;
        }
    }
    if(weight > weight_th) return;
    if(isnan(weight)) return;
    if(weightRec > weight_th or isnan(weightRec)){
        over_cnt++;
        cout<<Form("w = %g, wRec = %g, event count over threshold = %d", weight, weightRec, over_cnt)<<endl;
    }
#if SecondaryCorr 
    if(gf->Xiflag){
        for(auto ev:EventVars){
            for(auto sec:SecondaryCorrection){
                TString sec_name = sec.first;
                if(ev != sec_name) continue;
                TGraph* gr = sec.second;
                double var = event.GetEventVariableRec(sec_name);
                double eff_sec = gr->Eval(var);
                cout<<Form("Secondary correction for %s: var = %g, eff = %g", sec_name.Data(), var, eff_sec)<<endl;
                weightRec *= 1./eff_sec;
            }
        }
    }
#endif
    bool filled_weight = 0;
    for(auto p:particle){
        for(int iv=0;iv<variable.size();++iv){
            auto v = variable[iv];
            for(auto chk: CheckLists){
                if(DebugLog) cout<<Form("Filling histograms for particle %s, variable %s, check %s...", p.Data(), v.Data(), chk.Data())<<endl;
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
                if(w < weight_th and !isnan(w)){
                    hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1], DC[p][2], v), w);
                }
                for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                    auto v2 = variable[iv2];
                    if(pc == "P" or pc == "Pi1" or pc == "Pi2"){
                        if(v2 == "CosOpen") continue;
                    }
                    w = event.GetWeight(pc, num, den, v1, v2);
                    key = CorrectionHists(tgt, p, v, chk, pc, num, den, v1, v2);
                    if(w < weight_th and !isnan(w)){
                        hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1], DC[p][2], v), w);
                    }
                }
            }
        }
            }//chk
            if(SuffixCheck("XiPflag", gf)){
                if(DebugLog) cout<<Form("Filling histograms for particle %s, variable %s, XiPflag...", p.Data(), v.Data())<<endl;
                key = AcceptanceHistTitle1D(tgt, p, v, "XiPflag" + trig);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v));
                key = AcceptanceHistTitle1D(tgt, p, v, "XiPflagCor" + trig);
                hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v), weight);
                if(weightRec < weight_th and !isnan(weightRec)){
                    key = AcceptanceHistTitle1D(tgt, p, v, "XiPflagRec" + trig);
                    hMap[key]->Fill(GetVariable(DCRec[p][0], DCRec[p][1],DCRec[p][2], v));
                    key = AcceptanceHistTitle1D(tgt, p, v, "XiPflagRecCor" + trig);
                    hMap[key]->Fill(GetVariable(DCRec[p][0], DCRec[p][1],DCRec[p][2], v), weightRec);
                    if(!filled_weight){
                        hMap["h_weight"]->Fill(weightRec);
                        filled_weight = 1;
                    }
                }
            }
            if(SuffixCheck("GoodXiAndScatPGood", gf)){
                if(DebugLog) cout<<Form("Filling histograms for particle %s, variable %s, GoodXiAndScatPGood...", p.Data(), v.Data())<<endl;
                key = AcceptanceHistTitle1D(tgt, p, v, "GoodXiAndScatPGoodCor" + trig);
                if(weightRec < weight_th and !isnan(weightRec)){
                    hMap[key]->Fill(GetVariable(DC[p][0], DC[p][1],DC[p][2], v), weight);
                }
            }
        }//iv
    }//particle
    for(auto ev:EventVars){
        for(auto chk: CheckLists){
            if(!SuffixCheck(chk, gf)) continue;
            if(!TrigCheck(trig, gf)) continue;
            key = EventTitle(tgt, ev, chk);
            hMap[key]->Fill(event.GetEventVariable(ev));
            for(auto cp:CorrPars){
                TString num = cp.num;
                if(num != chk) continue;
                TString den = cp.den;
                TString pc = cp.p_cor;
                vector<TString> var_cor = cp.var_cor;
                TString v1 = var_cor[0];
                TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
                w = event.GetWeight(pc, num, den, v1, v2);
                key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1, v2);
                if(w < weight_th and !isnan(w)){
                    hMap[key]->Fill(event.GetEventVariable(ev), w);
                }
            }
        }
        if(SuffixCheck("XiPflag", gf)){
            key = EventTitle(tgt, ev, "XiPflag" + trig);
            hMap[key]->Fill(event.GetEventVariable(ev));
            key = EventTitle(tgt, ev, "XiPflagCor" + trig);
            hMap[key]->Fill(event.GetEventVariable(ev), weight);
            if(weightRec < weight_th and !isnan(weightRec)){
                key = EventTitle(tgt, ev, "XiPflagRec" + trig);
                hMap[key]->Fill(event.GetEventVariableRec(ev));
                key = EventTitle(tgt, ev, "XiPflagRecCor" + trig);
                hMap[key]->Fill(event.GetEventVariableRec(ev), weightRec);
            }
        }
        if(SuffixCheck("GoodXiAndScatPGood", gf)){
            key = EventTitle(tgt, ev, "GoodXiAndScatPGoodCor" + trig);
            if(weight < weight_th and !isnan(weight)){
                hMap[key]->Fill(event.GetEventVariable(ev), weight);
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
            if(maxi == 0 ){
                cout<<"Warning! Empty histogram for "<<key<<endl;
                maxi = 1;
                //Null results shuld also be recorded. 
            }
            key = AcceptanceHistTitle1D(tgt, p, v, "XiAcpt"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kBlue);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiRecAcpt"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kGreen+2);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiAcptCor"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kRed);
            hMap[key]->SetMarkerColor(kRed);
            hMap[key]->SetMarkerStyle(20);
            key = AcceptanceHistTitle1D(tgt, p, v, "XiRecAcptCor"+trig);
            hMap[key]->Scale(1./maxi);
            hMap[key]->SetLineColor(kOrange+7);
            hMap[key]->SetMarkerColor(kOrange+7);
            hMap[key]->SetMarkerStyle(20);

            key = AcceptanceHistTitle1D(tgt, p, v, "GoodXiAndScatPGoodCor"+trig);
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
            }//chklist
        }//iv
    }//particle
    //ctau hists
    for(auto ev:EventVars){
        double maxi = 0;
        key = EventTitle(tgt, ev, "Gen" + trig);
        maxi = hMap[key]->GetMaximum();
        if(maxi == 0 ){
            cout<<"Warning! Empty histogram for "<<key<<endl;
            maxi = 1;
            //Null results shuld also be recorded. 
        }
        key = EventTitle(tgt, ev, "XiAcpt"+trig);
        hMap[key]->Scale(1./maxi);
        hMap[key]->SetLineColor(kBlue);
        key = EventTitle(tgt, ev, "XiRecAcpt"+trig);
        hMap[key]->Scale(1./maxi);
        hMap[key]->SetLineColor(kGreen+2);
        key = EventTitle(tgt, ev, "XiAcptCor"+trig);
        hMap[key]->Scale(1./maxi);
        hMap[key]->SetLineColor(kRed);
        hMap[key]->SetMarkerColor(kRed);
        hMap[key]->SetMarkerStyle(20);
        key = EventTitle(tgt, ev, "XiRecAcptCor"+trig);
        hMap[key]->Scale(1./maxi);
        hMap[key]->SetLineColor(kOrange+7);
        hMap[key]->SetMarkerColor(kOrange+7);
        hMap[key]->SetMarkerStyle(20);
        key = EventTitle(tgt, ev, "GoodXiAndScatPGoodCor"+trig);
        hMap[key]->Scale(1./maxi);
        hMap[key]->SetLineColor(kGreen+2);
        hMap[key]->SetMarkerColor(kGreen+2);
        hMap[key]->SetMarkerStyle(20);
        for(auto chk: CheckLists){
            key = EventTitle(tgt, ev, chk);
            hMap[key]->Scale(1./maxi);
    for(auto cp:CorrPars){
        TString num = cp.num;
        if(num != chk) continue;
        TString den = cp.den;
        TString pc = cp.p_cor;
        for(int iv1 = 0; iv1 < variable.size();++iv1){
            auto v1 = variable[iv1];
            key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1);
            hMap[key]->Scale(1./maxi);
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                key = CorrectedEventTitle(tgt, ev, chk, pc, num, den, v1, v2);
                hMap[key]->Scale(1./maxi);
            }
        }
    }
        }//chklist
    }//ev

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
        if(y1 == 0) err = 0;
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

double CalcChi2(TH1* Gen, TH1* Cor, double trun = 0.0){
    //260730->trun 0.2 -> 0.
    //
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
        dels.push_back(abs(gen - cor));
        //dels.push_back((gen - cor)/err);
    }
    std::sort(dels.begin(), dels.end());
    int n = dels.size();
    int nc = 0;
    for(int i = n*trun; i < n ;++i){
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
            if(v1 == "DistT") continue;
            if(part_cor == "P" or part_cor == "Pi1" or part_cor == "Pi2"){
                if(v1 == "CosOpen") continue;
            }
            cor_chi2[{v1,""}] = CalcChi2(tgt, num, den, part_cor, v1);
            cout<<"chi2 for "<<part_cor<<" "<<num<<"/"<<den<<" with var "<<v1<<" is "<<cor_chi2[{v1,""}]<<endl;
            for(int iv2 = iv1+1; iv2 < variable.size();++iv2){
                auto v2 = variable[iv2];
                if(v2 == "DistT") continue;
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
            if(chi2_min < 0){
                chi2_min = chi2;
                best_var = var;
            }
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

void CheckAcceptance(vector<TFile*> files);
