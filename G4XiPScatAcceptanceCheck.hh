#include "G4XiPScatAcceptance.hh"
#include "TGraphErrors.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
extern TString figdir_base;
int ScatPRecIndex = 0; // Selected recoil candidate in XiResidualsMom_*; original convention.

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
        nbinx = 30; minx = mXi + mp; maxx = 3.2;
    }
    else if(var == "CThXiPScat"){
        axis = "cos#theta_{p}^{helicity} (#Xi p CM)";
        nbinx = 20; minx = -1; maxx = 1;
    }
    else if(var == "Collinearity"){
        axis = "Collinearity";
        nbinx = 50; minx = 0.5; maxx = 1;
    }
    return axis;
}
CorrPars = {
    {"GoodXiAndScatPGood", "GoodXi", "ScatP", {"CosTh", "Mom"}},
    {"GoodXi", "GoodLAndPi2Tracked", "Xi", {"CosTh","CosOpen"}},
    {"GoodLAndPi2Tracked", "GoodL", "Pi2", {"Mom", "CosPsi"}},
    {"GoodL", "PPi1Tracked", "L", {"CosTh", "Mom"}},
    {"PPi1Tracked", "PTracked", "Pi1", {"Mom", "CosPsi"}},
    {"PTracked", "Gen", "P", {"CosTh", "Mom"}}
};
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
    if(trig == "TrigB") return;
    trig = "TrigB";
    for(auto &chk:CheckLists) chk += trig;
    for(auto &cp:CorrPars){ cp.num += trig; cp.den += trig; }
    vector<pair<TString, int>> colors;
    for(auto c:colorMap) colors.push_back({c.first+trig, c.second});
    for(auto c:colors) colorMap[c.first] = c.second;
    vector<pair<TString, TString>> labels;
    for(auto label:Labels) labels.push_back({label.first+trig, label.second});
    for(auto label:labels) Labels[label.first] = label.second;
#endif
}
vector<vector<TString>> CorrectionVariables(TString part){
    vector<vector<TString>> vars;
    for(int i=0;i<variable.size();++i){
        TString v1 = variable[i];
        if(v1 == "DistT") continue; // Already represented by the original vertex slices.
        if(v1 == "CosOpen" and (part == "P" or part == "Pi1" or part == "Pi2")) continue;
        vars.push_back({v1, ""});
        for(int j=i+1;j<variable.size();++j){
            TString v2 = variable[j];
            if(v2 == "DistT") continue;
            if(v2 == "CosOpen" and (part == "P" or part == "Pi1" or part == "Pi2")) continue;
            vars.push_back({v1, v2});
        }
    }
    return vars;
}
void SetCorrectionCheckLists(){
    for(auto cp:CorrPars){
        if(find(CheckLists.begin(), CheckLists.end(), cp.num) == CheckLists.end()) CheckLists.push_back(cp.num);
        if(find(CheckLists.begin(), CheckLists.end(), cp.den) == CheckLists.end()) CheckLists.push_back(cp.den);
    }
}
vector<TString> FinalXiPChecks(){
    return {"XiPflag", "XiPflagCor", "XiPflagRec", "XiPflagRecCor",
            "GoodXiAndScatPGoodCor", "GoodXiAndScatPGoodRec", "GoodXiAndScatPGoodRecCor"};
}
TString CorrectionParameters(){
    TString text;
    for(auto cp:CorrPars){
        text += cp.num + "\t" + cp.den + "\t" + cp.p_cor + "\t" + cp.var_cor[0] + "\t";
        text += (cp.var_cor.size() > 1 ? cp.var_cor[1] : "-") + "\n";
    }
    return text;
}
void LoadBestCorrections(TString filename){
    ifstream input(filename.Data());
    if(!input) throw runtime_error(string("Cannot read best corrections: ") + filename.Data());
    vector<CorrectionParameter> selected;
    string line;
    getline(input, line); // TSV header.
    while(getline(input, line)){
        if(line.empty()) continue;
        istringstream row(line);
        string num, den, part, v1, v2, penalty;
        if(!(row >> num >> den >> part >> v1 >> v2 >> penalty)) throw runtime_error("Invalid BestCorrections.tsv row.");
        double score = strtod(penalty.c_str(), nullptr);
        if(!isfinite(score)) throw runtime_error("A correction stage has no usable candidate.");
        vector<TString> vars = {TString(v1)};
        if(v2 != "-") vars.push_back(TString(v2));
        selected.push_back({TString(num), TString(den), TString(part), vars});
    }
    if(selected.size() != CorrPars.size()) throw runtime_error("Best correction file has different stages.");
    for(int i=0;i<CorrPars.size();++i){
        if(selected[i].num != CorrPars[i].num or selected[i].den != CorrPars[i].den or selected[i].p_cor != CorrPars[i].p_cor)
            throw runtime_error("Best correction file has different stages or trigger settings.");
        vector<TString> vars = {selected[i].var_cor[0], selected[i].var_cor.size()>1 ? selected[i].var_cor[1] : TString("")};
        auto allowed = CorrectionVariables(selected[i].p_cor);
        if(find(allowed.begin(), allowed.end(), vars) == allowed.end()) throw runtime_error("Unknown correction variables.");
    }
    CorrPars = selected;
}
void BookCorrectionHistogram(TString key, TString title, int bins, double lo, double hi, int color = kBlack){
    hMap[key] = new TH1D(key, title, bins, lo, hi);
    hMap[key]->SetDirectory(nullptr);
    hMap[key]->Sumw2();
    hMap[key]->SetLineColor(color);
}
void InitializeCorrectionHistograms(TString tgt){
    SetCorrectionCheckLists();
    BookCorrectionHistogram("h_weight", ";Correction weight;Events", 1000, 0, 50);
    BookCorrectionHistogram("h_weightRec", ";Correction weightRec;Events", 1000, 0, 50);
    BookCorrectionHistogram("h_weightRatio", ";WeightRatio;Events", 1000, 0, 5);
    BookCorrectionHistogram("h_weightDiff", ";Correction weightDiff;Events", 1000, -15, 15);
    BookCorrectionHistogram("hgood_weight", ";Correction weight;Events", 1000, 0, 50);
    BookCorrectionHistogram("hgood_weightRec", ";Correction weightRec;Events", 1000, 0, 50);
    BookCorrectionHistogram("hgood_weightRatio", ";WeightRatio;Events", 1000, 0, 5);
    BookCorrectionHistogram("hgood_weightDiff", ";Correction weightDiff;Events", 1000, -15, 15);
    BookCorrectionHistogram("hbad_weight", ";Correction weight;Events", 1000, 0, 50);
    BookCorrectionHistogram("hbad_weightRec", ";Correction weightRec;Events", 1000, 0, 50);
    BookCorrectionHistogram("hbad_weightRatio", ";WeightRatio;Events", 1000, 0, 5);
    BookCorrectionHistogram("hbad_weightDiff", ";Correction weightDiff;Events", 1000, -15, 15);

    for(int ic = 0; ic < CorrPars.size(); ++ic){
        auto cp = CorrPars[ic];
        TString key = "h_weight_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weight;Events", 1000, 0, 5);
        key = "h_weightRec_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightRec;Events", 1000, 0, 5);
        key = "hgood_weight_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weight;Events", 1000, 0, 5);
        key = "hgood_weightRec_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightRec;Events", 1000, 0, 5);
        key = "hbad_weight_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weight;Events", 1000, 0, 5);
        key = "hbad_weightRec_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightRec;Events", 1000, 0, 5);

        key = "h_weightRatio_" + cp.p_cor;
        BookCorrectionHistogram(key, ";WeightRatio;Events", 1000, 0, 5);
        key = "h_weightDiff_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightDiff;Events", 1000, -15, 15);
        key = "hgood_weightRatio_" + cp.p_cor;
        BookCorrectionHistogram(key, ";WeightRatio;Events", 1000, 0, 5);
        key = "hgood_weightDiff_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightDiff;Events", 1000, -15, 15);
        key = "hbad_weightRatio_" + cp.p_cor;
        BookCorrectionHistogram(key, ";WeightRatio;Events", 1000, 0, 5);
        key = "hbad_weightDiff_" + cp.p_cor;
        BookCorrectionHistogram(key, ";Correction weightDiff;Events", 1000, -15, 15);
    }
    BookCorrectionHistogram("h_final_weights", ";Truth accepted/rejected, reco accepted/rejected;Events", 4, 0, 4);
    BookCorrectionHistogram("h_reco_available", ";XiPflag: valid/invalid reconstructed arrays;Events", 2, 0, 2);
    for(auto cp:CorrPars){
        for(auto vars:CorrectionVariables(cp.p_cor)){
            TString key = "WeightCheck" + Correction(cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
            BookCorrectionHistogram(key, ";Usable/rejected inverse-efficiency weights;Events", 2, 0, 2);
        }
    }
    int nbinx;
    double minx, maxx;
    for(auto chk:CheckLists){
        for(auto p:particle){
            for(auto v:variable){
                TString title = SetAxis(p, v, nbinx, minx, maxx);
                TString key = AcceptanceHistTitle1D(tgt, p, v, chk);
                BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx, colorMap[chk]);
                for(auto cp:CorrPars){
                    if(cp.num != chk) continue;
                    for(auto vars:CorrectionVariables(cp.p_cor)){
                        key = CorrectionHists(tgt, p, v, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
                        BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx, colorMap[chk]);
                    }
                }
            }
        }
        for(auto ev:EventVars){
            TString title = SetEventVarAxis(ev, nbinx, minx, maxx);
            TString key = EventTitle(tgt, ev, chk);
            BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx, colorMap[chk]);
            for(auto cp:CorrPars){
                if(cp.num != chk) continue;
                for(auto vars:CorrectionVariables(cp.p_cor)){
                    key = CorrectedEventTitle(tgt, ev, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
                    BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx, colorMap[chk]);
                }
            }
        }
    }
    for(auto conf:FinalXiPChecks()){
        for(auto p:particle){
            for(auto v:variable){
                TString title = SetAxis(p, v, nbinx, minx, maxx);
                TString key = AcceptanceHistTitle1D(tgt, p, v, conf+trig);
                BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx);
            }
        }
        for(auto ev:EventVars){
            TString title = SetEventVarAxis(ev, nbinx, minx, maxx);
            TString key = EventTitle(tgt, ev, conf+trig);
            BookCorrectionHistogram(key, key+";"+title, nbinx, minx, maxx);
        }
    }
}
void LoadOrAddHistogram(TString key, TFile* file, TString chk = ""){
    LoadOrAdd1DHistogram(key, file);
    if(chk != "") hMap[key]->SetLineColor(colorMap[chk]);
}
void LoadCorrectionHistograms(TFile* file,TString tgt){
    if(file == nullptr or file->IsZombie()) throw runtime_error("Cannot read correction worker file.");
    TNamed* parameters = dynamic_cast<TNamed*>(file->Get("CorrectionParameters"));
    if(parameters == nullptr or TString(parameters->GetTitle()) != CorrectionParameters())
        throw runtime_error("Correction worker files use different parameters. Rerun workers with matching CorrPars.");
    cout<<"Loading histograms from "<<file->GetName()<<endl;
    SetCorrectionCheckLists();
    for(auto key:{"h_weight", "h_final_weights", "h_reco_available",
         "h_weightRec", "hgood_weight", "hgood_weightRec", "hbad_weight", "hbad_weightRec",
        "h_weightRatio", "h_weightDiff", 
        "hgood_weightRatio", "hgood_weightDiff", "hbad_weightRatio", "hbad_weightDiff"
        }) LoadOrAddHistogram(key, file);
    for(auto cp:CorrPars){
        TString key = "h_weight_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "h_weightRec_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hgood_weight_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hgood_weightRec_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hbad_weight_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hbad_weightRec_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "h_weightRatio_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "h_weightDiff_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hgood_weightRatio_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hgood_weightDiff_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hbad_weightRatio_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
        key = "hbad_weightDiff_" + cp.p_cor;
        LoadOrAddHistogram(key, file);
    }
    for(auto cp:CorrPars){
        for(auto vars:CorrectionVariables(cp.p_cor))
            LoadOrAddHistogram("WeightCheck"+Correction(cp.p_cor, cp.num, cp.den, vars[0], vars[1]), file);
    }
    for(auto chk:CheckLists){
        for(auto p:particle){
            for(auto v:variable){
                LoadOrAddHistogram(AcceptanceHistTitle1D(tgt, p, v, chk), file, chk);
                for(auto cp:CorrPars){
                    if(cp.num != chk) continue;
                    for(auto vars:CorrectionVariables(cp.p_cor))
                        LoadOrAddHistogram(CorrectionHists(tgt, p, v, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1]), file, chk);
                }
            }
        }
        for(auto ev:EventVars){
            LoadOrAddHistogram(EventTitle(tgt, ev, chk), file, chk);
            for(auto cp:CorrPars){
                if(cp.num != chk) continue;
                for(auto vars:CorrectionVariables(cp.p_cor))
                    LoadOrAddHistogram(CorrectedEventTitle(tgt, ev, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1]), file, chk);
            }
        }
    }
    for(auto conf:FinalXiPChecks()){
        for(auto p:particle){
            for(auto v:variable) LoadOrAddHistogram(AcceptanceHistTitle1D(tgt, p, v, conf+trig), file);
        }
        for(auto ev:EventVars) LoadOrAddHistogram(EventTitle(tgt, ev, conf+trig), file);
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
        bool hasRec = false;
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
            // XiPflag is the only observed-event flag. No MC proton flag is consulted here.
            if(XiPOfflineFlag(gf)){
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
                    hasRec = true;
                    PScatRec.SetXYZ(gf->XipVtxScatPMom_x, gf->XipVtxScatPMom_y, gf->XipVtxScatPMom_z);
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
        bool HasReconstruction(){ return hasRec; }
        map<TString, vector<TVector3>> GetDataContainer(){
            return DC;
        }
        map<TString, vector<TVector3>> GetDataContainerRec(){
            return DCRec;
        }
        bool Good(){
            if(Xi.Mag() <1e-5) return false;
            if(isnan(Xi.Mag()) or isinf(Xi.Mag())) return false;
            #if TrigB
            if(gf->nhHtof < 2) return false;
            #endif
            return true;
        }
        double GetWeightDC(TString part, TString num, TString den, TString v1, TString v2 = "", map<TString, vector<TVector3>> DC_ = map<TString, vector<TVector3>>()){
            auto found = DC_.find(part);
            if(found == DC_.end() or found->second.size() < 3) return numeric_limits<double>::quiet_NaN();
            auto vectors = found->second;
            TString key = v2 == "" ? EffTitle1D(tgt, part, v1, num, den) : EffTitle2D(tgt, part, v1, v2, num, den);
            TString distkey;
#if PerDist
            if(v1 != "DistT" and v2 != "DistT" and part != "Xi" and part != "Pi2") distkey = GetDistKey(vectors[2]);
#endif
#if PerXiDist
            if(v1 != "DistT" and v2 != "DistT" and (part == "Xi" or part == "Pi2")) distkey = GetDistKey(vectors[2]);
#endif
            key += distkey;
            auto eff = effMap.find(key);
            if(eff == effMap.end() and v2 != ""){
                key = EffTitle2D(tgt, part, v2, v1, num, den) + distkey;
                eff = effMap.find(key);
                swap(v1, v2);
            }
            if(eff == effMap.end() or eff->second == nullptr)
                throw runtime_error(string("Efficiency histogram not found: ") + key.Data());
            double x = GetXiPVariable(vectors[0], vectors[1], vectors[2], part, v1);
            double y = v2 != "" ? GetXiPVariable(vectors[0], vectors[1], vectors[2], part, v2) : 0;
            if(!isfinite(x) or !isfinite(y)) return numeric_limits<double>::quiet_NaN();
            const TH1* total = eff->second->GetTotalHistogram();
            // Keep the standard ROOT map-bin lookup, including populated under/overflow bins.
            int bin = v2 == "" ? eff->second->FindFixBin(x) : eff->second->FindFixBin(x, y);
            if(total->GetBinContent(bin) <= 0) return numeric_limits<double>::quiet_NaN();
            double value = eff->second->GetEfficiency(bin);
            return value > 0 and isfinite(value) ? 1./value : numeric_limits<double>::quiet_NaN();
        }
        double GetWeight(TString part, TString num, TString den, TString v1, TString v2 = ""){
            return GetWeightDC(part, num, den, v1, v2, DC);
        }
        double GetWeightRec(TString part, TString num, TString den, TString v1, TString v2 = ""){
            if(!hasRec) return numeric_limits<double>::quiet_NaN();
            return GetWeightDC(part, num, den, v1, v2, DCRec);
        }
        double GetEventVariable(TString var){
            double val = numeric_limits<double>::quiet_NaN();
            if(var == "ctau_L"){
                double gb_L = L.Mag()/mLd;
                val = (VL - VXi).Mag() / gb_L;
            }
            else if(var == "ctau_Xi"){
                double gb_Xi = Xi.Mag()/mXi;
                val = (VXi - VXiProd).Mag() / gb_Xi;
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
                TVector3 TVXiProd = LVXiPScat.Vect(); // Pair lab direction: proton helicity axis.
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
            double val = numeric_limits<double>::quiet_NaN();
            if(var == "ctau_L"){
                double gb_L = LRec.Mag()/mLd;
                val = (VLRec - VXiRec).Mag() / gb_L;
            }
            else if(var == "ctau_Xi"){
                double gb_Xi = XiRec.Mag()/mXi;
                val = (VXiRec - VXiProdRec).Mag() / gb_Xi;
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
                TVector3 TVXiProd = LVXiPScat.Vect(); // Pair lab direction: proton helicity axis.
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
bool GoodCorrectionWeight(double w){ return isfinite(w) and w > 0 and w <= weight_th; }
void FillCorrectionHistogram(TString key, double value, double weight = 1.){
    if(isfinite(value) and isfinite(weight) and weight > 0) hMap[key]->Fill(value, weight);
}
void FillFinalXiP(Event &event, TString tgt, TString conf, double weight = 1., bool reconstructed = false){
    auto DC = reconstructed ? event.GetDataContainerRec() : event.GetDataContainer();
    for(auto p:particle){
        for(auto v:variable){
            double value = GetXiPVariable(DC[p][0], DC[p][1], DC[p][2], p, v);
            FillCorrectionHistogram(AcceptanceHistTitle1D(tgt, p, v, conf+trig), value, weight);
        }
    }
    for(auto ev:EventVars){
        double value = reconstructed ? event.GetEventVariableRec(ev) : event.GetEventVariable(ev);
        FillCorrectionHistogram(EventTitle(tgt, ev, conf+trig), value, weight);
    }
}
void FillHistograms(g4genfitcarbon* gf, TString tgt){
    Event event(gf, tgt);
    if(!event.Good()) return;
    auto DC = event.GetDataContainer();
    bool observed = gf->XiPflag;
    bool matched = gf->xigood and gf->p_scat_good;
    // References and uncorrected samples are independent of weight support or rejection.
    for(auto chk:CheckLists){
        if(!SuffixCheck(chk, gf) or !TrigCheck(trig, gf)) continue;
        for(auto p:particle){
            for(auto v:variable)
                FillCorrectionHistogram(AcceptanceHistTitle1D(tgt, p, v, chk), GetXiPVariable(DC[p][0], DC[p][1], DC[p][2], p, v));
        }
        for(auto ev:EventVars) FillCorrectionHistogram(EventTitle(tgt, ev, chk), event.GetEventVariable(ev));
    }
    if(observed){
        FillFinalXiP(event, tgt, "XiPflag");
        hMap["h_reco_available"]->Fill(event.HasReconstruction() ? .5 : 1.5);
        if(event.HasReconstruction()){
            FillFinalXiP(event, tgt, "XiPflagRec", 1., true);
            if(matched) FillFinalXiP(event, tgt, "GoodXiAndScatPGoodRec", 1., true);
        }
    }
    // Each candidate corrects its own numerator to its own denominator, as in XiAcceptanceCheck.
    for(auto cp:CorrPars){
        if(!SuffixCheck(cp.num, gf) or !TrigCheck(trig, gf)) continue;
        for(auto vars:CorrectionVariables(cp.p_cor)){
            double w = event.GetWeight(cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
            TString stats = "WeightCheck" + Correction(cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
            hMap[stats]->Fill(GoodCorrectionWeight(w) ? .5 : 1.5);
            if(!GoodCorrectionWeight(w)) continue;
            for(auto p:particle){
                for(auto v:variable){
                    TString key = CorrectionHists(tgt, p, v, cp.num, cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
                    FillCorrectionHistogram(key, GetXiPVariable(DC[p][0], DC[p][1], DC[p][2], p, v), w);
                }
            }
            for(auto ev:EventVars){
                TString key = CorrectedEventTitle(tgt, ev, cp.num, cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
                FillCorrectionHistogram(key, event.GetEventVariable(ev), w);
            }
        }
    }
    if(!observed and !matched) return;
    double weight = 1., weightRec = 1.;
    for(auto cp:CorrPars){
        TString v1 = cp.var_cor[0], v2 = cp.var_cor.size()>1 ? cp.var_cor[1] : TString("");
        double w_step = event.GetWeight(cp.p_cor, cp.num, cp.den, v1, v2);
        TString wkey = "h_weight_" + cp.p_cor;
        hMap[wkey]->Fill(w_step);
        weight *= w_step;
        // Do not require truth matching or an MC-only proton-tracking flag for reconstructed candidates.
        double w_stepRec = event.GetWeightRec(cp.p_cor, cp.num, cp.den, v1, v2);
        weightRec *= w_stepRec;
        wkey = "h_weightRec_" + cp.p_cor;
        hMap[wkey]->Fill(w_stepRec);
        wkey = "h_weightRatio_" + cp.p_cor;
        hMap[wkey]->Fill(w_stepRec/w_step);
        wkey = "h_weightDiff_" + cp.p_cor;
        hMap[wkey]->Fill(w_stepRec-w_step);
        if(matched){
            wkey = "hgood_weight_" + cp.p_cor;
            hMap[wkey]->Fill(w_step);
            wkey = "hgood_weightRec_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec);
           
            wkey = "hgood_weightRatio_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec/w_step);
            wkey = "hgood_weightDiff_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec-w_step);
        }
        else{
            wkey = "hbad_weight_" + cp.p_cor;
            hMap[wkey]->Fill(w_step);
            wkey = "hbad_weightRec_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec);
            
            wkey = "hbad_weightRatio_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec/w_step);
            wkey = "hbad_weightDiff_" + cp.p_cor;
            hMap[wkey]->Fill(w_stepRec-w_step);
        }
    }
    hMap["h_weight"]->Fill(weight);
    hMap["h_weightRec"]->Fill(weightRec);
    hMap["h_weightRatio"]->Fill(weightRec/weight);
    hMap["h_weightDiff"]->Fill(weightRec-weight);
    if(matched){
        hMap["hgood_weight"]->Fill(weight);
        hMap["hgood_weightRec"]->Fill(weightRec);
        hMap["hgood_weightRatio"]->Fill(weightRec/weight);
        hMap["hgood_weightDiff"]->Fill(weightRec-weight);
    }
    else{
        hMap["hbad_weight"]->Fill(weight);
        hMap["hbad_weightRec"]->Fill(weightRec);
        hMap["hbad_weightRatio"]->Fill(weightRec/weight);
        hMap["hbad_weightDiff"]->Fill(weightRec-weight);
    }
    if(GoodCorrectionWeight(weight)){
        hMap["h_final_weights"]->Fill(.5);
        if(observed){
            FillFinalXiP(event, tgt, "XiPflagCor", weight);
        }
        if(matched) FillFinalXiP(event, tgt, "GoodXiAndScatPGoodCor", weight);
    }
    else{ hMap["h_final_weights"]->Fill(1.5); ++over_cnt; }
    if(observed and event.HasReconstruction()){
        if(GoodCorrectionWeight(weightRec)){
            hMap["h_final_weights"]->Fill(2.5);
            FillFinalXiP(event, tgt, "XiPflagRecCor", weightRec, true);
            if(matched) FillFinalXiP(event, tgt, "GoodXiAndScatPGoodRecCor", weightRec, true);
        }
        else hMap["h_final_weights"]->Fill(3.5);
    }
}
void NormalizeHistograms(TString tgt){
    // Only for figure display; CheckAcceptance writes raw counts before calling this.
    for(auto p:particle){
        for(auto v:variable){
            double maxi = hMap[AcceptanceHistTitle1D(tgt, p, v, "Gen"+trig)]->GetMaximum();
            if(maxi <= 0) maxi = 1;
            for(auto chk:CheckLists){
                hMap[AcceptanceHistTitle1D(tgt, p, v, chk)]->Scale(1./maxi);
                for(auto cp:CorrPars){
                    if(cp.num != chk) continue;
                    for(auto vars:CorrectionVariables(cp.p_cor))
                        hMap[CorrectionHists(tgt, p, v, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1])]->Scale(1./maxi);
                }
            }
            for(auto conf:FinalXiPChecks()) hMap[AcceptanceHistTitle1D(tgt, p, v, conf+trig)]->Scale(1./maxi);
        }
    }
    for(auto ev:EventVars){
        double maxi = hMap[EventTitle(tgt, ev, "Gen"+trig)]->GetMaximum();
        if(maxi <= 0) maxi = 1;
        for(auto chk:CheckLists){
            hMap[EventTitle(tgt, ev, chk)]->Scale(1./maxi);
            for(auto cp:CorrPars){
                if(cp.num != chk) continue;
                for(auto vars:CorrectionVariables(cp.p_cor))
                    hMap[CorrectedEventTitle(tgt, ev, chk, cp.p_cor, cp.num, cp.den, vars[0], vars[1])]->Scale(1./maxi);
            }
        }
        for(auto conf:FinalXiPChecks()) hMap[EventTitle(tgt, ev, conf+trig)]->Scale(1./maxi);
    }
    vector<TString> styles = {"XiPflag", "XiPflagCor", "XiPflagRec", "XiPflagRecCor",
                             "GoodXiAndScatPGood", "GoodXiAndScatPGoodCor", "GoodXiAndScatPGoodRec", "GoodXiAndScatPGoodRecCor"};
    vector<int> colors = {kBlue, kRed, kBlue, kOrange+7, kGreen+2, kGreen+2, kViolet, kViolet};
    for(int i=0;i<styles.size();++i){
        for(auto p:particle){
            for(auto v:variable){
                TString key = AcceptanceHistTitle1D(tgt, p, v, styles[i]+trig);
                hMap[key]->SetLineColor(colors[i]); hMap[key]->SetMarkerColor(colors[i]); hMap[key]->SetMarkerStyle(20+i%3);
            }
        }
        for(auto ev:EventVars){
            TString key = EventTitle(tgt, ev, styles[i]+trig);
            hMap[key]->SetLineColor(colors[i]); hMap[key]->SetMarkerColor(colors[i]); hMap[key]->SetMarkerStyle(20+i%3);
        }
    }
}
void SetFinalXiPPlotRange(TH1* reference, TString base){
    // Include corrected yields and their weight errors, including mis-reconstructed candidates.
    double maximum = 0.;
    vector<TString> shown = {"Gen", "GoodXiAndScatPGood", "GoodXiAndScatPGoodCor"};
#if DrawGeneValue
    shown.push_back("XiPflag"); shown.push_back("XiPflagCor");
#endif
#if DrawResSmeared
    shown.push_back("XiPflagRec"); shown.push_back("XiPflagRecCor");
#endif
    for(auto conf:shown){
        TH1* hist = hMap[base+conf+trig];
        for(int bin=1;bin<=hist->GetNbinsX();++bin)
            maximum = max(maximum, hist->GetBinContent(bin)+hist->GetBinError(bin));
    }
    reference->GetYaxis()->SetRangeUser(0, maximum > 0 ? 1.2*maximum : 1.);
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
        double err = hypot(err1/y2, y1*err2/(y2*y2));
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

// D^2: absolute count differences, averaged over ordinary bins; no error weighting.
double CalcPenalty(TH1* Gen, TH1* Cor, double trun = 0.){
    vector<double> differences;
    for(int i=1;i<=Gen->GetNbinsX();++i) differences.push_back(abs(Gen->GetBinContent(i)-Cor->GetBinContent(i)));
    sort(differences.begin(), differences.end());
    double penalty = 0.;
    int first = int(differences.size()*trun), count = 0;
    for(int i=first;i<differences.size();++i){ penalty += differences[i]*differences[i]; ++count; }
    return count > 0 ? penalty/count : numeric_limits<double>::infinity();
}
double CalcPenalty(TString tgt, TString num, TString den, TString part_cor, TString v1, TString v2 = ""){
    double penalty = 0.;
    for(auto p:particle){
        for(auto v:variable){
            if(v == "CosOpen" and (p == "P" or p == "Pi1" or p == "Pi2")) continue;
            TString reference = AcceptanceHistTitle1D(tgt, p, v, den);
            TString corrected = CorrectionHists(tgt, p, v, num, part_cor, num, den, v1, v2);
            penalty += CalcPenalty(hMap[reference], hMap[corrected]);
        }
    }
    return penalty;
}
// Retain the original helper name for macro compatibility. This is not statistical chi-square.
double CalcChi2(TH1* Gen, TH1* Cor, double trun = 0.){ return CalcPenalty(Gen, Cor, trun); }
double CalcChi2(TString tgt, TString num, TString den, TString part_cor, TString v1, TString v2 = ""){
    return CalcPenalty(tgt, num, den, part_cor, v1, v2);
}
struct Chi2Map{
    TString num, den, part_cor;
    map<vector<TString>, double> cor_chi2;
};
vector<Chi2Map> chi2Maps;
vector<Chi2Map> chi2Maps_best;
void MakeChi2Map(TString tgt){
    chi2Maps.clear(); chi2Maps_best.clear();
    gSystem->mkdir(figdir_base, true);
    ofstream scores((figdir_base+"PenaltyScores.tsv").Data());
    ofstream best((figdir_base+"BestCorrections.tsv").Data());
    if(!scores or !best) throw runtime_error("Cannot write correction scores in figure directory.");
    scores << "num\tden\tparticle\tvar1\tvar2\tpenalty\tusable\trejected\n" << setprecision(17);
    best << "num\tden\tparticle\tvar1\tvar2\tpenalty\n" << setprecision(17);
    for(auto cp:CorrPars){
        map<vector<TString>, double> penalties;
        for(auto vars:CorrectionVariables(cp.p_cor)){
            TString stats = "WeightCheck"+Correction(cp.p_cor, cp.num, cp.den, vars[0], vars[1]);
            double usable = hMap[stats]->GetBinContent(1), rejected = hMap[stats]->GetBinContent(2);
            double penalty = usable > 0 ? CalcPenalty(tgt, cp.num, cp.den, cp.p_cor, vars[0], vars[1]) : numeric_limits<double>::infinity();
            penalties[vars] = penalty;
            scores << cp.num << '\t' << cp.den << '\t' << cp.p_cor << '\t' << vars[0] << '\t'
                   << (vars[1] == "" ? TString("-") : vars[1]) << '\t' << penalty << '\t' << usable << '\t' << rejected << '\n';
        }
        chi2Maps.push_back({cp.num, cp.den, cp.p_cor, penalties});
        double minimum = numeric_limits<double>::infinity();
        vector<TString> selected = penalties.begin()->first;
        for(auto result:penalties){
            if(result.second < minimum){ minimum = result.second; selected = result.first; }
        }
        chi2Maps_best.push_back({cp.num, cp.den, cp.p_cor, {{selected, minimum}}});
        best << cp.num << '\t' << cp.den << '\t' << cp.p_cor << '\t' << selected[0] << '\t'
             << (selected[1] == "" ? TString("-") : selected[1]) << '\t' << minimum << '\n';
        cout << "Best D^2 for " << cp.p_cor << ' ' << cp.num << '/' << cp.den << ": " << minimum
             << " with " << selected[0] << ' ' << selected[1] << endl;
    }
}

void CheckAcceptance(vector<TFile*> files, bool make_figures = true);
