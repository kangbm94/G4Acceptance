#define PerDist 1
#define PerXiDist 1
#define LVtxIsXi 1
#define TrigB 0
#define Recon 1
#define CH2 0
#define date 261001
#define PosShift 1
#define SecondaryCorr 0
#define DrawResSmeared 1
#define DrawGeneValue 0
#include "G4XiPScatAcceptanceCheck.hh"
#include <TCanvas.h>
TString WAcc,Target,Conf,LVtxConf,filename,figdir_base,file_dir, bufdir, fmap_dir,fmap_name;
TString tgt = "Carbon";
bool test_run = 0;
#if CH2
  tgt = "CH2";
#endif
double canv_size_x = 1500,canv_size_y = 1000;
double left_margin = 0.1, right_margin = 0.05, top_margin = 0.05, bottom_margin = 0.16;
double offset_x = 0.8, offset_y = 0.8;
double tsize = 0.09, lsize = 0.06;
void G4XiPScatAcceptanceCheck(){
  cout<<"Setting configurations..."<<endl;
  WAcc = "_WB";
  figdir_base = Form("figs_%d/%s/", date,tgt.Data());
  if(TrigB){
    figdir_base.ReplaceAll("figs","figs_TrigB");
  }
#if PerXiDist
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_PerXiDist").Data());
  bufdir.ReplaceAll(tgt.Data(),(tgt +"_PerXiDist").Data());
#endif
#if PerDist
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_PerDist").Data());
  bufdir.ReplaceAll(tgt.Data(),(tgt +"_PerDist").Data());
#endif
#if CH2
  file_dir = "rootfiles/Geant4CH2/6mmShift/wo_upstream_search/";
  Conf = "CH26mmShift";
  Target = "CH2";
#else
  file_dir = "rootfiles/Geant4Prod/XiPScat/";
  Conf = "ProdXiPScat6mmShift";
  Target = "Carbon";
#endif
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_6mmShift").Data());
  bufdir.ReplaceAll(tgt.Data(),(tgt +"_6mmShift").Data());
  fmap_dir = "./Maps/";
  bufdir = fmap_dir + Form("buf/%s/%d/",Conf.Data(),date);
  if(test_run){
    figdir_base.ReplaceAll(tgt,tgt + "_testrun");
  }
  LVtxConf = "_LVtxIsXi";
  figdir_base.ReplaceAll(tgt.Data(),(tgt + LVtxConf ).Data());
#if SecondaryCorr
  figdir_base.ReplaceAll(tgt.Data(),(tgt + "_SecondaryCorr").Data());
  bufdir.ReplaceAll(tgt.Data(),(tgt + "_SecondaryCorr").Data());
#endif
  cout<<"=======Usage======"<<endl;
  cout<<"G4XiPScatAcceptanceCheck(i,frac,ndiv) : Run the analysis for file index i, with fraction frac of ndiv parts"<<endl;
  cout<<"CheckAcceptanceAll() : Check the acceptance for all files"<<endl;
}

void G4XiPScatAcceptanceCheck(int i, int frac = 1, int ndiv = 1){
  if(file_dir == "") G4XiPScatAcceptanceCheck();
  SetStyle();
  TChain* tree = new TChain("tpc");
  cout<<"Loading files..."<<endl;
#if CH2
  filename = Form("XiReconCH2_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
#else
  filename = Form("XiPScatReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_XiPScat.root",i);
#endif
  if(FileCheck(file_dir+filename) == 0){
    cout<<"File "<<filename<<" not found. Exiting..."<<endl;
    return;
  }
  cout<<"Accpt file"<<endl;
#if CH2
  TFile* acpt_file;
  #if LVtxIsXi
    acpt_file = TFile::Open(Form("./Maps/CH2_ReconP0__LVtxIsXi_WB_6mmShifted_%d.root",260729));
    #if date >= 260817
      acpt_file = TFile::Open(Form("./Maps/CH2_ReconP0__LVtxIsXi_WB_6mmShifted_%d.root",260817));
    #endif
  #else
    //acpt_file = TFile::Open(Form("./Maps/%s_ReconPE42_%s%s_6mmShifted_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),260617));
  #endif
#else
  TFile* acpt_file;
  fmap_dir = "./Maps/";
  fmap_name = Form("%sXiPScat_Recon_P0_%s%s_6mmShifted_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),date);
  acpt_file = TFile::Open(fmap_dir + fmap_name);
#endif
  cout<<Form("Run %d, Target: %s, Date: %d", i,tgt.Data(), date)<<endl;
  cout<<"Map file :"<<acpt_file->GetName()<<endl;
  InitializeTriggerCondtions();
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  tree->Add(file_dir+filename);
  SetBranches(tree);
  InitializeCorrectionHistograms( tgt );
  LoadEfficiencies(acpt_file, tgt);
  auto ent = tree->GetEntries();
  cout<<"Entries : "<<ent<<endl;
  int start = 0,last = ent;
  if(ndiv != 1){
    int part = (ent / ndiv);
    start = part * (frac - 1);
    last = part * frac;
    if(frac == ndiv) last = ent;
  }
  for(auto i = start;i<last;i++){
    if(i%1000==0) cout << i << endl;
    Xi->GetEntry(i);
    if(i == start) cout<<"Start filling histograms..."<<endl;
    FillHistograms(Xi, tgt);
  }
  gSystem->mkdir(figdir_base, true);
  gSystem->mkdir(bufdir, true);
  TFile* out_file = TFile::Open(Form("%sCorrectedHistograms_%d_%d_%d.root", bufdir.Data(), i, frac, ndiv), "RECREATE");
  out_file->cd("");
  for(auto& [key, hist]: hMap){
    hist->Write();
  }
  out_file->Write();
  out_file->Close();

}
void CheckAcceptance(vector<TFile*> files){
  if(figdir_base == "") G4XiPScatAcceptanceCheck();
  SetStyle();
  gStyle->SetTitleSize(tsize,"XY");
  gStyle->SetLabelSize(lsize,"XY");
  gStyle->SetTitleOffset(offset_x,"X");
  gStyle->SetTitleOffset(offset_y,"Y");
  gStyle->SetPadLeftMargin(left_margin);
  gStyle->SetPadRightMargin(right_margin);
  gStyle->SetPadTopMargin(top_margin);
  gStyle->SetPadBottomMargin(bottom_margin);
  gROOT->ForceStyle();
  for(auto file: files){
    LoadCorrectionHistograms(file, tgt);
  }
  MakeChi2Map(tgt);
  NormalizeHistograms(tgt);
  TFile* out_file = TFile::Open(figdir_base + "CorrectionResults.root", "RECREATE");
  for(auto p:particle){
    //Each step by gen
    TString figdir = figdir_base  + "/StepByGen/";
    gSystem->mkdir(figdir, true);
    TString ct = Form("Canv_%s", p.Data());
    TCanvas* c = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    c->Divide(3,2,0.001,0.001);
		cout<<"StepByGen Particle: "<<p<<endl;
    for(int iv = 0; iv < variable.size();++iv){
      auto v1 = variable[iv];
      c->cd(iv+1);
      TString key = AcceptanceHistTitle1D(tgt, p, v1, "Gen"+trig);
			cout<<"Drawing "<<key<<endl;
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.6,0.75,0.9,0.95);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
      for(auto chk: CheckLists){
        if(chk == "Gen" + trig) continue;
        key = AcceptanceHistTitle1D(tgt, p, v1, chk);
        auto gr = MakeDivisionGraph(hMap[key], h_gen);
        cout<<"Divided by "<<key<<endl;
        gr->Draw("PE same");
        leg->AddEntry(gr, chk, "pe");
      }
      leg->Draw();
    }
    c->SaveAs(figdir + ct + ".pdf");
    
    //Step by step correction
    TString NumConf,DenConf; 
    for(auto cp:CorrPars){
      auto p_cor = cp.p_cor;
      if(p_cor == p){
        NumConf = cp.num;
        DenConf = cp.den;
				cout<<"StepByStep Particle: "<<p_cor<<" Num Den = "<<NumConf<<" Over "<<DenConf<<endl;
      }
    }
    figdir = figdir_base  + Form("/StepByStep/%s_Over_%s/", NumConf.Data(), DenConf.Data());
    gSystem->mkdir(figdir, true);
    for(auto parti: particle){
      gSystem->mkdir(figdir, true);
      ct = Form("CanvStepByStep_%s", parti.Data());
      TCanvas* c_sbs = new TCanvas(ct, ct, canv_size_x, canv_size_y);
      c_sbs->Divide(3,2,0.001,0.001);
      for(int iv = 0; iv < variable.size();++iv){
        auto v1 = variable[iv];
        c_sbs->cd(iv+1);
        TString key = AcceptanceHistTitle1D(tgt, parti, v1, DenConf);
        TH1* h_gen = (TH1*)hMap[key]->Clone();
        cout<<"Drawing "<<key<<endl;
        h_gen->GetYaxis()->SetRangeUser(0,1.5);
        h_gen->Draw("hist");
        cout<<"GenMaxi = "<<h_gen->GetMaximum()<<endl;
        h_gen->SetLineColor(kBlack);
        TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
        leg->SetTextSize(0.04);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(h_gen, DenConf + trig, "l");
        for(auto chk: CheckLists){
          if(chk == "Gen" + trig) continue;
          if(chk == DenConf) continue;
          key = AcceptanceHistTitle1D(tgt, parti, v1, chk);
          if(hMap[key]->Integral()>h_gen->Integral())continue;
          cout<<"Divided by "<<key<<"w Maxi "<<hMap[key]->GetMaximum()<<endl;
          auto gr = MakeDivisionGraph(hMap[key], h_gen);
          gr->Draw("PE same");
          leg->AddEntry(gr, chk, "pe");
        }
        leg->Draw();
      }
      c_sbs->SaveAs(figdir + ct + ".pdf");
    }

    //Recon with final correction
    TString figdir_rec = figdir_base + "/Recon/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr", p.Data());
    TCanvas* c_corr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    c_corr->Divide(3,2,0.001,0.001);
    for(int iv = 0; iv < variable.size();++iv){
      auto v1 = variable[iv];
      if(v1 == "CosOpen" and (p == "P" or p == "Pi1" or p == "Pi2")) continue;
      c_corr->cd(iv+1);
      TString key = AcceptanceHistTitle1D(tgt, p, v1, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0, h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, Labels["Gen" + trig], "l");
#if DrawGeneValue
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcpt"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, Labels["XiAcpt"+trig]+ "wo_resol", "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcptCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same"); 
      leg->AddEntry(h_xiacptcor, "Acceptance-corrected" + "wo_resol" , "pe");
#endif
#if DrawResSmeared
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiRecAcpt"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, Labels["XiAcpt"+trig], "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiRecAcptCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      leg->AddEntry(h_xirecacptcor, "Acceptance-corrected", "pe");
  #endif 
      //key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiCor"+trig);
      //auto h_goodxicor = (TH1*)hMap[key]->Clone();
      //h_goodxicor->Draw("pe same");
      //leg->AddEntry(h_goodxicor, Labels["GoodXiCor"], "pe");
      if(iv == 0) leg->Draw();
    }
    c_corr->SaveAs(figdir_rec + ct + ".pdf");
   
    //Recon with final correction and ratio
    figdir_rec = figdir_base + "/Recon_w_rat/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr_w_rat", p.Data());
    TCanvas* c_corr_rat = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    c_corr_rat->Divide(3,2,0.001,0.001);
    for(int iv = 0; iv < variable.size();++iv){
      auto v1 = variable[iv];
      if(v1 == "CosOpen" and (p == "P" or p == "Pi1" or p == "Pi2")) continue;
      c_corr_rat->cd(iv+1);
      TString key = AcceptanceHistTitle1D(tgt, p, v1, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0, h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, Labels["Gen" + trig], "l");
#if DrawGeneValue
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcpt"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, "XiAcpt"+trig, "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcptCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same"); 
      leg->AddEntry(h_xiacptcor, "XiAcptCor"+trig, "pe");
      auto gr = MakeDivisionGraph(h_xiacptcor, h_gen);
      gr->SetMarkerColor(kAzure);
      gr->Draw("PE same");
      leg->AddEntry(gr, "Ratio", "pe");
#endif
#if DrawResSmeared
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiRecAcpt"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, "XiRecAcpt"+trig, "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiRecAcptCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same"); 
      leg->AddEntry(h_xirecacptcor, "XiRecAcptCor"+trig, "pe");
      auto grrec = MakeDivisionGraph(h_xirecacptcor, h_gen);
      grrec->SetMarkerColor(kMagenta);
      grrec->Draw("PE same");
      leg->AddEntry(grrec, "Rec Ratio", "pe");
#endif
      //key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiCor"+trig);
      //auto h_goodxicor = (TH1*)hMap[key]->Clone();
      //h_goodxicor->Draw("pe same");
      //leg->AddEntry(h_goodxicor, "GoodXiCor", "pe");
      if(iv == 0) leg->Draw();
    }
    c_corr_rat->SaveAs(figdir_rec + ct + ".pdf");

    //Chi2 Evaluation and correction comparison
    vector<map<TString, double>> best_chi2_lists;
    for(int icm = 0; icm < chi2Maps_best.size();++icm){
      auto cp = chi2Maps_best[icm];
      double chi2 = cp.cor_chi2.begin()->second;
      TString pc = cp.part_cor;
      TString num = cp.num;
      TString den = cp.den;
      TString conf = pc + num + den;
      map<TString, double> best_chi2_map;
      best_chi2_map[pc+num+den] = chi2;
      best_chi2_lists.push_back(best_chi2_map);
      vector<TString> var_cor = cp.cor_chi2.begin()->first;
      TString v1 = var_cor[0];
      TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
      TString figdir = figdir_base + "BestChi2/" + p + "/";
      gSystem->mkdir(figdir, true); 
      ct = Form("Canv_%s_BestChi2Corr_%s_over_%s_for_%s_%s%s", p.Data(), num.Data(), den.Data(), pc.Data(), v1.Data(), v2.Data());
      TCanvas* c_corr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
      c_corr->Divide(3,2,0.001,0.001);
      for(int iv = 0; iv < variable.size();++iv){
        auto v = variable[iv];
        c_corr->cd(iv+1);
        TString key = AcceptanceHistTitle1D(tgt, p, v, den);
        TString key_num = AcceptanceHistTitle1D(tgt, p, v, num);
        TString key_cor = CorrectionHists(tgt, p, v, num, pc, num, den, v1, v2);
        TH1* h_den = (TH1*)hMap[key]->Clone();
        h_den->GetYaxis()->SetRangeUser(0, 1.5);
        h_den->SetLineColor(kBlack);
        h_den->Draw("hist");
        TH1* h_num = (TH1*)hMap[key_num]->Clone();
        h_num->SetLineColor(kBlue);
        h_num->Draw("hist same");
        TH1* h_cor = (TH1*)hMap[key_cor]->Clone();
        h_cor->SetLineColor(kRed);
        h_cor->Draw("hist same");
        auto gr = MakeDivisionGraph(h_cor, h_den);
        gr->Draw("PE same");
        TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
        leg->SetTextSize(0.04);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(h_den, Labels[den], "l");
        leg->AddEntry(h_num, Labels[num], "l");
        leg->AddEntry(h_cor, "Acceptance-corrected", "l");
        //leg->AddEntry(gr, Form("%s/%s cor w #chi^{2} %.4g", num.Data(), den.Data(), chi2), "pe");
        //leg->AddEntry(gr, Form("Acceptance-corrected yield; D^{2} %.4g", chi2), "pe");
        leg->AddEntry(gr, Form("Closure Ratio"), "pe");
        if(iv == 0)leg->Draw();
      }
      c_corr->SaveAs(figdir + ct + ".pdf");
    }
    for(int icm = 0; icm < chi2Maps.size(); ++icm){
      auto cm = chi2Maps[icm];
      TString pc = cm.part_cor;
      TString num = cm.num;
      TString den = cm.den;
      map<vector<TString>, double> cor_chi2 = cm.cor_chi2;
      for(auto c2 : cor_chi2){
        double chi2 = c2.second;
        vector<TString> var_cor = c2.first;
        TString conf = pc + num + den;
        bool skip = 0;
        for(auto best_chi2_map: best_chi2_lists){
          if(best_chi2_map.find(conf) != best_chi2_map.end()){
            if(chi2 == best_chi2_map[conf]) {
              skip = 1;
              break;
            }
          }
        }
        if(skip) continue;
        TString v1 = var_cor[0];
        TString v2 = var_cor.size() > 1 ? var_cor[1] : "";
        ct = Form("Canv_%s_Chi2Corr_%s_over_%s_for_%s_%s%s", p.Data(), num.Data(), den.Data(), pc.Data(), v1.Data(), v2.Data());
        TString figdir = figdir_base + "Chi2/" + conf + "/" ;
        cout<<"Making "<< ct << " in directory "<< figdir << " with chi2 " << chi2 << endl;
        gSystem->mkdir(figdir, true); 
        TCanvas* c_corr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
        c_corr->Divide(3,2,0.001,0.001);
        for(int iv = 0; iv < variable.size();++iv){
          auto v = variable[iv];
          c_corr->cd(iv+1);
          TString key = AcceptanceHistTitle1D(tgt, p, v, den);
          TString key_num = AcceptanceHistTitle1D(tgt, p, v, num);
          TString key_cor = CorrectionHists(tgt, p, v, num, pc, num, den, v1, v2);
          TH1* h_den = (TH1*)hMap[key]->Clone();
          h_den->GetYaxis()->SetRangeUser(0, 1.5);
          h_den->SetLineColor(kBlack);
          h_den->Draw("hist");
          TH1* h_num = (TH1*)hMap[key_num]->Clone();
          h_num->SetLineColor(kBlue);
          h_num->Draw("hist same");
          TH1* h_cor = (TH1*)hMap[key_cor]->Clone();
          h_cor->SetLineColor(kRed);
          h_cor->Draw("hist same");
          auto gr = MakeDivisionGraph(h_cor, h_den);
          gr->Draw("PE same");
          TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
          leg->SetTextSize(0.04);
          leg->SetFillStyle(0);
          leg->SetBorderSize(0);
          leg->AddEntry(h_den, Labels[den], "l");
          leg->AddEntry(h_num, Labels[num], "l");
          leg->AddEntry(h_cor, "Acceptance-corrected", "l");
          //leg->AddEntry(gr, Form("%s/%s cor w #chi^{2} %.4g", num.Data(), den.Data(), chi2), "pe");
          //leg->AddEntry(gr, Form("Acceptance-corrected yield; D^{2} %.4g", chi2), "pe");
          leg->AddEntry(gr, Form("Closure Ratio"), "pe");
          if(iv == 0)leg->Draw();
        }
        c_corr->SaveAs(figdir + ct + ".pdf");
      }
    }
  }//part
  figdir_base+= "/EventVars/";
  for(auto ev: EventVars){
    cout<<"EventVar: "<<ev<<endl;
    //Each step by gen
    TString figdir = figdir_base  + "StepByGen/";
    gSystem->mkdir(figdir, true);
    TString ct = Form("Canv_%s", ev.Data());
    TCanvas* c0 = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen"+trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
      for(auto chk: CheckLists){
        if(chk == "Gen" + trig) continue;
        key = EventTitle(tgt, ev, chk);
        auto gr = MakeDivisionGraph(hMap[key], h_gen);
        cout<<"Divided by "<<key<<endl;
        gr->Draw("PE same");
        leg->AddEntry(gr, chk, "pe");
      }
      leg->Draw();
    }
    c0->SaveAs(figdir + ct + ".pdf");

    //Step by step correction
    figdir = figdir_base  + "StepByStep/";
    gSystem->mkdir(figdir, true);
    ct = Form("CanvStepByStep_%s", ev.Data());
    TCanvas* c_sbs = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen"+trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
      for(auto chk: CheckLists){
        if(chk == "Gen" + trig) continue;
        key = EventTitle(tgt, ev, chk);
        if(hMap[key]->Integral()>h_gen->Integral())continue;
        cout<<"Divided by "<<key<<"w Maxi "<<hMap[key]->GetMaximum()<<endl;
        auto gr = MakeDivisionGraph(hMap[key], h_gen);
        gr->Draw("PE same");
        leg->AddEntry(gr, chk, "pe");
      }
    }
    c_sbs->SaveAs(figdir + ct + ".pdf");

    //Recon with final correction
    TString figdir_rec = figdir_base + "Recon/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr", ev.Data());
    TCanvas* c_corr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0, h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
      key = EventTitle(tgt, ev, "XiAcpt"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      key = EventTitle(tgt, ev, "XiRecAcpt"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      key = EventTitle(tgt, ev, "GoodXiCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      //h_goodxicor->Draw("pe same");
      key = EventTitle(tgt, ev, "XiAcptCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same");
      key = EventTitle(tgt, ev, "XiRecAcptCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      leg->AddEntry(h_xiacpt, "XiAcpt"+trig, "l");
      leg->AddEntry(h_xirecacpt, "XiRecAcpt"+trig, "l");
      //leg->AddEntry(h_goodxicor, "GoodXiCor", "pe");
      leg->AddEntry(h_xiacptcor, "XiAcptCor"+trig, "pe");
      leg->AddEntry(h_xirecacptcor, "XiRecAcptCor"+trig, "pe");
      leg->Draw();
    }
    c_corr->SaveAs(figdir_rec + ct + ".pdf");

    //Recon with final correction and ratio
    figdir_rec = figdir_base + "Recon_w_rat/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr_w_rat", ev.Data());
    TCanvas* c_corr_rat = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0, h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.75,0.9,0.95);
      leg->SetTextSize(0.04);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      key = EventTitle(tgt, ev, "XiAcpt"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      key = EventTitle(tgt, ev, "XiRecAcpt"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      key = EventTitle(tgt, ev, "GoodXiCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      //h_goodxicor->Draw("pe same");
      key = EventTitle(tgt, ev, "XiAcptCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same");
      key = EventTitle(tgt, ev, "XiRecAcptCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      auto gr = MakeDivisionGraph(h_xiacptcor, h_gen);
      gr->SetMarkerColor(kAzure);
      gr->Draw("PE same");
      auto grrec = MakeDivisionGraph(h_xirecacptcor, h_gen);
      grrec->SetMarkerColor(kMagenta);
      grrec->Draw("PE same");
      leg->AddEntry(h_gen, "Gen" + trig, "l");
      leg->AddEntry(h_xiacpt, "XiAcpt"+trig, "l");
      leg->AddEntry(h_xirecacpt, "XiRecAcpt"+trig, "l");
      //leg->AddEntry(h_goodxicor, "GoodXiCor", "pe");
      leg->AddEntry(h_xiacptcor, "XiAcptCor"+trig, "pe");
      leg->AddEntry(h_xirecacptcor, "XiRecAcptCor"+trig, "pe");
      leg->AddEntry(gr, "Ratio", "pe");
      leg->AddEntry(grrec, "Rec Ratio", "pe");
      leg->Draw();
      if(ev.Contains("Pol")){
        double mean_rat_gen =0,std_rat_gen=0, mean_rat_cor=0, std_rat_cor=0;
        int np = gr->GetN();
        for(int ip = 0; ip < np; ++ip){
          double x, y;
          gr->GetPoint(ip, x, y);
          mean_rat_gen += y;
          std_rat_gen += y*y;
          grrec->GetPoint(ip, x, y);
          mean_rat_cor += y;
          std_rat_cor += y*y;
        }
        if(np == 0){
          cout<<"No points in graph for "<<ev<<endl;
          np = 1;
          //Null results will still be recored.
        }
        mean_rat_gen /= np;
        std_rat_gen = sqrt(std_rat_gen/np - mean_rat_gen*mean_rat_gen);
        mean_rat_cor /= np;
        std_rat_cor = sqrt(std_rat_cor/np - mean_rat_cor*mean_rat_cor);
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->SetTextColor(kAzure);
        tex->DrawLatex(0.2, 0.3, Form("Gen Ratio: %.3f #pm %.5f", mean_rat_gen, std_rat_gen));
        tex->SetTextColor(kMagenta);
        tex->DrawLatex(0.2, 0.2, Form("Rec Ratio: %.3f #pm %.5f", mean_rat_cor, std_rat_cor));
      }
      // ROOT declares cd() with a null default path; pass an explicit empty
      // path to select this file while avoiding -Wnonnull from that default.
      out_file->cd("");
      gr->Write(Form("Graph_%s_Corr_w_rat", ev.Data()));
      grrec->Write(Form("Graph_%s_RecCorr_w_rat", ev.Data()));
      
    }
    c_corr_rat->SaveAs(figdir_rec + ct + ".pdf");
  }//EventVar
  out_file->cd("");
  for(auto& [key, hist]: hMap){
    hist->Write(); 
  }
  out_file->Write();
  out_file->Close();
}
void CheckAcceptanceAll(){
  if(figdir_base == "") G4XiPScatAcceptanceCheck();
  int nfile = 30;
  int ndiv = 5;
  int runs = ndiv;
  if(test_run){
    nfile = 1;
    runs = 1;
  }
  vector<TFile*> files;
  for(int i=0;i<nfile;++i){
    for(int frac = 1; frac <= runs; ++frac){
      TFile* file = TFile::Open(Form("%sCorrectedHistograms_%d_%d_%d.root", bufdir.Data(), i, frac, ndiv));
      files.push_back(file);
    }
  }
  CheckAcceptance(files);
}
void G4XiPScatAcceptanceCheckAll(){
  if(figdir_base == "")G4XiPScatAcceptanceCheck();
  SetStyle();
  double pxi = 0;
  int nfile = 30;
  bool test_run = 0;
  if(test_run) nfile = 1;
  TChain* tree = new TChain("tpc");
  cout<<"Loading files..."<<endl;
  bool file_check = 1;
  for(int i=0;i<nfile;++i){
#if CH2
    filename = Form("XiReconCH2_P_E42_%d_GenfitCarbonGeant4Ver5.root",i);
	#if PosShift
		filename = Form("XiReconCH2_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
	#endif
#else
    filename = Form("XiReconProd_P_E42_%d_GenfitCarbonGeant4Ver5.root",i);
	#if PosShift
		filename = Form("XiReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
	#endif
#endif
    cout<<"Loading "<<filename<<endl;
    if(FileCheck(file_dir+filename) == 0){
      file_check = 0;
      continue;
    }
    if(tree -> Add(file_dir+filename) != 1){
      cout<<"Error adding file "<<filename<<endl;
      file_check = 0;
    }
  }
  if(file_check == 0){
    cout<<"Error loading files. Exiting..."<<endl;
    return;
  }
  cout<<"Accpt file"<<endl;
#if CH2
  TFile* acpt_file;
  #if LVtxIsXi
    acpt_file = TFile::Open(Form("./Maps/CH2_ReconP0__LVtxIsXi_WB_6mmShifted_%d.root",260729));
  #else
    acpt_file = TFile::Open(Form("./Maps/%s_ReconPE42_%s%s_6mmShifted_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),260617));
  #endif
#else
  TFile* acpt_file;
  #if LVtxIsXi
    acpt_file = TFile::Open(Form("./Maps/Carbon_ReconP0__LVtxIsXi_WB_6mmShifted_%d.root",260726));
  #else
    acpt_file = TFile::Open(Form("./Maps/Carbon_ReconP0_WB_6mmShifted_%d.root",260726));
  #endif
#endif
  cout<<Form("Run Target: %s, Date: %d", tgt.Data(), date)<<endl;
  InitializeTriggerCondtions();
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
  InitializeCorrectionHistograms( tgt );
  LoadEfficiencies(acpt_file, tgt);
  auto ent = tree->GetEntries();
  cout<<"Entries : "<<ent<<endl;
  if(test_run) ent = ent / 10;
  for(auto i = 0;i<ent;i++){
    if(i%1000==0) cout << i << endl;
    Xi->GetEntry(i);
    FillHistograms(Xi, tgt);
  }
  TString figdir_base = Form("figs_%d/%s/", date,tgt.Data());
  if(TrigB){
    figdir_base.ReplaceAll("figs","figs_TrigB");
  }
#if PerXiDist
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_PerXiDist").Data());
#endif
#if PerDist
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_PerDist").Data());
#endif
#if PosShift
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_6mmShift").Data());
#endif
  if(test_run) figdir_base.ReplaceAll(tgt,tgt + "_testrun");
#if CH2
#else
  //figdir_base.ReplaceAll(tgt,tgt+ "ForcedPi2Ph");
#endif
#if LVtxIsXi
  figdir_base.ReplaceAll(tgt.Data(),(tgt +"_LVtxIsXi").Data());
#endif
  gSystem->mkdir(figdir_base, true);
  gSystem->mkdir(figdir_base + "rootfiles/", true);
  TFile* out_file = TFile::Open(bufdir + "rootfiles/AcceptanceCorrectionMaps.root", "RECREATE");
  out_file->cd("");
  for(auto& [key, hist]: hMap){
    hist->Write();
  }
  out_file->Write();
  vector<TFile*> files;
  files.push_back(out_file);
  CheckAcceptance(files);
}
