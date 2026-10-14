#define PerDist 1
#define PerXiDist 1
#define LVtxIsXi 1
#define TrigB 0
#define Recon 1
#define CH2 0
#define Debug 0
#define date 261006
#define date_map 261006
#define PosShift 1
#define SecondaryCorr 0
// Final XiPflag overlays: reconstructed or generated kinematics (both may be drawn).
#ifndef DrawResSmeared
#define DrawResSmeared 1
#endif
#ifndef DrawGeneValue
#define DrawGeneValue 0
#endif
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
  WAcc = "_WB";
  Conf = "ProdXiPScat6mmShift";
  Target = "Carbon";
  tgt = Target;
  LVtxConf = LVtxIsXi ? "_LVtxIsXi" : "";
  if(file_dir == "") file_dir = "rootfiles/Geant4Prod/XiPScat/";
  if(fmap_dir == "") fmap_dir = "./Maps/";
  if(!fmap_dir.EndsWith("/")) fmap_dir += "/";
  if(fmap_name == "") fmap_name = Form("%sXiPScat_Recon_P0_%s%s_6mmShifted_%d.root",Target.Data(),LVtxConf.Data(),WAcc.Data(),date_map);
  if(bufdir == "") bufdir = fmap_dir+Form("buf/%s/%d/Check/",Conf.Data(),date);
  if(!bufdir.EndsWith("/")) bufdir += "/";
  if(figdir_base == ""){
    figdir_base = Form("figs_%d/%sXiPScat/",date,tgt.Data());
#if TrigB
    figdir_base.ReplaceAll("figs", "figs_TrigB");
#endif
#if PerXiDist
    figdir_base.ReplaceAll(tgt, tgt+"_PerXiDist");
#endif
#if PerDist
    figdir_base.ReplaceAll(tgt, tgt+"_PerDist");
#endif
#if PosShift
    figdir_base.ReplaceAll(tgt, tgt+"_6mmShift");
#endif
#if LVtxIsXi
    figdir_base.ReplaceAll(tgt, tgt+"_LVtxIsXi");
#endif
    if(test_run) figdir_base.ReplaceAll(tgt, tgt+"_testrun");
  }
  if(!figdir_base.EndsWith("/")) figdir_base += "/";
  cout<<"G4XiPScatAcceptanceCheck(i,frac,ndiv): standard Xi-style closure worker"<<endl;
  cout<<"CheckAcceptanceAll(nfile,ndiv): merge and draw (default 10,10)"<<endl;
  cout<<"Optional fourth worker argument: BestCorrections.tsv"<<endl;
}
void PrepareXiPScatCheck(TString best_file = ""){
  if(file_dir == "" or fmap_name == "" or figdir_base == "") G4XiPScatAcceptanceCheck();
  if(!figdir_base.EndsWith("/")) figdir_base += "/";
  InitializeTriggerCondtions();
  if(best_file != "") LoadBestCorrections(best_file);
  ResetXiPScatHistograms();
  over_cnt = 0;
  InitializeCorrectionHistograms(tgt);
  TFile* acpt_file = TFile::Open(fmap_dir+fmap_name);
  cout<<"Xi-p efficiency map file: "<<fmap_dir+fmap_name<<endl;
  if(acpt_file == nullptr or acpt_file->IsZombie()) throw runtime_error("Cannot open Xi-p efficiency maps.");
  TNamed* format = dynamic_cast<TNamed*>(acpt_file->Get("XiPScatMapFormat"));
  if(format == nullptr or TString(format->GetTitle()) != "XiPScat_StandardXi_20261006")
    throw runtime_error("Regenerate Xi-p maps using the revised macros.");
  LoadEfficiencies(acpt_file, tgt);
  acpt_file->Close(); delete acpt_file;
}
TString XiPScatCheckInputFile(int i){
  TString name = Form("XiPScatReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_XiPScat.root",i);
  return file_dir+(file_dir.EndsWith("/") ? "" : "/")+name;
}
void WriteXiPScatCheckWorker(TString name){
  TFile* out_file = TFile::Open(name, "RECREATE");
  if(out_file == nullptr or out_file->IsZombie()) throw runtime_error("Cannot create correction worker file.");
  TNamed parameters("CorrectionParameters", CorrectionParameters().Data());
  out_file->WriteTObject(&parameters);
  for(auto hist:hMap) out_file->WriteTObject(hist.second);
  out_file->Write(); out_file->Close(); delete out_file;
}
void G4XiPScatAcceptanceCheck(int i, int frac = 1, int ndiv = 1, TString best_file = ""){
  if(ndiv < 1 or frac < 1 or frac > ndiv) throw runtime_error("Invalid worker fraction.");
  PrepareXiPScatCheck(best_file);
  filename = XiPScatCheckInputFile(i);
  if(!FileCheck(filename)) throw runtime_error(string("Cannot read input: ")+filename.Data());
  TChain* tree = new TChain("tpc");
  tree->Add(filename);
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
  Long64_t ent = tree->GetEntries(), part = ent/ndiv;
  Long64_t start = part*(frac-1), last = frac == ndiv ? ent : part*frac;
  for(Long64_t entry=start;entry<last;++entry){
    if(entry%1000 == 0) cout<<entry<<endl;
    Xi->GetEntry(entry);
    FillHistograms(Xi, tgt);
  }
  gSystem->mkdir(bufdir, true);
  WriteXiPScatCheckWorker(Form("%sCorrectedHistograms_%d_%d_%d.root",bufdir.Data(),i,frac,ndiv));
}
void CheckAcceptance(vector<TFile*> files, bool make_figures){
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
  if(files.empty()) throw runtime_error("No correction workers supplied.");
  InitializeTriggerCondtions();
  ResetXiPScatHistograms();
  for(auto file: files){
    LoadCorrectionHistograms(file, tgt);
  }
  MakeChi2Map(tgt);
  TFile* out_file = TFile::Open(figdir_base + "CorrectionResults.root", "RECREATE");
  if(out_file == nullptr or out_file->IsZombie()) throw runtime_error("Cannot create CorrectionResults.root.");
  TNamed parameters("CorrectionParameters", CorrectionParameters().Data());
  out_file->WriteTObject(&parameters);
  for(auto hist:hMap) out_file->WriteTObject(hist.second); // Original merged counts.
  if(!make_figures){ out_file->Write(); out_file->Close(); delete out_file; return; }
  NormalizeHistograms(tgt);
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
        h_gen->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_gen->GetMaximum()));
        h_gen->Draw("hist");
        cout<<"GenMaxi = "<<h_gen->GetMaximum()<<endl;
        h_gen->SetLineColor(kBlack);
        TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
        leg->SetTextSize(0.03);
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
      SetFinalXiPPlotRange(h_gen, AcceptanceHistTitle1D(tgt, p, v1, ""));
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, Labels["Gen" + trig], "l");
#if DrawGeneValue
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflag"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, Labels["XiPflag"+trig]+ "wo_resol", "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same"); 
      leg->AddEntry(h_xiacptcor, TString("Acceptance-corrected") + " (MC kinematics)" , "pe");
#endif
#if DrawResSmeared
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagRec"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, Labels["XiPflag"+trig], "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagRecCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      leg->AddEntry(h_xirecacptcor, "Acceptance-corrected", "pe");
  #endif 
      key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiAndScatPGood"+trig);
      auto h_goodxi = (TH1*)hMap[key]->Clone();
      h_goodxi->Draw("hist same");
      leg->AddEntry(h_goodxi, Labels["GoodXiAndScatPGood"], "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiAndScatPGoodCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      h_goodxicor->Draw("pe same");
      leg->AddEntry(h_goodxicor, Labels["GoodXiAndScatPGoodCor"], "pe");
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
      SetFinalXiPPlotRange(h_gen, AcceptanceHistTitle1D(tgt, p, v1, ""));
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, Labels["Gen" + trig], "l");
#if DrawGeneValue
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflag"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, "XiPflag"+trig, "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same"); 
      leg->AddEntry(h_xiacptcor, "XiPflagCor"+trig, "pe");
      auto gr = MakeDivisionGraph(h_xiacptcor, h_gen);
      gr->SetMarkerColor(kAzure);
      gr->Draw("PE same");
      leg->AddEntry(gr, "Ratio", "pe");
#endif
#if DrawResSmeared
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagRec"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, "XiPflagRec"+trig, "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiPflagRecCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same"); 
      leg->AddEntry(h_xirecacptcor, "XiPflagRecCor"+trig, "pe");
      auto grrec = MakeDivisionGraph(h_xirecacptcor, h_gen);
      grrec->SetMarkerColor(kMagenta);
      grrec->Draw("PE same");
      leg->AddEntry(grrec, "Rec Ratio", "pe");
#endif
      key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiAndScatPGood"+trig);
      auto h_goodxi = (TH1*)hMap[key]->Clone();
      h_goodxi->Draw("hist same");
      leg->AddEntry(h_goodxi, Labels["GoodXiAndScatPGood"], "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiAndScatPGoodCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      h_goodxicor->Draw("pe same");
      leg->AddEntry(h_goodxicor, Labels["GoodXiAndScatPGoodCor"], "pe");
      auto grgood = MakeDivisionGraph(h_goodxicor, h_gen);
      grgood->SetMarkerColor(kGreen+2);
      grgood->Draw("PE same");
      leg->AddEntry(grgood, "MC-matched closure ratio", "pe");
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
        h_den->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_den->GetMaximum()));
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
        TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
        leg->SetTextSize(0.03);
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
          h_den->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_den->GetMaximum()));
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
          TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
          leg->SetTextSize(0.03);
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
  TString event_figdir = figdir_base+"EventVars/";
  for(auto ev: EventVars){
    cout<<"EventVar: "<<ev<<endl;
    //Each step by gen
    TString figdir = event_figdir  + "StepByGen/";
    gSystem->mkdir(figdir, true);
    TString ct = Form("Canv_%s", ev.Data());
    TCanvas* c0 = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen"+trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
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
    figdir = event_figdir  + "StepByStep/";
    gSystem->mkdir(figdir, true);
    ct = Form("CanvStepByStep_%s", ev.Data());
    TCanvas* c_sbs = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen"+trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
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
    TString figdir_rec = event_figdir + "Recon/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr", ev.Data());
    TCanvas* c_corr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      SetFinalXiPPlotRange(h_gen, EventTitle(tgt, ev, ""));
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
#if DrawGeneValue
      key = EventTitle(tgt, ev, "XiPflag"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, "XiPflag"+trig, "l");
      key = EventTitle(tgt, ev, "XiPflagCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same");
      leg->AddEntry(h_xiacptcor, "XiPflagCor"+trig, "pe");
#endif
#if DrawResSmeared
      key = EventTitle(tgt, ev, "XiPflagRec"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, "XiPflagRec"+trig, "l");
      key = EventTitle(tgt, ev, "XiPflagRecCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      leg->AddEntry(h_xirecacptcor, "XiPflagRecCor"+trig, "pe");
#endif
      key = EventTitle(tgt, ev, "GoodXiAndScatPGood"+trig);
      auto h_goodxi = (TH1*)hMap[key]->Clone();
      h_goodxi->Draw("hist same");
      leg->AddEntry(h_goodxi, Labels["GoodXiAndScatPGood"], "l");
      key = EventTitle(tgt, ev, "GoodXiAndScatPGoodCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      h_goodxicor->Draw("pe same");
      leg->AddEntry(h_goodxicor, Labels["GoodXiAndScatPGoodCor"], "pe");
      auto grgood = MakeDivisionGraph(h_goodxicor, h_gen);
      grgood->SetMarkerColor(kGreen+2);
      grgood->Draw("PE same");
      leg->AddEntry(grgood, "MC-matched closure ratio", "pe");
      leg->Draw();
    }
    c_corr->SaveAs(figdir_rec + ct + ".pdf");

    //Recon with final correction and ratio
    figdir_rec = event_figdir + "Recon_w_rat/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr_w_rat", ev.Data());
    TCanvas* c_corr_rat = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    {
      TString key = EventTitle(tgt, ev, "Gen" + trig);
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      SetFinalXiPPlotRange(h_gen, EventTitle(tgt, ev, ""));
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
      leg->SetTextSize(0.03);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen" + trig, "l");
#if DrawGeneValue
      key = EventTitle(tgt, ev, "XiPflag"+trig);
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      leg->AddEntry(h_xiacpt, "XiPflag"+trig, "l");
      key = EventTitle(tgt, ev, "XiPflagCor"+trig);
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same");
      leg->AddEntry(h_xiacptcor, "XiPflagCor"+trig, "pe");
      auto gr = MakeDivisionGraph(h_xiacptcor, h_gen);
      gr->SetMarkerColor(kAzure);
      gr->Draw("PE same");
      leg->AddEntry(gr, "Ratio", "pe");
#endif
#if DrawResSmeared
      key = EventTitle(tgt, ev, "XiPflagRec"+trig);
      auto h_xirecacpt = (TH1*)hMap[key]->Clone();
      h_xirecacpt ->Draw("hist same");
      leg->AddEntry(h_xirecacpt, "XiPflagRec"+trig, "l");
      key = EventTitle(tgt, ev, "XiPflagRecCor"+trig);
      auto h_xirecacptcor = (TH1*)hMap[key]->Clone();
      h_xirecacptcor->Draw("pe same");
      leg->AddEntry(h_xirecacptcor, "XiPflagRecCor"+trig, "pe");
      auto grrec = MakeDivisionGraph(h_xirecacptcor, h_gen);
      grrec->SetMarkerColor(kMagenta);
      grrec->Draw("PE same");
      leg->AddEntry(grrec, "Rec Ratio", "pe");
#endif
      key = EventTitle(tgt, ev, "GoodXiAndScatPGood"+trig);
      auto h_goodxi = (TH1*)hMap[key]->Clone();
      h_goodxi->Draw("hist same");
      leg->AddEntry(h_goodxi, Labels["GoodXiAndScatPGood"], "l");
      key = EventTitle(tgt, ev, "GoodXiAndScatPGoodCor"+trig);
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      h_goodxicor->Draw("pe same");
      leg->AddEntry(h_goodxicor, Labels["GoodXiAndScatPGoodCor"], "pe");
      auto grgood = MakeDivisionGraph(h_goodxicor, h_gen);
      grgood->SetMarkerColor(kGreen+2);
      grgood->Draw("PE same");
      leg->AddEntry(grgood, "MC-matched closure ratio", "pe");
      leg->Draw();
#if DrawGeneValue || DrawResSmeared
      if(ev.Contains("Pol")){
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.04);
#if DrawGeneValue
        double mean_rat_gen = 0, std_rat_gen = 0;
        int np_gen = gr->GetN();
        for(int ip = 0; ip < np_gen; ++ip){
          double x, y;
          gr->GetPoint(ip, x, y);
          mean_rat_gen += y;
          std_rat_gen += y*y;
        }
        if(np_gen == 0){
          cout<<"No points in truth ratio for "<<ev<<endl;
          np_gen = 1;
        }
        mean_rat_gen /= np_gen;
        std_rat_gen = sqrt(max(0., std_rat_gen/np_gen - mean_rat_gen*mean_rat_gen));
        tex->SetTextColor(kAzure);
        tex->DrawLatex(0.2, 0.3, Form("Gen Ratio: %.3f #pm %.5f", mean_rat_gen, std_rat_gen));
#endif
#if DrawResSmeared
        double mean_rat_cor = 0, std_rat_cor = 0;
        int np_rec = grrec->GetN();
        for(int ip = 0; ip < np_rec; ++ip){
          double x, y;
          grrec->GetPoint(ip, x, y);
          mean_rat_cor += y;
          std_rat_cor += y*y;
        }
        if(np_rec == 0){
          cout<<"No points in reconstructed ratio for "<<ev<<endl;
          np_rec = 1;
        }
        mean_rat_cor /= np_rec;
        std_rat_cor = sqrt(max(0., std_rat_cor/np_rec - mean_rat_cor*mean_rat_cor));
        tex->SetTextColor(kMagenta);
        tex->DrawLatex(0.2, 0.2, Form("Rec Ratio: %.3f #pm %.5f", mean_rat_cor, std_rat_cor));
#endif
      }
#endif
      // ROOT declares cd() with a null default path; pass an explicit empty
      // path to select this file while avoiding -Wnonnull from that default.
      out_file->cd("");
#if DrawGeneValue
      gr->Write(Form("Graph_%s_Corr_w_rat", ev.Data()));
#endif
#if DrawResSmeared
      grrec->Write(Form("Graph_%s_RecCorr_w_rat", ev.Data()));
#endif
      
    }
    c_corr_rat->SaveAs(figdir_rec + ct + ".pdf");
  }//EventVar

  TString weight_figdir = figdir_base+"Weights/";
  gSystem->mkdir(weight_figdir, true);
  vector<TString> weight_config = {"total"};
  vector<TString> whist_config = {"h_weight", "h_weightRec", "hgood_weight", "hgood_weightRec","hbad_weight", "hbad_weightRec"};
  vector<TString> wRathist_config = {"h_weightRatio",  "hgood_weightRatio","hbad_weightRatio"};
  vector<TString> wDiffhist_config = {"h_weightDiff", "hgood_weightDiff","hbad_weightDiff"};
  for(auto cp:CorrPars){
    weight_config.push_back(cp.p_cor);
  }
  for(auto wconf: weight_config){
    TString wconf_dir = weight_figdir + "/Weights/";
    gSystem->mkdir(wconf_dir, true);
    vector<int> whist_colors ={
      kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange+7
    };
    TString ct = Form("Canv_Weights_%s", wconf.Data());
    TCanvas* c_w = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    TString tag = wconf == "total" ? "" : "_"+wconf;
    TLegend* leg = new TLegend(0.15,0.58,0.95,0.94);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    TLegend* leg_rat = new TLegend(0.15,0.58,0.95,0.94);
    leg_rat->SetTextSize(0.03);
    leg_rat->SetFillStyle(0);
    leg_rat->SetBorderSize(0);
    for(int iw = 0; iw < whist_config.size(); ++iw){
      TString whist = whist_config[iw];
      TString key = whist + tag;
      TH1* h_w = (TH1*)hMap[key]->Clone();
      h_w->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_w->GetMaximum()));
      h_w->Draw("hist same");
      h_w->SetLineColor(whist_colors[iw]);
      leg->AddEntry(h_w, whist, "l");
    }
    leg->Draw();
    c_w->SaveAs(wconf_dir + ct + ".pdf");
    wconf_dir = weight_figdir + "/WeightRatio/";
    gSystem->mkdir(wconf_dir, true);
    ct = Form("Canv_WeightRatio_%s", wconf.Data());
    TCanvas* c_wr = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    for(int iw = 0; iw < wRathist_config.size(); ++iw){
      TString whist = wRathist_config[iw];
      TString key = whist + tag;
      TH1* h_w = (TH1*)hMap[key]->Clone();
      leg_rat->AddEntry(h_w, whist, "l");
      h_w->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_w->GetMaximum()));
      h_w->Draw("hist same");
      h_w->SetLineColor(whist_colors[iw]);
    }
    leg_rat->Draw();
    c_wr->SaveAs(wconf_dir + ct + ".pdf");

    wconf_dir = weight_figdir + "/WeightDiff/";
    gSystem->mkdir(wconf_dir, true);
    ct = Form("Canv_WeightDiff_%s", wconf.Data());
    TCanvas* c_wd = new TCanvas(ct, ct, canv_size_x, canv_size_y);
    for(int iw = 0; iw < wDiffhist_config.size(); ++iw){
      TString whist = wDiffhist_config[iw];
      TString key = whist + tag;
      TH1* h_w = (TH1*)hMap[key]->Clone();
      h_w->GetYaxis()->SetRangeUser(0,max(1.5,1.5*h_w->GetMaximum()));
      h_w->Draw("hist same");
      h_w->SetLineColor(whist_colors[iw]);
    }
    leg_rat->Draw();
    c_wd->SaveAs(wconf_dir + ct + ".pdf");
  }


  TDirectory* display = out_file->mkdir("NormalizedPlots");
  for(auto& [key, hist]: hMap){
    display->WriteTObject(hist); 
  }
  out_file->Write();
  out_file->Close();
  delete out_file;
}
void CheckAcceptanceAll(int nfile = 10, int ndiv = 10, TString best_file = "", bool make_figures = true){
  if(figdir_base == "" or bufdir == "") G4XiPScatAcceptanceCheck();
  InitializeTriggerCondtions();
  if(best_file != "") LoadBestCorrections(best_file);
  if(test_run) nfile = 1;
  vector<TFile*> files;
  for(int i=0;i<nfile;++i){
    for(int frac=1;frac<=ndiv;++frac){
      TString name = Form("%sCorrectedHistograms_%d_%d_%d.root",bufdir.Data(),i,frac,ndiv);
      TFile* file = TFile::Open(name);
      if(file == nullptr or file->IsZombie()) throw runtime_error(string("Cannot read worker: ")+name.Data());
      files.push_back(file);
    }
  }
  CheckAcceptance(files, make_figures);
  for(auto file:files){ file->Close(); delete file; }
}
void G4XiPScatAcceptanceCheckAll(int nfile = 10, TString best_file = "", bool make_figures = true){
  PrepareXiPScatCheck(best_file);
  TChain* tree = new TChain("tpc");
  for(int i=0;i<nfile;++i){
    TString name = XiPScatCheckInputFile(i);
    if(!FileCheck(name)) throw runtime_error(string("Cannot read input: ")+name.Data());
    tree->Add(name);
  }
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
  Long64_t ent = tree->GetEntries();
  if(test_run) ent /= 10;
  for(Long64_t entry=0;entry<ent;++entry){ Xi->GetEntry(entry); FillHistograms(Xi, tgt); }
  gSystem->mkdir(bufdir, true);
  TString output = bufdir+"CorrectedHistograms_All.root";
  WriteXiPScatCheckWorker(output);
  TFile* file = TFile::Open(output);
  CheckAcceptance({file}, make_figures);
  file->Close(); delete file;
}
