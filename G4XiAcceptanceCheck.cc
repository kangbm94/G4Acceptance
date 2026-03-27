#define PerDist 1
#include "G4XiAcceptanceCheck.hh"
#define Recon 1
#define CH2 1
#include <TCanvas.h>
void G4XiAcceptanceCheck(){
  SetStyle();
  double pxi = 0;
	TString WAcc,Target,Conf;
  WAcc = "_WB";
  int nfile = 10;
  bool test_run = 0;
  if(test_run) nfile = 1;
#if CH2
  TString file_dir = "./rootfiles/Geant4CH2/";
	Conf = "CH2";
	Target = "CH2";
#else
  TString file_dir = "./rootfiles/Geant4Prod/";
	Conf = "Prod";
	Target = "Carbon";
#endif
  //file_dir += "WeightedFermi/";
  file_dir += "W_Acc/";
	int date = 260304;
  TString filename, figdir;
  TString tgt = "Carbon";
#if CH2
	tgt = "CH2";
#endif
  TChain* tree = new TChain("tpc");
	cout<<"Loading files..."<<endl;
	for(int j=0;j<3;++j){
		if(j==0) pxi=-1;
		if(j==1) pxi=0;
		if(j==2) pxi=1;
		if(j!=1) continue;
//    if(test_run and (pxi == -1 or pxi == 1)) continue;
	for(int i=0;i<nfile;++i){
#if CH2
		filename = Form("XiReconCH2_P_E42_%d_GenfitCarbonGeant4Ver5.root",i);
#else
		filename = Form("XiReconProd_P_E42_%d_GenfitCarbonGeant4Ver5.root",i);
#endif
		cout<<"Loading "<<filename<<endl;
		tree -> Add(file_dir+filename);
	}
	}
	cout<<"Accpt file"<<endl;
#if CH2
  TFile* acpt_file = TFile::Open("./Maps/CH2_ReconPE42__WB_260304.root");
#else
  TFile* acpt_file = TFile::Open("./Maps/Carbon_ReconPE42__WB_260227.root");
#endif
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
	InitializeCorrectionHistograms( tgt );
  LoadEfficiencies(acpt_file, tgt);
  int ent = tree->GetEntries();
	cout<<"Entries : "<<ent<<endl;
//  if(test_run) ent = ent / 10;
  for(int i = 0;i<ent;i++){
    if(i%1000==0) cout << i << endl;
    Xi->GetEntry(i);
    FillHistograms(Xi, tgt);
  }
  MakeChi2Map(tgt);
  NormalizeHistograms(tgt);
  TString figdir_base = Form("figs/%s/", tgt.Data());
	if(test_run) figdir_base.ReplaceAll(tgt,tgt + "_testrun");
  figdir_base.ReplaceAll(tgt,tgt+ "ForcedPi2Ph");
  gSystem->mkdir(figdir_base, true);
  TFile* out_file = TFile::Open(figdir_base + "AcceptanceCorrectionMaps.root", "RECREATE");
  out_file->cd();
  for(auto& [key, hist]: hMap){
    hist->Write();
  }
  out_file->Write();
  for(auto p:particle){
    TString figdir = figdir_base  + p + "/";
    gSystem->mkdir(figdir, true);
    TString ct = Form("Canv_%s", p.Data());
    TCanvas* c = new TCanvas(ct, ct, 1600, 1200);
    c->Divide(3,2);
    for(int iv = 0; iv < variable.size();++iv){
      auto v1 = variable[iv];
      c->cd(iv+1);
      TString key = AcceptanceHistTitle1D(tgt, p, v1, "Gen");
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0,h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen", "l");
      for(auto chk: CheckLists){
        if(chk == "Gen") continue;
        key = AcceptanceHistTitle1D(tgt, p, v1, chk);
        auto gr = MakeDivisionGraph(hMap[key], h_gen);
        gr->Draw("PE same");
        leg->AddEntry(gr, chk, "pe");
      }
      leg->Draw();
    }
    c->SaveAs(figdir + ct + ".pdf");

    TString figdir_rec = figdir_base + "/Recon/" + p + "/";
    gSystem->mkdir(figdir_rec, true);
    ct = Form("Canv_%s_FinalCorr", p.Data());
    TCanvas* c_corr = new TCanvas(ct, ct, 1600, 1200);
    c_corr->Divide(3,2);
    for(int iv = 0; iv < variable.size();++iv){
      auto v1 = variable[iv];
      c_corr->cd(iv+1);
      TString key = AcceptanceHistTitle1D(tgt, p, v1, "Gen");
      TH1* h_gen = (TH1*)hMap[key]->Clone();
      h_gen->GetYaxis()->SetRangeUser(0, h_gen->GetMaximum()*1.5);
      h_gen->Draw("hist");
      TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_gen, "Gen", "l");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcpt");
      auto h_xiacpt = (TH1*)hMap[key]->Clone();
      h_xiacpt ->Draw("hist same");
      key = AcceptanceHistTitle1D(tgt, p, v1, "GoodXiCor");
      auto h_goodxicor = (TH1*)hMap[key]->Clone();
      //h_goodxicor->Draw("pe same");
      key = AcceptanceHistTitle1D(tgt, p, v1, "XiAcptCor");
      auto h_xiacptcor = (TH1*)hMap[key]->Clone();
      h_xiacptcor->Draw("pe same"); 
      leg->AddEntry(h_xiacpt, "XiAcpt", "l");
      //leg->AddEntry(h_goodxicor, "GoodXiCor", "pe");
      leg->AddEntry(h_xiacptcor, "XiAcptCor", "pe");
      leg->Draw();
    }
    c_corr->SaveAs(figdir_rec + ct + ".pdf");

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
      TCanvas* c_corr = new TCanvas(ct, ct, 1600, 1200);
      c_corr->Divide(3,2);
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
        TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(h_den, den, "l");
        leg->AddEntry(h_num, num, "l");
        leg->AddEntry(gr, Form("%s/%s cor w #chi^{2} %.4g", num.Data(), den.Data(), chi2), "pe");
        leg->Draw();
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
        TCanvas* c_corr = new TCanvas(ct, ct, 1600, 1200);
        c_corr->Divide(3,2);
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
          TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
          leg->SetFillStyle(0);
          leg->SetBorderSize(0);
          leg->AddEntry(h_den, den, "l");
          leg->AddEntry(h_num, num, "l");
          leg->AddEntry(gr, Form("%s/%s cor w #chi^{2} %.4g", num.Data(), den.Data(), chi2), "pe");
          leg->Draw();
        }
        c_corr->SaveAs(figdir + ct + ".pdf");
      }
    }
  }
}
