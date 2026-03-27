#define LooseBin 0
#define Debug 0
#include "G4XiAcceptance.hh"
#define Recon 1
#define CH2 1
void G4XiAcceptance(int conf=0){

  gStyle -> SetOptStat(0);
//  double pxi = -0.3;
  double pxi = 0;
  int nfiles = 10;
	TString WAcc,Target,Conf;
  WAcc = "_WB";
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
  //file_dir += "W_Acc_Rot/";
	int date = 60327;
  TString filename, figdir;
  TChain* tree = new TChain("tpc");
//	filename = Form("%d_KpXiMM%s%s_P_%g_%d_GenfitCarbonGeant4Ver4.root",date,Conf.Data(),WAcc.Data(),pxi,0);
	for(int i=0;i<nfiles;++i){
		filename = Form("XiRecon%s_P_E42_%d_GenfitCarbonGeant4Ver5.root",Conf.Data(),i);
		cout<<"Loading "<<filename<<endl;
		tree -> Add(file_dir+filename);
	}
//	double ent_proton = tree->GetEntries();
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
	InitializeHistograms();
	double ent_tot = tree->GetEntries();
  #if Debug
    ent_tot = 1000;
  #endif
  for(int i = 0;i<ent_tot;i++){
    if(i%1000==0) cout << i << endl;
    Xi->GetEntry(i);
    FillHist(Xi);
  }
  MakeEfficiencies();

  int ic=0;
  TFile* fileOut;
  TString fout_dir = "./Maps/";
#if Recon
	figdir = Form("./figs/AcceptanceHists/ReconPol_%g_%s%s/", pxi,Target.Data(),WAcc.Data());
#else
  figdir = Form("./figs/AcceptanceHists/MissPol_%g_%s%s/", pxi,Target.Data(),WAcc.Data());
#endif
  TString fout_name = Form("%s_MM_%s_%d.root",Target.Data() ,WAcc.Data(),date);
#if Recon
	fout_name = Form("%s_ReconPE42_%s_%d.root",Target.Data() ,WAcc.Data(),date);
#else
#endif
#if LooseBin
  fout_name = fout_name.ReplaceAll(".root","_LooseBin.root");
	figdir.ReplaceAll("Hists","Hists_LooseBin");
#endif
	TString figdir2d = figdir + "2D/";
  gSystem->mkdir(figdir,1);
  gSystem->mkdir(figdir2d,1);
	pxi=0;
  fileOut = new TFile(fout_dir+fout_name,"recreate");
  fileOut->cd();
  for(auto hist:hMap){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    hist.second -> Draw("col");
    c -> SaveAs(figdir + hist.second -> GetName() + TString(".pdf"));
    ic++;
    hist.second -> Write();
  }
  for(auto hist:hMap2D){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    hist.second -> Draw("col");
    TString fig2d = "2D/";
    c -> SaveAs(figdir + fig2d + hist.second -> GetName() + TString(".pdf"));
    ic++;
    hist.second -> Write();
  }
  for(auto eff:effMap){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    eff.second -> Draw("colz");
    c -> SaveAs(figdir + eff.second -> GetName() + TString(".pdf"));
    ic++;
    eff.second -> Write();
  }
  for(auto eff:effMap){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    eff.second -> Draw("colz");
    TString fig2d = "2D/";
    c -> SaveAs(figdir + fig2d + eff.second -> GetName() + TString(".pdf"));
    ic++;
    eff.second -> Write();
  }
}
