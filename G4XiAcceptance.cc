#define LooseBin 0
#define Debug 0
#define date 260729
#define LVtxIsXi 1
#define SavePDF 0
#include "G4XiAcceptance.hh"
#define Recon 1
#define CH2 1
#define PosShift 1
void G4XiAcceptance(int conf=0){
	gErrorIgnoreLevel = kWarning;
  gStyle -> SetOptStat(0);
//  double pxi = -0.3;
  double pxi = 0;
  int nfiles = 30;
	TString WAcc,Target,Conf,LVtxConf;
#if LVtxIsXi
    LVtxConf = "_LVtxIsXi";
#endif
  WAcc = "_WB";
#if CH2
  TString file_dir = "./rootfiles/Geant4CH2/W_Acc/";
	Conf = "CH2";
	Target = "CH2";
#else
  TString file_dir = "./rootfiles/Geant4Prod/W_Acc/";
	Conf = "Prod";
	Target = "Carbon";
#endif
#if PosShift
#if CH2
  file_dir = "rootfiles/Geant4CH2/6mmShift/";
  Conf = "CH26mmShift";
#else
  file_dir = "rootfiles/Geant4Prod/6mmShift/";
  Conf = "Prod6mmShift";
#endif
#endif
  TString filename, figdir;
  TChain* tree = new TChain("tpc");
  bool file_check = 1;
	for(int i=0;i<nfiles;++i){
		filename = Form("XiRecon%s_P_E42_%d_GenfitCarbonGeant4Ver5.root",Conf.Data(),i);
#if PosShift
#if CH2
		filename = Form("XiReconCH2_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
#else
		filename = Form("XiReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
#endif
#endif
		cout<<"Loading "<<filename<<endl;
    if(FileCheck(file_dir+filename) == 0){
      continue;
      file_check = 0;
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
//	double ent_proton = tree->GetEntries();
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
	InitializeHistograms();
	auto ent_tot = tree->GetEntries();
  #if Debug
    ent_tot = 1000;
  #endif
  for(auto i = 0;i<ent_tot;i++){
    if(i%1000==0) cout << i << endl;
    Xi->GetEntry(i);
    FillHist(Xi);
  }

  cout<<"Closing files..."<<endl;
  cout<<"RunProceed. Making Eff Maps..."<<endl;
  MakeEfficiencies();
  cout<<"Maps made"<<endl;

  int ic=0;
  TFile* fileOut;
  TString fout_dir = "./Maps/";
#if Recon
	figdir = Form("./Maps_%d/AcceptanceHists/ReconPol_%s%s%s/", date,Conf.Data(),LVtxConf.Data(),WAcc.Data());
#else
  figdir = Form("./Maps_%d/AcceptanceHists/MissPol_%s%s/", date,Target.Data(),WAcc.Data());
#endif
  TString fout_name = Form("%s_MM_%s%s_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),date);
#if Recon
	fout_name = Form("%s_ReconPE42_%s%s_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),date);
#else
#endif
#if PosShift
  fout_name = Form("%s_ReconP0_%s%s_6mmShifted_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),date);
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
    if(SavePDF) c -> SaveAs(figdir + hist.second -> GetName() + TString(".pdf"));
    ic++;
    hist.second -> Write();
  }
  cout<<"1D hists done"<<endl;
  for(auto hist:hMap2D){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    hist.second -> Draw("col");
    TString fig2d = "2D/";
    if(SavePDF) c -> SaveAs(figdir + fig2d + hist.second -> GetName() + TString(".pdf"));
    ic++;
    hist.second -> Write();
  }
  cout<<"2D hists done"<<endl;
  for(auto eff:effMap){
    TCanvas* c = new TCanvas(Form("c%d", ic), Form("c%d", ic), 800, 800);
    TString fig2d = "2D/";
    eff.second -> Draw("colz");
    if(SavePDF){
      if(eff.second -> GetDimension() == 1){
        c -> SaveAs(figdir + eff.second -> GetName() + TString(".pdf"));
      }
      else{
        c -> SaveAs(figdir + fig2d + eff.second -> GetName() + TString(".pdf"));
      }
    }
    ic++;
    eff.second -> Write();
  }
  cout<<"All hists and effs written to file "<< fout_name << endl;
  fileOut->Write();
  cout<<"All done. Closing file."<<endl;
  fileOut->Close();
  exit(0);
}
