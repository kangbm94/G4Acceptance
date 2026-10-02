#define LooseBin 0
#define Debug 0
#define date 261001
#define LVtxIsXi 1
#define SavePDF 0
#include "G4XiPScatAcceptance.hh"
#define Recon 1
#define CH2 0
#define PosShift 1
gErrorIgnoreLevel = kWarning;
TString WAcc,Target,Conf,LVtxConf,filename,file_dir,fout_dir,fout_name,figdir2d;
gStyle -> SetOptStat(0);
bool file_check = 1;
void G4XiPScatAcceptance(){
  LVtxConf = "_LVtxIsXi";
  WAcc = "_WB";
#if CH2
  file_dir = "rootfiles/Geant4CH2/6mmShift/wo_upstream_search/";
  Conf = "CH26mmShift";
	Target = "CH2";
#else
  file_dir = "rootfiles/Geant4Prod/XiPScat/";
  Conf = "ProdXiPScat6mmShift";
	Target = "Carbon";
#endif
  fout_dir = "./Maps/";
  fout_name = Form("%sXiPScat_Recon_P0_%s%s_6mmShifted_%d.root",Target.Data(), LVtxConf.Data(),WAcc.Data(),date);
#if LooseBin
  fout_name = fout_name.ReplaceAll(".root","_LooseBin.root");
	figdir.ReplaceAll("Hists","Hists_LooseBin");
#endif
	figdir2d = figdir + "2D/";
  gSystem->mkdir(figdir,1);
  gSystem->mkdir(figdir2d,1);
  cout<<"=======Usage======"<<endl;
  cout<<"G4XiPScatAcceptance(i,frac,ndiv) : Run the analysis for file index i, with fraction frac of ndiv parts"<<endl;
  cout<<"G4MakeEfficiencies() : Make the efficiency maps from the histograms"<<endl;
}
void G4XiPScatAcceptance(int i, int frac = 1, int ndiv = 1){
  if(fout_dir == "") G4XiPScatAcceptance();
  TChain* tree = new TChain("tpc");
#if PosShift
#if CH2
  filename = Form("XiReconCH2_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_VtxFit.root",i);
#else
  filename = Form("XiPScatReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_XiPScat.root",i);
#endif
#endif
  cout<<"Loading "<<filename<<endl;
  if(FileCheck(file_dir+filename) == 0){
    cout<<"Error loading file "<<filename<<endl;
    return;
    file_check = 0;
  }
  if(tree -> Add(file_dir+filename) != 1){
    cout<<"Error adding file "<<filename<<endl;
    return;
    file_check = 0;
  }
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
	InitializeHistograms();
	auto ent = tree->GetEntries();
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
    FillHist(Xi);
  }
  TFile* fileOut;
  TString buf = fout_dir + Form("buf/%s/%d/",Conf.Data(),date);
  gSystem->mkdir(buf,1);
  TString file_parallel = fout_name;
  file_parallel.ReplaceAll(".root",Form("_%d_%d_%d.root",i,frac,ndiv));
  fileOut = new TFile(buf+file_parallel,"recreate");
  fileOut->cd();
  for(auto hist:hMap){
    hist.second -> Write();
  }
  cout<<"1D hists done"<<endl;
  for(auto hist:hMap2D){
    hist.second -> Write();
  }
  cout<<"2D hists done"<<endl;
  fileOut->Write();
  fileOut->Close();
}
void G4MakeEfficiencies(vector<TString> files){
  for(auto f:files){
    TFile* file = new TFile(f);
    LoadHistograms(file);
  }
  MakeEfficiencies();
  TFile* fileOut;
  fileOut = new TFile(fout_dir+fout_name,"recreate");
  fileOut->cd();
  for(auto eff:effMap){
    eff.second -> Write();  
  }
  cout<<"All effs written to file "<< fout_name << endl;
  fileOut->Write();
  cout<<"All done. Closing file."<<endl;
  fileOut->Close();
}
void G4MakeEfficiencies(){
  if(fout_dir == "") G4XiPScatAcceptance();
  TString buf = fout_dir + Form("buf/%s/%d/",Conf.Data(),date);
  gSystem->mkdir(buf,1);
  TString file_parallel = fout_name;
  int nfile = 10;
  int ndiv = 10;
  vector<TFile*> files;
  for(int i=0;i<nfile;++i){
    for(int frac = 1; frac <= ndiv; ++frac){
      file_parallel = fout_name;
      file_parallel.ReplaceAll(".root",Form("_%d_%d_%d.root",i,frac,ndiv));
      TFile* file = new TFile(buf+file_parallel);
      files.push_back(file);
    }
  }
  for(auto f:files){
    LoadHistograms(f);
  }
  MakeEfficiencies();
  TFile* fileOut;
  fileOut = new TFile(fout_dir+fout_name,"recreate");
  fileOut->cd();
  for(auto eff:effMap){
    eff.second -> Write();  
  }
  cout<<"All effs written to file "<< fout_name << endl;
  fileOut->Write();
  cout<<"All done. Closing file."<<endl;
  fileOut->Close();
}
void G4XiPScatAcceptanceAll(){
//  double pxi = -0.3;
  if(fout_dir == "") G4XiPScatAcceptance();
  double pxi = 0;
  int nfiles = 30;
  TChain* tree = new TChain("tpc");
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
  fileOut = new TFile(fout_dir+fout_name,"recreate");
	pxi=0;
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
