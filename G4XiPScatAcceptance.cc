#define LooseBin 0
#define Debug 0
#define date 261006
#define LVtxIsXi 1
#define SavePDF 0
#define Recon 1
#define CH2 0
#define PosShift 1
#include "G4XiPScatAcceptance.hh"
gErrorIgnoreLevel = kWarning;
TString WAcc,Target,Conf,LVtxConf,filename,file_dir,fout_dir,fout_name,figdir2d;
gStyle -> SetOptStat(0);
bool file_check = 1;
void G4XiPScatAcceptance(){
  LVtxConf = LVtxIsXi ? "_LVtxIsXi" : "";
  WAcc = "_WB";
  Conf = "ProdXiPScat6mmShift";
  Target = "Carbon";
  if(file_dir == "") file_dir = "rootfiles/Geant4Prod/XiPScat/";
  if(fout_dir == "") fout_dir = "./Maps/";
  if(!fout_dir.EndsWith("/")) fout_dir += "/";
  fout_name = Form("%sXiPScat_Recon_P0_%s%s_6mmShifted_%d.root",Target.Data(),LVtxConf.Data(),WAcc.Data(),date);
#if LooseBin
  fout_name.ReplaceAll(".root", "_LooseBin.root");
#endif
  figdir2d = figdir+"2D/";
  gSystem->mkdir(fout_dir, true);
  cout<<"G4XiPScatAcceptance(i,frac,ndiv): make worker histograms"<<endl;
  cout<<"G4MakeEfficiencies(nfile,ndiv): merge workers and make maps (default 10,10)"<<endl;
}
TString XiPScatInputFile(int i){
  TString name = Form("XiPScatReconProd_P_0_%d_6mmShifted_GenfitCarbonGeant4Ver16_XiPScat.root",i);
  return file_dir + (file_dir.EndsWith("/") ? "" : "/") + name;
}
void WriteXiPScatMaps(bool include_histograms = false){
  TFile* fileOut = TFile::Open(fout_dir+fout_name, "RECREATE");
  if(fileOut == nullptr or fileOut->IsZombie()) throw runtime_error("Cannot create Xi-p map file.");
  TNamed format("XiPScatMapFormat", "XiPScat_StandardXi_20261006");
  fileOut->WriteTObject(&format);
  if(include_histograms){
    for(auto hist:hMap) fileOut->WriteTObject(hist.second);
    for(auto hist:hMap2D) fileOut->WriteTObject(hist.second);
  }
  for(auto eff:effMap) fileOut->WriteTObject(eff.second);
  fileOut->Write(); fileOut->Close(); delete fileOut;
  cout<<"Efficiency maps written to "<<fout_dir+fout_name<<endl;
}
void G4XiPScatAcceptance(int i, int frac = 1, int ndiv = 1){
  if(fout_dir == "" or fout_name == "" or Conf == "") G4XiPScatAcceptance();
  if(ndiv < 1 or frac < 1 or frac > ndiv) throw runtime_error("Invalid worker fraction.");
  filename = XiPScatInputFile(i);
  if(!FileCheck(filename)) throw runtime_error(string("Cannot read input: ")+filename.Data());
  ResetXiPScatHistograms();
  TChain* tree = new TChain("tpc");
  tree->Add(filename);
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
  InitializeHistograms();
  Long64_t ent = tree->GetEntries(), part = ent/ndiv;
  Long64_t start = part*(frac-1), last = frac == ndiv ? ent : part*frac;
  for(Long64_t entry=start;entry<last;++entry){
    if(entry%1000 == 0) cout<<entry<<endl;
    Xi->GetEntry(entry);
    FillHist(Xi);
  }
  TString buf = fout_dir+Form("buf/%s/%d/",Conf.Data(),date);
  gSystem->mkdir(buf, true);
  TString output = fout_name;
  output.ReplaceAll(".root", Form("_%d_%d_%d.root",i,frac,ndiv));
  TFile* fileOut = TFile::Open(buf+output, "RECREATE");
  if(fileOut == nullptr or fileOut->IsZombie()) throw runtime_error("Cannot create histogram worker file.");
  TNamed format("XiPScatMapFormat", "XiPScat_StandardXi_20261006");
  fileOut->WriteTObject(&format);
  for(auto hist:hMap) fileOut->WriteTObject(hist.second);
  for(auto hist:hMap2D) fileOut->WriteTObject(hist.second);
  fileOut->Write(); fileOut->Close(); delete fileOut;
}
void G4MakeEfficiencies(vector<TString> files){
  if(fout_dir == "" or fout_name == "" or Conf == "") G4XiPScatAcceptance();
  if(files.empty()) throw runtime_error("No histogram workers supplied.");
  ResetXiPScatHistograms();
  for(auto name:files){
    TFile* file = TFile::Open(name);
    LoadHistograms(file);
    file->Close(); delete file;
  }
  MakeEfficiencies();
  WriteXiPScatMaps();
}
void G4MakeEfficiencies(int nfile = 10, int ndiv = 10){
  if(fout_dir == "" or fout_name == "" or Conf == "") G4XiPScatAcceptance();
  vector<TString> files;
  TString buf = fout_dir+Form("buf/%s/%d/",Conf.Data(),date);
  for(int i=0;i<nfile;++i){
    for(int frac=1;frac<=ndiv;++frac){
      TString name = fout_name;
      name.ReplaceAll(".root", Form("_%d_%d_%d.root",i,frac,ndiv));
      files.push_back(buf+name);
    }
  }
  G4MakeEfficiencies(files);
}
void G4MakeEfficiency(int nfile = 10, int ndiv = 10){ G4MakeEfficiencies(nfile, ndiv); }
void G4XiPScatAcceptanceAll(int nfile = 10){
  if(fout_dir == "" or fout_name == "" or Conf == "") G4XiPScatAcceptance();
  ResetXiPScatHistograms();
  TChain* tree = new TChain("tpc");
  for(int i=0;i<nfile;++i){
    TString name = XiPScatInputFile(i);
    if(!FileCheck(name)) throw runtime_error(string("Cannot read input: ")+name.Data());
    tree->Add(name);
  }
  g4genfitcarbon* Xi = new g4genfitcarbon(tree);
  SetBranches(tree);
  InitializeHistograms();
  for(Long64_t entry=0;entry<tree->GetEntries();++entry){ Xi->GetEntry(entry); FillHist(Xi); }
  MakeEfficiencies();
  WriteXiPScatMaps(true);
}
