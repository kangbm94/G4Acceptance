void SetStyle(){
  gStyle -> SetPalette(kRainBow);
//  gStyle -> SetPalette(kDarkBodyRadiator);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetMarkerStyle(24);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
 	gStyle->SetStatW(0.3);
 	gStyle->SetStatH(0.3);
	
	gStyle->SetTitleFontSize(0.06);
	gStyle->SetTitleFont(132,"x");//,“t”);
	gStyle->SetTitleFont(132,"t");//,“t”);
	gStyle->SetTitleFont(132,"y");//,“t”);
	gStyle->SetTitleOffset(0.9,"xy");
	gStyle->SetNdivisions(5);
	gStyle->SetNdivisions(5,"Y");
	gStyle->SetTitleSize(0.06,"X");
	gStyle->SetTitleSize(0.06,"Y");
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
//  gStyle -> SetLabelOffset(0.02,"xyz");
//	gStyle->SetTitleOffset(0.80,"X");
//	gStyle->
/*
	gStyle->SetStatBorderSize(0);	
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.1);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.15);
*/
  gStyle -> SetPadLeftMargin(0.15);
  gStyle -> SetPadBottomMargin(0.15);
	//gStyle -> SetPadRightMargin(0.15);
	gROOT->ForceStyle();
	/*
	gStyle->SetStatStyle(0);	
	gStyle->SetTitleOffset(0.8,"X");
	gStyle->cd();
  gStyle -> SetCanvasColor(0);
  gStyle -> SetFrameBorderMode(0);
  gStyle -> SetCanvasBorderMode(0);
  gStyle -> SetPadBorderMode(0);
  gStyle -> SetTitleBorderSize(0);
  gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadRightMargin(0.05);
  gStyle -> SetPadBottomMargin(0.12);
  gStyle -> SetTitleFillColor(0);
  gStyle -> SetPadTickX(1);
  gStyle -> SetPadTickY(1);
  gStyle -> SetNdivisions(510,"xyz");
*/
}
