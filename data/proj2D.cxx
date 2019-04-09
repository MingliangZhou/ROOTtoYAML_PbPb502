#include <iostream>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TLegend.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TView3D.h>
#include <TTree.h>
#include <TLatex.h>
#include <TROOT.h>
#include <stdio.h>
#include "TGraphErrors.h"

void proj2D()
{
  char name[300];
  TFile* fIn = new TFile("etnchPbPb.root");
  TFile* fOut = new TFile("etnchPbPb_proj.root", "RECREATE");
	TH1D* hTmp;
	TH2D* h2D;
  fOut->cd();
	h2D = (TH2D*)fIn->Get("hetnch_0_0");
	h2D->Rebin2D(6,4);
	h2D->Write();
	hTmp = (TH1D*)(((TH2D*)fIn->Get("hetnch_0_0")->Clone())->ProjectionX());
	hTmp->SetName("hEt");
	hTmp->Write();
	hTmp = (TH1D*)(((TH2D*)fIn->Get("hetnch_0_0")->Clone())->ProjectionY());
	hTmp->SetName("hNchRec");
	hTmp->Write();
  fOut->Close();
}
