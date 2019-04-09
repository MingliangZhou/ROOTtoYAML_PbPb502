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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"

char name[200];

void cutX(TGraphAsymmErrors*& gSts, double cuts_low, double cuts_up)
{
  int N = gSts->GetN();
  int i = 0;
  while(i<N)
  {
    double x; double y;
    gSts->GetPoint(i, x, y);
    if(x > cuts_up || x < cuts_low)
    {
      N --;
      gSts->RemovePoint(i);
    }
    else
    {
      i ++;
    }
  }
}

void extract(TFile* fOut, const char* fName, const char* cName, const char* hName, int nPad, bool modXbin, bool dummy=0)
{
	TFile* fIn = new TFile(fName);
	TCanvas* can = (TCanvas*)fIn->Get(cName);
	TGraphErrors* gSts;
	TGraphAsymmErrors* gSys;

	for(int iP=0; iP<nPad; iP++)
	{
		TPad* pad = (TPad*)can->GetListOfPrimitives()->At(iP);
		TObject* obj;
		TIter next(pad->GetListOfPrimitives());
		int iD = 0;
		while((obj=next()))
		{
			if(obj->InheritsFrom("TGraphErrors"))
			{
				fOut->cd();
				gSts = (TGraphErrors*)obj->Clone();
				if(dummy)
				{
					sprintf(name, "%s_pad%d_%s_Dummy%d",hName,iP+1,gSts->GetName(), iD);
					iD += 1;
				}
				else sprintf(name, "%s_pad%d_%s",hName,iP+1,gSts->GetName());
				gSts->SetName(name);
				if(modXbin)
				{
					for(int i=0; i<gSts->GetN(); i++)
					{
						double x;
						double y;
						gSts->GetPoint(i, x, y);
						gSts->SetPoint(i, 100*(1-x), y);
					}
				}
				gSts->Write();
			}
			else if(obj->InheritsFrom("TGraphAsymmErrors"))
			{
				fOut->cd();
				gSys = (TGraphAsymmErrors*)obj->Clone();
				if(dummy)
				{
					sprintf(name, "%s_pad%d_%s_Dummy%d",hName,iP+1,gSys->GetName(), iD);
					iD += 1;
				}
				else sprintf(name, "%s_pad%d_%s",hName,iP+1,gSys->GetName());
				gSys->SetName(name);
				if(modXbin)
				{
					for(int i=0; i<gSys->GetN(); i++)
					{
						double x;
						double y;
						gSys->GetPoint(i, x, y);
						gSys->SetPoint(i, 100*(1-x), y);
					}
				}
				if(nPad==2 && modXbin==1 && dummy==1) cutX(gSys, 0, 60);
				gSys->Write();
			}
			//cout<<"Reading: "<<obj->GetName()<<endl;
			obj->Draw();
		}
	}
}

void toPlot()
{
	TFile* fOut = new TFile("fromCanvas.root","RECREATE");
	extract(fOut, "5TeVPbPbsbys3.root", "comp_vn6vn4_har2_Cent", "fig07", 1, 1);
	extract(fOut, "5TeVPbPbsbys3.root", "comp_vn6vn4_har2_FCal", "fig18a", 1, 0);
	extract(fOut, "5TeVPbPbsbys3.root", "comp_vn6vn4_har2_Nchb", "fig18b", 1, 0);
	extract(fOut, "5TeVPbPbsbys_v1only.root", "papercompv1only_c4sub_Cent_har1", "fig08", 2, 1, 1);
	extract(fOut, "5TeVPbPbsbys_v1only.root", "papercompv1only_v4sub_Cent_har1", "fig09", 2, 1, 1);
	extract(fOut, "5TeVPbPbsbys3.root", "comp_vn2rat_FCal", "fig14a", 3, 0);
	extract(fOut, "5TeVPbPbsbys3.root", "comp_vn2rat_Nch", "fig14b", 3, 0);
	fOut->Close();
}







