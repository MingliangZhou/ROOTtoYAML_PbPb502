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

void grpow(TGraph*gr, double n)
{
	int N = gr->GetN();
	double *X = gr->GetX();
	double *Y = gr->GetY();
	double *Y1, *Y2;
	int type=0;
	if(gr->InheritsFrom("TGraphErrors")){
		Y1 = ((TGraphErrors*) gr)->GetEY();
	}else{
		type=1;
		Y1 = ((TGraphAsymmErrors*) gr)->GetEYlow();
		Y2 = ((TGraphAsymmErrors*) gr)->GetEYhigh();
	}
	for(int ip=0;ip<N;ip++){

		double val = Y[ip], val1=Y[ip]; int sign;
		sign=1; if(val<0) sign=-1;
		Y[ip] = pow(fabs(val),n)*sign;
		if(Y1[ip]<0) cout<<"HI"<<endl;
		val=val1-Y1[ip];      sign=1; if(val<0) sign=-1; Y1[ip] = Y[ip]-pow(fabs(val),n)*sign;

		if(type){
			val=val1+Y2[ip];    sign=1; if(val<0) sign=-1; Y2[ip] = pow(fabs(val),n)*sign-Y[ip];
		}
	}
}

TGraphAsymmErrors* toasygr(TGraphErrors*gr){
  int N = gr->GetN();
  double *X = gr->GetX();
  double *Y = gr->GetY();
  double *YE = gr->GetEY();
  TGraphAsymmErrors*gtmp = new TGraphAsymmErrors(N,X,Y,0,0,YE,YE);
  gtmp->SetName(gr->GetName());
  return gtmp;
}

void scalgr(TGraph*gr, double scal)
{

  int N = gr->GetN();
  double *X = gr->GetX();
  double *Y = gr->GetY();
  double *Ye1,*Ye2;
  double ty=1;
  if(gr->InheritsFrom("TGraphErrors")){
    Ye1 = ((TGraphErrors*) gr)->GetEY();
  }else{
    ty=0;
    Ye1 = ((TGraphAsymmErrors*) gr)->GetEYlow();
    Ye2 = ((TGraphAsymmErrors*) gr)->GetEYhigh();
  }
  double valtmp;
  for(int i=0;i<N;i++){
    if(scal>0){
      Y[i]*=scal;    Ye1[i]*=scal;    if(ty==0) Ye2[i]*=scal;
    }else{
      Y[i]*=scal;    Ye1[i]*=-scal;  if(ty==0) Ye2[i]*=-scal;
      if(ty==0){
  valtmp = Ye1[i]; Ye1[i]=Ye2[i]; Ye2[i]=valtmp;
      }
    }
  }
}

void cutX(TGraphErrors*& gSts, TGraphAsymmErrors*& gSys, double cuts_low, double cuts_up)
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
			gSys->RemovePoint(i);
		}
		else
		{
			i ++;
		}
	}
}

void cutX(TGraphAsymmErrors*& gSts, TGraphAsymmErrors*& gSys, double cuts_low, double cuts_up)
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
			gSys->RemovePoint(i);
		}
		else
		{
			i ++;
		}
	}
}

TGraph* grratio(TGraph*h1, TGraph*h2)
{
	char n[200];
	sprintf(n,"R%s",h1->GetName());
	TGraph *h3 = (TGraph*)h1->Clone(); h3->SetName(n);h3->SetTitle(n);

	double *x1 = h1->GetX();
	double *y1 = h1->GetY();  double *ty2 = h2->GetY();  double *y3 = h3->GetY();
	double *ey1u,*ey2u ,*ey3u,*ey1l,*ey2l ,*ey3l;
	double y2[1000];
	int type=0;
	if(h1->InheritsFrom("TGraphErrors")){
		ey1u = ((TGraphErrors*)h1)->GetEY();
		ey2u = ((TGraphErrors*)h2)->GetEY();
		ey3u =  ((TGraphErrors*)h3)->GetEY();
	}else{
		type=1;
		ey1u  = ((TGraphAsymmErrors*)h1)->GetEYhigh();
		ey2u = ((TGraphAsymmErrors*)h2)->GetEYhigh();//
		ey3u  = ((TGraphAsymmErrors*)h3)->GetEYhigh();

		ey1l  = ((TGraphAsymmErrors*)h1)->GetEYlow();
		ey2l = ((TGraphAsymmErrors*)h2)->GetEYlow(); //
		ey3l  = ((TGraphAsymmErrors*)h3)->GetEYlow();
	}

	int N= h1->GetN();
	for(int i=0;i<N;i++){
		y2[i] = h2->Eval(x1[i]);//
		y3[i] = ey3u[i]=0;
		if(y2[i]){
			y3[i] = y1[i]/y2[i];
			ey3u[i] = ey1u[i]/y2[i];
			if(type) ey3l[i] = ey1l[i]/y2[i];
		}
	}
	return h3;
}

void convertCent(TGraphErrors*& gSts, TGraphAsymmErrors*& gSys, TProfile* gCvt)
{
	for(int i=0; i<gSts->GetN(); i++)
	{
		double x; double y;
		gSts->GetPoint(i, x, y);
		x = gCvt->Interpolate(x);
		gSts->SetPoint(i, x, y);
		gSys->SetPoint(i, x, y);
	}
}

void toPaper()
{
	double cuts_c_low[3][5] = {{0, 0, 0, 0, 0}, {0, 0, 0*4.1, 0.05*4.1, 0.05*4.1}, {0, 0, 0*2800, 0.05*2800, 0.05*2800}};
	double cuts_c_up[3][5] = {{99, 45, 79, 60, 60}, {99, 99, 1.175*4.1, 1.175*4.1, 1.175*4.1}, {99, 99, 1.27*2800, 1.27*2800, 1.27*2800}};
	double cuts_sc_low[3][5] = {{0, 0, 0, 0, 0}, {0, 0, 0*4.1, 0.05*4.1, 0}, {0, 0, 0*2800, 0.05*2800, 0}};
	double cuts_sc_up[3][5] = {{99, 99, 70, 70, 99}, {99, 99, 1.18*4.1, 1.18*4.1, 99}, {99, 99, 1.27*2800, 1.27*2800, 99}};

	TFile* fCvt = new TFile("hist_cvt.root");
	TProfile* hCvt = (TProfile*)fCvt->Get("cvt_NchRec_FCal");
	hCvt->Rebin(10);

	char name[300];
	TFile* fIn;
	for(int iF=0; iF<3; iF++)
	{
		for(int iB=0; iB<5; iB++)
		{
			if(iF==0) sprintf(name,"hist_PbPb502_binCent_bin%d.root",iB);
			else if(iF==1) sprintf(name,"hist_PbPb502_binFCal_bin%d.root",iB);
			else sprintf(name,"hist_PbPb502_binNch_bin%d.root",iB);
			fIn = new TFile(name);
			if(iF==0) sprintf(name,"hist_PbPb502_binCent_bin%d_paper.root",iB);
			else if(iF==1) sprintf(name,"hist_PbPb502_binFCal_bin%d_paper.root",iB);
			else sprintf(name,"hist_PbPb502_binNch_bin%d_paper.root",iB);
			TFile* fOut = new TFile(name,"RECREATE");
			fOut->cd();
			for(int iV=1; iV<=4; iV++)
			{
				for(int iP=0; iP<=6; iP++)
				{
					TGraphErrors* gSts;
					TGraphAsymmErrors* gSys;

					// Figure 3
					TGraphAsymmErrors* gSts_v2;
					TGraphAsymmErrors* gSys_v2;
					sprintf(name, "sts_c2_3sub_Har%d_Pt%d", iV, iP);
					gSts_v2 = toasygr((TGraphErrors*)fIn->Get(name)->Clone());
					sprintf(name, "sts_v2_3sub_Har%d_Pt%d", iV, iP);
					gSts_v2->SetName(name);
					sprintf(name, "sys_c2_3sub_Har%d_Pt%d", iV, iP);
					gSys_v2 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_v2_3sub_Har%d_Pt%d", iV, iP);
					gSys_v2->SetName(name);

					grpow(gSts_v2, 1./2); grpow(gSys_v2, 1./2);
					gSts_v2->Write(); gSys_v2->Write();



					// Figure 4
					TGraphErrors* gSts_nc4;
					TGraphAsymmErrors* gSys_nc4;
					sprintf(name, "sts_nc4_1sub_Har%d_Pt%d", iV, iP);
					gSts_nc4 = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nc4_1sub_Har%d_Pt%d", iV, iP);
					gSys_nc4 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts_nc4, gSys_nc4, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_nc4->Write(); gSys_nc4->Write();

					sprintf(name, "sts_nc4_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nc4_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts, gSys, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts->Write(); gSys->Write();



					// Figure 5
					TGraphAsymmErrors* gSts_v4;
					TGraphAsymmErrors* gSys_v4;
					sprintf(name, "sts_c4_1sub_Har%d_Pt%d", iV, iP);
					gSts_v4 = toasygr((TGraphErrors*)fIn->Get(name)->Clone());
					sprintf(name, "sts_v4_1sub_Har%d_Pt%d", iV, iP);
					gSts_v4->SetName(name);
					sprintf(name, "sys_c4_1sub_Har%d_Pt%d", iV, iP);
					gSys_v4 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_v4_1sub_Har%d_Pt%d", iV, iP);
					gSys_v4->SetName(name);

					grpow(gSts_v4, 1./4); grpow(gSys_v4, 1./4);
					scalgr(gSts_v4, -1); scalgr(gSys_v4, -1);
					cutX(gSts_v4, gSys_v4, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_v4->Write(); gSys_v4->Write();

					gSts = (TGraphErrors*)grratio(gSts_v4, gSts_v2);
					sprintf(name, "sts_v42_1sub_Har%d_Pt%d", iV, iP);
					gSts->SetName(name);
					gSys = (TGraphAsymmErrors*)grratio(gSys_v4, gSys_v2);
					sprintf(name, "sys_v42_1sub_Har%d_Pt%d", iV, iP);
					gSys->SetName(name);

					cutX(gSts, gSys, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts->Write(); gSys->Write();



					// Figure 6
					TGraphErrors* gSts_nc6;
					TGraphAsymmErrors* gSys_nc6;
					sprintf(name, "sts_nc6_1sub_Har%d_Pt%d", iV, iP);
					gSts_nc6 = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nc6_1sub_Har%d_Pt%d", iV, iP);
					gSys_nc6 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts_nc6, gSys_nc6, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_nc6->Write(); gSys_nc6->Write();



					// Figure 7
					TGraphAsymmErrors* gSts_v6;
					TGraphAsymmErrors* gSys_v6;
					sprintf(name, "sts_c6_1sub_Har%d_Pt%d", iV, iP);
					gSts_v6 = toasygr((TGraphErrors*)fIn->Get(name)->Clone());
					sprintf(name, "sts_v6_1sub_Har%d_Pt%d", iV, iP);
					gSts_v6->SetName(name);
					sprintf(name, "sys_c6_1sub_Har%d_Pt%d", iV, iP);
					gSys_v6 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_v6_1sub_Har%d_Pt%d", iV, iP);
					gSys_v6->SetName(name);

					grpow(gSts_v6, 1./6); grpow(gSys_v6, 1./6);
					scalgr(gSts_v6, pow(4,-1./6)); scalgr(gSys_v6, pow(4,-1./6));
					cutX(gSts_v6, gSys_v6, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_v6->Write(); gSys_v6->Write();

					TGraphErrors* gSts_v64;
					TGraphAsymmErrors* gSys_v64;
					gSts_v64 = (TGraphErrors*)grratio(gSts_v6, gSts_v4);
					sprintf(name, "sts_v64_1sub_Har%d_Pt%d", iV, iP);
					gSts_v64->SetName(name);
					gSys_v64 = (TGraphAsymmErrors*)grratio(gSys_v6, gSys_v4);
					sprintf(name, "sys_v64_1sub_Har%d_Pt%d", iV, iP);
					gSys_v64->SetName(name);

					//cutX(gSts_v64, gSys_v64, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					cutX(gSts_v64, gSys_v64, 5, 70);
					gSts_v64->Write(); gSys_v64->Write();



					// Figure 8
					sprintf(name, "sts_c4_1sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_c4_1sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts, gSys, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts->Write(); gSys->Write();

					sprintf(name, "sts_c4_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_c4_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts, gSys, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts->Write(); gSys->Write();



					// Figure 9
					cutX(gSts_v4, gSys_v4, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_v4->Write(); gSys_v4->Write();

					sprintf(name, "sts_c4_3sub_Har%d_Pt%d", iV, iP);
					gSts_v4 = toasygr((TGraphErrors*)fIn->Get(name)->Clone());
					sprintf(name, "sts_v4_3sub_Har%d_Pt%d", iV, iP);
					gSts_v4->SetName(name);
					sprintf(name, "sys_c4_3sub_Har%d_Pt%d", iV, iP);
					gSys_v4 = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_v4_3sub_Har%d_Pt%d", iV, iP);
					gSys_v4->SetName(name);

					grpow(gSts_v4, 1./4); grpow(gSys_v4, 1./4);
					scalgr(gSts_v4, -1); scalgr(gSys_v4, -1);
					cutX(gSts_v4, gSys_v4, cuts_c_low[iF][iV], cuts_c_up[iF][iV]);
					gSts_v4->Write(); gSys_v4->Write();



					// Figure 10
					TGraphErrors* gSts_nsc;
					TGraphAsymmErrors* gSys_nsc;
					sprintf(name, "sts_nsc_1sub_Har%d_Pt%d", iV, iP);
					gSts_nsc = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nsc_1sub_Har%d_Pt%d", iV, iP);
					gSys_nsc = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts_nsc, gSys_nsc, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					gSts_nsc->Write(); gSys_nsc->Write();

					sprintf(name, "sts_nsc_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nsc_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					gSts->Write(); gSys->Write();



					// Figure 11



					// Figure 12
					TGraphErrors* gSts_nac;
					TGraphAsymmErrors* gSys_nac;
					sprintf(name, "sts_nac_1sub_Har%d_Pt%d", iV, iP);
					gSts_nac = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nac_1sub_Har%d_Pt%d", iV, iP);
					gSys_nac = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts_nac, gSys_nac, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					gSts_nac->Write(); gSys_nac->Write();

					sprintf(name, "sts_nac_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_nac_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					gSts->Write(); gSys->Write();



					// Figure 13



					// Figure 14


					
					// Figure 15



					// Figure 16
					if(iF==2)
					{
						convertCent(gSts_nc4, gSys_nc4, hCvt);
						sprintf(name, "cvt_sts_nc4_1sub_Har%d_Pt%d", iV, iP);
						gSts_nc4->SetName(name);
						sprintf(name, "cvt_sys_nc4_1sub_Har%d_Pt%d", iV, iP);
						gSys_nc4->SetName(name);
						gSts_nc4->Write(); gSys_nc4->Write();
					}

					

					// Figure 17
					if(iF==2)
					{
						convertCent(gSts_nc6, gSys_nc6, hCvt);
						sprintf(name, "cvt_sts_nc6_1sub_Har%d_Pt%d", iV, iP);
						gSts_nc6->SetName(name);
						sprintf(name, "cvt_sys_nc6_1sub_Har%d_Pt%d", iV, iP);
						gSys_nc6->SetName(name);
						gSts_nc6->Write(); gSys_nc6->Write();
					}



					// Figure 18
					if(iF==2)
					{
						convertCent(gSts_v64, gSys_v64, hCvt);
						sprintf(name, "cvt_sts_v64_1sub_Har%d_Pt%d", iV, iP);
						gSts_v64->SetName(name);
						sprintf(name, "cvt_sys_v64_1sub_Har%d_Pt%d", iV, iP);
						gSys_v64->SetName(name);
						gSts_v64->Write(); gSys_v64->Write();
					}



					// Figure 19



					// Figure 20
					if(iF==2)
					{
						convertCent(gSts_nsc, gSys_nsc, hCvt);
						sprintf(name, "cvt_sts_nsc_1sub_Har%d_Pt%d", iV, iP);
						gSts_nsc->SetName(name);
						sprintf(name, "cvt_sys_nsc_1sub_Har%d_Pt%d", iV, iP);
						gSys_nsc->SetName(name);
						gSts_nsc->Write(); gSys_nsc->Write();

						convertCent(gSts_nac, gSys_nac, hCvt);
						sprintf(name, "cvt_sts_nac_1sub_Har%d_Pt%d", iV, iP);
						gSts_nac->SetName(name);
						sprintf(name, "cvt_sys_nac_1sub_Har%d_Pt%d", iV, iP);
						gSys_nac->SetName(name);
						gSts_nac->Write(); gSys_nac->Write();
					}



					// Figure 28, 29
					sprintf(name, "sts_sc_1sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_sc_1sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					//cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					cutX(gSts, gSys, cuts_sc_low[iF][iV], 80);
					gSts->Write(); gSys->Write();

					sprintf(name, "sts_sc_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_sc_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					//cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					cutX(gSts, gSys, cuts_sc_low[iF][iV], 80);
					gSts->Write(); gSys->Write();


					// Figure 30
					sprintf(name, "sts_ac_1sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_ac_1sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					//cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					cutX(gSts, gSys, cuts_sc_low[iF][iV], 80);
					gSts->Write(); gSys->Write();

					sprintf(name, "sts_ac_3sub_Har%d_Pt%d", iV, iP);
					gSts = (TGraphErrors*)fIn->Get(name)->Clone();
					sprintf(name, "sys_ac_3sub_Har%d_Pt%d", iV, iP);
					gSys = (TGraphAsymmErrors*)fIn->Get(name)->Clone();
					//cutX(gSts, gSys, cuts_sc_low[iF][iV], cuts_sc_up[iF][iV]);
					cutX(gSts, gSys, cuts_sc_low[iF][iV], 80);
					gSts->Write(); gSys->Write();
				}
			}
			fOut->Close();
		}
	}
}
