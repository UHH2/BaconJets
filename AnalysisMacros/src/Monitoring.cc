#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLine.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include "../include/tdrstyle_mod15.h"
#include <TFile.h>
#include <iostream>
#include <assert.h>
#include "../include/CorrectionObject.h"
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include "../include/parameters.h"
#include <TProfile.h>
#include "../include/useful_functions.h"
#include <cmath>
#include <math.h>
#include <TMultiGraph.h>

using namespace std;

const int n_input = 3;
const int Fit_range[n_input+1] = {0, 12600, 19300, 35000};

const TString Name_range[n_input] = {"BCD", "EFearly", "FlateGH"};

void Save_2D_Plot(TH2D* hist, TString Method, TString Runnr, TString input_path){

  TLatex *tex_lumi = new TLatex();
  tex_lumi->SetNDC();
  tex_lumi->SetTextSize(0.045); 

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045); 

  TCanvas* p1 = new TCanvas("p1","p1",800,600);
  gStyle->SetOptStat(0000);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.1);
  hist->GetYaxis()->SetTitleOffset(1.5);

  if(Method.Contains("Param")){
    hist->GetZaxis()->SetRangeUser(0.8, 1.3);
  }
  else hist->SetMinimum(0);

  hist->Draw("colz");

  tex->DrawLatex(0.2,0.92,Method);
  tex_lumi->DrawLatex(0.53,0.92,"Run"+Runnr+" (13TeV)");
  p1->SaveAs(input_path+"/Monitoring/"+Method+"_"+Runnr+".pdf");
  delete p1;
}




void CorrectionObject::Monitoring(){
  cout << "--------------- Starting Monitoring() ---------------" << endl << endl;

  TFile* f_monitoring[n_input];
  int n_pt_ = max(n_pt,n_pt_HF);

  bool eta_cut_bool;

  TH2D *hist_A[n_input][n_eta_control-1][n_pt_-2];
  TProfile *pr_A[n_input][n_eta_control-1][n_pt_-2];


  // TGraphErrors *Npv_res[n_input][n_eta_control-1][n_pt-1];
  //  TH2D *hist_Npv[n_input][n_eta_control-1][n_pt-1];
  // TProfile *pr_Npv[n_input][n_eta_control-1][n_pt-1];

  TGraphErrors *rel_res[n_input][n_eta_control-1][n_pt_-2];
  TGraphErrors *mpf_res[n_input][n_eta_control-1][n_pt_-2];

  TGraphErrors *rel_ratio[n_input][n_eta_control-1][n_pt_-2];
  TGraphErrors *mpf_ratio[n_input][n_eta_control-1][n_pt_-2];

  TF1 *Fit_rel[n_input][n_eta_control-1][n_pt_-2];
  TF1 *Fit_mpf[n_input][n_eta_control-1][n_pt_-2];
 
  TH2D *hist_B[n_input][n_eta_control-1][n_pt_-2];
  TProfile *pr_B[n_input][n_eta_control-1][n_pt_-2];

  TH2D *chi2_rel_fit[n_input];
  TH2D *param_rel_fit[n_input];

  TH2D *chi2_mpf_fit[n_input];
  TH2D *param_mpf_fit[n_input];



  for(int i=0; i<n_input; i++){
    //    cout<<"Read in File: "<<endl;
    //  cout<<CorrectionObject::_input_path+"uhh2.AnalysisModuleRunner.DATA.DATA_Run"+Name_range[i]+"_AK4CHS.root"<<endl;

    f_monitoring[i] = new TFile(CorrectionObject::_input_path+"/uhh2.AnalysisModuleRunner.DATA.DATA_Run"+Name_range[i]+"_AK4CHS.root");
    //   cout<<"Create Hist"<<endl;
    for(int j=0; j<n_eta_control-1; j++){
     if(fabs(eta_bins_control[j])<fabs(eta_bins_control[j+1])){
       eta_cut_bool = fabs(eta_bins_control[j])>eta_cut;
      }
      else {
	eta_cut_bool = fabs(eta_bins_control[j+1])>eta_cut;
      }

      for(int k=1; k < (eta_cut_bool ?  n_pt_HF-2 : n_pt-2 ) ; k++ ){
	//	cout<<"j= "<<j<<" k="<<k<<" eta_range_control2[j]="<<eta_range_control2[j]<<" "<<eta_range_control2[j+1]<<" pt="<<(eta_cut_bool?pt_range_HF:pt_range)[k]<<" "<<(eta_cut_bool?pt_range_HF:pt_range)[k+1]<<endl;
        hist_A[i][j][k]    = (TH2D*)f_monitoring[i]    ->Get("Monitoring_Final/hist_data_A_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
	hist_A[i][j][k]->Print();
	pr_A[i][j][k]      = (TProfile*)hist_A[i][j][k] ->ProfileX(Form("prof_A_%i_%d_%d",i,j,k));
        
	hist_B[i][j][k]    = (TH2D*)f_monitoring[i]    ->Get("Monitoring_Final/hist_data_B_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
	hist_B[i][j][k]->Print();

	pr_B[i][j][k]      = (TProfile*)hist_B[i][j][k] ->ProfileX(Form("prof_B_%i_%d_%d",i,j,k));
	pr_B[i][j][k]->Print();
	/*
        hist_Npv[i][j][k]    = (TH2D*)f_monitoring[i]    ->Get("Monitoring_Final/hist_data_Npv_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
	hist_Npv[i][j][k]->Print();
	pr_Npv[i][j][k]      = (TProfile*)hist_A[i][j][k] ->ProfileX(Form("prof_Npv_%i_%d_%d",i,j,k));
	pr_Npv[i][j][k]->Print();
	*/
      }
    }
    
    chi2_mpf_fit[i]  = new TH2D("chi2_mpf_fit",""+Name_range[i],n_eta_control-1, eta_bins_control, n_pt_-2, pt_bins);
    param_mpf_fit[i] = new TH2D("param_mpf_fit",""+Name_range[i],n_eta_control-1, eta_bins_control, n_pt_-2, pt_bins);

    chi2_rel_fit[i]  = new TH2D("chi2_rel_fit",""+Name_range[i],n_eta_control-1, eta_bins_control, n_pt_-2, pt_bins);
    param_rel_fit[i] = new TH2D("param_rel_fit",""+Name_range[i],n_eta_control-1, eta_bins_control, n_pt_-2, pt_bins);

    //   cout<<"Finish Load Hists"<<endl;

    int n_lumi = 0;
    n_lumi = pr_B[i][1][1]->GetNbinsX();
    
    double bin_width = 0; 
    bin_width = pr_B[i][1][1]->GetXaxis()->GetBinWidth(1);
    
    double res_rel[n_eta_control-1][n_pt_-2][n_lumi-1];
    double res_mpf[n_eta_control-1][n_pt_-2][n_lumi-1];
    

    double err_res_rel[n_eta_control-1][n_pt_-2][n_lumi-1];
    double err_res_mpf[n_eta_control-1][n_pt_-2][n_lumi-1];

    /*
    double Npv_mean[n_eta_control-1][n_pt_-2][n_lumi-1];
    double err_Npv_mean[n_eta_control-1][n_pt_-2][n_lumi-1];
    */
    
    double xbin_tgraph[n_lumi],zero[n_lumi];
    double xbin = 0;
    for(int i=0;i<n_lumi;i++){
      xbin = xbin + bin_width;
      xbin_tgraph[i] = xbin;
      zero[i] = 0;
    }

    //  cout<<"After Lumi!"<<endl;

    for(int j=0; j<n_eta_control-1; j++){
     if(fabs(eta_bins_control[j])<fabs(eta_bins_control[j+1])){
       eta_cut_bool = fabs(eta_bins_control[j])>eta_cut;
     }
      else {
	eta_cut_bool = fabs(eta_bins_control[j+1])>eta_cut;
      }

      for(int k=1; k < (eta_cut_bool ?  n_pt_HF-2 : n_pt-2 ) ; k++ ){
	for(int l =0; l<n_lumi; l++){
	  
	  //	  cout<<"Fit_Rel_"+Name_range[i]+"_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]<<endl;

	  res_rel[j][k][l]=0;
	  err_res_rel[j][k][l]=0;
	  res_mpf[j][k][l]=0;
	  err_res_mpf[j][k][l]=0;

	  if(hist_A[i][j][k]->GetEntries()>100){
	    res_rel[j][k][l] = (1+pr_A[i][j][k]->GetBinContent(l+1))/(1-pr_A[i][j][k]->GetBinContent(l+1));
	    err_res_rel[j][k][l] = 2/(pow((1-pr_A[i][j][k]->GetBinContent(l+1)),2)) * pr_A[i][j][k]->GetBinError(l+1);
	  }	
	  if(hist_B[i][j][k]->GetEntries()>100){
	    res_mpf[j][k][l] = (1+pr_B[i][j][k]->GetBinContent(l+1))/(1-pr_B[i][j][k]->GetBinContent(l+1));
	    err_res_mpf[j][k][l] = 2/(pow((1-pr_B[i][j][k]->GetBinContent(l+1)),2)) * pr_B[i][j][k]->GetBinError(l+1);
	  }
	  /*
	  if(hist_Npv[i][j][k]->GetEntries()>100){
	    Npv_mean[j][k][l] = pr_Npv[i][j][k]->GetBinContent(l+1);
	    err_res_mpf[j][k][l] = pr_Npv[i][j][k]->GetBinError(l+1);
	  }
	  */
	}

	rel_res[i][j][k]= new TGraphErrors(n_lumi, xbin_tgraph, res_rel[j][k], zero, err_res_rel[j][k]);
	rel_res[i][j][k]= (TGraphErrors*)CleanEmptyPoints(rel_res[i][j][k]);

	mpf_res[i][j][k]= new TGraphErrors(n_lumi, xbin_tgraph, res_mpf[j][k], zero, err_res_mpf[j][k]);
	mpf_res[i][j][k]= (TGraphErrors*)CleanEmptyPoints(mpf_res[i][j][k]);

	/*
	Npv_res[i][j][k]= new TGraphErrors(n_lumi, xbin_tgraph, Npv_mean[j][k], zero, err_Npv_mean[j][k]);
	Npv_res[i][j][k]= (TGraphErrors*)CleanEmptyPoints(Npv_res[i][j][k]);
	*/

	Fit_rel[i][j][k] = new TF1("Fit_Rel_"+Name_range[i]+"_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1],"pol0", Fit_range[i] ,Fit_range[i+1]);
	Fit_mpf[i][j][k] = new TF1("Fit_MPF_"+Name_range[i]+"_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1],"pol0", Fit_range[i] ,Fit_range[i+1]);

	
	rel_res[i][j][k] ->Fit("Fit_Rel_"+Name_range[i]+"_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1],"","",Fit_range[i], Fit_range[i+1]);
	mpf_res[i][j][k] ->Fit("Fit_MPF_"+Name_range[i]+"_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1],"","",Fit_range[i], Fit_range[i+1]);

	rel_ratio[i][j][k] = (TGraphErrors*)BuildRatio(rel_res[i][j][k], Fit_rel[i][j][k]->GetParameter(0),Fit_rel[i][j][k]->GetParError(0));
	mpf_ratio[i][j][k] = (TGraphErrors*)BuildRatio(mpf_res[i][j][k], Fit_mpf[i][j][k]->GetParameter(0),Fit_mpf[i][j][k]->GetParError(0));
      
	double chi2_mpf = 0;
	if( Fit_mpf[i][j][k]->GetChisquare() / Fit_mpf[i][j][k]->GetNDF() ==  Fit_mpf[i][j][k]->GetChisquare() / Fit_mpf[i][j][k]->GetNDF()) chi2_mpf = Fit_mpf[i][j][k]->GetChisquare() / Fit_mpf[i][j][k]->GetNDF();

	double chi2_rel = 0;
	if( Fit_rel[i][j][k]->GetChisquare() / Fit_rel[i][j][k]->GetNDF() ==  Fit_rel[i][j][k]->GetChisquare() / Fit_rel[i][j][k]->GetNDF()) chi2_rel = Fit_rel[i][j][k]->GetChisquare() / Fit_rel[i][j][k]->GetNDF();
   
	chi2_mpf_fit[i]->SetBinContent(j+1, k+1, chi2_mpf);
	chi2_rel_fit[i]->SetBinContent(j+1, k+1, chi2_rel);

	param_mpf_fit[i]->SetBinContent(j+1, k+1, Fit_mpf[i][j][k]->GetParameter(0));
	param_rel_fit[i]->SetBinContent(j+1, k+1, Fit_rel[i][j][k]->GetParameter(0));	
      }
    }
    Save_2D_Plot(chi2_mpf_fit[i], "MPF_Chi2", Name_range[i], CorrectionObject::_input_path);
    Save_2D_Plot(chi2_rel_fit[i], "Rel_Chi2", Name_range[i], CorrectionObject::_input_path);

    Save_2D_Plot(param_mpf_fit[i], "MPF_Param", Name_range[i], CorrectionObject::_input_path);
    Save_2D_Plot(param_rel_fit[i], "Rel_Param", Name_range[i], CorrectionObject::_input_path);
  }

  //  cout<<"Setup canvas and draw several TGraphs"<<endl;


  TLine* ave = new TLine(0,1,35000,1);
  
  TLine* lineB = new TLine(5750, 0.82, 5750, 1.42);
  lineB -> SetLineStyle(2);
  TLine* lineC = new TLine(8330, 0.82, 8300, 1.42);
  lineC -> SetLineStyle(2);
  TLine* lineD = new TLine(12600, 0.82, 12600,1.42);
  lineD -> SetLineStyle(2);
  TLine* lineEFearly = new TLine(19300, 0.82, 19300,1.42);
  lineEFearly -> SetLineStyle(2);
  TLine* lineFlateG = new TLine(27300, 0.82, 27300,1.42);
  lineFlateG -> SetLineStyle(2);
  
  TLine* lineB_ratio = new TLine(5750, 0.86, 5750, 1.14);
  lineB_ratio -> SetLineStyle(2);
  TLine* lineC_ratio = new TLine(8330, 0.86, 8300, 1.14);
  lineC_ratio -> SetLineStyle(2);
  TLine* lineD_ratio = new TLine(12600, 0.86, 12600,1.14);
  lineD_ratio -> SetLineStyle(2);
  TLine* lineEFearly_ratio = new TLine(19300, 0.86, 19300,1.14);
  lineEFearly_ratio -> SetLineStyle(2);
  TLine* lineFlateG_ratio = new TLine(27300, 0.86, 27300,1.14);
  lineFlateG_ratio -> SetLineStyle(2);

  TLatex *tex1 = new TLatex();
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.036); 

  for(int j=0; j<n_eta_control-1; j++){
     if(fabs(eta_bins_control[j])<fabs(eta_bins_control[j+1])){
       eta_cut_bool = fabs(eta_bins_control[j])>eta_cut;
     }
      else {
	eta_cut_bool = fabs(eta_bins_control[j+1])>eta_cut;
      }

    for(int k=1; k<( eta_cut_bool ?  n_pt_HF-2 : n_pt-2 ); k++){
 
      //dummy for tdrCanvas 
      TH1D *h = new TH1D("h",";dummy;",5000,0,40000);
      h->SetMaximum(-0.4);
      h->SetMinimum(0.4);
      
      TH1D *d = new TH1D("d",";dummy;",5000,0,40000);
      d->SetMaximum(0.5);
      d->SetMinimum(1.5);
      
      TH1D *d2 = new TH1D("d2",";dummy;",5000,0,40000);
      d2->SetMaximum(0.5);
      d2->SetMinimum(1.5);
      



      TCanvas* c3 = tdrDiCanvas("c3",h,d,4,10,CorrectionObject::_lumitag);
      c3->cd(1);
      gStyle->SetOptFit(00000);
      
      h->GetXaxis()->SetTitle("Luminosity");
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.80);
      h->GetXaxis()->SetLimits(0,35000);
      h->GetXaxis()->SetLabelSize(0);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetRangeUser(0.82,1.42);
      h->GetYaxis()->SetTitle("R_{Rel}");

      TLegend *leg_rel_res;
      leg_rel_res = new TLegend(0.45,0.7,0.51,0.9,"","brNDC");//x+0.1
      leg_rel_res->SetBorderSize(0);
      leg_rel_res->SetTextSize(0.036);
      leg_rel_res->SetFillColor(10);
      leg_rel_res->SetFillStyle(1001);
      leg_rel_res->SetLineColor(1);
      leg_rel_res->SetTextFont(42);
      leg_rel_res->SetHeader("Rel Response, "+eta_range_control[j]+"#leq#eta<"+eta_range_control[j+1]); //+", #alpha<"+s_blpha_cut
 
      for(int i = 0; i<n_input; i++){
	rel_res[i][j][k]->SetMarkerStyle(2);
	if(i==0)     {  rel_res[i][j][k]->SetLineColor(kBlue);   rel_res[i][j][k]->SetMarkerColor(kBlue);}
	else if(i==1){  rel_res[i][j][k]->SetLineColor(kViolet); rel_res[i][j][k]->SetMarkerColor(kViolet);}
	else if(i==2){  rel_res[i][j][k]->SetLineColor(kGreen);  rel_res[i][j][k]->SetMarkerColor(kGreen);}
	else if(i==3){  rel_res[i][j][k]->SetLineColor(kRed);    rel_res[i][j][k]->SetMarkerColor(kRed);}
	leg_rel_res->AddEntry(rel_res[i][j][k], Name_range[i],"lep"); 
	
	rel_res[i][j][k]->Draw("PE SAME");
	c3->Modified();
      }

      lineB->Draw("SAME");
      lineC->Draw("SAME");
      lineD->Draw("SAME");
      lineEFearly->Draw("SAME");
      lineFlateG->Draw("SAME");
      
      leg_rel_res->Draw();
      
      tex1->DrawLatex(0.7,0.75, (eta_cut_bool?pt_range_HF:pt_range)[k]+"<p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
      c3->Modified();
      c3->Update();
      c3->cd(2);

      d->GetXaxis()->SetTitle("Luminosity");
      d->GetXaxis()->SetTitleSize(0.1);
      d->GetXaxis()->SetLabelSize(0.08);
      d->GetXaxis()->SetTitleOffset(1.2);
      d->GetXaxis()->SetLimits(0,35000);
      d->GetYaxis()->SetRangeUser(0.86,1.14);
      d->GetXaxis()->SetTitleOffset(1.2);
      d->GetYaxis()->SetLabelSize(0.08);
      d->GetYaxis()->SetTitle("Data/Fit");

      lineB_ratio->Draw("SAME");
      lineC_ratio->Draw("SAME");
      lineD_ratio->Draw("SAME");
      lineEFearly_ratio->Draw("SAME");
      lineFlateG_ratio->Draw("SAME");
      ave->Draw("SAME");
      
      for(int i = 0; i<n_input; i++){
	rel_ratio[i][j][k]->SetMarkerStyle(2);
	if(i==0)     {  rel_ratio[i][j][k]->SetLineColor(kBlue);   rel_ratio[i][j][k]->SetMarkerColor(kBlue);}
	else if(i==1){  rel_ratio[i][j][k]->SetLineColor(kViolet); rel_ratio[i][j][k]->SetMarkerColor(kViolet);}
	else if(i==2){  rel_ratio[i][j][k]->SetLineColor(kGreen);  rel_ratio[i][j][k]->SetMarkerColor(kGreen);}
	else if(i==3){  rel_ratio[i][j][k]->SetLineColor(kRed);    rel_ratio[i][j][k]->SetMarkerColor(kRed);}
	
	rel_ratio[i][j][k]->Draw("PE SAME");
	c3->Modified();
      }
      c3->Update();
      
      c3->SaveAs(CorrectionObject::_input_path+"/Monitoring/Rel_res_" + CorrectionObject::_generator_tag + "_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      delete c3;


      TCanvas* c4 = tdrDiCanvas("c4",h,d2,4,10,CorrectionObject::_lumitag);
      c4->cd(1);
      gStyle->SetOptFit(00000);
      
      h->GetXaxis()->SetTitle("Luminosity");
      h->GetXaxis()->SetTitleSize(0.1);
      h->GetXaxis()->SetTitleOffset(0.80);
      h->GetXaxis()->SetLimits(0,35000);
      h->GetXaxis()->SetLabelSize(0);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetRangeUser(0.82,1.42);
      h->GetYaxis()->SetTitle("R_{MPF}");

      TLegend *leg_mpf_res;
      leg_mpf_res = new TLegend(0.45,0.7,0.51,0.9,"","brNDC");//x+0.1
      leg_mpf_res->SetBorderSize(0);
      leg_mpf_res->SetTextSize(0.036);
      leg_mpf_res->SetFillColor(10);
      leg_mpf_res->SetFillStyle(1001);
      leg_mpf_res->SetLineColor(1);
      leg_mpf_res->SetTextFont(42);
      leg_mpf_res->SetHeader("MPF Response, "+eta_range_control[j]+"#leq#eta<"+eta_range_control[j+1]); //+", #alpha<"+s_blpha_cut
 
      for(int i = 0; i<n_input; i++){
	mpf_res[i][j][k]->SetMarkerStyle(2);
	if(i==0)     {  mpf_res[i][j][k]->SetLineColor(kBlue);   mpf_res[i][j][k]->SetMarkerColor(kBlue);}
	else if(i==1){  mpf_res[i][j][k]->SetLineColor(kViolet); mpf_res[i][j][k]->SetMarkerColor(kViolet);}
	else if(i==2){  mpf_res[i][j][k]->SetLineColor(kGreen);  mpf_res[i][j][k]->SetMarkerColor(kGreen);}
	else if(i==3){  mpf_res[i][j][k]->SetLineColor(kRed);    mpf_res[i][j][k]->SetMarkerColor(kRed);}
	leg_mpf_res->AddEntry(mpf_res[i][j][k], Name_range[i],"lep"); 
	
	mpf_res[i][j][k]->Draw("PE SAME");
	c4->Modified();
      }

      lineB->Draw("SAME");
      lineC->Draw("SAME");
      lineD->Draw("SAME");
      lineEFearly->Draw("SAME");
      lineFlateG->Draw("SAME");
      
      leg_mpf_res->Draw();
      
      tex1->DrawLatex(0.7,0.75, (eta_cut_bool?pt_range_HF:pt_range)[k]+"<p_{T}<"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]);
      c4->Modified();
      c4->Update();
      c4->cd(2);

      d2->GetXaxis()->SetTitle("Luminosity");
      d2->GetXaxis()->SetTitleSize(0.1);
      d2->GetXaxis()->SetTitleOffset(1.2);
      d2->GetXaxis()->SetLimits(0,35000);
      d2->GetYaxis()->SetRangeUser(0.86,1.14);
      d2->GetXaxis()->SetTitleOffset(1.2);
      d2->GetXaxis()->SetLabelSize(0.1);
      d2->GetYaxis()->SetLabelSize(0.1);
      d2->GetYaxis()->SetTitle("Data/Fit");

      lineB_ratio->Draw("SAME");
      lineC_ratio->Draw("SAME");
      lineD_ratio->Draw("SAME");
      lineEFearly_ratio->Draw("SAME");
      lineFlateG_ratio->Draw("SAME");
      ave->Draw("SAME");
      
      for(int i = 0; i<n_input; i++){
	mpf_ratio[i][j][k]->SetMarkerStyle(2);
	if(i==0)     {  mpf_ratio[i][j][k]->SetLineColor(kBlue);   mpf_ratio[i][j][k]->SetMarkerColor(kBlue);}
	else if(i==1){  mpf_ratio[i][j][k]->SetLineColor(kViolet); mpf_ratio[i][j][k]->SetMarkerColor(kViolet);}
	else if(i==2){  mpf_ratio[i][j][k]->SetLineColor(kGreen);  mpf_ratio[i][j][k]->SetMarkerColor(kGreen);}
	else if(i==3){  mpf_ratio[i][j][k]->SetLineColor(kRed);    mpf_ratio[i][j][k]->SetMarkerColor(kRed);}
	mpf_ratio[i][j][k]->Draw("PE SAME");
	c4->Modified();
      }
      c4->Update();
      
      c4->SaveAs(CorrectionObject::_input_path+"/Monitoring/MPF_res_" + CorrectionObject::_generator_tag + "_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+(eta_cut_bool?pt_range_HF:pt_range)[k]+"_"+(eta_cut_bool?pt_range_HF:pt_range)[k+1]+".pdf");
      delete c4;

      /*
      //Npv
      TCanvas* c5 = new TCanvas();
      tdrCanvas(c5,"c5",h,4,10,true,CorrectionObject::_lumitag);
      //   c5->cd(1);
      gStyle->SetOptFit(00000);
      
      h->GetXaxis()->SetTitle("Luminosity");
      h->GetXaxis()->SetTitleSize(0.036);
      h->GetXaxis()->SetTitleOffset(0.85);
      h->GetXaxis()->SetLimits(0,35000);
      h->GetXaxis()->SetLabelSize(0.036);
      h->GetYaxis()->SetRangeUser(0,20);
      h->GetYaxis()->SetTitle("Npv");

      TLegend *leg_Npv_res;
      leg_Npv_res = new TLegend(0.45,0.7,0.51,0.9,"","brNDC");//x+0.1
      leg_Npv_res->SetBorderSize(0);
      leg_Npv_res->SetTextSize(0.036);
      leg_Npv_res->SetFillColor(10);
      leg_Npv_res->SetFillStyle(1001);
      leg_Npv_res->SetLineColor(1);
      leg_Npv_res->SetTextFont(42);
      leg_Npv_res->SetHeader("NPV, "+eta_range_control[j]+"#leq#eta<"+eta_range_control[j+1]); //+", #alpha<"+s_blpha_cut
 
      for(int i = 0; i<n_input; i++){
	Npv_res[i][j][k]->SetMarkerStyle(2);
	if(i==0)     {  Npv_res[i][j][k]->SetLineColor(kBlue);   Npv_res[i][j][k]->SetMarkerColor(kBlue);}
	else if(i==1){  Npv_res[i][j][k]->SetLineColor(kViolet); Npv_res[i][j][k]->SetMarkerColor(kViolet);}
	else if(i==2){  Npv_res[i][j][k]->SetLineColor(kGreen);  Npv_res[i][j][k]->SetMarkerColor(kGreen);}
	else if(i==3){  Npv_res[i][j][k]->SetLineColor(kRed);    Npv_res[i][j][k]->SetMarkerColor(kRed);}
	leg_Npv_res->AddEntry(Npv_res[i][j][k], Name_range[i],"lep"); 
	
	Npv_res[i][j][k]->Draw("PE SAME");
	c5->Modified();
      }

      lineB->Draw("SAME");
      lineC->Draw("SAME");
      lineD->Draw("SAME");
      lineEFearly->Draw("SAME");
      lineFlateG->Draw("SAME");
      
      leg_Npv_res->Draw();
      
      tex1->DrawLatex(0.7,0.75, pt_range[k]+"<p_{T}<"+pt_range[k+1]);
      c5->Modified();
      c5->Update();
      */
      //   c5->cd(2);

      // d2->GetXaxis()->SetTitle("Luminosity");
      // d2->GetXaxis()->SetTitleSize(0.1);
      // d2->GetXaxis()->SetTitleOffset(1.2);
      // d2->GetXaxis()->SetLimits(0,35000);
      // d2->GetYaxis()->SetRangeUser(0.86,1.14);
      // d2->GetXaxis()->SetTitleOffset(1.2);
      // d2->GetYaxis()->SetTitle("Data/Fit");

      // lineB_ratio->Draw("SAME");
      // lineC_ratio->Draw("SAME");
      // lineD_ratio->Draw("SAME");
      // lineEFearly_ratio->Draw("SAME");
      // lineFlateG_ratio->Draw("SAME");
      // ave->Draw("SAME");
      
      // for(int i = 0; i<n_input; i++){
      // 	mpf_ratio[i][j][k]->SetMarkerStyle(2);
      // 	if(i==0)     {  mpf_ratio[i][j][k]->SetLineColor(kBlue);   mpf_ratio[i][j][k]->SetMarkerColor(kBlue);}
      // 	else if(i==1){  mpf_ratio[i][j][k]->SetLineColor(kViolet); mpf_ratio[i][j][k]->SetMarkerColor(kViolet);}
      // 	else if(i==2){  mpf_ratio[i][j][k]->SetLineColor(kGreen);  mpf_ratio[i][j][k]->SetMarkerColor(kGreen);}
      // 	else if(i==3){  mpf_ratio[i][j][k]->SetLineColor(kRed);    mpf_ratio[i][j][k]->SetMarkerColor(kRed);}
      // 	mpf_ratio[i][j][k]->Draw("PE SAME");
      // 	c5->Modified();
      // }
      // c5->Update();
      /*
      c5->SaveAs(CorrectionObject::_input_path+"/Monitoring/Npv_res_" + CorrectionObject::_generator_tag + "_eta_"+eta_range_control2[j]+"_"+eta_range_control2[j+1]+"_pT_"+pt_range[k]+"_"+pt_range[k+1]+".pdf");
      delete c5;
      */
    }//pt
  }//eta
}
