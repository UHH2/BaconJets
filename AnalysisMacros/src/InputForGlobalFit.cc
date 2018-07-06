// This script stores ratio between MC and DATA responses as function of pT in different eta bins and for different alpha values

#include "../include/parameters.h"
#include "../include/CorrectionObject.h"
#include "../include/useful_functions.h"

#include <iostream>
#include <cmath>
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"


using namespace std;


void CorrectionObject::InputForGlobalFit(){
  TStyle* m_gStyle = new TStyle();
  m_gStyle->SetOptFit(0);
  m_gStyle->SetOptStat(0);

  int n_pt_ = max(n_pt,n_pt_HF);
  bool eta_cut_bool;
  int n_pt_cutted;

  cout << "hello, it's macro to produce combination file for the Global Fit" << endl;
  //calculate ratio in MC to DATA responses

  double mc_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_mc_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double mc_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_mc_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double data_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_data_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double data_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_data_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; 

  double ratio_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_mpf[n_alpha_common][n_eta_common-1][n_pt_-1]; //error of ratio at pt,eta,alpha bins
  double ratio_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; //ratio at pt,eta,alpha bins
  double err_ratio_al_pTbal[n_alpha_common][n_eta_common-1][n_pt_-1]; //error of ratio at pt,eta,alpha bins

  TH1D *hdata_rel_r[n_pt_-1][n_eta_common-1][n_alpha_common];// pT-balance responce for data
  TH1D *hdata_mpf_r[n_pt_-1][n_eta_common-1][n_alpha_common];//MPF responce for data
  TH1D *hmc_rel_r[n_pt_-1][n_eta_common-1][n_alpha_common];// pT-balance responce for MC
  TH1D *hmc_mpf_r[n_pt_-1][n_eta_common-1][n_alpha_common];//MPF responce for MC

  TH1D *hmc_probejet_glu[n_alpha_common][n_eta_common-1];//   Math::Abs(pdgID)==21 --> glu
  TH1D *hmc_probejet_gluExt[n_alpha_common][n_eta_common-1];// 21||Undefined --> gluExt
  TH1D *hmc_probejet_b[n_alpha_common][n_eta_common-1];// 5--> b
  TH1D *hmc_probejet_c[n_alpha_common][n_eta_common-1];// 4-->c
  TH1D *hmc_probejet_s[n_alpha_common][n_eta_common-1];// 3--> s
  TH1D *hmc_probejet_ud[n_alpha_common][n_eta_common-1];// 1||2--> ud
  TH1D *hmc_probejet_undefined[n_alpha_common][n_eta_common-1];//   else --> Undefined
  TH1D *hmc_probejet_uds[n_alpha_common][n_eta_common-1];// 1||2||3 --> uds
  TH1D *hmc_probejet_total[n_alpha_common][n_eta_common-1];// total hist for normalisation

  TH1D *hmc_dijet_QQ[n_alpha_common][n_eta_common-1];// tag=quark(udcsb), probe=quark
  TH1D *hmc_dijet_GQ[n_alpha_common][n_eta_common-1];// tag=gluon, probe=quark(udcsb)
  TH1D *hmc_dijet_GG[n_alpha_common][n_eta_common-1];// tag=gluon, probe=gluon
  TH1D *hmc_dijet_QG[n_alpha_common][n_eta_common-1];// tag=quark(udcsb), probe=gluon
  TH1D *hmc_dijet_total[n_alpha_common][n_eta_common-1];// total hist for normalisation

  double flavorFrac_probejet_glu[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_glu[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_gluExt[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_gluExt[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_b[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_b[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_c[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_c[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_s[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_s[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_ud[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_ud[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_undefined[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_undefined[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_probejet_uds[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_probejet_uds[n_alpha_common][n_eta_common-1][n_pt_-1]; 

  double flavorFrac_dijet_QQ[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_dijet_QQ[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_dijet_GQ[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_dijet_GQ[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_dijet_GG[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_dijet_GG[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double flavorFrac_dijet_QG[n_alpha_common][n_eta_common-1][n_pt_-1]; 
  double err_flavorFrac_dijet_QG[n_alpha_common][n_eta_common-1][n_pt_-1]; 

  TH1D *hdata_num[n_alpha_common][n_eta_common-1];
  TH1D *hmc_num[n_alpha_common][n_eta_common-1];

  int count = 0;
  int count2= 0;
 
  TString name1 = "hist_data_rel_r_";
  TString name2 = "hist_data_mpf_r_";
  TString name3 = "hist_mc_rel_r_";
  TString name4 = "hist_mc_mpf_r_";
  TString name5 = "hist_data_num";
  TString name6 = "hist_mc_num";

  TString name11 = "hist_mc_probejet_glu";
  TString name12 = "hist_mc_probejet_gluExt";
  TString name13 = "hist_mc_probejet_b";
  TString name14 = "hist_mc_probejet_c";
  TString name15 = "hist_mc_probejet_s";
  TString name16 = "hist_mc_probejet_ud";
  TString name17 = "hist_mc_probejet_undefined";
  TString name18 = "hist_mc_probejet_uds";
  TString name19 = "hist_mc_probejet_total";

  TString name20 = "hist_mc_dijet_QQ";
  TString name21 = "hist_mc_dijet_GQ";
  TString name22 = "hist_mc_dijet_GG";
  TString name23 = "hist_mc_dijet_QG";
  TString name24 = "hist_mc_dijet_total";

  //Initialize with 0 values
  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      eta_cut_bool = fabs(eta_common_bins[j])>eta_cut;
      n_pt_cutted = ( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 );

      for(int k=0; k<n_pt_cutted; k++){

	ratio_al_mpf[i][j][k] = 0;
	err_ratio_al_mpf[i][j][k] = 0;
	ratio_al_pTbal[i][j][k] = 0;
	err_ratio_al_pTbal[i][j][k] = 0;
	mc_al_mpf[i][j][k] = 0;
	err_mc_al_mpf[i][j][k] = 0; 
	mc_al_pTbal[i][j][k] = 0;  
	err_mc_al_pTbal[i][j][k] = 0;
	data_al_mpf[i][j][k] = 0; 
	err_data_al_mpf[i][j][k] = 0; 
	data_al_pTbal[i][j][k] = 0; 
	err_data_al_pTbal[i][j][k] = 0; 

	flavorFrac_probejet_glu[i][j][k] = 0; 
	err_flavorFrac_probejet_glu[i][j][k] = 0; 
	flavorFrac_probejet_gluExt[i][j][k] = 0; 
	err_flavorFrac_probejet_gluExt[i][j][k] = 0; 
	flavorFrac_probejet_b[i][j][k] = 0; 
	err_flavorFrac_probejet_b[i][j][k] = 0; 
	flavorFrac_probejet_s[i][j][k] = 0; 
	err_flavorFrac_probejet_s[i][j][k] = 0; 
	flavorFrac_probejet_c[i][j][k] = 0; 
	err_flavorFrac_probejet_c[i][j][k] = 0; 
	flavorFrac_probejet_ud[i][j][k] = 0; 
	err_flavorFrac_probejet_ud[i][j][k] = 0; 
	flavorFrac_probejet_uds[i][j][k] = 0; 
	err_flavorFrac_probejet_uds[i][j][k] = 0; 
	flavorFrac_probejet_undefined[i][j][k] = 0; 
	err_flavorFrac_probejet_undefined[i][j][k] = 0; 
	flavorFrac_dijet_QQ[i][j][k] = 0; 
	err_flavorFrac_dijet_QQ[i][j][k] = 0; 
	flavorFrac_dijet_GQ[i][j][k] = 0; 
	err_flavorFrac_dijet_GQ[i][j][k] = 0; 
	flavorFrac_dijet_GG[i][j][k] = 0; 
	err_flavorFrac_dijet_GG[i][j][k] = 0; 
	flavorFrac_dijet_QG[i][j][k] = 0; 
	err_flavorFrac_dijet_QG[i][j][k] = 0; 

	TString name = name1; name+=count;
	hdata_rel_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name2;name+=count;
	hdata_mpf_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name3; name+=count;
	hmc_rel_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	name = name4; name+=count;
	hmc_mpf_r[k][j][i] = new TH1D(name,"",100, 0, 2.5);
	count++;

      }
      	TString name10 = name5; name10+=count2;
	hdata_num[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name6; name10+=count2;
	hmc_num[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name11; name10+=count2;
	hmc_probejet_glu[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name12; name10+=count2;
	hmc_probejet_gluExt[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name13; name10+=count2;
	hmc_probejet_b[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name14; name10+=count2;
	hmc_probejet_c[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name15; name10+=count2;
	hmc_probejet_s[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name16; name10+=count2;
	hmc_probejet_ud[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name17; name10+=count2;
	hmc_probejet_undefined[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name18; name10+=count2;
	hmc_probejet_uds[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name19; name10+=count2;
	hmc_probejet_total[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));

	name10 = name20; name10+=count2;
	hmc_dijet_QQ[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name21; name10+=count2;
	hmc_dijet_GQ[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name22; name10+=count2;
	hmc_dijet_GG[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name23; name10+=count2;
	hmc_dijet_QG[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));
	name10 = name24; name10+=count2;
	hmc_dijet_total[i][j] = new TH1D(name10,"",n_pt_cutted,(eta_cut_bool?pt_bins_HF:pt_bins));

	count2++;
    }
  }

 
 // Create the tree reader and its data containers
   TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
   TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
   TTreeReaderValue<Float_t> rel_r_data(myReader_DATA, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_data(myReader_DATA, "mpf_r");
   TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
   
   while (myReader_DATA.Next()) {
     for(int j=0; j<n_eta_common-1; j++){
       if(fabs(*probejet_eta_data)>eta_common_bins[j+1] || fabs(*probejet_eta_data)<eta_common_bins[j]) continue;
       eta_cut_bool = fabs(eta_common_bins[j])>eta_cut; 

       for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
   	   if(*pt_ave_data<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_data>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;

	     for(int i=0; i<n_alpha_common; i++){
	       if(*alpha_data>alpha_bins_common[i]) continue;
	       else{
		 hdata_rel_r[k][j][i]->Fill(*rel_r_data,*weight_data);
		 hdata_mpf_r[k][j][i]->Fill(*mpf_r_data,*weight_data);
	      }
	   }
	 }
     }
   }

   TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
   TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
   TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
   TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
   TTreeReaderValue<Float_t> rel_r_mc(myReader_MC, "rel_r");
   TTreeReaderValue<Float_t> mpf_r_mc(myReader_MC, "mpf_r");
   TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
   TTreeReaderValue<Int_t> flavorProbejet_mc(myReader_MC, "flavorProbejet"); 
   TTreeReaderValue<Int_t> flavorTagjet_mc(myReader_MC, "flavorBarreljet");

   while (myReader_MC.Next()) {
     for(int j=0; j<n_eta_common-1; j++){
       if(fabs(*probejet_eta_mc)>eta_common_bins[j+1] || fabs(*probejet_eta_mc)<eta_common_bins[j]) continue;
       eta_cut_bool = fabs(eta_common_bins[j])>eta_cut; 

       for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	 if(*pt_ave_mc<(eta_cut_bool?pt_bins_HF:pt_bins)[k] || *pt_ave_mc>(eta_cut_bool?pt_bins_HF:pt_bins)[k+1]) continue;
 
   	 for(int i=0; i<n_alpha_common; i++){
   	   if(*alpha_mc>alpha_bins_common[i]) continue;
   	   else{
   	     hmc_rel_r[k][j][i]->Fill(*rel_r_mc,*weight_mc);
   	     hmc_mpf_r[k][j][i]->Fill(*mpf_r_mc,*weight_mc);
	     if(*flavorProbejet_mc==21)
	       hmc_probejet_glu[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==21 || *flavorProbejet_mc<0)
	       hmc_probejet_gluExt[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==5)
	       hmc_probejet_b[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==4)
	       hmc_probejet_c[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==3)
	       hmc_probejet_s[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==1 || *flavorProbejet_mc==2)
	       hmc_probejet_ud[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc==1 || *flavorProbejet_mc==2 || *flavorProbejet_mc==3)
	       hmc_probejet_uds[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     if(*flavorProbejet_mc<0)
	       hmc_probejet_undefined[i][j]->Fill(*pt_ave_mc,*weight_mc);
	     hmc_probejet_total[i][j]->Fill(*pt_ave_mc,*weight_mc);

	     if(*flavorProbejet_mc>0 && *flavorProbejet_mc<6 && *flavorTagjet_mc>0 && *flavorTagjet_mc<6)
	       hmc_dijet_QQ[i][j]->Fill(*pt_ave_mc,*weight_mc); 
	      if(*flavorTagjet_mc==21 && *flavorProbejet_mc>0 && *flavorProbejet_mc<6 )
		hmc_dijet_GQ[i][j]->Fill(*pt_ave_mc,*weight_mc); 
	      if(*flavorTagjet_mc==21 && *flavorProbejet_mc==21)
		hmc_dijet_GG[i][j]->Fill(*pt_ave_mc,*weight_mc); 
	      if(*flavorTagjet_mc>0 && *flavorTagjet_mc<6 && *flavorProbejet_mc==21)
		hmc_dijet_QG[i][j]->Fill(*pt_ave_mc,*weight_mc); 
	      hmc_dijet_total[i][j]->Fill(*pt_ave_mc,*weight_mc); 

   	   }
   	 }
       }
     }
   }


  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      eta_cut_bool = fabs(eta_common_bins[j])>eta_cut; 
      for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	pair<double,double> res_mc_rel_r,res_data_rel_r;
	pair<double,double> res_mc_mpf_r,res_data_mpf_r;

	res_mc_rel_r = GetValueAndError(hmc_rel_r[k][j][i]);
	res_mc_mpf_r = GetValueAndError(hmc_mpf_r[k][j][i]);

	res_data_rel_r = GetValueAndError(hdata_rel_r[k][j][i]);
	res_data_mpf_r = GetValueAndError(hdata_mpf_r[k][j][i]);
	
	mc_al_mpf[i][j][k] = res_mc_mpf_r.first;
	err_mc_al_mpf[i][j][k] = res_mc_mpf_r.second;
	mc_al_pTbal[i][j][k] = res_mc_rel_r.first;
	err_mc_al_pTbal[i][j][k] = res_mc_rel_r.second;
	data_al_mpf[i][j][k] = res_data_mpf_r.first;
	err_data_al_mpf[i][j][k] = res_data_mpf_r.second;
	data_al_pTbal[i][j][k] = res_data_rel_r.first;
	err_data_al_pTbal[i][j][k] = res_data_rel_r.second;
	//	pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_mc_mpf_r,res_data_mpf_r);
	pair<double,double> ratio_res_mpf = Rmc_to_Rdata(res_data_mpf_r,res_mc_mpf_r); //Data/MC in sync with global fit framework
	ratio_al_mpf[i][j][k] = ratio_res_mpf.first;
	err_ratio_al_mpf[i][j][k] = ratio_res_mpf.second;
	//	pair<double,double> ratio_res_pTbal = Rmc_to_Rdata(res_mc_rel_r,res_data_rel_r);
	pair<double,double> ratio_res_pTbal = Rmc_to_Rdata(res_data_rel_r,res_mc_rel_r); //Data/MC in sync with global fit framework
	ratio_al_pTbal[i][j][k] = ratio_res_pTbal.first;
	err_ratio_al_pTbal[i][j][k] = ratio_res_pTbal.second;

	hmc_num[i][j]->SetBinContent(k+1, hmc_mpf_r[k][j][i]->GetEntries());
	hdata_num[i][j]->SetBinContent(k+1, hdata_mpf_r[k][j][i]->GetEntries());

      }
      hmc_probejet_glu[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_gluExt[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_b[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_c[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_s[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_ud[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_uds[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_probejet_undefined[i][j]->Divide(hmc_probejet_total[i][j]);
      hmc_dijet_QQ[i][j]->Divide(hmc_dijet_total[i][j]);
      hmc_dijet_GQ[i][j]->Divide(hmc_dijet_total[i][j]);
      hmc_dijet_GG[i][j]->Divide(hmc_dijet_total[i][j]);
      hmc_dijet_QG[i][j]->Divide(hmc_dijet_total[i][j]);

      for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	flavorFrac_probejet_glu[i][j][k] = hmc_probejet_glu[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_glu[i][j][k] = hmc_probejet_glu[i][j]->GetBinError(k);
	flavorFrac_probejet_gluExt[i][j][k] = hmc_probejet_gluExt[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_gluExt[i][j][k] = hmc_probejet_gluExt[i][j]->GetBinError(k);
	flavorFrac_probejet_b[i][j][k] = hmc_probejet_b[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_b[i][j][k] = hmc_probejet_b[i][j]->GetBinError(k);
	flavorFrac_probejet_c[i][j][k] = hmc_probejet_c[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_c[i][j][k] = hmc_probejet_c[i][j]->GetBinError(k);
	flavorFrac_probejet_s[i][j][k] = hmc_probejet_s[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_s[i][j][k] = hmc_probejet_s[i][j]->GetBinError(k);
	flavorFrac_probejet_ud[i][j][k] = hmc_probejet_ud[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_ud[i][j][k] = hmc_probejet_ud[i][j]->GetBinError(k);
	flavorFrac_probejet_uds[i][j][k] = hmc_probejet_uds[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_uds[i][j][k] = hmc_probejet_uds[i][j]->GetBinError(k);
	flavorFrac_probejet_undefined[i][j][k] = hmc_probejet_undefined[i][j]->GetBinContent(k);
	err_flavorFrac_probejet_undefined[i][j][k] = hmc_probejet_undefined[i][j]->GetBinError(k);

	flavorFrac_dijet_QQ[i][j][k] = hmc_dijet_QQ[i][j]->GetBinContent(k);
	err_flavorFrac_dijet_QQ[i][j][k] = hmc_dijet_QQ[i][j]->GetBinError(k);
	flavorFrac_dijet_GQ[i][j][k] = hmc_dijet_GQ[i][j]->GetBinContent(k);
	err_flavorFrac_dijet_GQ[i][j][k] = hmc_dijet_GQ[i][j]->GetBinError(k);
	flavorFrac_dijet_GG[i][j][k] = hmc_dijet_GG[i][j]->GetBinContent(k);
	err_flavorFrac_dijet_GG[i][j][k] = hmc_dijet_GG[i][j]->GetBinError(k);
	flavorFrac_dijet_QG[i][j][k] = hmc_dijet_QG[i][j]->GetBinContent(k);
	err_flavorFrac_dijet_QG[i][j][k] = hmc_dijet_QG[i][j]->GetBinError(k);

	//	cout<<"flavorFrac_probejet_glu[i][j][k] = "<<flavorFrac_probejet_glu[i][j][k];
	//	cout<<"err_flavorFrac_probejet_glu[i][j][k] = "<<err_flavorFrac_probejet_glu[i][j][k] <<endl;
      }  
    }
  }


  //Store results in TGraphErrors
  double res_xbin_tgraph[n_eta-1][n_pt_-1];// = {(pt_bins[0]+pt_bins[1])/2
  double res_zero[n_eta-1][n_pt_-1];
  for(int j=0; j<n_eta-1; j++){
    eta_cut_bool = fabs(eta_common_bins[j])>eta_cut; 
    for(int i=0;i<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 );i++){
      res_xbin_tgraph[j][i]=(( eta_cut_bool ?pt_bins_HF:pt_bins)[i]+( eta_cut_bool ?pt_bins_HF:pt_bins)[i+1])/2;
      res_zero[j][i]=(( eta_cut_bool ?pt_bins_HF:pt_bins)[i+1]-( eta_cut_bool ?pt_bins_HF:pt_bins)[i])/2 ;
    }
  }

  TGraphErrors *mpf_data[n_alpha_common][n_eta_common-1];
  TGraphErrors *rrel_data[n_alpha_common][n_eta_common-1];
  TGraphErrors *mpf_mc[n_alpha_common][n_eta_common-1];
  TGraphErrors *rrel_mc[n_alpha_common][n_eta_common-1];
  TGraphErrors *mpf_ratio[n_alpha_common][n_eta_common-1];
  TGraphErrors *rrel_ratio[n_alpha_common][n_eta_common-1];

  TGraphErrors *gr_flavorFrac_probejet_glu[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_gluExt[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_b[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_c[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_s[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_ud[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_uds[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_probejet_undefined[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_dijet_QQ[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_dijet_GQ[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_dijet_GG[n_alpha_common][n_eta_common-1];
  TGraphErrors *gr_flavorFrac_dijet_QG[n_alpha_common][n_eta_common-1];


  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],data_al_mpf[i][j],res_zero[j],err_data_al_mpf[i][j]);
      mpf_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],mc_al_mpf[i][j],res_zero[j],err_mc_al_mpf[i][j]);
      rrel_data[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],data_al_pTbal[i][j],res_zero[j],err_data_al_pTbal[i][j]);
      rrel_mc[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],mc_al_pTbal[i][j],res_zero[j],err_mc_al_pTbal[i][j]);
      mpf_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],ratio_al_mpf[i][j],res_zero[j],err_ratio_al_mpf[i][j]);
      rrel_ratio[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],ratio_al_pTbal[i][j],res_zero[j],err_ratio_al_pTbal[i][j]);

      gr_flavorFrac_probejet_glu[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_glu[i][j],res_zero[j],err_flavorFrac_probejet_glu[i][j]);
      gr_flavorFrac_probejet_gluExt[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_gluExt[i][j],res_zero[j],err_flavorFrac_probejet_gluExt[i][j]);
      gr_flavorFrac_probejet_b[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_b[i][j],res_zero[j],err_flavorFrac_probejet_b[i][j]);
      gr_flavorFrac_probejet_c[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_c[i][j],res_zero[j],err_flavorFrac_probejet_c[i][j]);
      gr_flavorFrac_probejet_s[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_s[i][j],res_zero[j],err_flavorFrac_probejet_s[i][j]);
      gr_flavorFrac_probejet_ud[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_ud[i][j],res_zero[j],err_flavorFrac_probejet_ud[i][j]);
      gr_flavorFrac_probejet_uds[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_uds[i][j],res_zero[j],err_flavorFrac_probejet_uds[i][j]);
      gr_flavorFrac_probejet_undefined[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_probejet_undefined[i][j],res_zero[j],err_flavorFrac_probejet_undefined[i][j]);

      gr_flavorFrac_dijet_QQ[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_dijet_QQ[i][j],res_zero[j],err_flavorFrac_dijet_QQ[i][j]);
      gr_flavorFrac_dijet_GQ[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_dijet_GQ[i][j],res_zero[j],err_flavorFrac_dijet_GQ[i][j]);
      gr_flavorFrac_dijet_GG[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_dijet_GG[i][j],res_zero[j],err_flavorFrac_dijet_GG[i][j]);
      gr_flavorFrac_dijet_QG[i][j] = new TGraphErrors(n_pt-1,res_xbin_tgraph[j],flavorFrac_dijet_QG[i][j],res_zero[j],err_flavorFrac_dijet_QG[i][j]);

    }
  }

  // Cleaning for empty points (with low statistic)
  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      mpf_data[i][j] = CleanEmptyPoints(mpf_data[i][j]);
      mpf_mc[i][j] = CleanEmptyPoints(mpf_mc[i][j]);
      rrel_data[i][j] = CleanEmptyPoints(rrel_data[i][j]);
      rrel_mc[i][j] = CleanEmptyPoints(rrel_mc[i][j]);
      mpf_ratio[i][j] = CleanEmptyPoints(mpf_ratio[i][j]);
      rrel_ratio[i][j] = CleanEmptyPoints(rrel_ratio[i][j]);

      hmc_num[i][j]->SetName("mc_RawNEvents_"+alpha_range_common[i]);
      hdata_num[i][j]->SetName("data_RawNEvents_"+alpha_range_common[i]);

      mpf_data[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      mpf_mc[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_data[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      rrel_mc[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);
      mpf_ratio[i][j]->SetName("mpfchs_dijet_"+alpha_range_common[i]);
      rrel_ratio[i][j]->SetName("ptchs_dijet_"+alpha_range_common[i]);

      gr_flavorFrac_probejet_glu[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_glu[i][j]);
      gr_flavorFrac_probejet_glu[i][j]->SetName("dijet_flavFraction_probejet_glu_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_gluExt[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_gluExt[i][j]);
      gr_flavorFrac_probejet_gluExt[i][j]->SetName("dijet_flavFraction_probejet_gluExt_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_b[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_b[i][j]);
      gr_flavorFrac_probejet_b[i][j]->SetName("dijet_flavFraction_probejet_b_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_c[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_c[i][j]);
      gr_flavorFrac_probejet_c[i][j]->SetName("dijet_flavFraction_probejet_c_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_s[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_s[i][j]);
      gr_flavorFrac_probejet_s[i][j]->SetName("dijet_flavFraction_probejet_s_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_ud[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_ud[i][j]);
      gr_flavorFrac_probejet_ud[i][j]->SetName("dijet_flavFraction_probejet_ud_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_uds[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_uds[i][j]);
      gr_flavorFrac_probejet_uds[i][j]->SetName("dijet_flavFraction_probejet_uds_"+alpha_range_common[i]);
      gr_flavorFrac_probejet_undefined[i][j] = CleanEmptyPoints(gr_flavorFrac_probejet_undefined[i][j]);
      gr_flavorFrac_probejet_undefined[i][j]->SetName("dijet_flavFraction_probejet_undefined_"+alpha_range_common[i]);

      gr_flavorFrac_dijet_QQ[i][j] = CleanEmptyPoints(gr_flavorFrac_dijet_QQ[i][j]);
      gr_flavorFrac_dijet_QQ[i][j]->SetName("dijet_flavFraction_dijet_QQ_"+alpha_range_common[i]);
      gr_flavorFrac_dijet_GQ[i][j] = CleanEmptyPoints(gr_flavorFrac_dijet_GQ[i][j]);
      gr_flavorFrac_dijet_GQ[i][j]->SetName("dijet_flavFraction_dijet_GQ_"+alpha_range_common[i]);
      gr_flavorFrac_dijet_GG[i][j] = CleanEmptyPoints(gr_flavorFrac_dijet_GG[i][j]);
      gr_flavorFrac_dijet_GG[i][j]->SetName("dijet_flavFraction_dijet_GG_"+alpha_range_common[i]);
      gr_flavorFrac_dijet_QG[i][j] = CleanEmptyPoints(gr_flavorFrac_dijet_QG[i][j]);
      gr_flavorFrac_dijet_QG[i][j]->SetName("dijet_flavFraction_dijet_QG_"+alpha_range_common[i]);

      }
    }

  //Save results in root file
  TFile* outputfile = new TFile(CorrectionObject::_outpath+"output/JEC_L2_Dijet_"+CorrectionObject::_jettag+"_"+CorrectionObject::_generator_tag+".root","RECREATE");
  outputfile->Print();

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      if(i==0){
	outputfile->mkdir("ratio/"+eta_output[j]);
	outputfile->mkdir("data/"+eta_output[j]);
	outputfile->mkdir("mc/"+eta_output[j]);
      }
      outputfile->cd("ratio/"+eta_output[j]);
      mpf_ratio[i][j]->Write();
      rrel_ratio[i][j]->Write();
    
      outputfile->cd("data/"+eta_output[j]);
      mpf_data[i][j]->Write();
      rrel_data[i][j]->Write();
      hdata_num[i][j]->Write();
       
      outputfile->cd("mc/"+eta_output[j]);
      mpf_mc[i][j]->Write();
      rrel_mc[i][j]->Write();
      hmc_num[i][j]->Write();
      gr_flavorFrac_probejet_glu[i][j]->Write();
      gr_flavorFrac_probejet_gluExt[i][j]->Write();
      gr_flavorFrac_probejet_b[i][j]->Write();
      gr_flavorFrac_probejet_c[i][j]->Write();
      gr_flavorFrac_probejet_s[i][j]->Write();
      gr_flavorFrac_probejet_ud[i][j]->Write();
      gr_flavorFrac_probejet_uds[i][j]->Write();
      gr_flavorFrac_probejet_undefined[i][j]->Write();
      gr_flavorFrac_dijet_QQ[i][j]->Write();
      gr_flavorFrac_dijet_GQ[i][j]->Write();
      gr_flavorFrac_dijet_GG[i][j]->Write();
      gr_flavorFrac_dijet_QG[i][j]->Write();
    }
 }

  cout<<"Draw result for alpha = "<<alpha_range_common[3]<<" eta = "<<eta_common_range[n_eta_common-2]<<" "<<eta_common_range[n_eta_common-1]<<endl;
  mpf_ratio[3][0]->Draw();

  outputfile->Write();
  outputfile->Close();




  //delete everything

  delete outputfile;

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      delete mpf_data[i][j];
      delete mpf_mc[i][j];
      delete rrel_data[i][j];
      delete rrel_mc[i][j];
      delete mpf_ratio[i][j];
      delete rrel_ratio[i][j];
    }
  }

  for(int i=0; i<n_alpha_common; i++){
    for(int j=0; j<n_eta_common-1; j++){
      eta_cut_bool = fabs(eta_common_bins[j])>eta_cut;
      for(int k=0; k<( eta_cut_bool ?  n_pt_HF-1 : n_pt-1 ); k++){
	delete hdata_rel_r[k][j][i];
	delete hdata_mpf_r[k][j][i];
	delete hmc_rel_r[k][j][i];
	delete hmc_mpf_r[k][j][i];
      }
    }
  }

  delete m_gStyle;
}
