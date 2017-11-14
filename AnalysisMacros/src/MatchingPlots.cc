#include "../include/parameters.h"
#include "../include/useful_functions.h"
#include "../include/CorrectionObject.h"
#include "../include/tdrstyle_mod15.h"

#include <TStyle.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TFile.h>


using namespace std;

void CorrectionObject::MatchingPlots(){
  cout << "--------------- Starting MatchingPlots() ---------------" << endl << endl;
  gStyle->SetOptStat(0);

  //Table with number of events in each pT- and eta-bin
  
  //Set up histos for ratios of responses
  TH1D *hmc_asymmetry[n_pt-1][n_eta-1];              // A for MC, all
  TH1D *hmc_asymmetry_jet_matched[n_pt-1][n_eta-1];  // A for MC, matched RECO and Gen Jet
  TH1D *hmc_asymmetry_matched[n_pt-1][n_eta-1];     // A for MC, matched RECO and Gen Jet + GenParticle

  TH1D *hdata_asymmetry[n_pt-1][n_eta-1];              // A for DATA, all

  TH1D *hmc_B[n_pt-1][n_eta-1];                     // B for MC, all
  TH1D *hmc_B_jet_matched[n_pt-1][n_eta-1];         // B for MC, matched RECO and Gen Jet
  TH1D *hmc_B_matched[n_pt-1][n_eta-1];             // B for MC, matched RECO and Gen Jet + GenParticle

  TH1D *hmc_jet3pt_matched[n_pt-1][n_eta-1];   // jet3 pt for MC, matched
  TH1D *hmc_jet3pt[n_pt-1][n_eta-1];   // jet3 pt for MC, matched
  TH1D *hdata_jet3pt[n_pt-1][n_eta-1];   // jet3 pt for DATA

  TH1D *hmc_zPV_matched[n_pt-1][n_eta-1];   // zPV for MC, matched
  TH1D *hmc_zPV[n_pt-1][n_eta-1];   // zPV for MC, matched
  TH1D *hdata_zPV[n_pt-1][n_eta-1];   // zPV for DATA

  TH1D *hmc_dRjet3_matched[n_pt-1][n_eta-1];   // dR between RECO jet3 and the closest jet for MC, matched in tag&probe
  TH1D *hmc_dRjet3[n_pt-1][n_eta-1];   // dR between RECO jet3 and the closest jet for MC,all

  int count = 0;
 
  TString name3 = "hist_mc_A_";
  TString name4 = "hist_mc_B_";
  TString name5 = "hist_mc_jet3pt_";
  TString name6 = "hist_data_jet3pt_";

  TString name7 = "hist_mc_zPV_";
  TString name8 = "hist_data_zPV_";

  TString name9 = "hist_mc_dRjet3_";

  TString name10 = "hist_data_A_";
 
  for(int j=0; j<n_eta-1; j++){
      TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
    for(int k=0; k<n_pt-1; k++){
      TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
      TString name = name3 + eta_name + "_" + pt_name;
      hmc_asymmetry[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name +="_matched_jet";
      hmc_asymmetry_jet_matched[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name3 + eta_name + "_" + pt_name + "matched_particle";
      hmc_asymmetry_matched[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);

      name = name4 + eta_name + "_" + pt_name;
      hmc_B[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name +="_matched_jet";
      hmc_B_jet_matched[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);
      name = name4 + eta_name + "_" + pt_name+ "matched_particle";
      hmc_B_matched[k][j] = new TH1D(name,"",nResponseBins, -1.2, 1.2);

      name = name5 + eta_name + "_" + pt_name;
      hmc_jet3pt[k][j] = new TH1D(name,"",50,0,pt_bins[k+1]);
      name +="_matched";
      hmc_jet3pt_matched[k][j] = new TH1D(name,"",50,0,pt_bins[k+1]);

      name = name6 + eta_name + "_" + pt_name;
      hdata_jet3pt[k][j] = new TH1D(name,"",50,0,pt_bins[k+1]);

      name = name7 + eta_name + "_" + pt_name;
      hmc_zPV[k][j] = new TH1D(name,"",100,-20,20);
      name +="_matched";
      hmc_zPV_matched[k][j] = new TH1D(name,"",100,-20,20);
      name = name8 + eta_name + "_" + pt_name;
      hdata_zPV[k][j] = new TH1D(name,"",100,-20,20);

      name = name9 + eta_name + "_" + pt_name;
      hmc_dRjet3[k][j] = new TH1D(name,"",100,0,10);
      name +="_matched";
      hmc_dRjet3_matched[k][j] = new TH1D(name,"",100,0,10);

      name = name10 + eta_name + "_" + pt_name;
      hdata_asymmetry[k][j]=new TH1D(name,"",nResponseBins, -1.2, 1.2);

      count++;
    }
  }



  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_mc(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_mc(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> asymmetry_mc(myReader_MC, "asymmetry");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> MET_mc(myReader_MC, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_mc(myReader_MC, "sum_jets_pt");
  TTreeReaderValue<Float_t> jet3_pt_mc(myReader_MC, "jet3_pt");
  TTreeReaderValue<Float_t> probejet_neutEmEF_mc(myReader_MC, "probejet_neutEmEF");
  TTreeReaderValue<Float_t> probejet_neutHadEF_mc(myReader_MC, "probejet_neutHadEF");
  TTreeReaderValue<Float_t> probejet_chEmEF_mc(myReader_MC, "probejet_chEmEF");
  TTreeReaderValue<Float_t> probejet_chHadEF_mc(myReader_MC, "probejet_chHadEF");
  TTreeReaderValue<Float_t> probejet_photonEF_mc(myReader_MC, "probejet_photonEF");
  TTreeReaderValue<Float_t> probejet_muonEF_mc(myReader_MC, "probejet_muonEF");
  TTreeReaderValue<Float_t> probejet_phi_mc(myReader_MC, "probejet_phi");
  TTreeReaderValue<Float_t> probejet_ptgen_mc(myReader_MC, "probejet_ptgen");
  TTreeReaderValue<Float_t> barreljet_ptgen_mc(myReader_MC, "barreljet_ptgen");
  TTreeReaderValue<Float_t> zPV_mc(myReader_MC, "Zpv");
  TTreeReaderValue<Int_t> flavorBarreljet_mc(myReader_MC, "flavorBarreljet");
  TTreeReaderValue<Int_t> flavorProbejet_mc(myReader_MC, "flavorProbejet");

  TTreeReaderValue<Float_t> dRjet3_mc(myReader_MC, "dR_jet3_RECO_GEN");

  

  while (myReader_MC.Next()) {
    if(*alpha_mc>alpha_cut) continue;
    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(*pt_ave_mc<pt_bins[k] || *pt_ave_mc>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	else{
	  if(*probejet_ptgen_mc>0 && *barreljet_ptgen_mc>0){    //RECO and Gen Jet Matched
	    hmc_asymmetry_jet_matched[k][j]->Fill(*asymmetry_mc,*weight_mc);
	    hmc_B_jet_matched[k][j]->Fill(*B_mc,*weight_mc);
	 }
	  if(*probejet_ptgen_mc>0 && *flavorProbejet_mc>0 && *barreljet_ptgen_mc>0 && *flavorBarreljet_mc>0 ){ //RECO and Gen Jet + GenParticle Matched
	    hmc_asymmetry_matched[k][j]->Fill(*asymmetry_mc,*weight_mc);
	    hmc_B_matched[k][j]->Fill(*B_mc,*weight_mc);
	    hmc_jet3pt_matched[k][j]->Fill(*jet3_pt_mc,*weight_mc);
	    hmc_zPV_matched[k][j]->Fill(*zPV_mc,*weight_mc);
	    hmc_dRjet3_matched[k][j]->Fill(*dRjet3_mc,*weight_mc);
	  }
	  hmc_asymmetry[k][j]->Fill(*asymmetry_mc,*weight_mc);
	  hmc_B[k][j]->Fill(*B_mc,*weight_mc);
	  hmc_jet3pt[k][j]->Fill(*jet3_pt_mc,*weight_mc);
	  hmc_zPV[k][j]->Fill(*zPV_mc,*weight_mc);
	  hmc_dRjet3[k][j]->Fill(*dRjet3_mc,*weight_mc);
	  //	  cout<<" dRjet3 = "<<*dRjet3_mc<<endl;
	}
      }
    }
  }
  


 //Get relevant information from DATA, loop over DATA events
  TTreeReader myReader_DATA("AnalysisTree", CorrectionObject::_DATAFile);
  TTreeReaderValue<Float_t> pt_ave_data(myReader_DATA, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_data(myReader_DATA, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_data(myReader_DATA, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_data(myReader_DATA, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_data(myReader_DATA, "alpha");
  TTreeReaderValue<Float_t> asymmetry_data(myReader_DATA, "asymmetry");
  TTreeReaderValue<Float_t> B_data(myReader_DATA, "B");
  TTreeReaderValue<Float_t> weight_data(myReader_DATA, "weight");
  TTreeReaderValue<Float_t> MET_data(myReader_DATA, "MET");
  TTreeReaderValue<Float_t> sum_jets_pt_data(myReader_DATA, "sum_jets_pt");
  TTreeReaderValue<Float_t> jet3_pt_data(myReader_DATA, "jet3_pt");
  // TTreeReaderValue<Float_t> zPV_data(myReader_DATA, "Zpv");

  while (myReader_DATA.Next()) {
    cout<<"Enter my Reader DATA!"<<endl;
    if(*alpha_data>alpha_cut) continue;
    //fill histos in bins of pt and eta
    for(int k=0; k<n_pt-1; k++){
      if(*pt_ave_data<pt_bins[k] || *pt_ave_data>pt_bins[k+1]) continue;
      for(int j=0; j<n_eta-1; j++){
	if(fabs(*probejet_eta_data)>eta_bins[j+1] || fabs(*probejet_eta_data)<eta_bins[j]) continue;
	else{
	  hdata_asymmetry[k][j]->Fill(*asymmetry_data, *weight_data);
	  hdata_jet3pt[k][j]->Fill(*jet3_pt_data,*weight_data);
	  // hdata_zPV[k][j]->Fill(*zPV_data,*weight_data);
	}
      }
    }
  }



 

  ofstream output;
  output.open(CorrectionObject::_outpath+"plots/control/Matched_Number_Events_Pt_Eta_bins_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".txt");
    

  output << "Number of matched events in each bin for MC" << endl;
  output << "|Eta|:          ";
  double n_tot_MC = 0;
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  for(int i=0; i<n_pt-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(j!=n_eta-2){
	if(hmc_B[i][j]->GetEntries()/1000 < 0.01)     output << hmc_B_matched[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 0.1) output << hmc_B_matched[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 1)   output << hmc_B_matched[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 <10)   output << hmc_B_matched[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_B[i][j]->GetEntries()/1000 <100)  output << hmc_B_matched[i][j]->GetEntries() << "  - ";
	else                                              output << hmc_B_matched[i][j]->GetEntries() << " - ";
      }
      else output << hmc_B_matched[i][j]->GetEntries() << endl;
      n_tot_MC+= hmc_B_matched[i][j]->GetEntries();
    }

  }
  output << endl << endl << "Total number of matched events in MC: " << n_tot_MC << endl;

  output << "Number of all events in each bin for MC" << endl;
  output << "|Eta|:          ";
  n_tot_MC = 0;
  for(int i=0; i<n_eta; i++) {
    if(i != n_eta-1) output << eta_range[i] << " -- ";
    else output << eta_range[i] << endl;
  }
  for(int i=0; i<n_pt-1; i++){
    if(i==0) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]  :    ";
    else if(i==1) output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "] :    ";
    else output << "pT = ["  << fixed << setprecision(0) << pt_bins[i] << "," << setprecision(0) << pt_bins[i+1] << "]:    ";

    for(int j=0; j<n_eta-1; j++){
      if(j!=n_eta-2){
	if(hmc_B[i][j]->GetEntries()/1000 < 0.01)     output << hmc_B[i][j]->GetEntries() << "      - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 0.1) output << hmc_B[i][j]->GetEntries() << "     - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 < 1)   output << hmc_B[i][j]->GetEntries() << "    - "; //<1000
	else if(hmc_B[i][j]->GetEntries()/1000 <10)   output << hmc_B[i][j]->GetEntries() << "   - "; //<10000
	else if(hmc_B[i][j]->GetEntries()/1000 <100)  output << hmc_B[i][j]->GetEntries() << "  - ";
	else                                              output << hmc_B[i][j]->GetEntries() << " - ";
      }
      else output << hmc_B[i][j]->GetEntries() << endl;
      n_tot_MC+= hmc_B[i][j]->GetEntries();
    }

  }
  output << endl << endl << "Total number of all events in MC: " << n_tot_MC << endl;
 

  // Dump 1-d distributions of A and B in bins of pT, eta

  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc_matched.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){     ///k=0 n_pt-1 
      hmc_B[k][j]->Write();
      hmc_B_jet_matched[k][j]->Write();
      hmc_B_matched[k][j]->Write();

      hmc_dRjet3[k][j]->Write();
      hmc_dRjet3_matched[k][j]->Write();
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;

 
  TFile* test_out_mc_A = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc_matched.root","RECREATE");
  for(int j=0; j<n_eta-1; j++){
    for(int k=0; k<n_pt-1; k++){
      hmc_asymmetry[k][j]->Write();
      hmc_asymmetry_jet_matched[k][j]->Write();
      hmc_asymmetry_matched[k][j]->Write();

      hmc_jet3pt[k][j]->Write();
      hmc_jet3pt_matched[k][j]->Write();

      hdata_jet3pt[k][j]->Write();
      hmc_zPV[k][j]->Write();
      hmc_zPV_matched[k][j]->Write();
      hdata_zPV[k][j]->Write();
    
      hdata_asymmetry[k][j]->Write();
    }
  }
  test_out_mc_A->Close();
  delete test_out_mc_A;

  

  //R_MC and R_DATA overlaid in the same plot as a function of pT, in bins of |eta|
  double val_rel_mc[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_mc[n_eta-1][n_pt-1]; //ratio at pt,eta
  double err_mpf_mc[n_eta-1][n_pt-1]; //error of ratio at pt,eta
 
  double val_rel_mc_matched[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_mc_matched[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_mc_matched[n_eta-1][n_pt-1]; //ratio at pt,eta
  double err_mpf_mc_matched[n_eta-1][n_pt-1]; //error of ratio at pt,eta

  double val_rel_mc_jet_matched[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_mc_jet_matched[n_eta-1][n_pt-1]; //error of ratio at pt,eta
  double val_mpf_mc_jet_matched[n_eta-1][n_pt-1]; //ratio at pt,eta
  double err_mpf_mc_jet_matched[n_eta-1][n_pt-1]; //error of ratio at pt,eta

  double val_rel_data[n_eta-1][n_pt-1]; //value at pt,eta
  double err_rel_data[n_eta-1][n_pt-1]; //error of ratio at pt,eta
 

  for(int i=0; i<n_eta-1; i++){
    for(int j=0; j<n_pt-1; j++){

      //get <A> and error on <A>
      pair <double,double> A_mc = GetValueAndError(hmc_asymmetry[j][i]);
      pair <double,double> B_mc = GetValueAndError(hmc_B[j][i]);
       if(i==(n_eta-2))
	 cout<<"j = "<<j<<" A_mc.first = "<<A_mc.first<<" A_mc.second = "<<A_mc.second<<endl;

      pair <double,double> A_data = GetValueAndError(hdata_asymmetry[j][i]);
 

      pair <double,double> A_mc_matched = GetValueAndError(hmc_asymmetry_matched[j][i]);
      pair <double,double> B_mc_matched = GetValueAndError(hmc_B_matched[j][i]);

      pair <double,double> A_mc_jet_matched = GetValueAndError(hmc_asymmetry_jet_matched[j][i]);
      pair <double,double> B_mc_jet_matched = GetValueAndError(hmc_B_jet_matched[j][i]);

      //build MPF and pt_bal and their errors


      pair<double,double> res_mc_rel_r,res_mc_rel_r_matched ,res_mc_rel_r_jet_matched;
      pair<double,double> res_mc_mpf_r,res_mc_mpf_r_matched, res_mc_mpf_r_jet_matched;
      pair<double,double> res_data_rel_r;
      res_mc_mpf_r.first = (1+B_mc.first)/(1-B_mc.first);
      res_mc_mpf_r.second = 2/(pow((1-B_mc.first),2)) * B_mc.second; 
      res_mc_rel_r.first = (1+A_mc.first)/(1-A_mc.first);
      res_mc_rel_r.second = 2/(pow((1-A_mc.first),2)) * A_mc.second;

      res_data_rel_r.first = (1+A_data.first)/(1-A_data.first);
      res_data_rel_r.second = 2/(pow((1-A_data.first),2)) * A_data.second;
    
      res_mc_mpf_r_matched.first = (1+B_mc_matched.first)/(1-B_mc_matched.first);
      res_mc_mpf_r_matched.second = 2/(pow((1-B_mc_matched.first),2)) * B_mc_matched.second; 
      res_mc_rel_r_matched.first = (1+A_mc_matched.first)/(1-A_mc_matched.first);
      res_mc_rel_r_matched.second = 2/(pow((1-A_mc_matched.first),2)) * A_mc_matched.second;

      res_mc_mpf_r_jet_matched.first = (1+B_mc_jet_matched.first)/(1-B_mc_jet_matched.first);
      res_mc_mpf_r_jet_matched.second = 2/(pow((1-B_mc_jet_matched.first),2)) * B_mc_jet_matched.second; 
      res_mc_rel_r_jet_matched.first = (1+A_mc_jet_matched.first)/(1-A_mc_jet_matched.first);
      res_mc_rel_r_jet_matched.second = 2/(pow((1-A_mc_jet_matched.first),2)) * A_mc_jet_matched.second;
    
      val_rel_mc[i][j] = res_mc_rel_r.first;
      err_rel_mc[i][j] = res_mc_rel_r.second;
      if(i==(n_eta-2))
       	cout<<"j = "<<j<<"  val_rel_mc[i][j] = "<< val_rel_mc[i][j]<<endl;

      val_mpf_mc[i][j] = res_mc_mpf_r.first;
      err_mpf_mc[i][j] = res_mc_mpf_r.second;

      val_rel_data[i][j] = res_data_rel_r.first;
      err_rel_data[i][j] = res_data_rel_r.second;
  
      val_rel_mc_matched[i][j] = res_mc_rel_r_matched.first;
      err_rel_mc_matched[i][j] = res_mc_rel_r_matched.second;
      val_mpf_mc_matched[i][j] = res_mc_mpf_r_matched.first;
      err_mpf_mc_matched[i][j] = res_mc_mpf_r_matched.second;

      val_rel_mc_jet_matched[i][j] = res_mc_rel_r_jet_matched.first;
      err_rel_mc_jet_matched[i][j] = res_mc_rel_r_jet_matched.second;
      val_mpf_mc_jet_matched[i][j] = res_mc_mpf_r_jet_matched.first;
      err_mpf_mc_jet_matched[i][j] = res_mc_mpf_r_jet_matched.second;
    }
  }

  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TCanvas* c_0 = new TCanvas();
  tdrCanvas(c_0,"c_0",h,4,10,true,CorrectionObject::_lumitag);

  
  for(int i=0; i<n_eta-1; i++){
    //Create and fill TGraphErrors
    double xbin_tgraph[n_pt-1];
    double zero[n_pt-1];
    for(int i=0;i<n_pt-1;i++){
      xbin_tgraph[i]=(pt_bins[i]+pt_bins[i+1])/2;
      zero[i]=(pt_bins[i+1]-pt_bins[i])/2 ;
    }


 

    TGraphErrors *graph_mpf_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_mc[i], zero, err_mpf_mc[i]);
    TGraphErrors *graph_mpf_mc_matched = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_mc_matched[i], zero, err_mpf_mc_matched[i]);
    TGraphErrors *graph_mpf_mc_jet_matched = new TGraphErrors(n_pt-1, xbin_tgraph, val_mpf_mc_jet_matched[i], zero, err_mpf_mc_jet_matched[i]);

    TGraphErrors *graph_rel_mc   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc[i], zero, err_rel_mc[i]);
    TGraphErrors *graph_rel_mc_matched = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc_matched[i], zero, err_rel_mc_matched[i]);
    TGraphErrors *graph_rel_mc_jet_matched = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_mc_jet_matched[i], zero, err_rel_mc_jet_matched[i]);

    TGraphErrors *graph_rel_data   = new TGraphErrors(n_pt-1, xbin_tgraph, val_rel_data[i], zero, err_rel_data[i]);

    graph_mpf_mc   = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc);
    graph_mpf_mc_matched = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc_matched);
    graph_mpf_mc_jet_matched = (TGraphErrors*)CleanEmptyPoints(graph_mpf_mc_jet_matched);

    graph_rel_mc   = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc);
    graph_rel_mc_matched = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc_matched);
    graph_rel_mc_jet_matched = (TGraphErrors*)CleanEmptyPoints(graph_rel_mc_jet_matched);

    graph_rel_data   = (TGraphErrors*)CleanEmptyPoints(graph_rel_data);


    graph_mpf_mc->SetTitle("");
    graph_mpf_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_mc->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_mpf_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_mc->SetMarkerColor(kRed);
    graph_mpf_mc->SetMarkerStyle(3);
    graph_mpf_mc->SetLineColor(kRed);

    graph_mpf_mc_matched->SetTitle("");
    graph_mpf_mc_matched->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_mc_matched->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_mc_matched->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_mc_matched->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_mpf_mc_matched->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_mc_matched->SetMarkerColor(kBlack);
    graph_mpf_mc_matched->SetMarkerStyle(26);
    graph_mpf_mc_matched->SetLineColor(kBlack);

    graph_mpf_mc_jet_matched->SetTitle("");
    graph_mpf_mc_jet_matched->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_mpf_mc_jet_matched->GetXaxis()->SetTitleSize(0.05);
    graph_mpf_mc_jet_matched->GetXaxis()->SetTitleOffset(0.80);
    graph_mpf_mc_jet_matched->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_mpf_mc_jet_matched->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_mpf_mc_jet_matched->SetMarkerColor(kBlue);
    graph_mpf_mc_jet_matched->SetMarkerStyle(32);
    graph_mpf_mc_jet_matched->SetLineColor(kBlue);

    graph_rel_mc->SetTitle("");
    graph_rel_mc->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_mc->GetXaxis()->SetTitleSize(0.05);
    graph_rel_mc->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_mc->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_mc->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc->SetMarkerColor(kRed);
    graph_rel_mc->SetMarkerStyle(3);
    graph_rel_mc->SetLineColor(kRed);

    graph_rel_mc_matched->SetTitle("");
    graph_rel_mc_matched->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_mc_matched->GetXaxis()->SetTitleSize(0.05);
    graph_rel_mc_matched->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_mc_matched->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_mc_matched->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc_matched->SetMarkerColor(kGreen);
    graph_rel_mc_matched->SetMarkerStyle(26);
    graph_rel_mc_matched->SetLineColor(kGreen);

    graph_rel_mc_jet_matched->SetTitle("");
    graph_rel_mc_jet_matched->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_mc_jet_matched->GetXaxis()->SetTitleSize(0.05);
    graph_rel_mc_jet_matched->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_mc_jet_matched->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_mc_jet_matched->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_mc_jet_matched->SetMarkerColor(kBlue);
    graph_rel_mc_jet_matched->SetMarkerStyle(32);
    graph_rel_mc_jet_matched->SetLineColor(kBlue);

    graph_rel_data->SetTitle("");
    graph_rel_data->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    graph_rel_data->GetXaxis()->SetTitleSize(0.05);
    graph_rel_data->GetXaxis()->SetTitleOffset(0.80);
    graph_rel_data->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    graph_rel_data->GetYaxis()->SetRangeUser(0.70,1.30);
    graph_rel_data->SetMarkerColor(kBlack);
    graph_rel_data->SetMarkerStyle(5);
    graph_rel_data->SetLineColor(kBlack);



    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    TString axistitle_mc = "R^{MC All}";
    TString axistitle_mc_matched = "R^{MC ParticleMatched}";
    TString axistitle_mc_jet_matched = "R^{MC JetMatched}";
    TString axistitle_data = "R^{Data}";

    /*
    axistitle_mc   += altitle;
    axistitle_mc_jet_matched += altitle;
    axistitle_mc_matched += altitle;
    axistitle_data += altitle;
    */

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 



    TCanvas* c_mpf = new TCanvas();
    tdrCanvas(c_mpf,"c_mpf",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_mpf_mc->Draw("P SAME");
    graph_mpf_mc_jet_matched->Draw("P SAME");
    graph_mpf_mc_matched->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_mpf;
    leg_mpf = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_mpf->SetBorderSize(0);
    leg_mpf->SetTextSize(0.036);
    leg_mpf->SetFillColor(10);
    leg_mpf->SetFillStyle(0);
    leg_mpf->SetLineColor(1);
    leg_mpf->SetTextFont(42);
    leg_mpf->SetHeader("MPF response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+", #alpha<"+s_alpha_cut);
    leg_mpf->AddEntry(graph_mpf_mc, axistitle_mc,"P");
    leg_mpf->AddEntry(graph_mpf_mc_jet_matched, axistitle_mc_jet_matched, "P");
    leg_mpf->AddEntry(graph_mpf_mc_matched, axistitle_mc_matched,"P");
    leg_mpf->Draw();

    //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");

    c_mpf->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_MPF_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("#bar{p}_{T} [GeV]");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.80);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.50,1.50);
    graph_rel_data->Draw("P SAME");
    graph_rel_mc->Draw("P SAME");
    graph_rel_mc_jet_matched->Draw("P SAME");
    graph_rel_mc_matched->Draw("P SAME");
    gPad->SetLogx();

    TLegend *leg_rel;
    leg_rel = new TLegend(0.35,0.72,0.51,0.92,"","brNDC");//x+0.1
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.036);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("p_{T}-balance response, "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]); 
    leg_rel->AddEntry(graph_rel_data, axistitle_data, "P");
    leg_rel->AddEntry(graph_rel_mc, axistitle_mc,"P");
    leg_rel->AddEntry(graph_rel_mc_jet_matched, axistitle_mc_jet_matched,"P");
    leg_rel->AddEntry(graph_rel_mc_matched, axistitle_mc_matched,"P");
    leg_rel->Draw();

    //tex->DrawLatex(0.53,0.91,CorrectionObject::_lumitag+"(13TeV)");

    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_Rel_Response_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] + ".pdf");


    //delete leg_rel;
    delete c_rel;
    //delete leg_mpf;
    delete c_mpf;
    delete tex;
    delete graph_rel_mc_matched;
    delete graph_rel_mc_jet_matched;
    delete graph_rel_mc;
    delete graph_rel_data;
    delete graph_mpf_mc_matched;
    delete graph_mpf_mc_jet_matched;
    delete graph_mpf_mc;
  }

  //Plot 1d jet3 distributions in a particular eta-bin for different pt-bins on different canvases



  // //********************************************************************  Plot all Control Hists ********************************************************************************


  //Get histo files
  // TFile* f_mpf_mc = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc.root","READ");
  // TFile* f_mpf_mc_matched = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc_matched.root","READ");
  TFile* f_rel_mc = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc_matched.root","READ");
  //  TFile* f_rel_mc_matched = new TFile(CorrectionObject::_outpath+"plots/control/A_1d_mc_matched.root","READ");
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    

    // TCanvas* c1 = new TCanvas();
    // tdrCanvas(c1,"c1",h,4,10,kSquare,"MC");
    // TLegend leg1 = tdrLeg(0.17,0.6,0.85,0.81);
    // leg1.SetNColumns(2);
   
    TH1D* htemp_mpf_mc,*htemp_mpf_mc_matched,*htemp_mpf_data;
    

    for(int j=0; j<n_pt-1; j++){   ///j=0 j<pt_n-1
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_mpf_mc = "hist_mc_jet3pt_"+eta_name+"_"+pt_name;
      htemp_mpf_mc = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      htemp_mpf_mc->Print();
      int n_ev = htemp_mpf_mc->GetEntries();
      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());
      h->GetXaxis()->SetTitle("jet3 pT");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,pt_bins[j+1]);
      h->SetMinimum(0.001);
      h->SetMaximum(1.0);
      // if(j<9) htemp_mpf_mc->SetLineColor(j+1);
      // else    htemp_mpf_mc->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc->SetMarkerColor(j+1);
      // else    htemp_mpf_mc->SetMarkerColor(j+31);
      htemp_mpf_mc->SetLineColor(kRed+j-4);
      htemp_mpf_mc->SetMarkerColor(kRed+j-4);
      htemp_mpf_mc->SetLineWidth(3);
      htemp_mpf_mc->SetMarkerStyle(20);
      //      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      //      leg1.AddEntry(htemp_mpf_mc, legname, "l");

      name_mpf_mc = "hist_mc_jet3pt_"+eta_name+"_"+pt_name+"_matched";
      htemp_mpf_mc_matched = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      n_ev = htemp_mpf_mc_matched->GetEntries();
      if(htemp_mpf_mc_matched->Integral() > 0)htemp_mpf_mc_matched->Scale(1/htemp_mpf_mc_matched->Integral());
      h->GetXaxis()->SetTitle("jet3 pT");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,pt_bins[j+1]);
      h->SetMinimum(0.001);
      h->SetMaximum(0.35);
      // if(j<9) htemp_mpf_mc_matched->SetLineColor(j+1);
      // else    htemp_mpf_mc_matched->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc_matched->SetMarkerColor(j+1);
      // else    htemp_mpf_mc_matched->SetMarkerColor(j+31);
      htemp_mpf_mc_matched->SetLineColor(kRed+j-4);
      htemp_mpf_mc_matched->SetMarkerColor(kRed+j-4);
      htemp_mpf_mc_matched->SetLineWidth(3);
      htemp_mpf_mc_matched->SetLineStyle(3);
      htemp_mpf_mc_matched->SetMarkerStyle(23);
      //      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");

      name_mpf_mc = "hist_data_jet3pt_"+eta_name+"_"+pt_name;
      htemp_mpf_data = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      n_ev = htemp_mpf_data->GetEntries();
      if(htemp_mpf_data->Integral() > 0) htemp_mpf_data->Scale(1/htemp_mpf_data->Integral());
      h->GetXaxis()->SetTitle("jet3 pT");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,pt_bins[j+1]);
      h->SetMinimum(0.001);
      h->SetMaximum(0.35);
      // if(j<9) htemp_mpf_mc_matched->SetLineColor(j+1);
      // else    htemp_mpf_mc_matched->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc_matched->SetMarkerColor(j+1);
      // else    htemp_mpf_mc_matched->SetMarkerColor(j+31);
      htemp_mpf_data->SetLineColor(kRed+j-3);
      htemp_mpf_data->SetMarkerColor(kRed+j-3);
      htemp_mpf_data->SetLineWidth(3);
      htemp_mpf_data->SetLineStyle(9);
      htemp_mpf_data->SetMarkerStyle(23);
      //      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");

      TCanvas* ctmp = new TCanvas();
      tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
      TLegend leg2 = tdrLeg(0.35,0.6,0.85,0.81);
      //      leg2.SetNColumns(2);
      leg2.AddEntry(htemp_mpf_mc, legname+" MC, all", "l");
      leg2.AddEntry(htemp_mpf_mc_matched, legname+" MC, matched", "l");
      leg2.AddEntry(htemp_mpf_data, legname+" DATA", "l");
      //      gPad->SetLogx();
      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");
      if(n_ev>0) htemp_mpf_data->Draw("HIST SAME");
      leg2.Draw();
      ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_Jet3_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] 
		   +"_" +pt_name +".pdf");
    }
  }


  //Plot zPV
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    

    // TCanvas* c1 = new TCanvas();
    // tdrCanvas(c1,"c1",h,4,10,kSquare,"MC");
    // TLegend leg1 = tdrLeg(0.17,0.6,0.85,0.81);
    // leg1.SetNColumns(2);
   
    TH1D* htemp_mpf_mc,*htemp_mpf_mc_matched,*htemp_mpf_data;
    

    for(int j=0; j<n_pt-1; j++){   ///j=0 j<pt_n-1
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_mpf_mc = "hist_mc_zPV_"+eta_name+"_"+pt_name;
      htemp_mpf_mc = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      htemp_mpf_mc->Print();
      int n_ev = htemp_mpf_mc->GetEntries();
      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());
      h->GetXaxis()->SetTitle("Z_{PV}, cm");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-20,20);
      h->SetMinimum(0.001);
      h->SetMaximum(1.0);
      // if(j<9) htemp_mpf_mc->SetLineColor(j+1);
      // else    htemp_mpf_mc->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc->SetMarkerColor(j+1);
      // else    htemp_mpf_mc->SetMarkerColor(j+31);
      htemp_mpf_mc->SetLineColor(kRed);
      htemp_mpf_mc->SetMarkerColor(kRed);
      htemp_mpf_mc->SetLineWidth(3);
      htemp_mpf_mc->SetMarkerStyle(20);
      //      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      //      leg1.AddEntry(htemp_mpf_mc, legname, "l");
      htemp_mpf_mc->SetFillColorAlpha(kRed,1);

      name_mpf_mc = "hist_mc_zPV_"+eta_name+"_"+pt_name+"_matched";
      htemp_mpf_mc_matched = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      n_ev = htemp_mpf_mc_matched->GetEntries();
      if(htemp_mpf_mc_matched->Integral() > 0)htemp_mpf_mc_matched->Scale(1/htemp_mpf_mc_matched->Integral());
      h->GetXaxis()->SetTitle("Z_{PV}, cm");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-20,20);
      h->SetMinimum(0.001);
      h->SetMaximum(0.35);
      // if(j<9) htemp_mpf_mc_matched->SetLineColor(j+1);
      // else    htemp_mpf_mc_matched->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc_matched->SetMarkerColor(j+1);
      // else    htemp_mpf_mc_matched->SetMarkerColor(j+31);
      htemp_mpf_mc_matched->SetLineColor(kRed-7);
      htemp_mpf_mc_matched->SetMarkerColor(kRed-7);
      htemp_mpf_mc_matched->SetLineWidth(3);
      htemp_mpf_mc_matched->SetLineStyle(3);
      //      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");
      htemp_mpf_mc_matched->SetFillColorAlpha(kRed-7,0.5);
      name_mpf_mc = "hist_data_zPV_"+eta_name+"_"+pt_name;
      htemp_mpf_data = (TH1D*)f_rel_mc->Get(name_mpf_mc);
      n_ev = htemp_mpf_data->GetEntries();
      if(htemp_mpf_data->Integral() > 0) htemp_mpf_data->Scale(1/htemp_mpf_data->Integral());
      h->GetXaxis()->SetTitle("Z_{PV}, cm");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(-20,20);
      h->SetMinimum(0.001);
      h->SetMaximum(0.1);
      // if(j<9) htemp_mpf_mc_matched->SetLineColor(j+1);
      // else    htemp_mpf_mc_matched->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc_matched->SetMarkerColor(j+1);
      // else    htemp_mpf_mc_matched->SetMarkerColor(j+31);
      htemp_mpf_data->SetLineColor(kBlack);
      htemp_mpf_data->SetMarkerColor(kBlack);
      // htemp_mpf_data->SetLineWidth(3);
      // htemp_mpf_data->SetLineStyle(9);
      htemp_mpf_data->SetMarkerColor(kBlack);
      htemp_mpf_data->SetMarkerStyle(20);
      
      //      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");

      TCanvas* ctmp = new TCanvas();
      tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
      TLegend leg2 = tdrLeg(0.35,0.6,0.85,0.81);
      //      leg2.SetNColumns(2);
      leg2.AddEntry(htemp_mpf_mc, legname+" MC, all", "f");
      leg2.AddEntry(htemp_mpf_mc_matched, legname+" MC, matched", "f");
      leg2.AddEntry(htemp_mpf_data, legname+" DATA", "lp");
      //      gPad->SetLogx();
      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");
      if(n_ev>0) htemp_mpf_data->Draw("SAME");
      leg2.Draw();
      ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_zPV_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] 
		   +"_" +pt_name +".pdf");
    }
  }
  delete f_rel_mc;

  //Plot dRjet3
  TFile* f_rel = new TFile(CorrectionObject::_outpath+"plots/control/B_1d_mc_matched.root","READ");
  f_rel->Print();
  for(int i=0; i<n_eta-1; i++){
    TString eta_name = "eta_"+eta_range2[i]+"_"+eta_range2[i+1];
    
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045); 
    TString text = eta_range[i] + " < |#eta| < " + eta_range[i+1];

    TLatex *tex_lumi = new TLatex();
    tex_lumi->SetNDC();
    tex_lumi->SetTextSize(0.045); 
    

    // TCanvas* c1 = new TCanvas();
    // tdrCanvas(c1,"c1",h,4,10,kSquare,"MC");
    // TLegend leg1 = tdrLeg(0.17,0.6,0.85,0.81);
    // leg1.SetNColumns(2);
   
    TH1D* htemp_mpf_mc,*htemp_mpf_mc_matched;
    

    for(int j=0; j<n_pt-1; j++){   ///j=0 j<pt_n-1
      TString pt_name = "pt_"+pt_range[j]+"_"+pt_range[j+1];
      TString legname = "p_{T} #in [" + pt_range[j] + "," + pt_range[j+1] + "]";
      TString name_mpf_mc = "hist_mc_dRjet3_"+eta_name+"_"+pt_name;
      htemp_mpf_mc = (TH1D*)f_rel->Get(name_mpf_mc);
      htemp_mpf_mc->Print();
      int n_ev = htemp_mpf_mc->GetEntries();
      if(htemp_mpf_mc->Integral() > 0)htemp_mpf_mc->Scale(1/htemp_mpf_mc->Integral());
      h->GetXaxis()->SetTitle("#Delta R_{(RECO,GEN)} jet3");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,10);
      h->SetMinimum(0.001);
      h->SetMaximum(1.0);
      // if(j<9) htemp_mpf_mc->SetLineColor(j+1);
      // else    htemp_mpf_mc->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc->SetMarkerColor(j+1);
      // else    htemp_mpf_mc->SetMarkerColor(j+31);
      htemp_mpf_mc->SetLineColor(kRed);
      htemp_mpf_mc->SetMarkerColor(kRed);
      htemp_mpf_mc->SetLineWidth(3);
      htemp_mpf_mc->SetMarkerStyle(20);
      //      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      //      leg1.AddEntry(htemp_mpf_mc, legname, "l");
      htemp_mpf_mc->SetFillColorAlpha(kRed,1);
      //      f_rel->Print();
      name_mpf_mc = "hist_mc_dRjet3_"+eta_name+"_"+pt_name+"_matched";
      htemp_mpf_mc_matched = (TH1D*)f_rel->Get(name_mpf_mc);
      htemp_mpf_mc_matched->Print();
      n_ev = htemp_mpf_mc_matched->GetEntries();
      if(htemp_mpf_mc_matched->Integral() > 0)htemp_mpf_mc_matched->Scale(1/htemp_mpf_mc_matched->Integral());
      h->GetXaxis()->SetTitle("#Delta R_{(RECO,GEN)} jet3");
      h->GetYaxis()->SetTitle("Normalized entries");
      h->GetYaxis()->SetTitleOffset(1.5);
      h->GetXaxis()->SetLimits(0,10);
      h->SetMinimum(0.001);
      h->SetMaximum(1.0);
      // if(j<9) htemp_mpf_mc_matched->SetLineColor(j+1);
      // else    htemp_mpf_mc_matched->SetLineColor(j+31);
      // if(j<9) htemp_mpf_mc_matched->SetMarkerColor(j+1);
      // else    htemp_mpf_mc_matched->SetMarkerColor(j+31);
      htemp_mpf_mc_matched->SetLineColor(kBlue);
      htemp_mpf_mc_matched->SetMarkerColor(kBlue);
      htemp_mpf_mc_matched->SetLineWidth(3);
      htemp_mpf_mc_matched->SetLineStyle(3);
      //      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");
      htemp_mpf_mc_matched->SetFillColorAlpha(kBlue,0.5);

      TCanvas* ctmp = new TCanvas();
      tdrCanvas(ctmp,"ctmp",h,4,10,kSquare,CorrectionObject::_lumitag);
      TLegend leg2 = tdrLeg(0.35,0.6,0.85,0.81);
      //      leg2.SetNColumns(2);
      leg2.AddEntry(htemp_mpf_mc, legname+" MC, all", "f");
      leg2.AddEntry(htemp_mpf_mc_matched, legname+" MC, matched", "f");
      //      gPad->SetLogx();
      if(n_ev>0) htemp_mpf_mc->Draw("HIST SAME");
      if(n_ev>0) htemp_mpf_mc_matched->Draw("HIST SAME");
      leg2.Draw();
      ctmp->SaveAs(CorrectionObject::_outpath+"plots/control/Matched_dRjet3_NormDistribution_MC_" + CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1]+"_" +pt_name +".pdf");
    }
  }

  delete f_rel;
}
