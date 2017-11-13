// Macro to plot responses with various pt binning on one plot
//------------------------------------------------------------

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

void CorrectionObject::GenResponsePlots_AllPtbinningOnePlot(TString flavor="All"){
  cout << "--------------- Starting GenResponsePlots_AllPtbinningOnePlot() ---------------" << endl << endl;
  gStyle->SetOptStat(0);
  TString flavorLabel = "";
  //Table with number of events in each pT- and eta-bin
  
  //Set up histos for ratios of responses
  //6 pt variables are considered for the binning
  //pT,ave,reco, pT,ave,gen, pT,probe,gen and pT, tag,gen; 
  //and pT,probe,reco and pT,tag,reco [would be quite useful for setting up the single-jet trigger L2Res analysis]

  TH1D *hmc_normjetpt[6][n_pt-1][n_eta-1];//  binning variable used for normalisation
  TH1D *hmc_probejetpt[6][n_pt-1][n_eta-1];// RECO probe jet pt devided to binning variable, e.g <pT,probe,RECO> = <pT,probe,RECO/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_probegenjetpt[6][n_pt-1][n_eta-1];// GEN probe jet pt devided to binning variable, e.g <pT,probe,GEN> = <pT,probe,GEN/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_tagjetpt[6][n_pt-1][n_eta-1];// RECO tag jet pt devided to binning variable, e.g <pT,tag,RECO> = <pT,tag,RECO/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_taggenjetpt[6][n_pt-1][n_eta-1];// GEN tag jet pt devided to binning variable, e.g <pT,tag,GEN> = <pT,tag,GEN/pT,ave,RECO> * <pT,ave,RECO>
  TH1D *hmc_B_RECO[6][n_pt-1][n_eta-1];// B_RECO
  TH1D *hmc_B_GEN[6][n_pt-1][n_eta-1];// B_GEN
  TString pt_binning_var_str[6] = {"#bar{p}^{RECO}_{T} [GeV]","p^{tag,GEN}_{T} [GeV]","p^{probe,GEN}_{T} [GeV]","#bar{p}^{GEN}_{T} [GeV]","p^{tag,RECO}_{T} [GeV]","p^{probe,RECO}_{T} [GeV]"};
  TString pt_binning_var_name[6] = {"__pT_ave_RECO__","__pT_tag_GEN__","__pT_probe_GEN__","__pT_ave_GEN__","__pT_tag_RECO__","__pT_probe_RECO__"};

  int count = 0;

  TString name99 = "hist_mc_normpt_";
  TString name100 = "hist_mc_probejet_norm_";
  TString name101 = "hist_mc_probegenjet_norm";
  TString name102 = "hist_mc_tagjet_norm";
  TString name103 = "hist_mc_taggenjet_norm";

  TString name104 = "hist_mc_B_RECO_";
  TString name105 = "hist_mc_B_GEN_";

  for(int bin=0; bin<6; bin++){
    //TString pt_binning_var_name
    for(int j=0; j<n_eta-1; j++){
      TString eta_name = "eta_"+eta_range2[j]+"_"+eta_range2[j+1];
      for(int k=0; k<n_pt-1; k++){
	TString pt_name = "pt_"+pt_range[k]+"_"+pt_range[k+1];
	TString name;
	name = name99 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_normjetpt[bin][k][j] = new TH1D(name,"",2000, 0, 2000);
	name = name100 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_probejetpt[bin][k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
	//	hmc_probejetpt[bin][k][j] = new TH1D(name,"",2000, 0.0, 2000.0);//TEST
	name = name101 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_probegenjetpt[bin][k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
	//hmc_probegenjetpt[bin][k][j] = new TH1D(name,"",2000, 0.0, 2000.0);//TEST
	name = name102 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_tagjetpt[bin][k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
	//hmc_tagjetpt[bin][k][j] = new TH1D(name,"",2000, 0.0, 2000.0);//TEST
	name = name103 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_taggenjetpt[bin][k][j] = new TH1D(name,"",nResponseBins, 0.0, 2.0);
	//hmc_taggenjetpt[bin][k][j] = new TH1D(name,"",2000, 0.0, 2000.0);//TEST

	name = name104 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_B_RECO[bin][k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);
	name = name105 + pt_binning_var_name[bin] + "_" + eta_name + "_" + pt_name;
	hmc_B_GEN[bin][k][j] = new TH1D(name,"",nResponseBins, -2.0, 2.0);

      }
    }
  }


  //Get relevant information from MC, loop over MC events 
  TTreeReader myReader_MC("AnalysisTree", CorrectionObject::_MCFile);
  TTreeReaderValue<Float_t> pt_ave_mc(myReader_MC, "pt_ave");
  TTreeReaderValue<Float_t> probejet_eta_mc(myReader_MC, "probejet_eta");
  TTreeReaderValue<Float_t> probejet_pt_mc(myReader_MC, "probejet_pt");
  TTreeReaderValue<Float_t> barreljet_pt_mc(myReader_MC, "barreljet_pt");
  TTreeReaderValue<Float_t> alpha_mc(myReader_MC, "alpha");
  TTreeReaderValue<Float_t> weight_mc(myReader_MC, "weight");
  TTreeReaderValue<Float_t> B_mc(myReader_MC, "B");
  TTreeReaderValue<Float_t> B_gen_mc(myReader_MC, "genB");

  TTreeReaderValue<Float_t> probejet_ptgen_mc(myReader_MC, "probejet_ptgen");
  TTreeReaderValue<Float_t> barreljet_ptgen_mc(myReader_MC, "barreljet_ptgen");

  TTreeReaderValue<Float_t> probejet_ptparton_mc(myReader_MC, "probejet_ptptcl");
  TTreeReaderValue<Float_t> barreljet_ptparton_mc(myReader_MC, "barreljet_ptptcl");

  TTreeReaderValue<Int_t> jet1_genID_mc(myReader_MC, "jet1_genID");
  TTreeReaderValue<Int_t> jet2_genID_mc(myReader_MC, "jet2_genID");
  TTreeReaderValue<Int_t> jet3_genID_mc(myReader_MC, "jet3_genID");

  TTreeReaderValue<Int_t> flavorProbejet_mc(myReader_MC, "flavorProbejet");
  TTreeReaderValue<Int_t> flavorTagjet_mc(myReader_MC, "flavorBarreljet");

  int icount=0;

  while (myReader_MC.Next()) {
  //  while (myReader_MC.Next() && icount<1e6) {
    double pt_binning_var = -100;
    for(int bin=0;bin<6;bin++){
      if(pt_binning_var_name[bin]=="__pT_ave_RECO__") pt_binning_var = *pt_ave_mc;//bin in pt_ave, RECO
      if(pt_binning_var_name[bin]=="__pT_tag_GEN__") pt_binning_var = *barreljet_ptgen_mc;//bin in pt_tag, GEN
      if(pt_binning_var_name[bin]=="__pT_probe_GEN__") pt_binning_var  = *probejet_ptgen_mc;//bin in pt_probe, GEN
      if(pt_binning_var_name[bin]=="__pT_ave_GEN__") pt_binning_var  =  0.5*(*barreljet_ptgen_mc+*probejet_ptgen_mc);//bin in pt_ave, GEN
      if(pt_binning_var_name[bin]=="__pT_tag_RECO__") pt_binning_var = *barreljet_pt_mc;//bin in pt_tag, RECO
      if(pt_binning_var_name[bin]=="__pT_probe_RECO__") pt_binning_var  = *probejet_pt_mc;//bin in pt_probe, RECO
      // if(fabs(*probejet_eta_mc)<0.26)      
      // 	cout<<"pt_binning_var_name[bin]: "<<pt_binning_var_name[bin]<<endl;
      if(*alpha_mc>alpha_cut) continue;
      //fill histos in bins of pt and eta
      for(int k=0; k<n_pt-1; k++){
	if(pt_binning_var<pt_bins[k] || pt_binning_var>pt_bins[k+1]) continue;
	for(int j=0; j<n_eta-1; j++){
	  if(fabs(*probejet_eta_mc)>eta_bins[j+1] || fabs(*probejet_eta_mc)<eta_bins[j]) continue;
	  else{
	    //	  bool matched;
	    if(*probejet_ptgen_mc<0 || *barreljet_ptgen_mc<0){ //not matched
	      //	    matched = false;
	    }
	    else{ //matched

	   // ///[BEGIN] Selection according to flavor of tag&probe jets---------------
	    bool flavor_sel=false;
	    if(flavor!="All"){ 
	      flavorLabel = ", "+flavor;

	      if(flavor=="QQ" && *flavorProbejet_mc>0 && *flavorProbejet_mc<6 && *flavorTagjet_mc>0 && *flavorTagjet_mc<6)
		flavor_sel=true;//QQ
	      if(flavor=="GQ" && *flavorTagjet_mc==21 && *flavorProbejet_mc>0 && *flavorProbejet_mc<6 )
		flavor_sel=true;//GQ
	      if(flavor=="GG" && *flavorTagjet_mc==21 && *flavorProbejet_mc==21)
		flavor_sel=true;//GG
	      if(flavor=="QG" && *flavorTagjet_mc>0 && *flavorTagjet_mc<6 && *flavorProbejet_mc==21)
		flavor_sel=true;//QG
	    }
	    else{
	      flavor_sel=true;
	      flavorLabel = "";
	    }
	    if(!flavor_sel) continue;
	    // ///[END] Selection according to flavor of tag&probe jets---------------

	    //    double pt_binning_var_norm = 1;//TEST without fancy normalisation
	    double pt_binning_var_norm = pt_binning_var;
	      hmc_normjetpt[bin][k][j]->Fill(pt_binning_var_norm,*weight_mc);
	      double probejetpt_norm = (*probejet_pt_mc)/pt_binning_var_norm;
	      hmc_probejetpt[bin][k][j]->Fill(probejetpt_norm,*weight_mc);
	      double probegenjetpt_norm = (*probejet_ptgen_mc)/pt_binning_var_norm;
	      hmc_probegenjetpt[bin][k][j]->Fill(probegenjetpt_norm,*weight_mc);
	      // double probepartonjetpt_norm = (*probejet_ptparton_mc)/pt_binning_var;
	      // hmc_probepartonjetpt[bin][k][j]->Fill(probepartonjetpt_norm,*weight_mc);

	      double tagjetpt_norm = (*barreljet_pt_mc)/pt_binning_var_norm;
	      hmc_tagjetpt[bin][k][j]->Fill(tagjetpt_norm,*weight_mc);
	      double taggenjetpt_norm = (*barreljet_ptgen_mc)/pt_binning_var_norm;
	      hmc_taggenjetpt[bin][k][j]->Fill(taggenjetpt_norm,*weight_mc);
	      // if(fabs(*probejet_eta_mc)<0.26){
	      // cout<<" barreljet_pt_mc = "<<*barreljet_pt_mc<<" barreljet_ptgen_mc = "<<*barreljet_ptgen_mc
	      // 	  <<" probejet_pt_mc = "<<*probejet_pt_mc<<" probejet_ptgen_mc = "<<*probejet_ptgen_mc<<endl;
	      // cout<<" tagjetpt_norm = "<<tagjetpt_norm <<" taggenjetpt_norm = "<<taggenjetpt_norm
	      // 	  <<" probejetpt_norm = "<<probejetpt_norm<<" probegenjetpt_norm  = "<<probegenjetpt_norm <<endl;
	      // cout<<""<<endl;
	      // }
	      // double tagpartonjetpt_norm = (*barreljet_ptparton_mc)/pt_binning_var;
	      // hmc_tagpartonjetpt[bin][k][j]->Fill(tagpartonjetpt_norm,*weight_mc);
	      hmc_B_RECO[bin][k][j]->Fill(*B_mc,*weight_mc);
	      hmc_B_GEN[bin][k][j]->Fill(*B_gen_mc,*weight_mc);
	    }
	  }
	}
	icount++;
      }
    }
  }

 
 
 

  // Dump 1-d distributions of A and B in bins of pT, eta

  TFile* test_out_mc_B = new TFile(CorrectionObject::_outpath+"plots/control/GenResponse_AllPtbinOnePlot_1d_mc_matched_"+flavor+".root","RECREATE");
  for(int bin=0;bin<6;bin++){
    for(int j=0; j<n_eta-1; j++){
      for(int k=0; k<n_pt-1; k++){     ///k=0 n_pt-1 
	hmc_normjetpt[bin][k][j]->Write();
	hmc_probejetpt[bin][k][j]->Write();
	hmc_probegenjetpt[bin][k][j]->Write();
	hmc_tagjetpt[bin][k][j]->Write();
	hmc_taggenjetpt[bin][k][j]->Write();
	hmc_B_RECO[bin][k][j]->Write();
	hmc_B_GEN[bin][k][j]->Write();
      }
    }
  }
  test_out_mc_B->Close();
  delete test_out_mc_B;

  double val_probejet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double err_probejet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double val_probegenjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double err_probegenjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  // double val_probepartonjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  // double err_probepartonjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double val_tagjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double err_tagjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double val_taggenjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta
  double err_taggenjet_pt[6][n_eta-1][n_pt-1]; //value at pt,eta

  double val_probeRECO_probeGEN[6][n_eta-1][n_pt-1];
  double err_probeRECO_probeGEN[6][n_eta-1][n_pt-1];
  double val_tagRECO_tagGEN[6][n_eta-1][n_pt-1];
  double err_tagRECO_tagGEN[6][n_eta-1][n_pt-1];

  double val_probeRECO_tagRECO[6][n_eta-1][n_pt-1];
  double err_probeRECO_tagRECO[6][n_eta-1][n_pt-1];

  double val_probeGEN_tagGEN[6][n_eta-1][n_pt-1];
  double err_probeGEN_tagGEN[6][n_eta-1][n_pt-1];

  double val_MPF_RECO[6][n_eta-1][n_pt-1];
  double err_MPF_RECO[6][n_eta-1][n_pt-1];

  double val_MPF_GEN[6][n_eta-1][n_pt-1];
  double err_MPF_GEN[6][n_eta-1][n_pt-1];

 
  for(int bin=0; bin<6; bin++){
    for(int i=0; i<n_eta-1; i++){
      for(int j=0; j<n_pt-1; j++){
	pair <double,double> normpt_mc = GetValueAndError(hmc_normjetpt[bin][j][i]); //<pt_bin> value used for normalisation

	if(i==0){
	  cout<<pt_binning_var_name[bin]<<endl;
	  cout<<" normpt_mc = "<<normpt_mc.first<<" +/- "<<normpt_mc.second<<endl;
	}
	pair <double,double> probejetpt_mc = GetValueAndError(hmc_probejetpt[bin][j][i]);
	val_probejet_pt[bin][i][j] = probejetpt_mc.first*normpt_mc.first;
	err_probejet_pt[bin][i][j] =  ErrorPropagation_AB(probejetpt_mc,normpt_mc);
	pair <double,double> probegenjetpt_mc = GetValueAndError(hmc_probegenjetpt[bin][j][i]);
	if(i==0)	
	  cout<<" probegenjetpt_mc = "<<probegenjetpt_mc.first<<" +/- "<<probegenjetpt_mc.second<<endl;
	val_probegenjet_pt[bin][i][j] = probegenjetpt_mc.first*normpt_mc.first;
	err_probegenjet_pt[bin][i][j] = ErrorPropagation_AB(probegenjetpt_mc,normpt_mc);
	// pair <double,double> probepartonjetpt_mc = GetValueAndError(hmc_probepartonjetpt[bin][j][i]);
	// val_probepartonjet_pt[bin][i][j] = probepartonjetpt_mc.first*normpt_mc.first;
	// err_probepartonjet_pt[bin][i][j] = ErrorPropagation_AB(probepartonjetpt_mc,normpt_mc);

	pair <double,double> tagjetpt_mc = GetValueAndError(hmc_tagjetpt[bin][j][i]);
	val_tagjet_pt[bin][i][j] = tagjetpt_mc.first*normpt_mc.first;
	err_tagjet_pt[bin][i][j] = ErrorPropagation_AB(tagjetpt_mc,normpt_mc);
	pair <double,double> taggenjetpt_mc = GetValueAndError(hmc_taggenjetpt[bin][j][i]);
	if(i==0)	
	  cout<<" taggenjetpt_mc = "<<taggenjetpt_mc.first<<" +/- "<<taggenjetpt_mc.second<<endl;

	val_taggenjet_pt[bin][i][j] = taggenjetpt_mc.first*normpt_mc.first;
	err_taggenjet_pt[bin][i][j] = ErrorPropagation_AB(taggenjetpt_mc,normpt_mc);
	if(i==0)	
	  cout<<" val_probegenjet_pt[bin][i][j] = "<<val_probegenjet_pt[bin][i][j]<<" val_taggenjet_pt[bin][i][j] = "<<val_taggenjet_pt[bin][i][j]<<endl;
	if(i==0)
	  cout<<" "<<endl;
	// pair <double,double> tagpartonjetpt_mc = GetValueAndError(hmc_tagpartonjetpt[bin][j][i]);
	// val_tagpartonjet_pt[bin][i][j] = tagpartonjetpt_mc.first*normpt_mc.first;
	// err_tagpartonjet_pt[bin][i][j] = ErrorPropagation_AB(tagpartonjetpt_mc,normpt_mc);


	pair <double,double> B_RECO_mc = GetValueAndError(hmc_B_RECO[bin][j][i]);
	val_MPF_RECO[bin][i][j] = (1+B_RECO_mc.first)/(1-B_RECO_mc.first);
	err_MPF_RECO[bin][i][j] = 2/(pow((1-B_RECO_mc.first),2)) * B_RECO_mc.second;
	pair <double,double> B_GEN_mc = GetValueAndError(hmc_B_GEN[bin][j][i]);
	val_MPF_GEN[bin][i][j] = (1+B_GEN_mc.first)/(1-B_GEN_mc.first);
	err_MPF_GEN[bin][i][j] = 2/(pow((1-B_GEN_mc.first),2)) * B_GEN_mc.second;


	if(val_probegenjet_pt[bin][i][j]>0){
	  val_probeRECO_probeGEN[bin][i][j] = val_probejet_pt[bin][i][j]/val_probegenjet_pt[bin][i][j];
	  pair<double,double> tmp1; tmp1.first = val_probejet_pt[bin][i][j]; tmp1.second = err_probejet_pt[bin][i][j];
	  pair<double,double> tmp2; tmp2.first = val_probegenjet_pt[bin][i][j]; tmp2.second = err_probegenjet_pt[bin][i][j];
	  err_probeRECO_probeGEN[bin][i][j] = ErrorPropagation_AoverB(tmp1,tmp2);
	}
	else{
	  val_probeRECO_probeGEN[bin][i][j] = 0;  err_probeRECO_probeGEN[bin][i][j] =0;
	}
	if(val_taggenjet_pt[bin][i][j]>0){
	  val_tagRECO_tagGEN[bin][i][j] = val_tagjet_pt[bin][i][j]/val_taggenjet_pt[bin][i][j];
	  pair<double,double> tmp1; tmp1.first = val_tagjet_pt[bin][i][j]; tmp1.second = err_tagjet_pt[bin][i][j];
	  pair<double,double> tmp2; tmp2.first = val_taggenjet_pt[bin][i][j]; tmp2.second = err_taggenjet_pt[bin][i][j];
	  err_tagRECO_tagGEN[bin][i][j] = ErrorPropagation_AoverB(tmp1,tmp2);

	}
	else{
	  val_tagRECO_tagGEN[bin][i][j] = 0;  err_tagRECO_tagGEN[bin][i][j] =0;
	}
	if(val_tagjet_pt[bin][i][j]>0){
	  val_probeRECO_tagRECO[bin][i][j] = val_probejet_pt[bin][i][j]/val_tagjet_pt[bin][i][j];
	  pair<double,double> tmp1; tmp1.first = val_probejet_pt[bin][i][j]; tmp1.second = err_probejet_pt[bin][i][j];
	  pair<double,double> tmp2; tmp2.first = val_tagjet_pt[bin][i][j]; tmp2.second = err_tagjet_pt[bin][i][j];
	  err_probeRECO_tagRECO[bin][i][j] = ErrorPropagation_AoverB(tmp1,tmp2);

	}
	else{
	  val_probeRECO_tagRECO[bin][i][j] = 0;  err_probeRECO_tagRECO[bin][i][j] =0;
	}
	if(val_taggenjet_pt[bin][i][j]>0){
	  val_probeGEN_tagGEN[bin][i][j] = val_probegenjet_pt[bin][i][j]/val_taggenjet_pt[bin][i][j];
	  pair<double,double> tmp1; tmp1.first = val_probegenjet_pt[bin][i][j]; tmp1.second = err_probegenjet_pt[bin][i][j];
	  pair<double,double> tmp2; tmp2.first = val_taggenjet_pt[bin][i][j]; tmp2.second = err_taggenjet_pt[bin][i][j];
	  err_probeGEN_tagGEN[bin][i][j] = ErrorPropagation_AoverB(tmp1,tmp2);
	  if(i==0) 
	    cout<<"val_probeGEN_tagGEN[bin][i][j] = "<<val_probeGEN_tagGEN[bin][i][j]<<endl;
	}
	else{
	  val_probeGEN_tagGEN[bin][i][j] = 0;  err_probeGEN_tagGEN[bin][i][j] =0;
	}
      }
    }
  }


  //dummy for tdrCanvas
  TH1D *h = new TH1D("h",";dummy;",41,0,5.191);
  h->SetMaximum(1.2);
  h->SetMinimum(0.8);

  TH1D *hEF = new TH1D("hEF",";dummy;",1000,0,5.191);

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

    TString alVal;
    alVal.Form("%0.2f\n",alpha_cut);
    TString altitle = "{#alpha<"+alVal+"}";
    TLegend *leg_rel;
    //    leg_rel = new TLegend(0.45,0.15,0.91,0.49,"","brNDC");//x+0.1
    leg_rel = new TLegend(0.36,0.65,0.91,0.90,"","brNDC");//x+0.1
    leg_rel->SetNColumns(2);
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.030);
    leg_rel->SetFillColor(10);
    leg_rel->SetFillStyle(0);
    leg_rel->SetLineColor(1);
    leg_rel->SetTextFont(42);
    leg_rel->SetHeader("(<p^{probe,GEN}_{T}>/<p^{tag,GEN}_{T}>)_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+" "+flavorLabel); 
    TCanvas* c_rel = new TCanvas();
    tdrCanvas(c_rel,"c_rel",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("p_{T} [GeV] (see label)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);

    TGraphErrors *graph_probeGEN_tagGEN[6];
    for(int bin=0; bin<6; bin++){
      graph_probeGEN_tagGEN[bin]   = new TGraphErrors(n_pt-1, xbin_tgraph, val_probeGEN_tagGEN[bin][i], zero, err_probeGEN_tagGEN[bin][i]);
      graph_probeGEN_tagGEN[bin]   = (TGraphErrors*)CleanEmptyPoints(graph_probeGEN_tagGEN[bin]);
      graph_probeGEN_tagGEN[bin]->SetTitle("");
      graph_probeGEN_tagGEN[bin]->SetMarkerColor(kGreen-1+bin);
      graph_probeGEN_tagGEN[bin]->SetMarkerStyle(20+bin);
      graph_probeGEN_tagGEN[bin]->SetMarkerSize(1.8);
      graph_probeGEN_tagGEN[bin]->SetLineColor(kGreen-1+bin);
      leg_rel->AddEntry(graph_probeGEN_tagGEN[bin],pt_binning_var_str[bin],"P");
      graph_probeGEN_tagGEN[bin]->Draw("P SAME");
    }
    gPad->SetLogx();
    leg_rel->Draw();
    c_rel->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_probeGEN_tagGEN__allPtbinsAtOnePlot_"+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +"_"+ flavor + ".pdf");

    TLegend *leg_reco;
    //    leg_reco = new TLegend(0.45,0.15,0.91,0.49,"","brNDC");//x+0.1
    leg_reco = new TLegend(0.36,0.65,0.91,0.90,"","brNDC");//x+0.1
    leg_reco->SetNColumns(2);
    leg_reco->SetBorderSize(0);
    leg_reco->SetTextSize(0.030);
    leg_reco->SetFillColor(10);
    leg_reco->SetFillStyle(0);
    leg_reco->SetLineColor(1);
    leg_reco->SetTextFont(42);
    leg_reco->SetHeader("(<p^{probe,RECO}_{T}>/<p^{tag,RECO}_{T}>)_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+" "+flavorLabel); 
    TCanvas* c_reco = new TCanvas();
    tdrCanvas(c_reco,"c_reco",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("p_{T} [GeV] (see label)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);

    TGraphErrors *graph_probeRECO_tagRECO[6];
    for(int bin=0; bin<6; bin++){
      graph_probeRECO_tagRECO[bin]   = new TGraphErrors(n_pt-1, xbin_tgraph, val_probeRECO_tagRECO[bin][i], zero, err_probeRECO_tagRECO[bin][i]);
      graph_probeRECO_tagRECO[bin]   = (TGraphErrors*)CleanEmptyPoints(graph_probeRECO_tagRECO[bin]);
      graph_probeRECO_tagRECO[bin]->SetTitle("");
      graph_probeRECO_tagRECO[bin]->SetMarkerColor(kRed-1+bin);
      graph_probeRECO_tagRECO[bin]->SetMarkerStyle(20+bin);
      graph_probeRECO_tagRECO[bin]->SetMarkerSize(1.8);
      graph_probeRECO_tagRECO[bin]->SetLineColor(kRed-1+bin);
      leg_reco->AddEntry(graph_probeRECO_tagRECO[bin],pt_binning_var_str[bin],"P");
      graph_probeRECO_tagRECO[bin]->Draw("P SAME");
    }
    gPad->SetLogx();
    leg_reco->Draw();

    c_reco->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_probeRECO_tagRECO__allPtbinsAtOnePlot_"+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +"_"+ flavor + ".pdf");


    TLegend *leg_mpf;
    //    leg_mpf = new TLegend(0.45,0.15,0.91,0.49,"","brNDC");//x+0.1
    leg_mpf = new TLegend(0.36,0.65,0.91,0.90,"","brNDC");//x+0.1
    leg_mpf->SetNColumns(2);
    leg_mpf->SetBorderSize(0);
    leg_mpf->SetTextSize(0.030);
    leg_mpf->SetFillColor(10);
    leg_mpf->SetFillStyle(0);
    leg_mpf->SetLineColor(1);
    leg_mpf->SetTextFont(42);
    leg_mpf->SetHeader("((1+<B_{RECO}>)/(1-<B_{RECO}>))_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+" "+flavorLabel); 
    TCanvas* c_MPF_reco = new TCanvas();
    tdrCanvas(c_MPF_reco,"c_MPF_reco",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("p_{T} [GeV] (see label)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);

    TGraphErrors *graph_MPF_RECO[6];
    for(int bin=0; bin<6; bin++){
      graph_MPF_RECO[bin]   = new TGraphErrors(n_pt-1, xbin_tgraph, val_MPF_RECO[bin][i], zero, err_MPF_RECO[bin][i]);
      graph_MPF_RECO[bin]   = (TGraphErrors*)CleanEmptyPoints(graph_MPF_RECO[bin]);
      graph_MPF_RECO[bin]->SetTitle("");
      graph_MPF_RECO[bin]->SetMarkerColor(kPink-1+bin);
      graph_MPF_RECO[bin]->SetMarkerStyle(20+bin);
      graph_MPF_RECO[bin]->SetMarkerSize(1.8);
      graph_MPF_RECO[bin]->SetLineColor(kPink-1+bin);
      leg_mpf->AddEntry(graph_MPF_RECO[bin],pt_binning_var_str[bin],"P");
      graph_MPF_RECO[bin]->Draw("P SAME");
    }
    gPad->SetLogx();
    leg_mpf->Draw();

    c_MPF_reco->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_MPF_RECO__allPtbinsAtOnePlot_"+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +"_"+ flavor + ".pdf");


    TLegend *leg_mpf_gen;
    //    leg_mpf_gen = new TLegend(0.45,0.15,0.91,0.49,"","brNDC");//x+0.1
    leg_mpf_gen = new TLegend(0.36,0.65,0.91,0.90,"","brNDC");//x+0.1
    leg_mpf_gen->SetNColumns(2);
    leg_mpf_gen->SetBorderSize(0);
    leg_mpf_gen->SetTextSize(0.030);
    leg_mpf_gen->SetFillColor(10);
    leg_mpf_gen->SetFillStyle(0);
    leg_mpf_gen->SetLineColor(1);
    leg_mpf_gen->SetTextFont(42);
    leg_mpf_gen->SetHeader("((1+<B_{GEN}>)/(1-<B_{GEN}>))_"+altitle+", "+eta_range[i]+"#leq|#eta|<"+eta_range[i+1]+" "+flavorLabel); 
    TCanvas* c_MPF_gen = new TCanvas();
    tdrCanvas(c_MPF_gen,"c_MPF_gen",h,4,10,true,CorrectionObject::_lumitag);
    h->GetXaxis()->SetTitle("p_{T} [GeV] (see label)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLimits(30,pt_bins[n_pt-1]+100);
    h->GetYaxis()->SetRangeUser(0.70,1.30);

    TGraphErrors *graph_MPF_GEN[6];
    for(int bin=0; bin<6; bin++){
      graph_MPF_GEN[bin]   = new TGraphErrors(n_pt-1, xbin_tgraph, val_MPF_GEN[bin][i], zero, err_MPF_GEN[bin][i]);
      graph_MPF_GEN[bin]   = (TGraphErrors*)CleanEmptyPoints(graph_MPF_GEN[bin]);
      graph_MPF_GEN[bin]->SetTitle("");
      graph_MPF_GEN[bin]->SetMarkerColor(kBlue-1+bin);
      graph_MPF_GEN[bin]->SetMarkerStyle(20+bin);
      graph_MPF_GEN[bin]->SetMarkerSize(1.8);
      graph_MPF_GEN[bin]->SetLineColor(kBlue-1+bin);
      leg_mpf_gen->AddEntry(graph_MPF_GEN[bin],pt_binning_var_str[bin],"P");
      graph_MPF_GEN[bin]->Draw("P SAME");
    }
    gPad->SetLogx();
    leg_mpf_gen->Draw();

    c_MPF_gen->SaveAs(CorrectionObject::_outpath+"plots/control/GenResponse_MPF_GEN__allPtbinsAtOnePlot_"+ CorrectionObject::_generator_tag + "_eta_" + eta_range2[i] + "_" + eta_range2[i+1] +"_"+ flavor + ".pdf");



  }
}
