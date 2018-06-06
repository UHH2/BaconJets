#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include "../include/useful_functions.h"

#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>

#include "../include/tdrstyle_mod15.h"






using namespace std;

void CorrectionObject::L2ResAllRuns(){
  cout << "--------------- Starting L2ResAllRuns() ---------------" << endl << endl;

  TCanvas* c1 = new TCanvas();
  TH1D *h = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h->SetMaximum(1.2); //1.2
  h->SetMinimum(0.8); //0.8
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  tdrCanvas(c1,"c1",h,4,10,true,CorrectionObject::_lumitag);

  TLegend leg1 = tdrLeg(0.20,0.19,0.65,0.4);

  TLine *line = new TLine(0,1,5.191,1);              
                  
  int n_runnr = 3;
  TString runnr_v[n_runnr]={"BCD","EFearly","FlateGH"};
  TString legend_v[n_runnr] ={"Run2016BCD","Run2016EF","Run2016GH"};

  TCanvas* c2 = new TCanvas;
  tdrCanvas(c2,"c2",h,4,10,true,"All Runs");
  TFile* f_Res_rel;
  TH1D* res_logpt_rel_kfsrfit;
  for(int j=0;j<n_runnr;j++){
    TString runnr=runnr_v[j];
    TString legend=legend_v[j];
    TString path =CorrectionObject::_input_path+"Run" + runnr + "/" + "Histo_Res_DiJet_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root";
    f_Res_rel = new TFile(path,"READ");   
    res_logpt_rel_kfsrfit = (TH1D*)f_Res_rel->Get("res_logpt_dijet");
    res_logpt_rel_kfsrfit->SetMarkerStyle(1);
    res_logpt_rel_kfsrfit->SetLineWidth(3);
    res_logpt_rel_kfsrfit->SetLineColor(1+j); 
    if(j == 2) res_logpt_rel_kfsrfit->SetLineColor(kGreen-2); 
    if(j == 4) res_logpt_rel_kfsrfit->SetLineColor(kCyan-3); 
    leg1.AddEntry(res_logpt_rel_kfsrfit,"Run2016"+runnr,"L");
    res_logpt_rel_kfsrfit->Draw("E SAME");  

  }

  line->Draw("SAME");
  leg1.Draw();

  TString save_as = CorrectionObject::_input_path+"/L2Res_Rel_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_ALL_runs_BCD_EFearly_FlateG_H.pdf";
  cout << "save-as = " << save_as << endl;
  c2->SaveAs(save_as); 





  //delete everything
  delete res_logpt_rel_kfsrfit;
  delete f_Res_rel;
  delete c2;
  delete line;
  //delete leg1;
  delete h;
  delete c1;



  TCanvas* c3 = new TCanvas();
  TH1D *h1 = new TH1D("h",";|#eta|;Relative correction",41,0,5.191);
  h1->SetMaximum(1.2); //1.2
  h1->SetMinimum(0.8); //0.8
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleSize(0.05);
  tdrCanvas(c3,"c3",h1,4,10,true,CorrectionObject::_lumitag);

  TLegend leg2 = tdrLeg(0.20,0.17,0.45,0.36);
  leg2 . SetHeader("MPF");      

  TString JetDescrib;                                                                    
  if (CorrectionObject::_collection=="AK4CHS") JetDescrib = "Anti-k_{t} R = 0.4, PF+CHS";
  if (CorrectionObject::_collection=="AK4Puppi") JetDescrib = "Anti-k_{t} R = 0.4, PF+PUPPI";
  if (CorrectionObject::_collection=="AK8CHS") JetDescrib = "Anti-k_{t} R = 0.8, PF+CHS";
  if (CorrectionObject::_collection=="AK8Puppi") JetDescrib = "Anti-k_{t} R = 0.8, PF+PUPPI"; 


  TLine *line2 = new TLine(0,1,5.191,1);              
                            
  //  TString runnr_v[5]={"BCD","EFearly","FlateG","H", "BCDEFGH"};
  TCanvas* c4 = new TCanvas;
  tdrCanvas(c4,"c4",h1,4,10,true,"Run2016  36.5 fb^{-1}");
  TFile* f_Res_mpf;
  TH1D* res_logpt_mpf_kfsrfit;
  for(int j=0;j<n_runnr;j++){
    TString runnr=runnr_v[j];
    TString legend=legend_v[j];
    TString path =CorrectionObject::_input_path+"Run" + runnr + "/" + "Histo_Res_MPF_L1_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+".root";
    f_Res_mpf = new TFile(path,"READ");   
    res_logpt_mpf_kfsrfit = (TH1D*)f_Res_mpf->Get("res_logpt_mpf");
    //   res_logpt_mpf_kfsrfit->SetMarkerStyle(20);
    res_logpt_mpf_kfsrfit->SetMarkerStyle(1);
    res_logpt_mpf_kfsrfit->SetLineWidth(3.5);
    /*
    res_logpt_mpf_kfsrfit->GetYaxis()->SetTitleSize(1.8);
    res_logpt_mpf_kfsrfit->GetXaxis()->SetTitleSize(1.8);
    res_logpt_mpf_kfsrfit->GetYaxis()->SetLabelSize(1.8);
    res_logpt_mpf_kfsrfit->GetXaxis()->SetLabelSize(1.8);
    */
    if(j == 0) res_logpt_mpf_kfsrfit->SetLineColor(kBlack);
    // if(j == 0)  res_logpt_mpf_kfsrfit->SetMarkerColor(kBlack);
    if(j == 1) res_logpt_mpf_kfsrfit->SetLineColor(kRed); 
    // if(j == 1)  res_logpt_mpf_kfsrfit->SetMarkerColor(kRed);
    if(j == 2) res_logpt_mpf_kfsrfit->SetLineColor(kBlue); 
    // if(j == 2)  res_logpt_mpf_kfsrfit->SetMarkerColor(kBlue);
    leg2.AddEntry(res_logpt_mpf_kfsrfit,legend,"LP");
    res_logpt_mpf_kfsrfit->Draw("E SAME");  

  }

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);  
  tex->DrawLatex(0.45,0.87,JetDescrib);

  line2->Draw("SAME");
  leg2.Draw("SAME");

  save_as = CorrectionObject::_input_path+"/L2Res_MPF_LOGLIN_"+CorrectionObject::_generator_tag+"_"+CorrectionObject::_jettag+"_nominal_ALL_runs_BCD_EFearly_FlateG_H.pdf";
  cout << "save-as = " << save_as << endl;
  c4->SaveAs(save_as); 





  //delete everything
  delete res_logpt_mpf_kfsrfit;
  delete f_Res_mpf;
  //delete c2;
  delete line2;
  //delete leg1;
  delete h1;
  delete c3;




}
