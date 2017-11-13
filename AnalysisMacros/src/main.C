#include <cmath>
#include <iostream>
#include "../include/CorrectionObject.h"
#include "../include/parameters.h"
#include <TString.h>


using namespace std;

int main(){

  //************************************************************
  //
  // Declare file directories 
  // Coose objects witch will be processed (Different Run eras) 
  // 
  // Coose analysis macros 
  // Use *_eta() to take into account negative eta dependency
  // To execute all, choose FullCycle_CorrectFormulae
  //
  //************************************************************

  cout << "Hello from main(). What am I going to do?" << endl << endl;


  TString generator    = "pythia";
  bool    closure_test    = false;
  bool    trigger_fwd     = true;     //Use for Weight Calc
  bool    trigger_central = true;     //Use for Weight Calc
  TString collection    = "AK4CHS";
  //  TString collection    = "AK4PUPPI";
  //  TString collection    = "AK8CHS";
  //   TString input_path   = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup/";
  //   TString weight_path  = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup_ForWeights/";

  // TString input_path   = "/nfs/dust/cms/user/karavdia/L2Res_fromJens/AK4CHS/MC_Reweighted_chsMET_NewSF/"; //copy of Jen's files
  // TString weight_path  = "/nfs/dust/cms/user/karavdia/L2Res_fromJens/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup_ForWeights/"; //copy of Jen's files

  //  TString input_path   = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ReWeighted/";
  //  TString input_path   = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ReWeighted_MC/";
  //  TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/"; 
   TString weight_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/"; 

  //  TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS/"; //AK4CHS
  //  TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv/"; //AK4CHS
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight/"; //AK4CHS, re-miniAOD
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_SFs_Legacy07Aug_woEtaPhiCleaning/"; //AK4CHS, legacy

   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning/"; //AK4CHS, re-miniAOD
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_MatchedOnly/"; //AK4CHS, re-miniAOD

   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_SFs_Legacy07Aug/"; //AK4CHS, legacy
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__wo_JER_SFs/"; //AK4CHS
  //  TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_jet3above25GeV/"; //AK4CHS

   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_PUPPI/";  //AK4PUPPI
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_AK8CHS/";  //AK8CHS
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_AK8CHS_woJERSFs/";  //AK8CHS

   // TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_MEonlyGENparticle/";//AK4CHS, all GEN particles
   // TString input_path = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_MEonlyGENparticle_addPartonsGenRecoIDs_dRpartonConeSize/"; 
  // TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_MEonlyGENparticle_AK8CHS/";//AK4CHS, all GEN particles
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_allGENparticle_AK8CHS/";//AK8CHS, all GEN particles

   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_NotWeighted_CHS_Zpv_PUreweight_3rdJetFlavor__w_JER_woEtaPhiCleaning_MEonlyGENparticle_AK8CHS_dRpartonConeSize/";//AK8CHS, all GEN particles
      TString input_path  = "/nfs/dust/cms/user/karavdia/JERC_flavorStudies/JEC_Summer16_03Feb2017_V4_AK4CHS//";//AK4CHS, cone size (gen jet, parton) study
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JERC_flavorStudies/JEC_Summer16_03Feb2017_V4_AK4CHS_JERup/";//AK4CHS, JER up
   //   TString input_path  = "/nfs/dust/cms/user/karavdia/JERC_flavorStudies/JEC_Summer16_03Feb2017_V4_AK4CHS_JERdown/";//AK4CHS, JER down

  //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
    vector<CorrectionObject> Objects;
  
    Objects.emplace_back(CorrectionObject("BCDEFGH", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("FlateG", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("H", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
 
 
    cout << "testobject is " << Objects[0] << endl;

    //Weight Calcualtion for QCD pT binned and no trigger splitting 
    //    for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();


      if(trigger_central && !trigger_fwd){
	for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(true);  //Central Triggers
      }
      else if(!trigger_central && trigger_fwd){
	for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(false);  //FWD Triggers 
      }
      else {
	  cout<<"No Weight Calculation??"<<endl;
      }


      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();
      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae();
      // // // // // //       //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae_eta();  //extended eta range to negative Values 

      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(true);   //MPF method
      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(false);  //pT bal method
// // // //       //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(true); //extended eta range to negative Values 
// // // //       //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(false); //extended eta range to negative Values 

      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();
// // // //       // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput_eta();
// // // //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit(); //Mikkos Macro 
// // // //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_eta_0_13(); //Mikkos Macro
      
      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();
     //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].MatchingPlots();
     // ////// GenResponsePlots with params: flavor={"All","QQ","QG","GQ","GQ"}, ptbinning={"pT_ave_RECO","pT_tag_GEN","pT_probe_GEN","pT_ave_GEN"}
     //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("All","pT_ave_RECO");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("All","pT_tag_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("All","pT_probe_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("All","pT_ave_GEN");

      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QQ","pT_ave_RECO");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QQ","pT_tag_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QQ","pT_probe_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QQ","pT_ave_GEN");

      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QG","pT_ave_RECO");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QG","pT_tag_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QG","pT_probe_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("QG","pT_ave_GEN");

      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GQ","pT_ave_RECO");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GQ","pT_tag_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GQ","pT_probe_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GQ","pT_ave_GEN");

      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GG","pT_ave_RECO");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GG","pT_tag_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GG","pT_probe_GEN");
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots("GG","pT_ave_GEN");

      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots_AllPtbinningOnePlot("All"); //All, QQ, QG, GG, GQ
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots_AllPtbinningOnePlot("QQ"); //All, QQ, QG, GG, GQ
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots_AllPtbinningOnePlot("QG"); //All, QQ, QG, GG, GQ
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots_AllPtbinningOnePlot("GQ"); //All, QQ, QG, GG, GQ
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots_AllPtbinningOnePlot("GG"); //All, QQ, QG, GG, GQ

// //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae_eta(); //extended eta range to negative Values
 
// // // // //Run all macros to calculate L2Res corrections 
// // // //       // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
// // // //       // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae_eta();  //For Closure Test
 
// // // // //Macros to compare different Runs 
// // // //    // Objects[0].L2ResAllRuns();
// // // //    // Objects[0].L2ResOverlay(true);
// // // //    // Objects[0].L2ResOverlay(false);

// // // // //Compare up/nominal/down Variations of JER
// // // //    // Objects[0].L2Res_JEC();

  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
