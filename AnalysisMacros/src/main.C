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

  cout << "Hello from main(). What am I going to do?\nWill it involve dead regions and high jetiness?\n13371n6 1n 7h3 57uff" << endl << endl;
  


  TString generator    = "pythia";
  bool    closure_test    = false;
  bool    trigger_fwd     = false;     //Use for Weight Calc
  bool    trigger_central = true;     //Use for Weight Calc
  TString collection    = "AK4CHS";

//   TString input_path   = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup/";
  //   TString weight_path  = "/nfs/dust/cms/user/multh/JEC/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup_ForWeights/";

  // TString input_path   = "/nfs/dust/cms/user/karavdia/L2Res_fromJens/AK4CHS/MC_Reweighted_chsMET_NewSF/"; //copy of Jen's files
  // TString weight_path  = "/nfs/dust/cms/user/karavdia/L2Res_fromJens/AK4CHS/MC_Reweighted_chsMET_NewTriggerSetup_ForWeights/"; //copy of Jen's files

  //  TString input_path   = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ReWeighted/";
  //  TString input_path   = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ReWeighted_MC/";
  // TString input_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/";
  TString input_path  = "/nfs/dust/cms/user/garbersc/forBaconJets/2016ReReco/Residuals/Summer16_03Feb2017_V3/AK4CHS/MC_Reweighted_chsMET_NewSF/DATA_RunB_promptReco_2017.root"; 
  TString weight_path  = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/"; 
  TString input_path_MC = "/nfs/dust/cms/user/karavdia/JEC_Summer16_V8_ForWeights/uhh2.AnalysisModuleRunner.MC.QCDPt50toInf_pythia8_AK4CHS.root";

  //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
    vector<CorrectionObject> Objects;
  
    // Objects.emplace_back(CorrectionObject("BCDEFGH", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //      Objects.emplace_back(CorrectionObject("FlateG", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    //    Objects.emplace_back(CorrectionObject("H", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
    Objects.emplace_back(CorrectionObject("B", generator,collection, input_path, input_path_MC, weight_path, closure_test, trigger_fwd, trigger_central));
 
    cout << "testobject is " << Objects[0] << endl;

    //Weight Calcualtion for QCD pT binned and no trigger splitting 
    //    for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();


      // if(trigger_central && !trigger_fwd){
      // 	for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(true);  //Central Triggers
      // }
      // else if(!trigger_central && trigger_fwd){
      // 	for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights_TriggerThresholds(false);  //FWD Triggers 
      // }
      // else {
      // 	  cout<<"No Weight Calculation??"<<endl;
      // }


//       for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();
//       for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae();
// // // //       //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae_eta();  //extended eta range to negative Values 

// std::cout<<"\nStarting Pt_Extrapolation_Alternative_CorrectFormulae(true)\n"<<std::endl;
    //   for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(true);   //MPF method
    // std::cout<<"\nStarting Pt_Extrapolation_Alternative_CorrectFormulae(false)\n"<<std::endl;
    //   for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(false);  //pT bal method

    // // // //       //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(true); //extended eta range to negative Values 
// // // //       //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(false); //extended eta range to negative Values
      
  // std::cout<<"\nStarting L2ResOutput()\n"<<std::endl;
      // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();

// //       // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput_eta();
// //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit(); //Mikkos Macro 
// //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_eta_0_13(); //Mikkos Macro

  std::cout<<"\nStarting FinalControlPlots_CorrectFormulae()\n"<<std::endl;
      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();
      
//      for(unsigned int i=0; i<Objects.size(); i++) Objects[i].MatchingPlots();

// //    // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae_eta(); //extended eta range to negative Values

     
// // // //Run all macros to calculate L2Res corrections 
   // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
// // //       // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae_eta();  //For Closure Test
 
// // // //Macros to compare different Runs 
// // //    // Objects[0].L2ResAllRuns();
// // //    // Objects[0].L2ResOverlay(true);
// // //    // Objects[0].L2ResOverlay(false);

// // // //Compare up/nominal/down Variations of JER
// // //    // Objects[0].L2Res_JEC();

  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  return 0;
}
