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
  bool    closure_test    = true;
  bool    trigger_fwd     = true;     //Use for Weight Calc
  bool    trigger_central = true;     //Use for Weight Calc
  TString collection    = "AK4CHS";

/*  
  TString input_path    ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth/";
  TString weight_path   ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth/";
*/

/*  
  TString input_path    ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth_JERUp/";
  TString weight_path   ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth_JERUp/";
*/
 
  TString input_path    ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth_JERDown/";
  TString weight_path   ="/nfs/dust/cms/user/multh/JEC/2016Legacy/ClosureTest/Summer16_07Aug2017_V18/AK4CHS/MC_NoReweighted_CHS_newMCTruth_JERDown/";
 

  //eine Klasse: enthaelt Info ueber runnr, Generator, collection, Strings zu MC/DATA-files, memberfunctions: controlPlots, kFSR etc.
  vector<CorrectionObject> Objects;
  
  Objects.emplace_back(CorrectionObject("BCDEFGH", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  Objects.emplace_back(CorrectionObject("BCD", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  Objects.emplace_back(CorrectionObject("EFearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  Objects.emplace_back(CorrectionObject("FlateGH", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  
  // Objects.emplace_back(CorrectionObject("B", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("BCDEFearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
 
  //Objects.emplace_back(CorrectionObject("B", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
   // Objects.emplace_back(CorrectionObject("C", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("D", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("E", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("Fearly", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("Flate", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  // Objects.emplace_back(CorrectionObject("G", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central)); 
  // Objects.emplace_back(CorrectionObject("H", generator,collection, input_path, weight_path, closure_test, trigger_fwd, trigger_central));
  
  
  cout << "testobject is " << Objects[0] << endl;
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Loop_MC_Data();

  //Weight Calcualtion for QCD pT binned and no trigger splitting 
  //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CalculateMCWeights();
  
// for(unsigned int i=0; i<Objects.size(); i++) Objects[i].ControlPlots();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectExtrapolation();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae("0","7");
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae("0","10");
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae("0","14");
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae("0","16");
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].kFSR_CorrectFormulae_eta();  //extended eta range to negative Values 
  
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(true);   //MPF method
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae(false);  //pT bal method
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(true); //extended eta range to negative Values 
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Pt_Extrapolation_Alternative_CorrectFormulae_eta(false); //extended eta range to negative Values 
  
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput();
 //for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput_CorrectExtrapolation();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].L2ResOutput_eta(); 
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit(); //Mikkos Macro 
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_eta_0_13(); //Mikkos Macro
 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].InputForGlobalFit_MC(); //Finer pT binning for MC

 // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].MatchingPlots();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].GenResponsePlots();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FinalControlPlots_CorrectFormulae_eta(); //extended eta range to negative Values
  
  //Run all macros to calculate L2Res corrections 
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae();
  // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FullCycle_CorrectFormulae_eta();  //For Closure Test
  
  //Monitoring Macro
  // for(unsigned int i=0; i<1; i++) Objects[i].Monitoring();
  
  //Derive tirgger thresholds
  //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Derive_Thresholds_alternativeWay();
  //  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].Derive_Thresholds_SingleMuonCrossCheck("HLT_Mu27");
  
  // // // // //Macros to compare different Runs 
 //  Objects[0].L2ResAllRuns();
//  Objects[0].L2ResOverlay(true);
  Objects[0].L2ResOverlay_JEC();
  
  // // // // //Compare up/nominal/down Variations of JER
  // Objects[0].L2Res_JEC();    Objects.size()

  // // for(unsigned int i=0; i<Objects.size(); i++) Objects[i].FlavorCorrection_TTree();


  
  cout << endl << "Closing MC and DATA files." << endl;
  for(unsigned int i=0; i<Objects.size(); i++) Objects[i].CloseFiles();
  cout << "Going to return 0 now, cya." << endl << endl;
  
  return 0;
}
