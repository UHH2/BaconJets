#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/EventHelper.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/L1Jet.h"
#include "UHH2/core/include/Utils.h"
//#include "../include/JECAnalysisHists.h"
#include "../include/JECCrossCheckHists.h"
//#include "../include/JECRunnumberHists.h"
//#include "../include/JECAnalysisFinalStateHadronsHists.h"
//#include "../include/JECAnalysisRecoGenMatchedHists.h" //uses PF and GEN candidates
#include "../include/JECAnalysisRecoGenMatchedHistsFractions.h" //uses PF and GEN fractions only
#include "../include/JECAnalysisPUjetsHists.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/LumiSelection.h> //includes also LuminosityHists.h
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/CleaningModules.h"
#include <UHH2/common/include/NSelections.h> 
#include <UHH2/common/include/MuonIds.h> 
#include <UHH2/common/include/ElectronIds.h>
#include "UHH2/common/include/PrintingModules.h"

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/constants.h"

//#include "UHH2/BaconJets/include/LumiHists.h"


#include "TClonesArray.h"
#include "TString.h"
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2Poly.h"


using namespace std;
using namespace uhh2;

class AnalysisModule_JetGenRecoMatch_RelValsTest: public uhh2::AnalysisModule {

  public:
  explicit AnalysisModule_JetGenRecoMatch_RelValsTest(uhh2::Context&);
    virtual bool process(uhh2::Event&) override;
    ~AnalysisModule_JetGenRecoMatch_RelValsTest();
    void init_JEC(uhh2::Context& ctx);
    void declare_output(uhh2::Context& ctx);
    void init_hists(uhh2::Context& ctx);

protected:

    //useful booleans
  bool isMC, ispuppi;
  bool is2016v2, is2016v3, is2017, is2018;
  // correctors
  std::unique_ptr<JetCorrector> jet_corrector_A,jet_corrector_B, jet_corrector_C, jet_corrector_D, jet_corrector_E, jet_corrector_F, 
    jet_corrector_G,jet_corrector_H, jet_corrector_MC;
  std::unique_ptr<GenericJetResolutionSmearer> jetER_smearer; 
    

  // cleaners
   std::unique_ptr<JetCleaner> jetID;
   std::unique_ptr<JetCleaner> jetcleaner;
   std::unique_ptr<GenJetCleaner> genjetcleaner;
   std::unique_ptr<MuonCleaner>     muoSR_cleaner;   
   std::unique_ptr<ElectronCleaner> eleSR_cleaner;    

    // selections
    std::unique_ptr<uhh2::Selection> lumi_sel;
    std::unique_ptr<uhh2::AndSelection> metfilters_sel;

    std::unique_ptr<uhh2::Selection> triggerSiMu_sel;

    std::unique_ptr<uhh2::Selection> trigger40_sel;
    std::unique_ptr<uhh2::Selection> trigger60_sel;
    std::unique_ptr<uhh2::Selection> trigger80_sel;
    std::unique_ptr<uhh2::Selection> trigger140_sel;
    std::unique_ptr<uhh2::Selection> trigger200_sel;
    std::unique_ptr<uhh2::Selection> trigger260_sel;
    std::unique_ptr<uhh2::Selection> trigger320_sel;
    std::unique_ptr<uhh2::Selection> trigger400_sel;
    std::unique_ptr<uhh2::Selection> trigger500_sel;

    std::unique_ptr<uhh2::Selection> trigger60_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger80_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger100_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger160_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger220_HFJEC_sel;
    std::unique_ptr<uhh2::Selection> trigger300_HFJEC_sel;  
  
  //    std::unique_ptr<LumiHists> h_monitoring_final;
  
    //// Data/MC scale factors
    std::unique_ptr<uhh2::AnalysisModule> pileupSF;
  unique_ptr<AnalysisModule>  Jet_printer;
  unique_ptr<AnalysisModule>  GenJet_printer;
  unique_ptr<AnalysisModule> GenParticles_printer;
  unique_ptr<AnalysisModule> LumiWeight_module;

  //    Event::Handle<float> tt_jet1_ptGen;  Event::Handle<float> tt_jet2_ptGen;  Event::Handle<float> tt_jet3_ptGen;
  Event::Handle<float> tt_gen_pthat; Event::Handle<float> tt_gen_weight;  Event::Handle<float> tt_gen_PUpthat;
  Event::Handle<int> tt_nPU;
  Event::Handle<float> tt_ev_weight;
  Event::Handle<int> tt_jet_n;//number of jets
  Event::Handle<int> tt_matchedjet_n;//number of RECO jets matched to GEN
  Event::Handle<int> tt_genjet_n;//number of gen jets
  Event::Handle<int> tt_genptcl_n;//number of partons with pt>0
  Event::Handle<int> tt_nvertices;
  Event::Handle<float> tt_rho;//event energy density
  Event::Handle<int> tt_nGoodvertices;

  Event::Handle<int> tt_run; 
  Event::Handle<int> tt_evID;
  Event::Handle<int> tt_lumiSec;

  Event::Handle<float> tt_jet1_pt;  Event::Handle<float> tt_jet2_pt;   Event::Handle<float> tt_jet3_pt;
  Event::Handle<float> tt_jet1_eta;  Event::Handle<float> tt_jet2_eta;   Event::Handle<float> tt_jet3_eta;
  Event::Handle<float> tt_jet1_phi;  Event::Handle<float> tt_jet2_phi;   Event::Handle<float> tt_jet3_phi;

  Event::Handle<float> tt_jet1_chEM;  Event::Handle<float> tt_jet2_chEM;   Event::Handle<float> tt_jet3_chEM;
  Event::Handle<float> tt_jet1_neutEM;  Event::Handle<float> tt_jet2_neutEM;   Event::Handle<float> tt_jet3_neutEM;
  Event::Handle<float> tt_jet1_chHAD;  Event::Handle<float> tt_jet2_chHAD;   Event::Handle<float> tt_jet3_chHAD;
  Event::Handle<float> tt_jet1_neutHAD;  Event::Handle<float> tt_jet2_neutHAD;   Event::Handle<float> tt_jet3_neutHAD;

  //  std::unique_ptr<JECAnalysisRecoGenMatchedHists> h_matched_all; //uses PF and GEN candidates
  //  std::unique_ptr<JECAnalysisRecoGenMatchedHists> h_matched_pt[n_pt-1];//uses PF and GEN candidates


  std::unique_ptr<JECAnalysisRecoGenMatchedHistsFractions> h_matched_all;//uses PF and GEN fractions only
  std::unique_ptr<JECAnalysisRecoGenMatchedHistsFractions> h_matched_pt[n_pt-1];//uses PF and GEN fractions only


  std::unique_ptr<JECAnalysisPUjetsHists> h_goodRECO_all;
  std::unique_ptr<JECAnalysisPUjetsHists> h_badRECO_all;
  std::unique_ptr<JECAnalysisPUjetsHists> h_goodRECO_pt[n_pt-1];
  std::unique_ptr<JECAnalysisPUjetsHists> h_badRECO_pt[n_pt-1];

  // std::unique_ptr<JECAnalysisHists> h_nocuts, h_sel, h_dijet, h_match, h_final;
  // std::unique_ptr<JECAnalysisHists> h_trg40, h_trg60, h_trg80, h_trg140, h_trg200,h_trg260,h_trg320,h_trg400,h_trg500;
  // std::unique_ptr<JECAnalysisHists> h_trgHF60, h_trgHF80,h_trgHF100, h_trgHF160,h_trgHF220, h_trgHF300;   

  std::unique_ptr<JECCrossCheckHists> h_input,h_beforeCleaner,h_afterCleaner,h_afternVts;
  uhh2bacon::Selection sel;
  
  bool debug;
  bool ispythia8;
    bool no_genp;
  bool split_JEC_DATA, split_JEC_MC, ClosureTest, apply_weights, apply_lumiweights, apply_unflattening, apply_smear, apply_METoverPt_cut, apply_EtaPhi_cut, trigger_central, trigger_fwd, ts, onlyBtB, apply_L1seed_from_bx1_filter;
  //    double lumiweight;
    string jetLabel;
   TString dataset_version;
  // TString JEC_Version;
    JetId Jet_PFID;
    int n_evt;
    std::unique_ptr<TFile> f_weights;
    
    // std::map<run_lumi, double> rl2lumi;
    // std::map<run_lumi, double> rl2intlumi;
    // TBranch * brun ;
    // TBranch * blumiblock;
    // TBranch * bilumi;
    // double integrated_lumi;
    // vector<run_lumi> upper_binborders_runnrs;
    // vector<double> lumi_in_bins;

    double L1METptThresh;
    double minJetPt;
    double minGenJetPt;
    double qScale_max;
    double eta_thresh_low;
    double eta_thresh_high;

  //  std::vector<uhh2bacon::run_lumi_ev>  unprefirableSummary;
  //    bool useUnprefireable;

  uhh2::GenericEvent::Handle<std::vector<L1Jet>> handle_l1jet_seeds;
  

};

void AnalysisModule_JetGenRecoMatch_RelValsTest::init_JEC(uhh2::Context& ctx){
  //FixME: switch from AK4 to AK8, ClosureTest
  // Jet energy corrections
  std::vector<std::string> JEC_AK4CHS_MC, JEC_AK8CHS_MC, JEC_AK4Puppi_MC, JEC_AK8Puppi_MC;
  std::vector<std::string> JEC_AK4CHS_A, JEC_AK4CHS_B, JEC_AK4CHS_C, JEC_AK4CHS_D, JEC_AK4CHS_E, JEC_AK4CHS_F, JEC_AK4CHS_G, JEC_AK4CHS_H;
  std::vector<std::string> JEC_AK4Puppi_A, JEC_AK4Puppi_B, JEC_AK4Puppi_C, JEC_AK4Puppi_D, JEC_AK4Puppi_E, JEC_AK4Puppi_F, JEC_AK4Puppi_G, JEC_AK4Puppi_H;
  std::vector<std::string> JEC_AK8CHS_A, JEC_AK8CHS_B, JEC_AK8CHS_C, JEC_AK8CHS_D, JEC_AK8CHS_E, JEC_AK8CHS_F, JEC_AK8CHS_G, JEC_AK8CHS_H;
  std::vector<std::string> JEC_AK8Puppi_A, JEC_AK8Puppi_B, JEC_AK8Puppi_C, JEC_AK8Puppi_D, JEC_AK8Puppi_E, JEC_AK8Puppi_F, JEC_AK8Puppi_G, JEC_AK8Puppi_H;

  std::vector<std::string> JEC_AK4CHS_MC_L1RC, JEC_AK8CHS_MC_L1RC, JEC_AK4Puppi_MC_L1RC, JEC_AK8Puppi_MC_L1RC;
  std::vector<std::string> JEC_AK4CHS_A_L1RC, JEC_AK4CHS_B_L1RC, JEC_AK4CHS_C_L1RC, JEC_AK4CHS_D_L1RC, JEC_AK4CHS_E_L1RC, JEC_AK4CHS_F_L1RC, JEC_AK4CHS_G_L1RC, JEC_AK4CHS_H_L1RC;
  std::vector<std::string> JEC_AK4Puppi_A_L1RC, JEC_AK4Puppi_B_L1RC, JEC_AK4Puppi_C_L1RC, JEC_AK4Puppi_D_L1RC, JEC_AK4Puppi_E_L1RC, JEC_AK4Puppi_F_L1RC, JEC_AK4Puppi_G_L1RC,JEC_AK4Puppi_H_L1RC;
  std::vector<std::string> JEC_AK8CHS_A_L1RC, JEC_AK8CHS_B_L1RC, JEC_AK8CHS_C_L1RC, JEC_AK8CHS_D_L1RC, JEC_AK8CHS_E_L1RC, JEC_AK8CHS_F_L1RC, JEC_AK8CHS_G_L1RC, JEC_AK8CHS_H_L1RC;
  std::vector<std::string> JEC_AK8Puppi_A_L1RC, JEC_AK8Puppi_B_L1RC, JEC_AK8Puppi_C_L1RC, JEC_AK8Puppi_D_L1RC, JEC_AK8Puppi_E_L1RC, JEC_AK8Puppi_F_L1RC, JEC_AK8Puppi_G_L1RC, JEC_AK8Puppi_H_L1RC;
  
  if(!is2017 && !is2018 && !is2016v2 && !is2016v3) cout<<"Hm, some is wrong. Did you mention year in Version parameter of xml file?"<<endl;
  if(is2016v2){
    cout<<"AnalysisModule_JetGenRecoMatch_RelValsTest uses JEC for 2016 data/MC"<<endl;
    JEC_AK4CHS_B       = JERFiles::Summer16_07Aug2017_V11_B_L123_AK4PFchs_DATA;
    JEC_AK4CHS_C       = JERFiles::Summer16_07Aug2017_V11_C_L123_AK4PFchs_DATA;
    JEC_AK4CHS_D       = JERFiles::Summer16_07Aug2017_V11_D_L123_AK4PFchs_DATA;
    JEC_AK4CHS_E       = JERFiles::Summer16_07Aug2017_V11_E_L123_AK4PFchs_DATA;
    JEC_AK4CHS_F       = JERFiles::Summer16_07Aug2017_V11_F_L123_AK4PFchs_DATA;
    JEC_AK4CHS_G       = JERFiles::Summer16_07Aug2017_V11_G_L123_AK4PFchs_DATA;
    JEC_AK4CHS_H       = JERFiles::Summer16_07Aug2017_V11_H_L123_AK4PFchs_DATA;
    JEC_AK4CHS_MC       = JERFiles::Summer16_07Aug2017_V11_L123_AK4PFchs_MC;

    JEC_AK8CHS_B       = JERFiles::Summer16_07Aug2017_V11_B_L123_AK8PFchs_DATA;
    JEC_AK8CHS_C       = JERFiles::Summer16_07Aug2017_V11_C_L123_AK8PFchs_DATA;
    JEC_AK8CHS_D       = JERFiles::Summer16_07Aug2017_V11_D_L123_AK8PFchs_DATA;
    JEC_AK8CHS_E       = JERFiles::Summer16_07Aug2017_V11_E_L123_AK8PFchs_DATA;
    JEC_AK8CHS_F       = JERFiles::Summer16_07Aug2017_V11_F_L123_AK8PFchs_DATA;
    JEC_AK8CHS_G       = JERFiles::Summer16_07Aug2017_V11_G_L123_AK8PFchs_DATA;
    JEC_AK8CHS_H       = JERFiles::Summer16_07Aug2017_V11_H_L123_AK8PFchs_DATA;
    JEC_AK8CHS_MC       = JERFiles::Summer16_07Aug2017_V11_L123_AK8PFchs_MC;

    JEC_AK4Puppi_B = JERFiles::Summer16_07Aug2017_V11_B_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C = JERFiles::Summer16_07Aug2017_V11_C_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D = JERFiles::Summer16_07Aug2017_V11_D_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_E = JERFiles::Summer16_07Aug2017_V11_E_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_F = JERFiles::Summer16_07Aug2017_V11_F_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_G = JERFiles::Summer16_07Aug2017_V11_G_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_H = JERFiles::Summer16_07Aug2017_V11_H_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK4PFPuppi_MC;

    JEC_AK8Puppi_B = JERFiles::Summer16_07Aug2017_V11_B_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C = JERFiles::Summer16_07Aug2017_V11_C_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D = JERFiles::Summer16_07Aug2017_V11_D_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_E = JERFiles::Summer16_07Aug2017_V11_E_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_F = JERFiles::Summer16_07Aug2017_V11_F_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_G = JERFiles::Summer16_07Aug2017_V11_G_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_H = JERFiles::Summer16_07Aug2017_V11_H_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC = JERFiles::Summer16_07Aug2017_V11_L123_AK8PFPuppi_MC;

    JEC_AK4CHS_B_L1RC       = JERFiles::Summer16_07Aug2017_V11_B_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_C_L1RC       = JERFiles::Summer16_07Aug2017_V11_C_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_D_L1RC       = JERFiles::Summer16_07Aug2017_V11_D_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_E_L1RC       = JERFiles::Summer16_07Aug2017_V11_E_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_F_L1RC       = JERFiles::Summer16_07Aug2017_V11_F_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_G_L1RC       = JERFiles::Summer16_07Aug2017_V11_G_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_H_L1RC       = JERFiles::Summer16_07Aug2017_V11_H_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_MC_L1RC       = JERFiles::Summer16_07Aug2017_V11_L1RC_AK4PFchs_MC;

    JEC_AK8CHS_B_L1RC       = JERFiles::Summer16_07Aug2017_V11_B_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_C_L1RC       = JERFiles::Summer16_07Aug2017_V11_C_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_D_L1RC       = JERFiles::Summer16_07Aug2017_V11_D_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_E_L1RC       = JERFiles::Summer16_07Aug2017_V11_E_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_F_L1RC       = JERFiles::Summer16_07Aug2017_V11_F_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_G_L1RC      = JERFiles::Summer16_07Aug2017_V11_G_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_H_L1RC       = JERFiles::Summer16_07Aug2017_V11_H_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_MC_L1RC       = JERFiles::Summer16_07Aug2017_V11_L1RC_AK8PFchs_MC;

    JEC_AK4Puppi_B_L1RC = JERFiles::Summer16_07Aug2017_V11_B_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C_L1RC = JERFiles::Summer16_07Aug2017_V11_C_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D_L1RC = JERFiles::Summer16_07Aug2017_V11_D_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_E_L1RC = JERFiles::Summer16_07Aug2017_V11_E_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_F_L1RC = JERFiles::Summer16_07Aug2017_V11_F_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_G_L1RC = JERFiles::Summer16_07Aug2017_V11_G_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_H_L1RC = JERFiles::Summer16_07Aug2017_V11_H_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC_L1RC = JERFiles::Summer16_07Aug2017_V11_L1RC_AK4PFPuppi_MC;

    JEC_AK8Puppi_B_L1RC = JERFiles::Summer16_07Aug2017_V11_B_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C_L1RC = JERFiles::Summer16_07Aug2017_V11_C_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D_L1RC = JERFiles::Summer16_07Aug2017_V11_D_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_E_L1RC = JERFiles::Summer16_07Aug2017_V11_E_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_F_L1RC = JERFiles::Summer16_07Aug2017_V11_F_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_G_L1RC = JERFiles::Summer16_07Aug2017_V11_G_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_H_L1RC = JERFiles::Summer16_07Aug2017_V11_H_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC_L1RC = JERFiles::Summer16_07Aug2017_V11_L1RC_AK8PFPuppi_MC;
  }

  if(is2017){
    cout<<"AnalysisModule_JetGenRecoMatch_RelValsTest uses JEC for 2017 data/MC"<<endl;
    JEC_AK4CHS_B       = JERFiles::Fall17_17Nov2017_V32_B_L123_AK4PFchs_DATA;
    JEC_AK4CHS_C       = JERFiles::Fall17_17Nov2017_V32_C_L123_AK4PFchs_DATA;
    JEC_AK4CHS_D       = JERFiles::Fall17_17Nov2017_V32_D_L123_AK4PFchs_DATA;
    JEC_AK4CHS_E       = JERFiles::Fall17_17Nov2017_V32_E_L123_AK4PFchs_DATA;
    JEC_AK4CHS_F       = JERFiles::Fall17_17Nov2017_V32_F_L123_AK4PFchs_DATA;
    JEC_AK4CHS_MC       = JERFiles::Fall17_17Nov2017_V32_L123_AK4PFchs_MC;

    JEC_AK8CHS_B       = JERFiles::Fall17_17Nov2017_V32_B_L123_AK8PFchs_DATA;
    JEC_AK8CHS_C       = JERFiles::Fall17_17Nov2017_V32_C_L123_AK8PFchs_DATA;
    JEC_AK8CHS_D       = JERFiles::Fall17_17Nov2017_V32_D_L123_AK8PFchs_DATA;
    JEC_AK8CHS_E       = JERFiles::Fall17_17Nov2017_V32_E_L123_AK8PFchs_DATA;
    JEC_AK8CHS_F       = JERFiles::Fall17_17Nov2017_V32_F_L123_AK8PFchs_DATA;
    JEC_AK8CHS_MC       = JERFiles::Fall17_17Nov2017_V32_L123_AK8PFchs_MC;

    JEC_AK4Puppi_B = JERFiles::Fall17_17Nov2017_V32_B_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C = JERFiles::Fall17_17Nov2017_V32_C_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D = JERFiles::Fall17_17Nov2017_V32_D_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_E = JERFiles::Fall17_17Nov2017_V32_E_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_F = JERFiles::Fall17_17Nov2017_V32_F_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC = JERFiles::Fall17_17Nov2017_V32_L123_AK4PFPuppi_MC;

    JEC_AK8Puppi_B = JERFiles::Fall17_17Nov2017_V32_B_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C = JERFiles::Fall17_17Nov2017_V32_C_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D = JERFiles::Fall17_17Nov2017_V32_D_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_E = JERFiles::Fall17_17Nov2017_V32_E_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_F = JERFiles::Fall17_17Nov2017_V32_F_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC = JERFiles::Fall17_17Nov2017_V32_L123_AK8PFPuppi_MC;

    JEC_AK4CHS_B_L1RC       = JERFiles::Fall17_17Nov2017_V32_B_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_C_L1RC       = JERFiles::Fall17_17Nov2017_V32_C_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_D_L1RC       = JERFiles::Fall17_17Nov2017_V32_D_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_E_L1RC       = JERFiles::Fall17_17Nov2017_V32_E_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_F_L1RC       = JERFiles::Fall17_17Nov2017_V32_F_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_MC_L1RC       = JERFiles::Fall17_17Nov2017_V32_L1RC_AK4PFchs_MC;

    JEC_AK8CHS_B_L1RC       = JERFiles::Fall17_17Nov2017_V32_B_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_C_L1RC       = JERFiles::Fall17_17Nov2017_V32_C_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_D_L1RC       = JERFiles::Fall17_17Nov2017_V32_D_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_E_L1RC       = JERFiles::Fall17_17Nov2017_V32_E_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_F_L1RC       = JERFiles::Fall17_17Nov2017_V32_F_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_MC_L1RC       = JERFiles::Fall17_17Nov2017_V32_L1RC_AK8PFchs_MC;

    JEC_AK4Puppi_B_L1RC = JERFiles::Fall17_17Nov2017_V32_B_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C_L1RC = JERFiles::Fall17_17Nov2017_V32_C_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D_L1RC = JERFiles::Fall17_17Nov2017_V32_D_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_E_L1RC = JERFiles::Fall17_17Nov2017_V32_E_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_F_L1RC = JERFiles::Fall17_17Nov2017_V32_F_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC_L1RC = JERFiles::Fall17_17Nov2017_V32_L1RC_AK4PFPuppi_MC;

    JEC_AK8Puppi_B_L1RC = JERFiles::Fall17_17Nov2017_V32_B_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C_L1RC = JERFiles::Fall17_17Nov2017_V32_C_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D_L1RC = JERFiles::Fall17_17Nov2017_V32_D_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_E_L1RC = JERFiles::Fall17_17Nov2017_V32_E_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_F_L1RC = JERFiles::Fall17_17Nov2017_V32_F_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC_L1RC = JERFiles::Fall17_17Nov2017_V32_L1RC_AK8PFPuppi_MC;
  }
  if(is2018){
    cout<<"AnalysisModule_JetGenRecoMatch_RelValsTest uses JEC for 2018 data/MC"<<endl;
    JEC_AK4CHS_A       = JERFiles::Autumn18_V8_A_L123_AK4PFchs_DATA;
    JEC_AK4CHS_B       = JERFiles::Autumn18_V8_B_L123_AK4PFchs_DATA;
    JEC_AK4CHS_C       = JERFiles::Autumn18_V8_C_L123_AK4PFchs_DATA;
    JEC_AK4CHS_D       = JERFiles::Autumn18_V8_D_L123_AK4PFchs_DATA;
    JEC_AK4CHS_MC       = JERFiles::Autumn18_V8_L123_AK4PFchs_MC;

    JEC_AK8CHS_A       = JERFiles::Autumn18_V8_A_L123_AK8PFchs_DATA;
    JEC_AK8CHS_B       = JERFiles::Autumn18_V8_B_L123_AK8PFchs_DATA;
    JEC_AK8CHS_C       = JERFiles::Autumn18_V8_C_L123_AK8PFchs_DATA;
    JEC_AK8CHS_D       = JERFiles::Autumn18_V8_D_L123_AK8PFchs_DATA;
    JEC_AK8CHS_MC      = JERFiles::Autumn18_V8_L123_AK8PFchs_MC;

    JEC_AK4Puppi_A = JERFiles::Autumn18_V8_A_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_B = JERFiles::Autumn18_V8_B_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C = JERFiles::Autumn18_V8_C_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D = JERFiles::Autumn18_V8_D_L123_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC = JERFiles::Autumn18_V8_L123_AK4PFPuppi_MC;

    JEC_AK8Puppi_A = JERFiles::Autumn18_V8_A_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_B = JERFiles::Autumn18_V8_B_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C = JERFiles::Autumn18_V8_C_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D = JERFiles::Autumn18_V8_D_L123_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC = JERFiles::Autumn18_V8_L123_AK8PFPuppi_MC;

    JEC_AK4CHS_A_L1RC       = JERFiles::Autumn18_V8_A_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_B_L1RC      = JERFiles::Autumn18_V8_B_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_C_L1RC       = JERFiles::Autumn18_V8_C_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_D_L1RC       = JERFiles::Autumn18_V8_D_L1RC_AK4PFchs_DATA;
    JEC_AK4CHS_MC_L1RC       = JERFiles::Autumn18_V8_L1RC_AK4PFchs_MC;

    JEC_AK8CHS_A_L1RC       = JERFiles::Autumn18_V8_A_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_B_L1RC       = JERFiles::Autumn18_V8_B_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_C_L1RC       = JERFiles::Autumn18_V8_C_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_D_L1RC       = JERFiles::Autumn18_V8_D_L1RC_AK8PFchs_DATA;
    JEC_AK8CHS_MC_L1RC       = JERFiles::Autumn18_V8_L1RC_AK8PFchs_MC;

    JEC_AK4Puppi_A_L1RC = JERFiles::Autumn18_V8_A_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_B_L1RC = JERFiles::Autumn18_V8_B_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_C_L1RC = JERFiles::Autumn18_V8_C_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_D_L1RC = JERFiles::Autumn18_V8_D_L1RC_AK4PFPuppi_DATA;
    JEC_AK4Puppi_MC_L1RC = JERFiles::Autumn18_V8_L1RC_AK4PFPuppi_MC;

    JEC_AK8Puppi_A_L1RC = JERFiles::Autumn18_V8_A_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_B_L1RC = JERFiles::Autumn18_V8_B_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_C_L1RC = JERFiles::Autumn18_V8_C_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_D_L1RC = JERFiles::Autumn18_V8_D_L1RC_AK8PFPuppi_DATA;
    JEC_AK8Puppi_MC_L1RC = JERFiles::Autumn18_V8_L1RC_AK8PFPuppi_MC;
  }

  //jet_corrector_B.reset(new JetCorrector(ctx, JEC_corr_B, JEC_corr_B_L1RC));
  if(jetLabel == "AK4CHS" || jetLabel == "AK4PUPPI"){
    if(jetLabel == "AK4CHS"){
      cout<<"JECs: AK4CHS"<<endl;
      if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK4CHS_A,JEC_AK4CHS_A_L1RC));
      jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK4CHS_B,JEC_AK4CHS_B_L1RC));
      jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK4CHS_C,JEC_AK4CHS_C_L1RC));
      jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4CHS_D,JEC_AK4CHS_D_L1RC));
      if(is2017 or is2016v2 or is2016v3){
	//  jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4CHS_D,JEC_AK4CHS_D_L1RC));//tmp
	jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK4CHS_E,JEC_AK4CHS_E_L1RC));
	jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK4CHS_F,JEC_AK4CHS_F_L1RC));
	if(is2016v2 or is2016v3){
	  jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK4CHS_G,JEC_AK4CHS_G_L1RC));
	  jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK4CHS_H,JEC_AK4CHS_H_L1RC));
	}
      }
      jet_corrector_MC.reset(new JetCorrector(ctx,            JEC_AK4CHS_MC,JEC_AK4CHS_MC_L1RC));    
    }
    else{
      cout<<"JECs: AK4PUPPI"<<endl;
      if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK4Puppi_A));
      jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK4Puppi_B));
      jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK4Puppi_C));
      jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK4Puppi_D));
      if(is2017 or is2016v2 or is2016v3){
	jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK4Puppi_E));
	jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK4Puppi_F));
	if(is2016v2 or is2016v3){
	  jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK4Puppi_G));
	  jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK4Puppi_H));
	}
      }
      jet_corrector_MC.reset(new JetCorrector(ctx, JEC_AK4Puppi_MC));    
    }
  }
  else if (jetLabel == "AK8CHS" || jetLabel == "AK8PUPPI"){
    if(jetLabel == "AK8CHS"){
      cout<<"JECs: AK8CHS"<<endl;
      if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK8CHS_A,JEC_AK8CHS_A_L1RC));
      jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK8CHS_B,JEC_AK8CHS_B_L1RC));
      jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK8CHS_C,JEC_AK8CHS_C_L1RC));
      jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK8CHS_D,JEC_AK8CHS_D_L1RC));
      if(is2017 or is2016v2 or is2016v3){
	//  jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK8CHS_D,JEC_AK8CHS_D_L1RC));//tmp
	jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK8CHS_E,JEC_AK8CHS_E_L1RC));
	jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK8CHS_F,JEC_AK8CHS_F_L1RC));
	if(is2016v2 or is2016v3){
	  jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK8CHS_G,JEC_AK8CHS_G_L1RC));
	  jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK8CHS_H,JEC_AK8CHS_H_L1RC));
	}
      }
      jet_corrector_MC.reset(new JetCorrector(ctx,            JEC_AK8CHS_MC,JEC_AK8CHS_MC_L1RC));    
      //      cout<<"use AK8CHS JECs!"<<endl;
    }
    else{
      cout<<"JECs: AK8PUPPI"<<endl;
      if(is2018) jet_corrector_A.reset(new JetCorrector(ctx, JEC_AK8Puppi_A));
      jet_corrector_B.reset(new JetCorrector(ctx, JEC_AK8Puppi_B));
      jet_corrector_C.reset(new JetCorrector(ctx, JEC_AK8Puppi_C));
      jet_corrector_D.reset(new JetCorrector(ctx, JEC_AK8Puppi_D));
      if(is2017 or is2016v2 or is2016v3){
	jet_corrector_E.reset(new JetCorrector(ctx, JEC_AK8Puppi_E));
	jet_corrector_F.reset(new JetCorrector(ctx, JEC_AK8Puppi_F));
	if(is2016v2 or is2016v3){
	  jet_corrector_G.reset(new JetCorrector(ctx, JEC_AK8Puppi_G));
	  jet_corrector_H.reset(new JetCorrector(ctx, JEC_AK8Puppi_H));
	}
      }
      cout<<"use AK8PUPPI JECs!"<<endl;
      jet_corrector_MC.reset(new JetCorrector(ctx, JEC_AK8Puppi_MC));    
    }
  }
  else
    cout<<"Have no idea about JEC for jetLabel= "<<jetLabel<<endl;
}

void AnalysisModule_JetGenRecoMatch_RelValsTest::declare_output(uhh2::Context& ctx){
 //Store only vars useful for debuging later
    tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
    tt_gen_PUpthat = ctx.declare_event_output<float>("gen_PUpthat");
    tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
    tt_run = ctx.declare_event_output<int>("run");
    tt_evID = ctx.declare_event_output<int>("eventID");
    tt_lumiSec = ctx.declare_event_output<int>("lumi_sec");
    tt_nvertices = ctx.declare_event_output<int>("nvertices");
    tt_nGoodvertices = ctx.declare_event_output<int>("nGoodvertices");
    tt_nPU = ctx.declare_event_output<int>("nPU");
    tt_ev_weight = ctx.declare_event_output<float>("weight");
    tt_jet_n= ctx.declare_event_output<int>("Njet");
    tt_genjet_n= ctx.declare_event_output<int>("Ngenjet");
    tt_matchedjet_n= ctx.declare_event_output<int>("Nmatchedrecojet");
    tt_rho = ctx.declare_event_output<float>("rho");

    tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
    tt_jet1_eta = ctx.declare_event_output<float>("jet1_eta");
    tt_jet1_phi = ctx.declare_event_output<float>("jet1_phi");
    tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
    tt_jet2_eta = ctx.declare_event_output<float>("jet2_eta");
    tt_jet2_phi = ctx.declare_event_output<float>("jet2_phi");
    tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
    tt_jet3_eta = ctx.declare_event_output<float>("jet3_eta");
    tt_jet3_phi = ctx.declare_event_output<float>("jet3_phi");

    tt_jet1_chEM = ctx.declare_event_output<float>("jet1_chEM");
    tt_jet2_chEM = ctx.declare_event_output<float>("jet2_chEM");
    tt_jet3_chEM = ctx.declare_event_output<float>("jet3_chEM");
    tt_jet1_neutEM = ctx.declare_event_output<float>("jet1_neutEM");
    tt_jet2_neutEM = ctx.declare_event_output<float>("jet2_neutEM");
    tt_jet3_neutEM = ctx.declare_event_output<float>("jet3_neutEM");
    tt_jet1_chHAD = ctx.declare_event_output<float>("jet1_chHAD");
    tt_jet2_chHAD = ctx.declare_event_output<float>("jet2_chHAD");
    tt_jet3_chHAD = ctx.declare_event_output<float>("jet3_chHAD");
    tt_jet1_neutHAD = ctx.declare_event_output<float>("jet1_neutHAD");
    tt_jet2_neutHAD = ctx.declare_event_output<float>("jet2_neutHAD");
    tt_jet3_neutHAD = ctx.declare_event_output<float>("jet3_neutHAD");


}

void AnalysisModule_JetGenRecoMatch_RelValsTest::init_hists(uhh2::Context& ctx){
    h_input.reset(new JECCrossCheckHists(ctx,"CrossCheck_input"));
    //    h_lumisel.reset(new JECCrossCheckHists(ctx,"CrossCheck_lumisel"));
    h_beforeCleaner.reset(new JECCrossCheckHists(ctx,"CrossCheck_beforeCleaner"));
    h_afterCleaner.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterCleaner"));
   
    h_afternVts.reset(new JECCrossCheckHists(ctx,"CrossCheck_afternVts"));

    h_matched_all.reset(new JECAnalysisRecoGenMatchedHistsFractions(ctx,"MatchedRecoGen_all"));
    for(int i=0;i<n_pt-1;i++){
      TString ptname = pt_range[i]; ptname +="_"; ptname +=pt_range[i+1];
      TString histname = "MatchedRecoGen_"; histname+=ptname;
      cout<<histname.Data()<<endl;
      h_matched_pt[i].reset(new JECAnalysisRecoGenMatchedHistsFractions(ctx,histname.Data()));
    }

    h_goodRECO_all.reset(new JECAnalysisPUjetsHists(ctx,"goodRECO_all"));
    h_badRECO_all.reset(new JECAnalysisPUjetsHists(ctx,"badRECO_all"));

    for(int i=0;i<n_pt-1;i++){
      TString ptname = pt_range[i]; ptname +="_"; ptname +=pt_range[i+1];
      TString histname = "goodRECO_"; histname+=ptname;
      h_goodRECO_pt[i].reset(new JECAnalysisPUjetsHists(ctx,histname.Data()));
      histname = "badRECO_"; histname+=ptname;
      h_badRECO_pt[i].reset(new JECAnalysisPUjetsHists(ctx,histname.Data()));
    }



    // h_nocuts.reset(new JECAnalysisHists(ctx,"NoCuts"));
    // h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    // h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    // h_sel.reset(new JECAnalysisHists(ctx,"Selection"));
    // h_final.reset(new JECAnalysisHists(ctx,"Final"));
    
    // h_trg40.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve40"));
    // h_trg60.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve60"));
    // h_trg80.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve80"));
    // h_trg140.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve140"));
    // h_trg200.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve200"));
    // h_trg260.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve260"));
    // h_trg320.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve320"));
    // h_trg400.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve400"));
    // h_trg500.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve500"));

    // h_trgHF60.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve60_HFJEC"));
    // h_trgHF80.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve80_HFJEC"));
    // h_trgHF100.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve100_HFJEC"));
    // h_trgHF160.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve160_HFJEC"));
    // h_trgHF220.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve220_HFJEC"));
    // h_trgHF300.reset(new JECAnalysisHists(ctx,"HLT_DiPFJetAve300_HFJEC"));
    
    // h_lumi_nocuts.reset(new LuminosityHists(ctx,"Lumi_noCuts"));  
    // h_lumi_sel.reset(new LuminosityHists(ctx,"Lumi_Selection"));
    // h_lumi_dijet.reset(new LuminosityHists(ctx,"Lumi_diJet"));
    // h_lumi_match.reset(new LuminosityHists(ctx,"Lumi_JetMatching"));
    // h_lumi_final.reset(new LuminosityHists(ctx,"Lumi_Final"));
    // h_lumi_Trig40.reset(new LuminosityHists(ctx,"Lumi_Trig40"));  
    // h_lumi_Trig60.reset(new LuminosityHists(ctx,"Lumi_Trig60")); 
    // h_lumi_Trig80.reset(new LuminosityHists(ctx,"Lumi_Trig80")); 
    // h_lumi_Trig140.reset(new LuminosityHists(ctx,"Lumi_Trig140")); 
    // h_lumi_Trig200.reset(new LuminosityHists(ctx,"Lumi_Trig200")); 
    // h_lumi_Trig260.reset(new LuminosityHists(ctx,"Lumi_Trig260")); 
    // h_lumi_Trig320.reset(new LuminosityHists(ctx,"Lumi_Trig320")); 
    // h_lumi_Trig400.reset(new LuminosityHists(ctx,"Lumi_Trig400")); 
    // h_lumi_Trig500.reset(new LuminosityHists(ctx,"Lumi_Trig500"));

    // h_lumi_TrigHF60.reset(new LuminosityHists(ctx,"Lumi_TrigHF60")); 
    // h_lumi_TrigHF80.reset(new LuminosityHists(ctx,"Lumi_TrigHF80")); 
    // h_lumi_TrigHF100.reset(new LuminosityHists(ctx,"Lumi_TrigHF100")); 
    // h_lumi_TrigHF160.reset(new LuminosityHists(ctx,"Lumi_TrigHF160")); 
    // h_lumi_TrigHF220.reset(new LuminosityHists(ctx,"Lumi_TrigHF220")); 
    // h_lumi_TrigHF300.reset(new LuminosityHists(ctx,"Lumi_TrigHF300"));
    // h_monitoring_final.reset(new LumiHists(ctx, "Monitoring_Final"));
    
}

 AnalysisModule_JetGenRecoMatch_RelValsTest::AnalysisModule_JetGenRecoMatch_RelValsTest(uhh2::Context & ctx) :
    sel(ctx)
  {
    // no_genp = false;
    //    no_genp=false;
    no_genp = true;
    //    cout << " no_genp = " << no_genp << endl;
    try{
      debug = ctx.get("Debug") == "true";
    }
    catch(const runtime_error& error){
      debug = false;
    }

    try{
      ispythia8 = ctx.get("MCgenerator") == "Pythia8";
    }
    catch(const runtime_error& error){
      ispythia8 = false;
    }

    //    const bool ispuppi = (ctx.get("is_puppi") == "true");
    ispuppi = false; //TEST
    cout << "Is this running on puppi: " << ispuppi << endl;
    is2016v2 = (ctx.get("dataset_version").find("2016v2") != std::string::npos);
    is2016v3 = (ctx.get("dataset_version").find("2016v3") != std::string::npos);
    is2017 = (ctx.get("dataset_version").find("2017") != std::string::npos);
    is2018 = (ctx.get("dataset_version").find("2018") != std::string::npos);
 
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    cout << "start" << endl;
    
   
    
    
    isMC = (ctx.get("dataset_type") == "MC");
    //    cout<<" isMC = "<< isMC<<" no_genp = "<<no_genp<<endl;
    if(isMC && !no_genp) cout<<"!!! WARNING, no genparticle are used! !!!"<<endl;
    //// COMMON MODULES

    L1METptThresh = stod(ctx.get("L1METptThresh"));
    eta_thresh_low = stod(ctx.get("eta_thresh_low"));
      
    if(!isMC) lumi_sel.reset(new LumiSelection(ctx));
    /* MET filters */
    if(!isMC){
      metfilters_sel.reset(new uhh2::AndSelection(ctx, "metfilters"));
      metfilters_sel->add<TriggerSelection>("1-good-vtx", "Flag_goodVertices"); 
      metfilters_sel->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter"); 
      metfilters_sel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");        
      metfilters_sel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
      metfilters_sel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter"); 
      metfilters_sel->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
      metfilters_sel->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter");
      metfilters_sel->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
      //    metfilters_sel->add<TriggerSelection>("ecalBadCalibFilter","Flag_ecalBadCalibFilter");
      metfilters_sel->add<TriggerSelection>("ecalBadCalibFilter","ecalBadCalibReducedMINIAODFilter");
    }
    jetLabel = ctx.get("JetLabel");
    init_JEC(ctx);//set up JECs
    // if(ispuppi)
    //   Jet_PFID = JetPFID(JetPFID::WP_TIGHT_PUPPI);
    // else
    //   Jet_PFID = JetPFID(JetPFID::WP_TIGHT_CHS);
    // jetID.reset(new JetCleaner(ctx, Jet_PFID));

    //remove low pt jets
    minJetPt = stod(ctx.get("minJetPt"));
    minGenJetPt = stod(ctx.get("minGenJetPt"));
    qScale_max = stod(ctx.get("qScale"));
    jetcleaner.reset(new JetCleaner(ctx, minJetPt, 5.2));
    genjetcleaner.reset(new GenJetCleaner(ctx, minGenJetPt, 5.2));

    //Lepton cleaner
    const     MuonId muoSR(AndId<Muon>    (MuonID(Muon::CutBasedIdTight),PtEtaCut  (15, 2.4)));
    const ElectronId eleSR(AndId<Electron>(ElectronID_MVA_Fall17_loose_noIso , PtEtaSCCut(15, 2.4)));  
    muoSR_cleaner.reset(new     MuonCleaner(muoSR));
    eleSR_cleaner.reset(new ElectronCleaner(eleSR)); 

    
    //new

    dataset_version = ctx.get("dataset_version");
    ClosureTest = (ctx.get("ClosureTest") == "true");
    apply_METoverPt_cut = (ctx.get("METoverPt_cut") == "true");
    apply_EtaPhi_cut = (ctx.get("EtaPhi_cut") == "true");
    //    JEC_Version = ctx.get("JEC_Version");

    apply_L1seed_from_bx1_filter =  (ctx.get("Apply_L1Seed_From_BX1_Filter") == "true" && !isMC);
    apply_smear = (ctx.get("Apply_MC_Smear")=="true" && isMC);

    split_JEC_MC   = false; //Different MC corrections only existed for Spring16_25ns_V8* 
    split_JEC_DATA = true; //TODO check the JEC!!!

    if(debug) std::cout<<"isMC: "<<isMC<<"  split_JEC_MC: "<<split_JEC_MC<<"  split_JEC_DATA: "<<split_JEC_DATA <<"   ClosureTest: "<<ClosureTest<<std::endl;
    
    //JER Smearing for corresponding JEC-Version
    if(isMC){
      if(is2018){ 
	cout<<"2018 JER smearer ..."<<endl;
	jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", JERSmearing::SF_13TeV_Autumn18_V1,"2018/Autumn18_V1_MC_PtResolution_AK4PFchs.txt"));//JER SFs obtained with Pythia8 Flat
	cout<<"JER smearer is ready!"<<endl;
      }
      if(is2017){
	cout<<"2017 JER smearer ..."<<endl;
	jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", JERSmearing::SF_13TeV_Fall17_V3_RunBCDEF_Madgraph,"2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));//JER SFs obtained with Madgraph
      }
      if(is2016v2 || is2016v3){ 
	cout<<"2016 JER smearer ..."<<endl;
	jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets",  JERSmearing::SF_13TeV_Summer16_25nsV1,"2016/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"));
      }
    }
    if(debug) cout<<"clean output"<<endl;

    // //output
    ctx.undeclare_all_event_output();//clean
    declare_output(ctx);//store only varaibles useful for dijet analysis

    init_hists(ctx);
    if(debug) cout<<"reset hists \n";
    
   
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    GenJet_printer.reset(new GenJetPrinter("GenJet-Printer", 0));
    
    if(!no_genp) 
      GenParticles_printer.reset(new GenParticlesPrinter(ctx));
     

    LumiWeight_module.reset(new MCLumiWeight(ctx));

    cout<<"end of AnalyseModule Constructor\n";
    
  };

 AnalysisModule_JetGenRecoMatch_RelValsTest::~AnalysisModule_JetGenRecoMatch_RelValsTest() {
  }

 bool AnalysisModule_JetGenRecoMatch_RelValsTest::process(Event & event) {
 //###############################################################
    //
    //Selection Module for Basic jet checks 
    // Veto events with leptons
    // Matches RECO and GEN jets
    // Fills histograms and stored variables for later analysis
    //###############################################################
   
   n_evt++;
   if(debug){
     cout << endl << "++++++++++ NEW EVENT +++++++++" << endl << endl;
     cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
     //      cout << "pfparticles.size() = " <<event.pfparticles->size()<<endl;
   }
   // Weight modules
   LumiWeight_module->process(event);

   //Dump Input
   h_input->fill(event);
   sel.SetEvent(event);

   if(debug) cout<<"after first fill\n";
   
   //LEPTON selection
   muoSR_cleaner->process(event);
   sort_by_pt<Muon>(*event.muons); 
   if(debug )std::cout<<"#muons = "<<event.muons->size()<<std::endl;
   
   eleSR_cleaner->process(event);
   sort_by_pt<Electron>(*event.electrons);
   if(debug) std::cout<<"#electrons = "<<event.electrons->size()<<std::endl;
   
   if (event.electrons->size()>0 || event.muons->size()>0) return false; //veto events with leptons
   if(debug) cout<<"No leptons in the event";    

   h_beforeCleaner->fill(event);
   
   //################################ JEC #######################################
   jet_corrector_MC->process(event);//Attention: no proper JEC for RelVals

   //Apply JER to all jet collections
   if(apply_smear){  if(jetER_smearer.get()){
       if(debug) cout<<"jet smearing will be done\n";
       jetER_smearer->process(event);
       if(debug) cout<<"after jet smearing\n";
     }
   }

   //################################ Cleaning low pt jets  #######################################    
   jetcleaner->process(event);
   genjetcleaner->process(event);
  
   if(event.jets->size()<1) return false;

   //################################ Cleaning for Flat QCD MC #######################################    

    float gen_pthat = 0; //pt hat (from QCD simulation)
    float gen_weight = 0;
    if(!event.isRealData){
      //      event.weight = 1.0; //TEST: Flat MC, set weight to 1
      gen_weight = event.weight;
      gen_pthat = event.genInfo->qScale();
      if(debug) cout<<"gen_pthat = "<<gen_pthat<<endl;
    }
    float nvertices = event.pvs->size(); 
    event.set(tt_nvertices, nvertices);
    float nPU = 0 ;
    if(!event.isRealData) 
      nPU = event.genInfo->pileup_TrueNumInteractions();
    double pu_pthat = -1;
    if(!event.isRealData) pu_pthat = event.genInfo->PU_pT_hat_max();
    if(debug) std::cout<<"pu_pthat = "<<pu_pthat<<" gen_pthat = "<<gen_pthat<<std::endl;
    if(!event.isRealData && !no_genp){
      if(!sel.PUpthat()) return false; //Pu_pt_hat/pt_hat Selection
      //      if(!sel.PtaveVsQScale(qScale_max)) return false;//Ptave/pt_hat selection; MadGraph 1.5, 1.2 Pythia8 1.2, Herwig++ 1.2
    }
    h_afterCleaner->fill(event);
    if(debug) cout<<"Passsed PUpthat and Ptave/pt_hat selection! "<<endl;
    event.set(tt_gen_pthat,gen_pthat);
    event.set(tt_gen_PUpthat,pu_pthat);
    event.set(tt_gen_weight,gen_weight);
    event.set(tt_run,event.run);
    event.set(tt_evID,event.event);
    event.set(tt_lumiSec,event.luminosityBlock);
    event.set(tt_nPU,nPU);
    event.set(tt_ev_weight,event.weight);
    event.set(tt_rho,event.rho);  


    //################################ Good vertex #######################################
    //FixME: vertex selection still makes sense??
    //good primary vertex
    int nGoodVts = sel.goodPVertex();
    if(nGoodVts<=0) return false;
    event.set(tt_nGoodvertices, nGoodVts);
    h_afternVts->fill(event);

    
    //    h_afterJEC->fill(event);



    //################################ Fill vars for jets ######################################
    event.set(tt_jet1_pt,-1); event.set(tt_jet1_eta,-1); event.set(tt_jet1_phi,-1);  
    event.set(tt_jet2_pt,-1); event.set(tt_jet2_eta,-1); event.set(tt_jet2_phi,-1);  
    event.set(tt_jet3_pt,-1); event.set(tt_jet3_eta,-1); event.set(tt_jet3_phi,-1);  
    event.set(tt_jet1_chEM,-1); event.set(tt_jet2_chEM,-1);    event.set(tt_jet3_chEM,-1);
    event.set(tt_jet1_neutEM,-1); event.set(tt_jet2_neutEM,-1);    event.set(tt_jet3_neutEM,-1);
    event.set(tt_jet1_chHAD,-1); event.set(tt_jet2_chHAD,-1);    event.set(tt_jet3_chHAD,-1);
    event.set(tt_jet1_neutHAD,-1); event.set(tt_jet2_neutHAD,-1);    event.set(tt_jet3_neutHAD,-1);

    for(unsigned int j=0; j<event.jets->size(); j++){
      if(j>3) break;
      Jet recoj = event.jets->at(j);
      if(j<1){
	event.set(tt_jet1_pt,recoj.pt());  
	event.set(tt_jet1_eta,recoj.eta());  
	event.set(tt_jet1_phi,recoj.phi());  
	event.set(tt_jet1_chEM,recoj.chargedEmEnergyFraction());  
	event.set(tt_jet1_neutEM,recoj.neutralEmEnergyFraction());
	event.set(tt_jet1_chHAD,recoj.chargedHadronEnergyFraction());    
	event.set(tt_jet1_neutHAD,recoj.neutralHadronEnergyFraction());    
      }
      if(j>0 && j<2){
	event.set(tt_jet2_pt,recoj.pt());  
	event.set(tt_jet2_eta,recoj.eta());  
	event.set(tt_jet2_phi,recoj.phi());  
	event.set(tt_jet2_chEM,recoj.chargedEmEnergyFraction());  
	event.set(tt_jet2_neutEM,recoj.neutralEmEnergyFraction());
	event.set(tt_jet2_chHAD,recoj.chargedHadronEnergyFraction());    
	event.set(tt_jet2_neutHAD,recoj.neutralHadronEnergyFraction());    
      }
      if(j>1 && j<3){
	event.set(tt_jet3_pt,recoj.pt());  
	event.set(tt_jet3_eta,recoj.eta());  
	event.set(tt_jet3_phi,recoj.phi()); 
	event.set(tt_jet3_chEM,recoj.chargedEmEnergyFraction());  
	event.set(tt_jet3_neutEM,recoj.neutralEmEnergyFraction());
	event.set(tt_jet3_chHAD,recoj.chargedHadronEnergyFraction());    
	event.set(tt_jet3_neutHAD,recoj.neutralHadronEnergyFraction());    
      }
    }

    //################################ Match GEN and RECO jet  #######################################    
    if(debug){
      Jet_printer->process(event);
    }
    if(debug && isMC){
      GenJet_printer->process(event);
      //GenParticles_printer->process(event);
    }
    const unsigned int jet_n = event.jets->size();
    const unsigned int genjet_n = event.genjets->size();
    event.set(tt_jet_n,jet_n);
    event.set(tt_genjet_n,genjet_n);
    event.set(tt_matchedjet_n,-1);

    //arrays contain one idx for each jet in the event. If -1: unmatched, else: idx of the closest jet with dR<=dr_cut
    int idx_jet_matching_GenRecojets[genjet_n];//for each GEN jet store RECO matched id
    int idx_matched_RecoGenjets[jet_n]; //for each RECO store matched GEN jet id

    if(isMC){    
      double dr_cut = 0;
      if(jetLabel == "AK4CHS" || jetLabel == "AK4PUPPI") dr_cut = 0.2;
      //      if(jetLabel == "AK4CHS" || jetLabel == "AK4PUPPI") dr_cut = 0.1; //TEST
      else if (jetLabel == "AK8CHS" || jetLabel == "AK8PUPPI") dr_cut = 0.4;
      else throw runtime_error("AnalysisModule_JetGenRecoMatch_RelValsTest: Invalid jet-label specified.");
      for(unsigned int i=0; i<event.genjets->size(); i++){
 	double dR_min = 99999; int idx_matching_jet = -1;
 	for(unsigned int j=0; j<event.jets->size(); j++){
 	  double dR = deltaR(event.jets->at(j), event.genjets->at(i));
	  if(dR<dR_min){
 	    dR_min = dR; 
 	    if(dR_min<dr_cut) idx_matching_jet = j;
 	  }
 	}
 	idx_jet_matching_GenRecojets[i] = idx_matching_jet;
	idx_matched_RecoGenjets[idx_matching_jet] = i;
 	if(debug) cout << "the jet matching the genjet no. " << i << " is jet no. " << idx_matching_jet <<" with dR_min="<<dR_min<<" (for reco jet #"<<idx_matching_jet<<")"<< endl;
	if(debug && idx_matching_jet>-1)
	  cout<<" with pt,rec/pt,gen = "<<event.jets->at(idx_matching_jet).pt()/event.genjets->at(i).pt()<< endl;
      }
    
      //fill histograms
      int n_matched_jets = 0;
      for(unsigned int i=0; i<event.genjets->size(); i++){
	if(idx_jet_matching_GenRecojets[i]>-1){
	  n_matched_jets++;
	  if(debug) cout<< "Fill plots the jet matching the genjet no. " << i << " is jet no. " << idx_jet_matching_GenRecojets[i] << endl;
	  h_matched_all->fill(event,i,idx_jet_matching_GenRecojets[i]);//event, gen_id, reco_id
	  bool isStop=false;
	  for(int j=0;j<n_pt-1;j++){
	    // if(event.genjets->at(i).pt()>120) 
	    //   cout<<"AAAA! event.genjets->at(i).pt() = "<<event.genjets->at(i).pt()<<" pt_bins[j] = "<<pt_bins[j]<<" pt_bins[j+1] = "<<pt_bins[j+1]<<endl;
	    if((event.genjets->at(i).pt())>pt_bins[j] && ((event.genjets->at(i).pt())<=pt_bins[j+1])){
	      // if(event.genjets->at(i).pt()>120) 
	      // 	cout<<" Picked up! "<<event.genjets->at(i).pt()<<" "<<event.jets->at(idx_jet_matching_GenRecojets[i]).pt()<<" pt_bins[j] = "<<pt_bins[j]<<" pt_bins[j+1] = "<<pt_bins[j+1]<<endl;
	      h_matched_pt[j]->fill(event,i,idx_jet_matching_GenRecojets[i]);//event, gen_id, reco_id
	      isStop = true;
	    }
	    if(isStop) break;
	  }
	}
      }
      event.set(tt_matchedjet_n,n_matched_jets);
      for(unsigned int j=0; j<event.jets->size(); j++){
	if(idx_matched_RecoGenjets[j]>-1){//RECO jet is matched
	  h_goodRECO_all->fill(event,j);//event, reco_id
	  bool isStop=false;
	  for(int k=0;k<n_pt-1;k++){
	    if((event.jets->at(j).pt())>pt_bins[k] && ((event.jets->at(j).pt())<=pt_bins[k+1])){
	      h_goodRECO_pt[k]->fill(event,j);//event, reco_id
	      isStop = true;
	    }
	    if(isStop) break;
	  }
	}
	else{//RECO jet is not matched, most probably PU jet
	  h_badRECO_all->fill(event,j);//event, reco_id
	  bool isStop=false;
	  for(int k=0;k<n_pt-1;k++){
	    if((event.jets->at(j).pt())>pt_bins[k] && ((event.jets->at(j).pt())<=pt_bins[k+1])){
	      h_badRECO_pt[k]->fill(event,j);//event, reco_id
	      isStop = true;
	    }
	    if(isStop) break;
	  }
	}
      }
    }
    //################################ Raw RECO jet  #######################################    

    //################################ End of story  #######################################    
    if(debug) cout<<"THIS IS THE END"<<endl;
    return true;
 }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(AnalysisModule_JetGenRecoMatch_RelValsTest)
