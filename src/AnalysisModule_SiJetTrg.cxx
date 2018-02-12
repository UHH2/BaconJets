#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/EventHelper.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Selection.h"
#include "../include/JECAnalysisHists.h"
#include "../include/JECCrossCheckHists.h"
#include "../include/JECRunnumberHists.h"

#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/LumiSelection.h> //includes also LuminosityHists.h
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/CleaningModules.h"
#include <UHH2/common/include/NSelections.h> 
#include <UHH2/common/include/MuonIds.h> 
#include <UHH2/common/include/ElectronIds.h>
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/constants.h"

#include "TClonesArray.h"
#include "TString.h"
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"

#include <stdexcept>

using namespace std;
using namespace uhh2;

class AnalysisModule_SiJetTrg: public uhh2::AnalysisModule {

  public:
    explicit AnalysisModule_SiJetTrg(uhh2::Context&);
    virtual bool process(uhh2::Event&) override;
    ~AnalysisModule_SiJetTrg();
    
  protected:
    // correctors
    std::unique_ptr<JetCorrector> jet_corrector, jet_corrector_BCD, jet_corrector_EFearly, jet_corrector_FlateG, jet_corrector_H;
    std::unique_ptr<GenericJetResolutionSmearer> jetER_smearer; 

  // cleaners
   std::unique_ptr<JetLeptonCleaner> jetleptoncleaner, JLC_BCD, JLC_EFearly, JLC_FlateG, JLC_H;
   std::unique_ptr<JetCleaner> jetcleaner;
   // std::unique_ptr<AnalysisModule>     muoSR_cleaner;   
   std::unique_ptr<AnalysisModule> eleSR_cleaner;

  std::vector<std::unique_ptr<AnalysisModule>> lep_cleaner;  

    // selections
    std::unique_ptr<uhh2::Selection> lumi_sel;
    std::unique_ptr<uhh2::AndSelection> metfilters_sel;     
    //    std::unique_ptr<uhh2::AndSelection> trigger_sel;

    std::unique_ptr<uhh2::Selection> triggerSiMu_sel;

    std::unique_ptr<uhh2::Selection> minBias_sel;
    std::unique_ptr<uhh2::Selection> trigger40_sel;
    std::unique_ptr<uhh2::Selection> trigger60_sel;
    std::unique_ptr<uhh2::Selection> trigger80_sel;
    std::unique_ptr<uhh2::Selection> trigger140_sel;
    std::unique_ptr<uhh2::Selection> trigger200_sel;
    std::unique_ptr<uhh2::Selection> trigger260_sel;
    std::unique_ptr<uhh2::Selection> trigger320_sel;
    std::unique_ptr<uhh2::Selection> trigger400_sel;
    std::unique_ptr<uhh2::Selection> trigger450_sel;
    std::unique_ptr<uhh2::Selection> trigger500_sel;
   
    //// Data/MC scale factors
    std::unique_ptr<uhh2::AnalysisModule> pileupSF;
    unique_ptr<AnalysisModule>  Jet_printer, GenParticles_printer;

    Event::Handle<float> tt_jet1_ptGen;  Event::Handle<float> tt_jet2_ptGen;  Event::Handle<float> tt_jet3_ptGen;
    Event::Handle<float> tt_gen_pthat; Event::Handle<float> tt_gen_weight;  Event::Handle<float> tt_gen_PUpthat;

 //Different energy fractions in jets
    Event::Handle<float> tt_probejet_neutEmEF; Event::Handle<float> tt_probejet_neutHadEF; 
    Event::Handle<float> tt_probejet_chEmEF; Event::Handle<float> tt_probejet_chHadEF; 
    Event::Handle<float> tt_probejet_photonEF; Event::Handle<float> tt_probejet_muonEF; 

    Event::Handle<float> tt_jet1_pt;     Event::Handle<float> tt_jet2_pt;     Event::Handle<float> tt_jet3_pt;
    Event::Handle<float> tt_jet1_ptRaw;  Event::Handle<float> tt_jet2_ptRaw;  Event::Handle<float> tt_jet3_ptRaw;
  Event::Handle<float> tt_jet1_pt_onoff_Resp;     Event::Handle<float> tt_jet2_pt_onoff_Resp;
    Event::Handle<int> tt_nvertices;
    Event::Handle<float> tt_probejet_eta;  Event::Handle<float> tt_probejet_phi; Event::Handle<float> tt_probejet_pt; Event::Handle<float> tt_probejet_ptRaw;
    Event::Handle<float> tt_barreljet_eta;  Event::Handle<float> tt_barreljet_phi; Event::Handle<float> tt_barreljet_pt; Event::Handle<float> tt_barreljet_ptRaw;
    Event::Handle<float> tt_probejet_ptgen;         Event::Handle<float> tt_barreljet_ptgen;     
    Event::Handle<float> tt_pt_ave;
    Event::Handle<float> tt_alpha;
    Event::Handle<float> tt_rel_r; Event::Handle<float> tt_mpf_r; 
    Event::Handle<float> tt_asymmetry;
    Event::Handle<float> tt_B;
    Event::Handle<float> tt_MET;
    Event::Handle<int> tt_nPU;
    Event::Handle<int> tt_matchJetId_0;
    Event::Handle<int> tt_matchJetId_1;		
    Event::Handle<float> tt_ev_weight;
    Event::Handle<float> tt_jets_pt;//sum of jets pT
    Event::Handle<int> tt_jet_n;//number of jets
    Event::Handle<float> tt_rho;//event energy density
    Event::Handle<int> tt_nGoodvertices;
    Event::Handle<int> tt_partonFlavor; //only MC
    Event::Handle<int> tt_flavorBarreljet, tt_flavorProbejet, tt_flavorLeadingjet, tt_flavorSubleadingjet; //only MC
    Event::Handle<float> tt_responseBarreljet, tt_responseProbejet; //only MC
    Event::Handle<float> tt_response_leadingjet;
    Event::Handle<float> tt_had_n_Efrac, tt_had_ch_Efrac, tt_mu_Efrac, tt_ph_Efrac;
    Event::Handle<float> tt_inst_lumi, tt_integrated_lumi_in_bin, tt_integrated_lumi;
    Event::Handle<int> tt_lumibin;
    Event::Handle<int> tt_Nmuon; Event::Handle<float> tt_muon_pt;
    Event::Handle<int> tt_Nele; Event::Handle<float> tt_ele_pt;
    Event::Handle<float> tt_dR_jet3_barreljet;  Event::Handle<float> tt_dR_jet3_probejet;
    Event::Handle<int> tt_triggerSiMu;
    Event::Handle<int> tt_minBias;  
    Event::Handle<int> tt_trigger40;
    Event::Handle<int> tt_trigger60;
    Event::Handle<int> tt_trigger80;
    Event::Handle<int> tt_trigger140;
    Event::Handle<int> tt_trigger200;
    Event::Handle<int> tt_trigger260;
    Event::Handle<int> tt_trigger320;
    Event::Handle<int> tt_trigger400;
    Event::Handle<int> tt_trigger450;
    Event::Handle<int> tt_trigger500;
  
    std::unique_ptr<JECAnalysisHists> h_nocuts, h_sel, h_dijet, h_match, h_final;
  std::unique_ptr<JECAnalysisHists> h_trgSiMu, h_minBias ,h_trg40, h_trg60, h_trg80, h_trg140, h_trg200,h_trg260,h_trg320,h_trg400,h_trg450,h_trg500;
    std::unique_ptr<LuminosityHists> h_lumi_nocuts, h_lumi_sel, h_lumi_dijet, h_lumi_match, h_lumi_final;    
  std::unique_ptr<LuminosityHists> h_lumi_TrigSiMu, h_lumi_minBias ,h_lumi_Trig40, h_lumi_Trig60, h_lumi_Trig80, h_lumi_Trig140, h_lumi_Trig200, h_lumi_Trig260, h_lumi_Trig320, h_lumi_Trig400,h_lumi_Trig450, h_lumi_Trig500;
    std::unique_ptr<JECRunnumberHists> h_runnr_input;
    std::unique_ptr<JECCrossCheckHists> h_input,h_lumisel, h_beforeCleaner,h_afterCleaner,h_2jets,h_beforeJEC,h_afterJEC,h_afterJER,h_afterMET,h_beforeTriggerData,h_afterTriggerData,h_beforeFlatFwd,h_afterFlatFwd,h_afterPtEtaReweight,h_afterLumiReweight,h_afterUnflat,h_afternVts;
    uhh2bacon::Selection sel;

    bool debug;
  bool isMC, split_JEC_DATA, split_JEC_MC, ClosureTest, apply_weights, apply_lumiweights, apply_unflattening, apply_METoverPt_cut, apply_EtaPhi_cut, trigger_central, trigger_fwd, ts, onlyBtB;
    double lumiweight;
    string jetLabel;
    TString dataset_version, JEC_Version;
    JetId Jet_PFID;
    int n_evt;
    std::unique_ptr<TFile> f_weights;
    
    std::map<run_lumi, double> rl2lumi;
    std::map<run_lumi, double> rl2intlumi;
    TBranch * brun ;
    TBranch * blumiblock;
    TBranch * bilumi;
    double integrated_lumi;
    vector<run_lumi> upper_binborders_runnrs;
    vector<double> lumi_in_bins;

  };

  AnalysisModule_SiJetTrg::AnalysisModule_SiJetTrg(uhh2::Context & ctx) :
		  sel(ctx)
  {
    try{
      debug = ctx.get("Debug") == "true";
    }
    catch(const runtime_error& error){
      debug = false;
    }
    
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    cout << "start" << endl;
    isMC = (ctx.get("dataset_type") == "MC");
    //// COMMON MODULES
    if(!isMC) lumi_sel.reset(new LumiSelection(ctx));
    /* MET filters */ 
    metfilters_sel.reset(new uhh2::AndSelection(ctx, "metfilters")); 
    metfilters_sel->add<TriggerSelection>("1-good-vtx", "Flag_goodVertices"); 
    metfilters_sel->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter"); 
    metfilters_sel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");        
    metfilters_sel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
    metfilters_sel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter"); 
    metfilters_sel->add<TriggerSelection>("CSCTightHalo2016Filter", "Flag_CSCTightHalo2016Filter"); 
   
    Jet_PFID = JetPFID(JetPFID::WP_LOOSE);
    //Jet_PFID = JetPFID(JetPFID::WP_TIGHT);
    jetcleaner.reset(new JetCleaner(ctx, Jet_PFID));
    
//Lepton cleaner
    MuonId muoSR;
    muoSR = AndId<Muon>(MuonID(Muon::CutBasedIdTight),PtEtaCut  (15.0, 2.4));
    const ElectronId eleSR(AndId<Electron>(ElectronID_PHYS14_25ns_tight , PtEtaSCCut(15.0, 2.4)));
   // muoSR_cleaner.reset(new     MuonCleaner(muoSR));
    lep_cleaner.emplace_back(new     MuonCleaner(muoSR));
    eleSR_cleaner.reset(new ElectronCleaner(eleSR)); 

//#############################################  Trigger  #########################################################
    trigger_central = (ctx.get("Trigger_Central") == "true");
    trigger_fwd     = (ctx.get("Trigger_FWD") == "true");

    ts  = (ctx.get("Trigger_Single") == "true"); //if true use single jet trigger, if false di jet trigger TODO collapse the SiJet and DiJEt AnalysisModules into one?
    // ts = true;
    onlyBtB = (ctx.get("Only_BtB") == "true");
    if(debug) cout<<"onlyBtb is "<<onlyBtB<<endl;
        
    if(!isMC){
    const std::string& triggerSiMu = ctx.get("triggerSiMu", "NULL");

    const std::string& minBias = ctx.get("minBias", "NULL");
    
    const std::string& trigger40 = ctx.get("trigger40", "NULL");
    const std::string& trigger60 = ctx.get("trigger60", "NULL");
    const std::string& trigger80 = ctx.get("trigger80", "NULL");
    const std::string& trigger140 = ctx.get("trigger140", "NULL");
    const std::string& trigger200 = ctx.get("trigger200", "NULL");
    const std::string& trigger260 = ctx.get("trigger260", "NULL");
    const std::string& trigger320 = ctx.get("trigger320", "NULL");
    const std::string& trigger400 = ctx.get("trigger400", "NULL");
    const std::string& trigger450 = ctx.get("trigger450", "NULL");
    const std::string& trigger500 = ctx.get("trigger500", "NULL");

      if(triggerSiMu != "NULL") triggerSiMu_sel.reset(new TriggerSelection(triggerSiMu));
      else triggerSiMu_sel.reset(new uhh2::AndSelection(ctx));

      if(minBias != "NULL") minBias_sel.reset(new TriggerSelection(minBias));
      else minBias_sel.reset(new uhh2::AndSelection(ctx));
      
      if(trigger40 != "NULL" and ts) trigger40_sel.reset(new TriggerSelection(trigger40));
      else trigger40_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger60 != "NULL" and ts) trigger60_sel.reset(new TriggerSelection(trigger60));
      else trigger60_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger80 != "NULL" and ts) trigger80_sel.reset(new TriggerSelection(trigger80));
      else trigger80_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger140 != "NULL" and ts) trigger140_sel.reset(new TriggerSelection(trigger140));
      else trigger140_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger200 != "NULL" and ts) trigger200_sel.reset(new TriggerSelection(trigger200));
      else trigger200_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger260 != "NULL" and ts) trigger260_sel.reset(new TriggerSelection(trigger260));
      else trigger260_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger320 != "NULL" and ts) trigger320_sel.reset(new TriggerSelection(trigger320));
      else trigger320_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger400 != "NULL" and ts) trigger400_sel.reset(new TriggerSelection(trigger400));
      else trigger400_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger450 != "NULL" and ts) trigger450_sel.reset(new TriggerSelection(trigger450));
      else trigger450_sel.reset(new uhh2::AndSelection(ctx));
      if(trigger500 != "NULL" and ts) trigger500_sel.reset(new TriggerSelection(trigger500));
      else trigger500_sel.reset(new uhh2::AndSelection(ctx));
    }

    //new
    jetLabel = ctx.get("JetLabel");
    dataset_version = ctx.get("dataset_version");
    ClosureTest = (ctx.get("ClosureTest") == "true");
    apply_METoverPt_cut = (ctx.get("METoverPt_cut") == "true");
    apply_EtaPhi_cut = (ctx.get("EtaPhi_cut") == "true");
    JEC_Version = ctx.get("JEC_Version");

    split_JEC_MC   = false; //Different MC corrections only existed for Spring16_25ns_V8* 
    split_JEC_DATA = false; //TODO check the JEC!!!
    
    std::vector<std::string> JEC_corr,       JEC_corr_BCD,       JEC_corr_EFearly,       JEC_corr_FlateG,       JEC_corr_H,      JEC_corr_MC_FlateGH;
    std::vector<std::string> JEC_corr_L1RC,  JEC_corr_BCD_L1RC,  JEC_corr_EFearly_L1RC,  JEC_corr_FlateG_L1RC,  JEC_corr_H_L1RC, JEC_corr_MC_FlateGH_L1RC;
    if(isMC){
      //for MC
      if(jetLabel == "AK4CHS"){
	if(!ClosureTest){
	  //residuals
	  if(JEC_Version == "Summer16_23Sep2016_V4"){
	    JEC_corr              = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;           //noRes only for DATA ;), only one version for MC for deriving Summer16_23Sep2016
	    JEC_corr_L1RC         = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;           
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;           //ReReco Data + Moriond17 MC
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }
	  //JEC-Version:V7 
	  else if(JEC_Version == "Summer16_03Feb2017_V3"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;           //noRes only for DATA ;), only one version for MC for deriving Summer16_23Sep2016
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;           
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;          
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }
	  //JEC-Version MC: V6 
	  else if(JEC_Version == "Summer16_03Feb2017_V4" || JEC_Version == "Summer16_03Feb2017_V5" || JEC_Version == "Summer16_03Feb2017_V6"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;           //noRes only for DATA ;), only one version for MC for deriving Summer16_23Sep2016
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;           
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;          
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }

	  else throw runtime_error("In AnalysisModule_SiJetTrg.cxx: Invalid JEC_Version for deriving residuals on AK4CHS, MC specified ("+JEC_Version+") ");
	}
	//closure
	else{
	 if(JEC_Version == "Summer16_23Sep2016_V4"){
	    JEC_corr              = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;               //ReReco Data + Summer16 MC
	    JEC_corr_L1RC         = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }
	 //JEC-Version MC: V1
	  else if(JEC_Version == "Summer16_03Feb2017_V3"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;               //ReReco Data + Summer16 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V1_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V1_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }
	 //JEC-Version MC: V6
	  else if(JEC_Version == "Summer16_03Feb2017_V4" || JEC_Version == "Summer16_03Feb2017_V5" || JEC_Version == "Summer16_03Feb2017_V6"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;               //ReReco Data + Summer16 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    //dummies, in this version, MC is not split
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V6_L123_AK4PFchs_MC;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V6_L1RC_AK4PFchs_MC;
	    cout << "This is MC, JECs used are: ";
	    for(unsigned int i=0; i<JEC_corr.size(); i++) cout << JEC_corr[i] << ", ";
	    cout << endl;
	  }

	 else throw runtime_error("In AnalysisModule_SiJetTrg.cxx: Invalid JEC_Version for closure test on AK4CHS, MC specified.");
	}
      }
    }
    else { 
      //for DATA
      if(jetLabel == "AK4CHS"){
	if(!ClosureTest){
	  //residuals
	 if(JEC_Version == "Summer16_23Sep2016_V4"){
	    JEC_corr              = JERFiles::Summer16_23Sep2016_V4_H_L123_noRes_AK4PFchs_DATA;  //ReReco Data + Moriond17 MC
	    JEC_corr_L1RC         = JERFiles::Summer16_23Sep2016_V4_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_23Sep2016_V4_BCD_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_23Sep2016_V4_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_23Sep2016_V4_EF_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_23Sep2016_V4_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_23Sep2016_V4_G_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_23Sep2016_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_23Sep2016_V4_H_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_23Sep2016_V4_H_L1RC_AK4PFchs_DATA;
	  }
	 //JEC-Version: V7
	  else if(JEC_Version == "Summer16_03Feb2017_V3"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V3_H_L123_noRes_AK4PFchs_DATA;  //ReReco Data + Moriond17 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V3_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V3_BCD_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V3_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V3_EF_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V3_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V3_G_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V3_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V3_H_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V3_H_L1RC_AK4PFchs_DATA;
	  }
	 //JEC-Version Data: V4
	  else if(JEC_Version == "Summer16_03Feb2017_V4"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V4_H_L123_noRes_AK4PFchs_DATA;  //ReReco Data + Moriond17 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V4_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V4_BCD_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V4_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V4_EF_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V4_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V4_G_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V4_H_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V4_H_L1RC_AK4PFchs_DATA;
	  }
	 //JEC-Version Data: V5
	  else if(JEC_Version == "Summer16_03Feb2017_V5"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V5_H_L123_noRes_AK4PFchs_DATA;  //ReReco Data + Moriond17 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V5_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V5_BCD_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V5_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V5_EF_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V5_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V5_G_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V5_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V5_H_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V5_H_L1RC_AK4PFchs_DATA;
	  }
	 //JEC-Version Data: V6
	  else if(JEC_Version == "Summer16_03Feb2017_V6"){
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V6_H_L123_noRes_AK4PFchs_DATA;  //ReReco Data + Moriond17 MC V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V6_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V6_BCD_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V6_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V6_EF_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V6_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V6_G_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V6_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V6_H_L123_noRes_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V6_H_L1RC_AK4PFchs_DATA;
	  }
	  else throw runtime_error("In AnalysisModule_SiJetTrg.cxx: Invalid JEC_Version for deriving residuals on AK4CHS, DATA specified.");
	}
	else{
	 if(JEC_Version == "Summer16_23Sep2016_V4"){
	    //closure
	    JEC_corr              = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;              //ReReco Data + Summer16 MC
	    JEC_corr_L1RC         = JERFiles::Summer16_23Sep2016_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_23Sep2016_V4_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_23Sep2016_V4_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_23Sep2016_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_23Sep2016_V4_H_L1RC_AK4PFchs_DATA;
	    cout << "JEC for DATA: Summer16_23Sep2016_V4_BCD/EFearly/FlateG/H_L123_AK4PFchs_DATA;" << endl;
	  }
	 //JEC-Version: V7
	  else if(JEC_Version == "Summer16_03Feb2017_V3"){
	    //closure
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V3_H_L123_AK4PFchs_DATA;              //ReReco Data + Summer16 V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V3_H_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V3_BCD_L123_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V3_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V3_EF_L123_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V3_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V3_G_L123_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V3_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V3_H_L123_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V3_H_L1RC_AK4PFchs_DATA;
	    cout << "JEC for DATA: Summer16_03Feb2017_V3_BCD/EFearly/FlateG/H_L123_AK4PFchs_DATA;" << endl;
	  }
	 //JEC-Version Data: V4
	  else if(JEC_Version == "Summer16_03Feb2017_V4"){
	    //closure
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V4_H_L123_AK4PFchs_DATA;              //ReReco Data + Summer16 V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V4_BCD_L123_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V4_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V4_EF_L123_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V4_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V4_G_L123_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V4_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V4_H_L123_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V4_H_L1RC_AK4PFchs_DATA;
	    cout << "JEC for DATA: Summer16_03Feb2017_V4_BCD/EFearly/FlateG/H_L123_AK4PFchs_DATA;" << endl;
	  }
	  else if(JEC_Version == "Summer16_03Feb2017_V5"){
	    //closure
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V5_G_L123_AK4PFchs_DATA;              //ReReco Data + Summer16 V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V5_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V5_BCD_L123_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V5_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V5_EF_L123_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V5_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V5_G_L123_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V5_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V5_H_L123_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V5_H_L1RC_AK4PFchs_DATA;
	    cout << "JEC for DATA: Summer16_03Feb2017_V5_BCD/EFearly/FlateG/H_L123_AK4PFchs_DATA;" << endl;
	  }
	  else if(JEC_Version == "Summer16_03Feb2017_V6"){
	    //closure
	    JEC_corr              = JERFiles::Summer16_03Feb2017_V6_H_L123_AK4PFchs_DATA;              //ReReco Data + Summer16 V7
	    JEC_corr_L1RC         = JERFiles::Summer16_03Feb2017_V6_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_BCD          = JERFiles::Summer16_03Feb2017_V6_BCD_L123_AK4PFchs_DATA;
	    JEC_corr_BCD_L1RC     = JERFiles::Summer16_03Feb2017_V6_BCD_L1RC_AK4PFchs_DATA;
	    JEC_corr_EFearly      = JERFiles::Summer16_03Feb2017_V6_EF_L123_AK4PFchs_DATA;
	    JEC_corr_EFearly_L1RC = JERFiles::Summer16_03Feb2017_V6_EF_L1RC_AK4PFchs_DATA;
	    JEC_corr_FlateG       = JERFiles::Summer16_03Feb2017_V6_G_L123_AK4PFchs_DATA;
	    JEC_corr_FlateG_L1RC  = JERFiles::Summer16_03Feb2017_V6_G_L1RC_AK4PFchs_DATA;
	    JEC_corr_H            = JERFiles::Summer16_03Feb2017_V6_H_L123_AK4PFchs_DATA;
	    JEC_corr_H_L1RC       = JERFiles::Summer16_03Feb2017_V6_H_L1RC_AK4PFchs_DATA;
	    cout << "JEC for DATA: Summer16_03Feb2017_V6_BCD/EFearly/FlateG/H_L123_AK4PFchs_DATA;" << endl;
	  }

	 else throw runtime_error("In AnalysisModule_SiJetTrg.cxx: Invalid JEC_Version for closure test on AK4CHS, DATA specified.");
	}
      }
    }
      
      //Only for closure test
      if(ClosureTest){
	//DATA
	if(!isMC){
	  if(split_JEC_DATA){ //these only exist for DATA 
	  jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_corr_BCD, JEC_corr_BCD_L1RC));
	  jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_corr_EFearly,  JEC_corr_EFearly_L1RC));
	  jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_corr_FlateG,  JEC_corr_FlateG_L1RC));
	  jet_corrector_H.reset(new JetCorrector(ctx, JEC_corr_H, JEC_corr_H_L1RC));
	  JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_BCD));
	  JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_corr_EFearly));
	  JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_corr_FlateG));
	  JLC_H.reset(new JetLeptonCleaner(ctx, JEC_corr_H));
	  }
	  else{
	    jet_corrector.reset(new JetCorrector(ctx, JEC_corr, JEC_corr_L1RC));	
	    jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
	  }
	}

	//MC
	 //For MC: only one version of JECs exists
	else if(isMC){
	  if(split_JEC_MC){
	    jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_corr_BCD, JEC_corr_BCD_L1RC));
	    jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_corr_EFearly, JEC_corr_EFearly_L1RC));
	    jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_corr_FlateG, JEC_corr_FlateG_L1RC));
	    jet_corrector_H.reset(new JetCorrector(ctx, JEC_corr_H, JEC_corr_H_L1RC));
	    JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_BCD));
	    JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_corr_EFearly));
	    JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_corr_FlateG));
	    JLC_H.reset(new JetLeptonCleaner(ctx, JEC_corr_H));
	  }
	  else{
	    jet_corrector.reset(new JetCorrector(ctx, JEC_corr, JEC_corr_L1RC));
	    jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
	  }
	}
      }
      //for residuals
      else{
	//DATA
	if(!isMC){
	  if(split_JEC_DATA){
	    jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_corr_BCD, JEC_corr_BCD_L1RC));
	    jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_corr_EFearly, JEC_corr_EFearly_L1RC));
	    jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_corr_FlateG, JEC_corr_FlateG_L1RC));
	    jet_corrector_H.reset(new JetCorrector(ctx, JEC_corr_H, JEC_corr_H_L1RC));
	    JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_BCD));
	    JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_corr_EFearly));
	    JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_corr_FlateG));
	    JLC_H.reset(new JetLeptonCleaner(ctx, JEC_corr_H));
	  }
	  else{
	    jet_corrector.reset(new JetCorrector(ctx, JEC_corr, JEC_corr_L1RC));
	    jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
	  }
	}
	//MC
	
	else if(isMC){
	  if(split_JEC_MC){
	    jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_corr_BCD, JEC_corr_BCD_L1RC));
	    jet_corrector_EFearly.reset(new JetCorrector(ctx, JEC_corr_EFearly, JEC_corr_EFearly_L1RC));
	    jet_corrector_FlateG.reset(new JetCorrector(ctx, JEC_corr_FlateG, JEC_corr_FlateG_L1RC));
	    jet_corrector_H.reset(new JetCorrector(ctx, JEC_corr_H, JEC_corr_H_L1RC));
	    JLC_BCD.reset(new JetLeptonCleaner(ctx, JEC_corr_BCD));
	    JLC_EFearly.reset(new JetLeptonCleaner(ctx, JEC_corr_EFearly));
	    JLC_FlateG.reset(new JetLeptonCleaner(ctx, JEC_corr_FlateG));
	    JLC_H.reset(new JetLeptonCleaner(ctx, JEC_corr_H));
	  }
	  else{
	    jet_corrector.reset(new JetCorrector(ctx, JEC_corr, JEC_corr_L1RC));
	    jetleptoncleaner.reset(new JetLeptonCleaner(ctx, JEC_corr));
	    cout << "setting up jet_corrector and JLC for MC, non-split JEC." << endl;
	  }
	}
      }
    
//JER Smearing for corresponding JEC-Version
      // if(isMC){
      // 	if(JEC_Version == "Summer16_23Sep2016_V4") jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017)); 
      // 	else if(JEC_Version == "Summer16_03Feb2017_V3") jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017)); 
      // 	else if(JEC_Version == "Summer16_03Feb2017_V4") jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017));
      // 	else if(JEC_Version == "Summer16_03Feb2017_V5") jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017));
      // 	else if(JEC_Version == "Summer16_03Feb2017_V6") jetER_smearer.reset(new GenericJetResolutionSmearer(ctx, "jets", "genjets", true, JERSmearing::SF_13TeV_2016_03Feb2017));
      // 	else throw runtime_error("In AnalysisModule_SiJetTrg.cxx: When setting up JER smearer, invalid 'JEC_Version' was specified.");
      // }
     
    //output
    ctx.undeclare_all_event_output();   
    // //pileup (define it after undeclaring all other variables to keep the weights in the output)
    // pileupSF.reset(new MCPileupReweight(ctx));

    //    tt_dijet_event = ctx.declare_event_output<dijet_event>("dijet");
    //Store only vars needed for the dijet analysis
    tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
    tt_gen_PUpthat = ctx.declare_event_output<float>("gen_PUpthat");
    tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
    tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
    tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
    tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
    tt_jet1_ptRaw = ctx.declare_event_output<float>("jet1_ptRaw");
    tt_jet2_ptRaw = ctx.declare_event_output<float>("jet2_ptRaw");
    tt_jet3_ptRaw = ctx.declare_event_output<float>("jet3_ptRaw");
    tt_jet1_ptGen = ctx.declare_event_output<float>("jet1_ptGen");
    tt_jet2_ptGen = ctx.declare_event_output<float>("jet2_ptGen");
    tt_jet3_ptGen = ctx.declare_event_output<float>("jet3_ptGen");

    tt_jet1_pt_onoff_Resp = ctx.declare_event_output<float>("jet1_pt_onoff_Resp");
    tt_jet2_pt_onoff_Resp = ctx.declare_event_output<float>("jet2_pt_onoff_Resp");

    tt_probejet_neutEmEF = ctx.declare_event_output<float>("probejet_neutEmEF");
    tt_probejet_neutHadEF = ctx.declare_event_output<float>("probejet_neutHadEF");
    tt_probejet_chEmEF = ctx.declare_event_output<float>("probejet_chEmEF");
    tt_probejet_chHadEF = ctx.declare_event_output<float>("probejet_chHadEF");
    tt_probejet_photonEF = ctx.declare_event_output<float>("probejet_photonEF");
    tt_probejet_muonEF = ctx.declare_event_output<float>("probejet_muonEF");

    tt_nvertices = ctx.declare_event_output<int>("nvertices");
    tt_nGoodvertices = ctx.declare_event_output<int>("nGoodvertices");
    tt_probejet_eta = ctx.declare_event_output<float>("probejet_eta");
    tt_probejet_phi = ctx.declare_event_output<float>("probejet_phi");
    tt_probejet_pt = ctx.declare_event_output<float>("probejet_pt");
    tt_probejet_ptgen = ctx.declare_event_output<float>("probejet_ptgen");
    tt_probejet_ptRaw = ctx.declare_event_output<float>("probejet_ptRaw");
    tt_barreljet_eta = ctx.declare_event_output<float>("barreljet_eta");
    tt_barreljet_phi = ctx.declare_event_output<float>("barreljet_phi");
    tt_barreljet_pt = ctx.declare_event_output<float>("barreljet_pt");
    tt_barreljet_ptgen = ctx.declare_event_output<float>("barreljet_ptgen");
    tt_barreljet_ptRaw = ctx.declare_event_output<float>("barreljet_ptRaw");
    tt_pt_ave = ctx.declare_event_output<float>("pt_ave");
    tt_alpha = ctx.declare_event_output<float>("alpha");
    tt_rel_r = ctx.declare_event_output<float>("rel_r");
    tt_mpf_r = ctx.declare_event_output<float>("mpf_r");
    tt_asymmetry = ctx.declare_event_output<float>("asymmetry");
    tt_B = ctx.declare_event_output<float>("B");
    tt_MET = ctx.declare_event_output<float>("MET");
    tt_nPU = ctx.declare_event_output<int>("nPU");
    tt_matchJetId_0 = ctx.declare_event_output<int>("matchJetId_0");
    tt_matchJetId_1 = ctx.declare_event_output<int>("matchJetId_1");   
    tt_ev_weight = ctx.declare_event_output<float>("weight");
    tt_jets_pt= ctx.declare_event_output<float>("sum_jets_pt");
    tt_jet_n= ctx.declare_event_output<int>("Njet");
    tt_rho = ctx.declare_event_output<float>("rho");
    tt_partonFlavor = ctx.declare_event_output<int>("partonFlavor");
    tt_flavorBarreljet = ctx.declare_event_output<int>("flavorBarreljet");
    tt_responseBarreljet = ctx.declare_event_output<float>("responseBarreljet_genp");
    tt_flavorProbejet = ctx.declare_event_output<int>("flavorProbejet");
    tt_responseProbejet = ctx.declare_event_output<float>("responseProbejet_genp");
    tt_flavorLeadingjet = ctx.declare_event_output<int>("flavorLeadingjet");
    tt_flavorSubleadingjet = ctx.declare_event_output<int>("flavorSubleadingjet");
    tt_response_leadingjet = ctx.declare_event_output<float>("leadingjet_response");
    tt_inst_lumi = ctx.declare_event_output<float>("instantaneous_lumi");
    tt_integrated_lumi_in_bin = ctx.declare_event_output<float>("integrated_lumi_in_bin");
    tt_lumibin = ctx.declare_event_output<int>("lumibin");
    tt_Nmuon = ctx.declare_event_output<int>("Nmuon");
    tt_muon_pt = ctx.declare_event_output<float>("muon_pt");
    tt_ele_pt = ctx.declare_event_output<float>("electron_pt");
    tt_Nele = ctx.declare_event_output<int>("Nelectron");
    tt_integrated_lumi = ctx.declare_event_output<float>("integrated_lumi");

    tt_triggerSiMu = ctx.declare_event_output<int>("triggerSiMu");

    tt_minBias = ctx.declare_event_output<int>("minBias");
    tt_trigger40 = ctx.declare_event_output<int>("trigger40");
    tt_trigger60 = ctx.declare_event_output<int>("trigger60");
    tt_trigger80 = ctx.declare_event_output<int>("trigger80");
    tt_trigger140 = ctx.declare_event_output<int>("trigger140");
    tt_trigger200 = ctx.declare_event_output<int>("trigger200");
    tt_trigger260 = ctx.declare_event_output<int>("trigger260");
    tt_trigger320 = ctx.declare_event_output<int>("trigger320");
    tt_trigger400 = ctx.declare_event_output<int>("trigger400");
    tt_trigger450 = ctx.declare_event_output<int>("trigger450");
    tt_trigger500 = ctx.declare_event_output<int>("trigger500");
    
    tt_dR_jet3_barreljet = ctx.declare_event_output<float>("dR_jet3_barreljet");
    tt_dR_jet3_probejet = ctx.declare_event_output<float>("dR_jet3_probejet");

    h_runnr_input.reset(new JECRunnumberHists(ctx,"Runnr_input"));
    h_input.reset(new JECCrossCheckHists(ctx,"CrossCheck_input"));
    h_lumisel.reset(new JECCrossCheckHists(ctx,"CrossCheck_lumisel"));
    h_beforeCleaner.reset(new JECCrossCheckHists(ctx,"CrossCheck_beforeCleaner"));
    h_afterCleaner.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterCleaner"));
    h_2jets.reset(new JECCrossCheckHists(ctx,"CrossCheck_2jets"));
    h_beforeJEC.reset(new JECCrossCheckHists(ctx,"CrossCheck_beforeJEC"));
    h_afterJEC.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterJEC"));
    h_afterJER.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterJER"));
    h_afterMET.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterMET"));
    h_beforeTriggerData.reset(new JECCrossCheckHists(ctx,"CrossCheck_beforeTriggerData"));
    h_afterTriggerData.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterTriggerData"));
    h_beforeFlatFwd.reset(new JECCrossCheckHists(ctx,"CrossCheck_beforeFlatFwd"));
    h_afterFlatFwd.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterFlatFwd"));
    h_afterPtEtaReweight.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterPtEtaReweight"));
    h_afterLumiReweight.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterLumiReweight"));
    h_afterUnflat.reset(new JECCrossCheckHists(ctx,"CrossCheck_afterUnflat"));
    h_afternVts.reset(new JECCrossCheckHists(ctx,"CrossCheck_afternVts"));
    h_nocuts.reset(new JECAnalysisHists(ctx,"NoCuts"));
    h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    h_sel.reset(new JECAnalysisHists(ctx,"Selection"));
    h_final.reset(new JECAnalysisHists(ctx,"Final"));
   
    //Single Jet Trigger Studies
    
    h_trgSiMu.reset(new JECAnalysisHists(ctx,ctx.get("triggerSiMu_DirName", "NULL")));
    
    h_minBias.reset(new JECAnalysisHists(ctx,"HLT_minBias"));
    
    h_trg40.reset(new JECAnalysisHists(ctx,"HLT_PFJet40"));
    h_trg60.reset(new JECAnalysisHists(ctx,"HLT_PFJet60"));
    h_trg80.reset(new JECAnalysisHists(ctx,"HLT_PFJet80"));
    h_trg140.reset(new JECAnalysisHists(ctx,"HLT_PFJet140"));
    h_trg200.reset(new JECAnalysisHists(ctx,"HLT_PFJet200"));
    h_trg260.reset(new JECAnalysisHists(ctx,"HLT_PFJet260"));
    h_trg320.reset(new JECAnalysisHists(ctx,"HLT_PFJet320"));
    h_trg400.reset(new JECAnalysisHists(ctx,"HLT_PFJet400"));
    h_trg450.reset(new JECAnalysisHists(ctx,"HLT_PFJet450"));
    h_trg500.reset(new JECAnalysisHists(ctx,"HLT_PFJet500"));
    
    h_lumi_nocuts.reset(new LuminosityHists(ctx,"Lumi_noCuts"));  
    h_lumi_sel.reset(new LuminosityHists(ctx,"Lumi_Selection"));
    h_lumi_dijet.reset(new LuminosityHists(ctx,"Lumi_diJet"));
    h_lumi_match.reset(new LuminosityHists(ctx,"Lumi_JetMatching"));
    h_lumi_final.reset(new LuminosityHists(ctx,"Lumi_Final"));

    h_lumi_TrigSiMu.reset(new LuminosityHists(ctx,"Lumi_TrigSiMu"));   

    h_lumi_Trig40.reset(new LuminosityHists(ctx,"Lumi_minBias"));      
    h_lumi_Trig40.reset(new LuminosityHists(ctx,"Lumi_Trig40"));  
    h_lumi_Trig60.reset(new LuminosityHists(ctx,"Lumi_Trig60")); 
    h_lumi_Trig80.reset(new LuminosityHists(ctx,"Lumi_Trig80")); 
    h_lumi_Trig140.reset(new LuminosityHists(ctx,"Lumi_Trig140")); 
    h_lumi_Trig200.reset(new LuminosityHists(ctx,"Lumi_Trig200")); 
    h_lumi_Trig260.reset(new LuminosityHists(ctx,"Lumi_Trig260")); 
    h_lumi_Trig320.reset(new LuminosityHists(ctx,"Lumi_Trig320")); 
    h_lumi_Trig400.reset(new LuminosityHists(ctx,"Lumi_Trig400")); 
    h_lumi_Trig450.reset(new LuminosityHists(ctx,"Lumi_Trig450"));
    h_lumi_Trig500.reset(new LuminosityHists(ctx,"Lumi_Trig500"));
     
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));
 
    n_evt = 0;
    TString name_weights = ctx.get("MC_Weights_Path");
    apply_weights = (ctx.get("Apply_Weights") == "true" && isMC);
    //cout<<"Apply Weights: "<<apply_weights<<endl;

    apply_lumiweights = (ctx.get("Apply_Lumiweights") == "true" && isMC);
    apply_unflattening = (ctx.get("Apply_Unflattening") == "true" && isMC);
    if(apply_weights && apply_lumiweights) throw runtime_error("In AnalysisModule_SiJetTrg.cxx: 'apply_weights' and 'apply_lumiweights' are set 'true' simultaneously. This won't work, please decide on one");

    
    string lumifile = ctx.get("lumi_file");
    std::unique_ptr<TFile> file(TFile::Open(lumifile.c_str(), "read"));
    TTree * tree = dynamic_cast<TTree*>(file->Get("AnalysisTree"));
    if(!tree){
      throw runtime_error("LuminosityHists: Did not find TTree 'AnalysisTree' in file ;" + lumifile + "'");
    }
    // only fetch branches we really need:
    brun = tree->GetBranch("run");
    blumiblock = tree->GetBranch("luminosityBlock");
    bilumi = tree->GetBranch("intgRecLumi");

    run_lumi rl;
    double ilumi;
    double intlumi_pb = 0;
    brun->SetAddress(&rl.run);
    blumiblock->SetAddress(&rl.lumiblock);
    bilumi->SetAddress(&ilumi);

    //loop over all lumiblocks to save the map between run/lumiblock and stored lumi of the lumiblock (to be divided by 23s)
    auto ientries = tree->GetEntries();
    for(auto ientry = 0l; ientry < ientries; ientry++){
      for(auto b : {brun, blumiblock, bilumi}){
	b->GetEntry(ientry);
      }
      double ilumi_pb = ilumi * 1e-6; // convert units in file (microbarn) to pb.
      intlumi_pb += ilumi_pb;
      rl2lumi.insert(make_pair(rl, ilumi_pb));
      rl2intlumi.insert(make_pair(rl, intlumi_pb));
    }
    
    double ilumi_current_bin = 0.0;
    run_lumi last_entry;
    for(const auto & rl : rl2lumi){
      ilumi_current_bin += rl.second;
      if(ilumi_current_bin >= 2000){
	upper_binborders_runnrs.push_back(rl.first);
	lumi_in_bins.push_back(ilumi_current_bin - rl.second);
	ilumi_current_bin = ilumi_current_bin - 2000;
      }
      last_entry = rl.first;
    }
    upper_binborders_runnrs.push_back(last_entry); //this is not exactly an UPPER limit because it is equal to the highest possible entry, not greater than it...created exception for this case.
    lumi_in_bins.push_back(ilumi_current_bin);
  };


  AnalysisModule_SiJetTrg::~AnalysisModule_SiJetTrg() {}

  bool AnalysisModule_SiJetTrg::process(Event & event) {
 //###############################################################
    //
    //Selection Module for L2 Residual Calculation  
    //
    //Select Di-Jet Events 
    //Define Barrel and Probe Jet
    //Use possible third Jet to estimate alpha
    //Apply MC-Weights for Reweighting (Second Iteration)
    //
    //###############################################################

    n_evt++;
    if(debug){
      cout << endl << "++++++++++ NEW EVENT +++++++++" << endl << endl;
      cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
    }
    if(debug) cout<<"before first fill\n";   
    //Dump Input
    h_input->fill(event);

    //LEPTON selection

    if(debug) cout<<"before muon cleaning\n";
    lep_cleaner.back()->process(event);
    // muoSR_cleaner->process(event);
    if(debug) cout<<"after muon cleaning\n";
        
    sort_by_pt<Muon>(*event.muons); 
    if(debug) std::cout<<"#muons = "<<event.muons->size()<<std::endl;

    eleSR_cleaner->process(event);
    sort_by_pt<Electron>(*event.electrons);
    if(debug)  std::cout<<"#electrons = "<<event.electrons->size()<<std::endl;

    if (event.electrons->size()>0 || event.muons->size()>0) return false; 

    h_runnr_input->fill(event);

    /* CMS-certified luminosity sections */
    if(event.isRealData){
      if(!lumi_sel->passes(event)){
	if(debug) std::cout<<"Event not in lumi selection!"<<std::endl;
	return false;
      }
      else  h_lumisel->fill(event);
    }
 
    // MET filters   
    //    if(!metfilters_sel->passes(event)) return false;   

    int event_in_lumibin = -1;
    double fill_event_integrated_lumi = 0;
    double inst_lumi = -1;
    double int_lumi_event = -1;
    if(event.isRealData){
      run_lumi rl_event{event.run, event.luminosityBlock};
      double lumiblock_lumi = rl2lumi[rl_event];
      inst_lumi = lumiblock_lumi/23; //FIXME why /23 ???
      int_lumi_event = rl2intlumi[rl_event];

      vector<run_lumi>::iterator it;
      if(!(rl_event < upper_binborders_runnrs.back())){
	if(upper_binborders_runnrs.back() < rl_event) throw runtime_error("AnalysisModule_SiJetTrg: run_lumi of event greater than largest bin-border.");
	else it = prev(upper_binborders_runnrs.end()); //force the entries with the highest run_lumi to enter the last bin instead of overflow.
      }
      else it = upper_bound(upper_binborders_runnrs.begin(), upper_binborders_runnrs.end(), rl_event); //find the first entry in the vector of binborders that is greater than rl_event
      
      event_in_lumibin = distance(upper_binborders_runnrs.begin(), it); //find how many elements of the vector of binborders are smaller than 'it', this is the bin to be filled
      fill_event_integrated_lumi = lumi_in_bins.at(event_in_lumibin);
    }
    
    h_beforeCleaner->fill(event);

//############### Jet Cleaner and First Selection (N_Jets >=2) ##############################
    unsigned int n_jets_beforeCleaner = event.jets->size();
    //JetID
    if(jetLabel == "AK4CHS" || jetLabel == "AK8CHS") jetcleaner->process(event);
    unsigned int n_jets_afterCleaner = event.jets->size();
    //discard events if not all jets fulfill JetID instead of just discarding single jets
    if (n_jets_beforeCleaner != n_jets_afterCleaner) return false;
    // sort_by_pt<Jet>(*event.jets);
    //h_cleaner->fill(event);
    //DEBUG: pt order check!!!
    // float _pt=10000.;
    // for(unsigned int i = 0;i<n_jets_afterCleaner ;i++){
    //   if(_pt<event.jets->at(i).pt()){
    // 	cout<<"!!!! wrong pt ordering !!!!!!";
    // 	throw runtime_error(" wrong pt ordering");
    // 	return false;
    //   }
    //   _pt = event.jets->at(i).pt();
    // }
    // cout<<"event jet was pt ordered\n";
      
    h_afterCleaner->fill(event);

    const int jet_n = event.jets->size();
    if(jet_n<2) return false;
    h_2jets->fill(event);
//###########################################################################################
  
//####################  Select and Apply proper JEC-Versions for every Run ##################

    bool apply_global = true;
    
    h_beforeJEC->fill(event);
 
    if(debug) std::cout <<" before jetleptoncleaner  "<<std::endl;
    if (event.electrons->size()!=0 || event.muons->size()!=0)    jetleptoncleaner->process(event);
    if(debug) std::cout <<" after jetleptoncleaner  "<<std::endl;
    jet_corrector->process(event);
    if(debug) std::cout <<" after jet_corrector  "<<std::endl;
    
    h_afterJEC->fill(event);

//##############################################################################################
    
//################################  Apply JER and MET  #########################################

    //Apply JER to all jet collections
    if(jetER_smearer.get()) jetER_smearer->process(event);
    h_afterJER->fill(event); 

    //correct MET only AFTER smearing the jets
    if(apply_global){
      jet_corrector->correct_met(event);
    }
    h_afterMET->fill(event); 
//##############################################################################################

    int jetid_0 = 0;
    int jetid_1 = 1;
    int jetid_2 = 2;

    int jetid_0_last = 0;
    int jetid_1_last = 1;
    int jetid_2_last = 2;   
  
 //Calculate pt_ave
   Jet* jet1 = &event.jets->at(0);// leading jet
   Jet* jet2 = &event.jets->at(1);// sub-leading jet
   float jet1_pt = jet1->pt();
   float jet2_pt = jet2->pt();
   float pt_ave = (jet1_pt + jet2_pt)/2.;    
    
////###############################################  Trigger  ##################################
 
    int triggerSiMu  = 0;

    int minBias  = 0;   
    int trigger40  = 0;
    int trigger60  = 0;
    int trigger80  = 0;
    int trigger140 = 0;
    int trigger200 = 0;
    int trigger260 = 0;
    int trigger320 = 0;
    int trigger400 = 0;
    int trigger450 = 0;
    int trigger500 = 0;

    int triggerDi40  = 0;
    int triggerDi60  = 0;
    int triggerDi80  = 0;
    int triggerDi140 = 0;
    int triggerDi200 = 0;
    int triggerDi260 = 0;
    int triggerDi320 = 0;
    int triggerDi400 = 0;
    int triggerDi500 = 0;
    
    bool pass_triggerSiMu=false;
    bool pass_minBias=false;   
    bool pass_trigger40=false; bool pass_trigger60=false; bool pass_trigger80=false;
    bool pass_trigger140=false; bool pass_trigger200=false; bool pass_trigger260=false;
    bool pass_trigger320=false; bool pass_trigger400=false; bool pass_trigger450=false;  bool pass_trigger500=false;

    bool pass_triggerDi40=false; bool pass_triggerDi60=false; bool pass_triggerDi80=false;
    bool pass_triggerDi140=false; bool pass_triggerDi200=false; bool pass_triggerDi260=false;
    bool pass_triggerDi320=false; bool pass_triggerDi400=false; bool pass_triggerDi500=false;

    double trg_thresh[10] = {s_Pt_Ave40_cut,s_Pt_Ave60_cut,s_Pt_Ave80_cut,s_Pt_Ave140_cut,s_Pt_Ave200_cut,s_Pt_Ave260_cut,s_Pt_Ave320_cut,s_Pt_Ave400_cut, s_Pt_Ave450_cut, s_Pt_Ave500_cut};

    int n_trig = 0;
    bool pass_trigger;
    
    if(event.isRealData){
        pass_minBias = (minBias_sel->passes(event));

      if(ts){
	pass_trigger40 = (pass_trigger40 && pt_ave>trg_thresh[0]   && pt_ave<trg_thresh[1]);
	pass_trigger60 = (pass_trigger60 && pt_ave>trg_thresh[1]   && pt_ave<trg_thresh[2]);
	pass_trigger80 = (pass_trigger80 && pt_ave>trg_thresh[2]   && pt_ave<trg_thresh[3]); 
	pass_trigger140 = (pass_trigger140 && pt_ave>trg_thresh[3] && pt_ave<trg_thresh[4]);
	pass_trigger200 = (pass_trigger200 && pt_ave>trg_thresh[4] && pt_ave<trg_thresh[5]);
	pass_trigger260 = (pass_trigger260 && pt_ave>trg_thresh[5] && pt_ave<trg_thresh[6]); 
	pass_trigger320 = (pass_trigger320 && pt_ave>trg_thresh[6] && pt_ave<trg_thresh[7]); 
	pass_trigger400 = (pass_trigger400 && pt_ave>trg_thresh[7] && pt_ave<trg_thresh[8]);
	pass_trigger450 = (pass_trigger450 && pt_ave>trg_thresh[8] && pt_ave<trg_thresh[9]);
	pass_trigger500 = (pass_trigger500 && pt_ave>trg_thresh[9]);
      }

      pass_trigger = pass_minBias || pass_trigger40 || pass_trigger60 || pass_trigger80 || pass_trigger140 || pass_trigger200  || pass_trigger260 || pass_trigger320 || pass_trigger400 || pass_trigger450 || pass_trigger500;
    
      if(debug){
	cout << "before triggers: " << endl;
	cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
      }
      h_beforeTriggerData->fill(event);

      if(!pass_trigger){
	return false;
      }
        //Count Events passed Trigger
      if(pass_minBias){ n_trig++; minBias = 1;}
      if(pass_trigger40){ n_trig++; trigger40 = 1;}
      if(pass_trigger60){ n_trig++; trigger60 = 1;}
      if(pass_trigger80){ n_trig++; trigger80 = 1;}
      if(pass_trigger140){ n_trig++; trigger140 = 1;}
      if(pass_trigger200){ n_trig++; trigger200 = 1;}
      if(pass_trigger260){ n_trig++; trigger260 = 1;}
      if(pass_trigger320){ n_trig++; trigger320 = 1;}
      if(pass_trigger400){ n_trig++; trigger400 = 1;}
      if(pass_trigger500){ n_trig++; trigger500 = 1;}
      if(pass_trigger450){ n_trig++; trigger450 = 1;}
      
    }

    h_afterTriggerData->fill(event);

    if(onlyBtB){
     //turn jet2 around and check dR to jet1
     float eta1 = jet1->eta();
     float eta2 = -1.*jet2->eta();
     float phi1 = jet1->phi();
     float phi2 = TVector2::Phi_mpi_pi(jet2->phi()+3.14159265358979f);
     
     float deta = eta1 - eta2;
     float dphi =  TVector2::Phi_mpi_pi(phi1 - phi2);
     float dR = TMath::Sqrt( deta*deta + dphi*dphi );

     if(debug) cout<<"in BtB if with dR "<<dR<<endl;
     
     if(dR>0.4) return false;
   }

//###############################  Declare Probe and Barrel Jet  ###############################
    
    Jet* jet_probe = jet1;
    Jet* jet_barrel = jet2;
    if ((fabs(jet1->eta())<s_eta_barr)&&(fabs(jet2->eta())<s_eta_barr)) {
      int ran = rand();
      int numb = ran % 2;
      if(numb==0){
	jet_probe = jet2;
	jet_barrel = jet1;
      }
    } 
    else if ((fabs(jet1->eta())<s_eta_barr)||(fabs(jet2->eta())<s_eta_barr)){
      if(fabs(jet1->eta())<s_eta_barr){
	jet_probe = jet2;
	jet_barrel = jet1;
      }
      else{
	jet_probe = jet1;
	jet_barrel = jet2;
      }
    }

//##############################################################################################

    //read or calculated values for dijet events
    float gen_pthat = 0; //pt hat (from QCD simulation)
    float gen_weight = 0;
    if(!event.isRealData){
      gen_weight = event.weight;
      gen_pthat = event.genInfo->binningValues()[0];// only for pythia8 samples //todo: for herwig, madgraph
    }
    float nvertices = event.pvs->size(); 
    float nPU = 0 ;//todo for data?
    if(!event.isRealData) nPU = event.genInfo->pileup_TrueNumInteractions();

    float matchJetId_0 = -10.;
    float matchJetId_1 = -10.;
    
    float genjet1_pt = 0;
    float genjet2_pt = 0;
    float genjet3_pt = 0;
    if(isMC){
      if(event.genjets->size()>0)genjet1_pt = event.genjets->at(0).pt();
      if(event.genjets->size()>1)genjet2_pt = event.genjets->at(1).pt();
      if(event.genjets->size()>2)genjet3_pt = event.genjets->at(2).pt();
    }

    auto factor_raw1 = jet1->JEC_factor_raw();
    auto factor_raw2 = jet2->JEC_factor_raw();
    float jet1_ptRaw = jet1_pt*factor_raw1;
    float jet2_ptRaw = jet2_pt*factor_raw2;
    float probejet_eta = jet_probe->eta(); 
    float  probejet_phi = jet_probe->phi(); 
    float  probejet_pt = jet_probe->pt(); 
    auto factor_raw_probe = jet_probe->JEC_factor_raw();
    float probejet_ptRaw = probejet_pt*factor_raw_probe;
    float barreljet_eta = jet_barrel->eta(); 
    float  barreljet_phi = jet_barrel->phi(); 
    float  barreljet_pt = jet_barrel->pt(); 
    auto factor_raw_barrel = jet_barrel->JEC_factor_raw();
    float barreljet_ptRaw = barreljet_pt*factor_raw_barrel;

//##########################  Get third Jet for alpha, asymmetry calculation  ##################
	      
    double dR_jet3_barreljet = -1;
    double dR_jet3_probejet = -1;

   float jet3_pt = 0; float jet3_ptRaw = 0;
    if(jet_n>2){
      Jet* jet3 = &event.jets->at(2);
      jet3_pt = jet3->pt();
      auto factor_raw3 = jet3->JEC_factor_raw();
      jet3_ptRaw = jet3_pt*factor_raw3;
      dR_jet3_barreljet = deltaR(*jet3, *jet_barrel);
      dR_jet3_probejet = deltaR(*jet3, *jet_probe);
    }
    float alpha = jet3_pt/pt_ave;
    float asymmetry = (probejet_pt - barreljet_pt)/(probejet_pt + barreljet_pt);
    float rel_r = (1+asymmetry)/(1-asymmetry);
    TVector2 pt, met;
    met.Set(event.met->pt() * cos(event.met->phi()),event.met->pt() * sin(event.met->phi()));
    pt.Set(barreljet_pt * cos(barreljet_phi),barreljet_pt* sin(barreljet_phi));
    float mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
    float B = (met.Px()*pt.Px() + met.Py()*pt.Py())/((probejet_pt + barreljet_pt) * sqrt(pt.Px()*pt.Px() + pt.Py()*pt.Py())); //vec_MET*vec_ptbarr/(2ptave*ptbarr)

    float jets_pt = 0;
    for(int i=2;i<jet_n;i++){
      jets_pt += ((Jet*)&event.jets->at(i))->pt();
    }
    
    
//##############################################################################################

   h_afterUnflat->fill(event);
    
    int flavor = 0;
    float onoffDummy =0.;
    
 //fill that containers!
    double pu_pthat = -1;
    if(!event.isRealData) pu_pthat = event.genInfo->PU_pT_hat_max();
    event.set(tt_gen_pthat,gen_pthat);
    event.set(tt_gen_PUpthat,pu_pthat);
    event.set(tt_gen_weight,gen_weight);
    event.set(tt_jet1_pt,jet1_pt);
    event.set(tt_jet2_pt,jet2_pt);
    event.set(tt_jet3_pt,jet3_pt);
    event.set(tt_jet1_ptRaw,jet1_ptRaw);
    event.set(tt_jet2_ptRaw,jet2_ptRaw);
    event.set(tt_jet3_ptRaw,jet3_ptRaw);
    event.set(tt_jet1_ptGen,genjet1_pt);
    event.set(tt_jet2_ptGen,genjet2_pt);
    event.set(tt_jet3_ptGen,genjet3_pt);
    event.set(tt_jet1_pt_onoff_Resp,onoffDummy);
    event.set(tt_jet2_pt_onoff_Resp,onoffDummy);
    event.set(tt_nvertices,nvertices);

    event.set(tt_probejet_neutEmEF,jet_probe->neutralEmEnergyFraction());
    event.set(tt_probejet_neutHadEF,jet_probe->neutralHadronEnergyFraction());
    event.set(tt_probejet_chEmEF,jet_probe->chargedEmEnergyFraction());
    event.set(tt_probejet_chHadEF,jet_probe->chargedHadronEnergyFraction());
    event.set(tt_probejet_photonEF,jet_probe->photonEnergyFraction());
    event.set(tt_probejet_muonEF,jet_probe->muonEnergyFraction());

    event.set(tt_probejet_eta,probejet_eta);
    event.set(tt_probejet_phi,probejet_phi);
    event.set(tt_probejet_pt,probejet_pt);
    event.set(tt_probejet_ptRaw,probejet_ptRaw);
    event.set(tt_barreljet_eta,barreljet_eta);
    event.set(tt_barreljet_phi,barreljet_phi);
    event.set(tt_barreljet_pt,barreljet_pt);
    event.set(tt_barreljet_ptRaw,barreljet_ptRaw);
    event.set(tt_pt_ave,pt_ave);
    event.set(tt_alpha,alpha);
    event.set(tt_asymmetry,asymmetry);
    event.set(tt_B,B);
    event.set(tt_MET,event.met->pt());
    event.set(tt_rel_r,rel_r);
    event.set(tt_mpf_r,mpf_r);
    event.set(tt_nPU,nPU);
    event.set(tt_ev_weight,event.weight);
    event.set(tt_jets_pt,jets_pt);
    event.set(tt_jet_n,jet_n);
    event.set(tt_rho,event.rho);    
    event.set(tt_partonFlavor,flavor);
   
    event.set(tt_dR_jet3_barreljet,dR_jet3_barreljet);
    event.set(tt_dR_jet3_probejet,dR_jet3_probejet);

    event.set(tt_inst_lumi,inst_lumi);
    event.set(tt_integrated_lumi_in_bin,fill_event_integrated_lumi);
    event.set(tt_lumibin,event_in_lumibin);
    event.set(tt_integrated_lumi,int_lumi_event);

    event.set(tt_triggerSiMu, triggerSiMu);    
    event.set(tt_minBias, minBias);
    event.set(tt_trigger40, trigger40);
    event.set(tt_trigger60, trigger60);
    event.set(tt_trigger80, trigger80);
    event.set(tt_trigger140, trigger140);
    event.set(tt_trigger200, trigger200);
    event.set(tt_trigger260, trigger260);
    event.set(tt_trigger320, trigger320);
    event.set(tt_trigger400, trigger400);
    event.set(tt_trigger450, trigger450);
    event.set(tt_trigger500, trigger500);

    sel.SetEvent(event);
    
//##################################################   Advanced Selections   ###################
    //good primary vertex
    int nGoodVts = sel.goodPVertex();

    if(debug){
      cout << "debug is on" << endl;
      cout << "before good vertex selection : " << endl;
      cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
    }

    if(nGoodVts<=0) return false;
    if(debug) cout << "after good vertex selection " << endl;
    event.set(tt_nGoodvertices, nGoodVts);
    if(debug) cout << "after good vertex set " << endl;      
    h_afternVts->fill(event);
    if(debug) cout << "after good vertex fill " << endl;
    
    matchJetId_0 = -10;
    matchJetId_1 = -10;
    event.set(tt_matchJetId_0, matchJetId_0);
    event.set(tt_matchJetId_1, matchJetId_1);

    //DiJet-Events      
    if(!sel.DiJet()) return false;
    if(debug) cout << "after diJet sel" << endl;
    h_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

//Pu_pt_hat/pt_hat Selection
    if(!event.isRealData){
      if(!sel.PUpthat(event)) return false;
    }
    h_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    //MET/pt - Cut
    if(apply_METoverPt_cut && event.get(tt_MET)/(event.get(tt_jets_pt)+event.get(tt_barreljet_pt)+event.get(tt_probejet_pt))>0.2) return false; //skip events with large MET contribution  

    //PhiEta Region cleaning
    if(apply_EtaPhi_cut && !sel.EtaPhiCleaning(event)) return false; 

    if(debug){
     cout << "before 'dijet advanced selection' : " << endl;
     cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
    }

//Advanced Selection: DiJet Events
    if(!sel.DiJetAdvanced(event)) return false;   
    h_dijet->fill(event);
    h_lumi_dijet->fill(event);
    h_match->fill(event);
    h_lumi_match->fill(event);

   if(!event.isRealData){
     if(!sel.PtMC(event)) return false; // For MC only one Pt threshold
   }
	      
    //fill single trigger for efficiency plots before the dijet selection
    if(event.isRealData){
    if(debug) cout << "in is real data" << endl;

    if(pass_minBias) {
      h_minBias->fill(event);
      // h_lumi_minBias->fill(event);
    }     
    if(debug) cout << "after min Bias fill" << endl;
                
      if(debug) cout<<"AnalysisModule_SiJetTrg matchJetId0: "<<matchJetId_0<<endl;
      
    float matchJetId_0_last = -10.;
    float matchJetId_1_last = -10.;

    bool passes_Di[9] = {pass_trigger40,pass_trigger60,pass_trigger80,pass_trigger140,pass_trigger200,pass_trigger260,pass_trigger320,pass_trigger400,pass_trigger500};
    unique_ptr<JECAnalysisHists>* h_trgSi[9] = {&h_trg40, &h_trg60, &h_trg80, &h_trg140, &h_trg200,&h_trg260,&h_trg320,&h_trg400,&h_trg500};
    // unique_ptr< LuminosityHists>* h_lumiSi[10] =  {&h_lumi_Trig40, &h_lumi_Trig60, &h_lumi_Trig80, &h_lumi_Trig140, &h_lumi_Trig200, &h_lumi_Trig260, &h_lumi_Trig320, &h_lumi_Trig400,&h_lumi_Trig450, &h_lumi_Trig500};

    for(int i = 0; i<9; i++){
      if(passes_Di[i]){
	matchJetId_0_last = matchJetId_0;
	matchJetId_1_last = matchJetId_1;
	matchJetId_0 = sel.FindMatchingJet(0,trg_vals_Di[i]);
	matchJetId_1 = sel.FindMatchingJet(1,trg_vals_Di[i]);
	event.set(tt_matchJetId_0, matchJetId_0);
	event.set(tt_matchJetId_1, matchJetId_1);
	if(debug) cout<<"AnalysisModule_noPtBinning matchJetId0: "<<matchJetId_0<<endl;
	(*(h_trgSi[i]))->fill(event);
	// (*(h_lumiSi[i]))->fill(event);
	if(matchJetId_0_last!= -10. || matchJetId_1_last!= -10.){
	  if(matchJetId_0 != matchJetId_0_last || matchJetId_1 != matchJetId_1_last){
	    cout<<"new jet id differed for different trg.  jet id 0 was matched to "<<matchJetId_0<<" instead of "<<matchJetId_0_last<<", jet id 1 was matched to "<<matchJetId_1<<" instead of "<<matchJetId_1_last<<endl;
	    if(matchJetId_0_last < matchJetId_0 && matchJetId_0_last >= 0) matchJetId_0 = matchJetId_0_last;
	    if( ( matchJetId_1_last != matchJetId_0 && matchJetId_1 != matchJetId_0 && matchJetId_1_last < matchJetId_1 && matchJetId_1_last >= 0 ) || ( matchJetId_1 == matchJetId_0 ) ) matchJetId_1 = matchJetId_1_last;
	  }
	}
      }
    }

    event.set(tt_matchJetId_0, matchJetId_0);
    event.set(tt_matchJetId_1, matchJetId_1);
    sel.SetEvent(event);

    }//Data
    else{    
      if(debug){
	cout << "before Pt selection (MC only) : " << endl;
	cout << " Evt# "<<event.event<<" Run: "<<event.run<<" " << endl;
      }
    }
    
//###############################################################################################

    if (event.get(tt_alpha) < 0.3) {
      h_sel->fill(event);
      // h_lumi_sel->fill(event);
    }
    h_final->fill(event);
    // h_lumi_final->fill(event);

    event.set(tt_Nmuon,event.muons->size());
    if(event.muons->size()>0)  
      event.set(tt_muon_pt,event.muons->at(0).pt());
    else 
      event.set(tt_muon_pt,0);
    event.set(tt_Nele,event.electrons->size());
    if(event.electrons->size()>0)  
      event.set(tt_ele_pt,event.electrons->at(0).pt());
    else 
      event.set(tt_ele_pt,0);
     
    if(isMC){    
      double flavor_barreljet = 0;
      double response_barreljet = 0;
      double flavor_probejet = 0;
      double response_probejet = 0;
      double flavor_leadingjet = 0;
      double flavor_subleadingjet = 0;
      const unsigned int genjets_n = event.genjets->size();
      int idx_jet_matching_genjet[genjets_n];
      double probejet_ptgen = -1; 
      double barreljet_ptgen = -1; 

      //match genp to gen-jets
      int idx_j=0;
      int idx_genp_min = -1;
      //this array contains one idx for each jet in the event. If -1: unmatched, else: idx of the closest genpart with dR<=0.2
      int idx_matched_jets[jet_n];
      for(int i=0; i<jet_n; i++){
    	idx_matched_jets[i] = -1;
      }

      //matching gen- and reco-jets
      for(unsigned int i=0; i<event.genjets->size(); i++){
    	double dR_min = 99999; int idx_matching_jet = -1;
    	for(unsigned int j=0; j<event.jets->size(); j++){
    	  double dR = deltaR(event.jets->at(j), event.genjets->at(i));
    	  if(debug) cout << "dR between GenJet " << i << " and RecoJet " << j << ": " << dR << endl;
    	  if(dR<dR_min){
    	    dR_min = dR; 
    	    if(dR_min<0.1) idx_matching_jet = j;
    	  }
    	}
    	idx_jet_matching_genjet[i] = idx_matching_jet;
    	if(debug) cout << "the jet matching the genjet no. " << i << " is jet no. " << idx_matching_jet << endl;
      }
      /////////////////////

      for(Particle & genj : *event.genjets){
    	double dr_min = 99999;
    	double dr_cut = 0;
    	if(jetLabel == "AK4CHS" || jetLabel == "AK4PUPPI") dr_cut = 0.2;
    	else if (jetLabel == "AK8CHS" || jetLabel == "AK8PUPPI")dr_cut = 0.4;
    	else throw runtime_error("AnalysisModule_SiJetTrg.cxx: Invalid jet-label specified.");

    	int idx_g = 0;
    	for(GenParticle & genp: *event.genparticles){
    	  double dr = deltaR(genj,genp);
    	  if(dr < dr_min){
    	    dr_min = dr;
    	    idx_genp_min = idx_g;	
    	  }	
    	  if(debug) cout << "dr between genjet " << idx_j << " and genp (flavor: " << genp.flavor() << ") " << idx_g << "= " << dr << endl;
    	  idx_g++;
    	}
    	if(dr_min <= dr_cut) {
    	  if(debug) cout << "genjet " << idx_j << " is matched to genparticle " << idx_genp_min << " of flavor " << event.genparticles->at(idx_genp_min).flavor() << " within dR = " << dr_min << ". " <<  endl; 
    	  if(idx_jet_matching_genjet[idx_j] >= 0) idx_matched_jets[idx_jet_matching_genjet[idx_j]] = idx_genp_min;
    	}
    	idx_j++;
      }
      //only consider jets that could be matched to a genparticle, these shall take the partons flavor by definition

      if(debug && event.genjets->size() <2) cout << "WARNING: GENjets size < 2" << endl;

      //only consider the barreljet, is it leading or sub-leading jet?
      int idx_barreljet = -1;
      if(fabs(jet1->pt() - jet_barrel->pt()) < 0.001) idx_barreljet = 0;
      else if (fabs(jet2->pt() - jet_barrel->pt()) < 0.001) idx_barreljet = 1;
      else throw runtime_error("first two jets are not the barrel jets, how could this happen?");
    
      //obtain flavor of the barreljet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[idx_barreljet] != -1){	
    	flavor_barreljet = fabs(event.genparticles->at(idx_matched_jets[idx_barreljet]).flavor());
    	response_barreljet = jet_barrel->pt() / event.genparticles->at(idx_matched_jets[idx_barreljet]).pt();
        barreljet_ptgen = event.genparticles->at(idx_matched_jets[idx_barreljet]).pt(); 
      }
      else{
    	flavor_barreljet = -1;
    	response_barreljet = -1;
        barreljet_ptgen = -1;
      }
      if(debug) cout << "barreljet is jet no. " << idx_barreljet << ", alpha = " << event.get(tt_alpha) << ", flavor of barreljet = " << flavor_barreljet << endl;
    

      //also for probe jets
      int idx_probejet = fabs(idx_barreljet - 1);
      //obtain flavor of the probejet
      //-1: unmatched,  >0: flavor of matching genparticle 
      if(idx_matched_jets[idx_probejet] != -1){
    	flavor_probejet = fabs(event.genparticles->at(idx_matched_jets[idx_probejet]).flavor());
    	response_probejet = jet_probe->pt() / event.genparticles->at(idx_matched_jets[idx_probejet]).pt();
    	probejet_ptgen = event.genparticles->at(idx_matched_jets[idx_probejet]).pt();
      }
      else{ 
    	flavor_probejet = -1;
    	response_probejet = -1;
        probejet_ptgen = -1;
      }
      if(debug) cout << "probejet is jet no. " << idx_probejet << ", alpha = " << event.get(tt_alpha) << ", flavor of probejet = " << flavor_probejet << endl;
      
      
      //same for leading jet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[0] != -1) flavor_leadingjet = fabs(event.genparticles->at(idx_matched_jets[0]).flavor());
      else flavor_leadingjet = -1;
      if(debug) cout << "leadingjet is jet no. " << 0 << ", alpha = " << event.get(tt_alpha) << ", flavor of leadingjet = " << flavor_leadingjet << endl;
      

      //same for subleading jet
      //-1: unmatched, 0: alpha too large, >0: flavor of matching genparticle 
      if(idx_matched_jets[1] != -1) flavor_subleadingjet = fabs(event.genparticles->at(idx_matched_jets[1]).flavor());
      else flavor_subleadingjet = -1;
      if(debug) cout << "subleadingjet is jet no. " << 1 << ", alpha = " << event.get(tt_alpha) << ", flavor of subleadingjet = " << flavor_subleadingjet << endl;

      event.set(tt_flavorBarreljet,flavor_barreljet);   
      event.set(tt_responseBarreljet,response_barreljet);
      event.set(tt_barreljet_ptgen,barreljet_ptgen); 
      event.set(tt_flavorProbejet,flavor_probejet);  
      event.set(tt_responseProbejet,response_probejet);   
      event.set(tt_flavorLeadingjet,flavor_leadingjet);  
      event.set(tt_flavorSubleadingjet,flavor_subleadingjet);
      event.set(tt_probejet_ptgen,probejet_ptgen);   

      //response of leading jet
      //find corresponding genjet
      int idx_corresponding_genjet = -1;
      for(unsigned int i=0; i<event.genjets->size(); i++){
    	if(debug) cout << idx_jet_matching_genjet[i] << endl;
    	if(idx_jet_matching_genjet[i] == 0) idx_corresponding_genjet = i;
      }
      double response_jet1 = -1;
      if(idx_corresponding_genjet != -1) response_jet1 = event.jets->at(0).pt() / event.genjets->at(idx_corresponding_genjet).pt();
      event.set(tt_response_leadingjet,response_jet1);  

    } //isMC
    else{
      event.set(tt_flavorBarreljet,-1);   
      event.set(tt_responseBarreljet,-1);   
      event.set(tt_flavorProbejet,-1);  
      event.set(tt_responseProbejet,-1);  
      event.set(tt_flavorLeadingjet,-1);  
      event.set(tt_flavorSubleadingjet,-1);  
      event.set(tt_response_leadingjet,-1.); 
      event.set(tt_probejet_ptgen,-1.);  
      event.set(tt_barreljet_ptgen,-1);   
    }

   return true;
  }
	     
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(AnalysisModule_SiJetTrg)
// XD XD XD XD XD XD XD XD XD
// // https://xkcd.com/530/
    
