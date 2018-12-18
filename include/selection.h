//#include <iostream>

#include "TClonesArray.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
//#include "UHH2/BaconTrans/baconheaders/TEventInfo.hh"

namespace uhh2bacon {

class Selection {

    private:
    uhh2::Context& context;
    uhh2::Event* event;

    /* uhh2::Event::Handle<TClonesArray> h_jets; */
    /* uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo; */
    /* uhh2::Event::Handle<TClonesArray> h_pv; */

    //Additional vars in Event, specific for dijet
    uhh2::Event::Handle<float> tt_gen_pthat; uhh2::Event::Handle<float> tt_gen_weight;
    uhh2::Event::Handle<float> tt_jet1_pt;     uhh2::Event::Handle<float> tt_jet2_pt;     uhh2::Event::Handle<float> tt_jet3_pt;
    uhh2::Event::Handle<float> tt_jet1_ptRaw;  uhh2::Event::Handle<float> tt_jet2_ptRaw;  uhh2::Event::Handle<float> tt_jet3_ptRaw;
    uhh2::Event::Handle<int> tt_nvertices;
    uhh2::Event::Handle<float> tt_probejet_eta;  uhh2::Event::Handle<float> tt_probejet_phi; uhh2::Event::Handle<float> tt_probejet_pt; uhh2::Event::Handle<float> tt_probejet_ptRaw;
    uhh2::Event::Handle<float> tt_barreljet_eta;  uhh2::Event::Handle<float> tt_barreljet_phi; uhh2::Event::Handle<float> tt_barreljet_pt; uhh2::Event::Handle<float> tt_barreljet_ptRaw;
    uhh2::Event::Handle<float> tt_pt_ave;
    uhh2::Event::Handle<float> tt_alpha;
    uhh2::Event::Handle<float> tt_rel_r; uhh2::Event::Handle<float> tt_mpf_r; uhh2::Event::Handle<float> tt_asymmetry; uhh2::Event::Handle<int> tt_nPU;

    TString Cut_Dir;
    TString dataset_version;

    TFile* cut_map;
    TH2D* h_map;

    public:
    Selection(uhh2::Context & ctx);
    ~Selection();

    void SetEvent(uhh2::Event& evt);
    // bool Trigger(uhh2::Event& evt);
    bool PtMC(uhh2::Event& evt); //apply lowest Pt cut on MC
    bool DiJet();
    bool DiJetAdvanced(uhh2::Event& evt);
    bool ThirdJetSelection(uhh2::Event& evt);
    int goodPVertex();
    bool triggerFired(float bin1, float bin2);

    //  bool FullSelection();
    //Cut events with pthat in PU higher than pthat
    bool PUpthat(uhh2::Event& evt);
    bool EtaPhi_HCAL(uhh2::Event& evt);
    bool EtaPhiCleaning(uhh2::Event& evt);
    bool EnergyEtaCut();

};

}
