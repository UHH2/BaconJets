//#include <iostream>

#include "TClonesArray.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"


#include "UHH2/bacondataformats/interface/TEventInfo.h"

namespace uhh2bacon {

class Selection {

    private:
    uhh2::Context& context;
    uhh2::Event* event;

    uhh2::Event::Handle<TClonesArray> h_jets;
    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;
    uhh2::Event::Handle<TClonesArray> h_pv;
   // float csv_threshold;
    public:
    Selection(uhh2::Context & ctx);
    ~Selection();

    void SetEvent(uhh2::Event& evt);
    bool Trigger();
    bool DiJet();
    bool DiJetAdvanced();
    bool AlphaCut();
    bool goodPVertex();
    bool triggerFired(float bin1, float bin2);

//     TString WP_LOOSE;
//     TString WP_MEDIUM;
//     TString WP_TIGHT;
//     explicit jetIds(wp working_point);
    bool jetIds(float csv_threshold);

    bool FullSelection();


};

}
