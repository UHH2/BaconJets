#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TClonesArray.h"
#include "UHH2/bacondataformats/interface/TEventInfo.h"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.h"
using namespace baconhep;
/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class JECAnalysisHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    JECAnalysisHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    void fill(const uhh2::Event & ev, const int rand);
    virtual ~JECAnalysisHists();

private:
    uhh2::Event::Handle<TClonesArray> h_jets;
    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;
    uhh2::Event::Handle<TClonesArray> h_pv;
};
