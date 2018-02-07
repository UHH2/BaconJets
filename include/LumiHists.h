#include "UHH2/core/include/Hists.h"
#include <TString.h>

#pragma once

using namespace std;

/** \brief Create the "lumi plot", i.e. event yield vs. time in bins of equal integrated luminosity
 * 
 * Configuration from context:
 *  - "lumi_file": path to root file with luminosity information
 *  - "lumihists_lumi_per_bin": integrated luminosity per bin in the histogram (optional; default: 50.0)
 */


class LumiHists: public uhh2::Hists {
public:
    LumiHists(uhh2::Context & ctx,
                    const std::string & dirname,
                    const std::string & triggername = "",
                    bool do_inst_lumi_hist = false);

    virtual void fill(const uhh2::Event & ev) override;
    
private:
    
    // save the upper bin borders of those run/lumi numbers to
    // still include in the bin. Has size = nbins - 1, where
    // nbins is the number of bins in the lumi histogram
    std::vector<run_lumi> upper_binborders;

    std::map<run_lumi, double> rl2lumi;

    int n_pt = 12;
    // const int n_eta = 19;
    int n_eta = 37;

    TH1D * hlumi;
    TH2D * hAsymLumi[37][12];
    TH2D * hBsymLumi[37][12];
    TProfile * pr_AsymLumi[37][12];
    TH1D * hinstlumi;
    TH1D * hinstlumi_ref;
   double lumi_per_bin;
   const std::string triggername_;

 const TString pt_range[12]= {"51", "73", "95", "129", "163", "230", "299", "365", "453", "566", "1000", "2000"};
 const double pt_bins[12]={51, 73, 95, 129, 163, 230, 299, 365, 453, 566, 1000, 2000};

 /*
   TString eta_range[37] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};
   double eta_bins[37]    = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
 */


const TString eta_range[37] = {"-5191","-3839","-3489","-3139","-2964","-2853", "-265", "-25", "-2322", "-2172", "-193", "-1653", "-1479", "-1305", "-1044", "-0783", "-0522", "-0261","00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};
const double eta_bins[37]     = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

   bool do_inst_lumi_hist_;
};
