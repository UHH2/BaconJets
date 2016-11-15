#pragma once

#include <TString.h>

using namespace std;


const double alpha_cut = 0.3;



const int n_alpha = 9;
const int n_eta = 19;
const int n_pt = 9;
const int nResponseBins = 100;
const int n_etabarr=5; // needed for the normalization to 1 in the barrel

const TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045"};
const TString pt_range[n_pt]= {"56.000", "78.000", "100.000", "168.000", "232.000", "300.000", "366.000", "453.000", "562.000"};
const TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
const TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};

const double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450};
const double pt_bins[n_pt]       = {56, 78,  100, 168, 232, 300, 366, 453, 562};
const double eta_bins[n_eta]     = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};



