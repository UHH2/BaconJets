#pragma once

#include <TString.h>

using namespace std;

const double eta_cut = 2.8; //Central/FWD triggers
// const double eta_cut = 100; //Central triggers only

const int nResponseBins = 100;// 100
const int n_etabarr=5; // needed for the normalization to 1 in the barrel

//Test: Energy Fractions of Jets
const double neutEMEF_threshold = 0.025;
const double chHadEF_threshold  = 0.025;

//Alpha: 
const double alpha_cut = 0.3;
const TString s_alpha_cut = "0.3";


const int n_alpha = 9;
const int n_alpha_cut = 5;

const double al_cut[n_alpha_cut]={0.3, 0.25, 0.2, 0.15, 0.1};
const TString al_cut_string[n_alpha_cut]={"03","025", "02", "015", "01"};

const TString alpha_range[n_alpha] = {"a005", "a010", "a015", "a020", "a025", "a030", "a035", "a040", "a045"};
const double alpha_bins[n_alpha] = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350,  0.400, 0.450};

const int n_trigger = 9;
const int triggerVal[n_trigger] = {40, 60, 80, 140, 200, 260, 320, 400, 500};

const int n_trigger_HF = 6;
const int triggerVal_HF[n_trigger_HF] = {60, 80, 100, 160, 220, 300};


//Pt bins:
const int n_pt = 11;
const TString pt_range[n_pt]= {
  "51",
  "74",
  "96",
  "165",
  "232",
  "300",
  "366",
  "456",
  "569",
  "1000",
  "2000"};

const double pt_bins[n_pt] = {
  51,
  74, 
  96,
  165,
  232,
  300,
  366,
  456,
  569,
  1000,
  2000};
/*
const int n_pt = 11;
const TString pt_range[n_pt]= {
  "51",
  "74",
  "96",
  "163",
  "230",
  "299",
  "365",
  "453",
  "566",
  "1000",
  "2000"};

const double pt_bins[n_pt] = {
  51,
  74, 
  96,
  163,
  230,
  299,
  365,
  453,
  566,
  1000,
  2000};
*/
const int n_pt_HF = 10;
const TString pt_range_HF[n_pt_HF]={
  "51",
  "72",
  "95", 
  "118",
  "188",
  "257", 
  "354",
  "450",
  "1000",
  "2000"};

const double pt_bins_HF[n_pt_HF]{
  51,
  72,
  95,
  118,
  188,
  257,
  354,
  450,
  1000,
  2000};


//Eta bins:
//Abs eta range:
const int n_eta = 19;

const TString eta_range[n_eta] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839","5.191"};
const TString eta_range2[n_eta] = {"00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139","3489", "3839", "5191"}; 
const double eta_bins[n_eta]     = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139,3.489, 3.839, 5.191};


//Negative Eta Range for Final_ControlPlots_CorrectFormulae.cc

const int n_eta_control = 37;

const TString eta_range_control[n_eta_control] = {"-5.191","-3.839","-3.489","-3.139","-2.964","-2.853", "-2.65", "-2.5", "-2.322", "-2.172", "-1.93", "-1.653", "-1.479", "-1.305", "-1.044", "-0.783", "-0.522", "-0.261"," 0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};
const TString eta_range_control2[n_eta_control] = {"-5191","-3839","-3489","-3139","-2964","-2853", "-265", "-25", "-2322", "-2172", "-193", "-1653", "-1479", "-1305", "-1044", "-0783", "-0522", "-0261","00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"};
const double eta_bins_control[n_eta_control]     = {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};



const double eta_bins2[19]     = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

const int n_alpha_common = 4;
const TString alpha_range_common[n_alpha_common] = {"a10","a15", "a20", "a30"};//for Global fit we produce result only in few alpha values
const double alpha_bins_common[n_alpha_common] = {0.100, 0.150, 0.200, 0.300};

const int n_eta_common = 19;
const double eta_common_bins[n_eta_common] ={0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
const TString eta_common_range[n_eta_common] = {"0.000", "0.261", "0.522", "0.783", "1.044", "1.305", "1.479", "1.653", "1.930", "2.172", "2.322", "2.500", "2.650", "2.853", "2.964", "3.139", "3.489", "3.839", "5.191"};

const TString eta_output[n_eta_common-1] = {"eta_00_03", "eta_03_05","eta_05_08","eta_08_10","eta_10_13","eta_13_15","eta_15_17", "eta_17_19", "eta_19_22", "eta_22_23", "eta_23_25", "eta_25_27", "eta_27_29", "eta_29_30", "eta_30_31", "eta_31_35", "eta_35_38", "eta_38_52"};  //Anastasia comparison macro


// ************************************************************************************
//Wide binning
const int n_eta_common_2 = 7;
const double eta_common_bins_2[n_eta_common] ={0, 1.305, 1.93, 2.5, 2.964, 3.139, 5.191};
const TString eta_common_range_2[n_eta_common] = {"0.000", "1.305",  "1.930",  "2.500", "2.964", "3.139", "5.191"};

const TString eta_output_2[n_eta_common-1] = {"eta_00_13", "eta_13_19", "eta_19_25", "eta_25_30", "eta_30_32", "eta_32_52"};  


// ************************************************************************************
