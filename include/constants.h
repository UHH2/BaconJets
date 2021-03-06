#ifndef  CONSTANTS_H
#define  CONSTANTS_H

/** \brief Binning **/
static std::vector<double>   eta_range  = {0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.853, 2.964, 3.139, 5.232};
static std::vector<double>   pt_range   = {66, 107, 191, 240, 306, 379, 468, 900};
static std::vector<double>   alpha_range= {0.00, 0.10, 0.20, 0.30, 0.40};

/** \brief Dijet event selection **/
// barrel region (|eta| < 1.3)
static float s_eta_barr = 1.3;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.7)
static float s_delta_phi = 2.7;
// cut on the asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.70
static float s_asymm = 0.7;
// relative third jet fraction pt_rel = 2*j3_pt/(j1_pt + j2_pt) < 0.2
static float s_pt_rel = 0.2;

/** \brief good Primary Vertex reconstruction **/
// more than four tracks
static float s_n_PvTracks = 4;
// PV is located within 24cm in z vertex
static float s_n_Pv_z = 24.0;
// PV is located within 2cm in xy direction from the nominal interaction point
 static float s_n_Pv_xy = 2.0;


/** \brief The CSV V2 b-tag as jet id **/
/// Uses the preliminary thresholds for that tagger recommended here:
///  https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
/// 
/// Note that this will certainly need updates once 13TeV recommendations are available.
static float s_WP_LOOSE = 0.423f;
static float s_WP_MEDIUM = 0.814f;
static float s_WP_TIGHT = 0.941f;

static float s_working_point_csv_threshold = s_WP_MEDIUM;

/** \brief The trigger thresholds of pt_ave **/
/// Used for the 2012 data
static float s_Pt_Ave40_cut   = 66; //Hennings: 62
static float s_Pt_Ave80_cut   = 107;//Hennings: 107
static float s_Pt_Ave140_cut  = 191;//Hennings: 175
static float s_Pt_Ave200_cut  = 240;//Hennings: 242
static float s_Pt_Ave260_cut  = 306;//Hennings: 310
static float s_Pt_Ave320_cut  = 379;//Hennings: 379
static float s_Pt_Ave400_cut  = 468;//Hennings: 467


/** \brief Jet Resolution Smearering **/
// doing the matching from GEN to RECO
static float s_delta_R   = 0.3; 
//constant numberstaken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
// from 8TeV JER measurement.
constexpr const size_t n = 7;
static float eta_hi[n]    = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
static float c_nominal[n] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
static float c_up[n]      = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
static float c_down[n]    = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};


/** \brief Dijet event weighting **/
// in each pt-region defined by the trigger thresholds # of events in MC and in data should equal
static float scale_factor1[] = {2.14563, 4.27422,7.57198,  10.7035,17.899, 15.6095, 4.64629};
static float scale_factor2[] = {242.617, 3514.77, 49549.8, 180540, 693177, 1.83449e+06, 1.63241e+07};



#endif
