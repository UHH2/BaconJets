#include "UHH2/BaconJets/include/selection.h"

#include <iostream>
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/PrimaryVertex.h"
// #include "UHH2/BaconTrans/baconheaders/TJet.hh"
//#include "UHH2/BaconTrans/baconheaders/TVertex.hh"
#include "UHH2/BaconJets/include/constants.h"

#include "TVector2.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace std;
namespace uhh2bacon {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  // auto jetCollection = ctx.get("jetCollection");
  // h_jets = ctx.declare_event_input<TClonesArray>(jetCollection);
  // //    h_jets = context.declare_event_input<TClonesArray>("AK4PFCHS");
  // // h_jets = context.declare_event_input<TClonesArray>("AK4PFPUPPI");
  // h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
  // h_pv = context.declare_event_input<TClonesArray>("PV");

  tt_gen_pthat = ctx.declare_event_output<float>("gen_pthat");
  tt_gen_weight = ctx.declare_event_output<float>("gen_weight");
  tt_jet1_pt = ctx.declare_event_output<float>("jet1_pt");
  tt_jet2_pt = ctx.declare_event_output<float>("jet2_pt");
  tt_jet3_pt = ctx.declare_event_output<float>("jet3_pt");
  tt_jet1_ptRaw = ctx.declare_event_output<float>("jet1_ptRaw");
  tt_jet2_ptRaw = ctx.declare_event_output<float>("jet2_ptRaw");
  tt_jet3_ptRaw = ctx.declare_event_output<float>("jet3_ptRaw");
  tt_nvertices = ctx.declare_event_output<int>("nvertices");
  tt_probejet_eta = ctx.declare_event_output<float>("probejet_eta");
  tt_probejet_phi = ctx.declare_event_output<float>("probejet_phi");
  tt_probejet_pt = ctx.declare_event_output<float>("probejet_pt");
  tt_probejet_ptRaw = ctx.declare_event_output<float>("probejet_ptRaw");
  tt_barreljet_eta = ctx.declare_event_output<float>("barreljet_eta");
  tt_barreljet_phi = ctx.declare_event_output<float>("barreljet_phi");
  tt_barreljet_pt = ctx.declare_event_output<float>("barreljet_pt");
  tt_barreljet_ptRaw = ctx.declare_event_output<float>("barreljet_ptRaw");
  tt_pt_ave = ctx.declare_event_output<float>("pt_ave");
  tt_alpha = ctx.declare_event_output<float>("alpha");
  tt_rel_r = ctx.declare_event_output<float>("rel_r");
  tt_mpf_r = ctx.declare_event_output<float>("mpf_r");
  tt_asymmetry = ctx.declare_event_output<float>("asymmetry");
  tt_nPU = ctx.declare_event_output<int>("nPU");

  Cut_Dir = ctx.get("Cut_dir");
  dataset_version = ctx.get("dataset_version");
  bool isMC = (ctx.get("dataset_type") == "MC");

  if(!isMC){
    if(dataset_version.Contains("RunB")){
      cut_map = new TFile(Cut_Dir+"hotjets-runB.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunC")){
      cut_map = new TFile(Cut_Dir+"hotjets-runC.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunD")){
      cut_map = new TFile(Cut_Dir+"hotjets-runD.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunE")){
      cut_map = new TFile(Cut_Dir+"hotjets-runE.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunFe")){
      cut_map = new TFile(Cut_Dir+"hotjets-runFe.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunFl")){
      cut_map = new TFile(Cut_Dir+"hotjets-runFl.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunG")){
      cut_map = new TFile(Cut_Dir+"hotjets-runG.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else if(dataset_version.Contains("RunH")){
      cut_map = new TFile(Cut_Dir+"hotjets-runH.root","READ");
      h_map = (TH2D*) cut_map->Get("h2hotfilter");
      h_map->SetDirectory(0);
      cut_map->Close();
    }
    else{
      throw std::runtime_error("In File selection.cxx: No cleaning map for selected Run!");
    }
  }
 


}

void Selection::SetEvent(uhh2::Event& evt)
{
   event = &evt;
   assert(event);
}

bool Selection::PtMC(uhh2::Event& evt)
{
  assert(event);
  //  std::cout<<"evt.get(tt_pt_ave) = "<<evt.get(tt_pt_ave)<<" s_Pt_Ave40_cut = "<<s_Pt_Ave40_cut<<std::endl;
  if (evt.get(tt_pt_ave) < s_Pt_AveMC_cut) 
    return false;
  return true;
}


bool Selection::DiJet()
{
    assert(event);
    const int njets = event->jets->size();
    if (njets>=2) return true;


   //  const TClonesArray & js = event->get(h_jets);
//     const baconhep::TEventInfo & info = event->get(h_eventInfo);
//     baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
//     assert(eventInfo);

//     int njets = js.GetEntries();
// //     std::cout << "hallo"<<std::endl;

//     // njets >= 2
//     if (njets>=2) return true;

    return false;
}
bool Selection::DiJetAdvanced(uhh2::Event& evt)
{
    assert(event);


  //   const TClonesArray & js = event->get(h_jets);
//     const baconhep::TEventInfo & info = event->get(h_eventInfo);
// //   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
// //   assert(jet);

//     baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
//     assert(eventInfo);
// //   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
// //   assert(jet);

//     int njets = js.GetEntries();

    const int njets = event->jets->size();
    if (njets < 2) return false;
    Jet* jet1 = &event->jets->at(0);// leading jet
    Jet* jet2 = &event->jets->at(1);// sub-leading jet

    // at least one barrel jet
    if((fabs(jet1->eta()) >= s_eta_barr) && (fabs(jet2->eta()) >= s_eta_barr)) return false; 

    // delta phi > 2.7
    double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jet1->phi() - jet2->phi()));
    if (deltaPhi < s_delta_phi) return false;

    // |asymm| < 0.7
    if (fabs((event->get(tt_jet2_pt) - event->get(tt_jet1_pt)) / (event->get(tt_jet2_pt) + event->get(tt_jet1_pt))) > s_asymm) return false;

    //(pTgen1 < 1.5*pThat || pTreco1 < 1.5* pTgen1)
    if(!event->isRealData){
      if(event->genjets->size() < 1) return false;
      if(!(event->genjets->at(0).pt() < 1.5*event->genInfo->binningValues()[0] || event->jets->at(0).pt() < 1.5*event->genjets->at(0).pt())) return false;
    }
    return true;
}
  bool Selection::ThirdJetSelection(uhh2::Event& evt)
  {
    assert(event);

    const int njets = event->jets->size();
    if(njets < 2) return false;
    
    if(njets > 2) {
      Jet* jet3 = &event->jets->at(2);
      // cout<<"jet 3 eta: "<<fabs(jet3->eta())<<"  < 1.3 (Barrel)"<<endl;
      if(fabs(jet3->eta()) > 1.3) return false;
    }
    
    return true;
  }


  int Selection::goodPVertex()
  {
    assert(event);
    Int_t nvertices = event->pvs->size();
    // require in the event that there is at least one reconstructed vertex
    if(nvertices<=0) return 0;//false;
    float nPrVer = 0;
    // pick the first (i.e. highest sum pt) verte
    for (int i=0; i<nvertices; i++){
      PrimaryVertex* vertices = &event->pvs->at(i);
      // require that the vertex meets certain criteria
      //if(vertices->nTracksFit)
      //std::cout<<" vertices->nTracksFit = "<<vertices->nTracksFit<<" "<<fabs(vertices->z)<<" cut at "<<s_n_PvTracks<<std::endl;
      //        if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
      //SHOULD BE LIKE
      //cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.")
      // std::cout<<fabs(vertices->z())<<" "<<fabs(vertices->rho())<<" "<<vertices->ndof()<<" "<<vertices->chi2()<<" "<<std::endl;
      if((fabs(vertices->z()) < 24.) && (fabs(vertices->rho()) < 2.) && (vertices->ndof() >= 4) 
	 && (vertices->chi2()) > 0){
	nPrVer++;
      }
    }
    return nPrVer;
    //    std::cout<<" nPrVer = "<<nPrVer<<" all vtxs = "<<nvertices<<std::endl;
    //    event->set(tt_nGoodvertices,nPrVer); 
    // goodVtx = nPrVer; 
    // if(nPrVer<=0) return false;
    // else
    //   return true;
  }

// bool Selection::goodPVertex()
// {
//     assert(event);

//     const TClonesArray & pvs = event->get(h_pv);

//     Int_t nvertices = pvs.GetEntries();
//     // require in the event that there is at least one reconstructed vertex
//     if(nvertices<=0) return false;
//     float nPrVer = 0;
//     // pick the first (i.e. highest sum pt) verte
//     for (int i=0; i<nvertices; i++){
//         baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
//         // require that the vertex meets certain criteria

//         if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
//             nPrVer++;
//         }
//     }

//  return true;
// }



// bool Selection::FullSelection()
// {
//     return DiJet()&&goodPVertex();

// }


  bool Selection::PUpthat(uhh2::Event& evt)
  {
    assert(event);
  
   double  pt_hat = event->genInfo->binningValues()[0];
   double  PU_pt_hat = event->genInfo->PU_pT_hat_max();
  
   double Ratio = PU_pt_hat/pt_hat;

    if(Ratio < 1) return true;

    return false;
  }

  bool Selection::EtaPhi_HCAL(uhh2::Event& evt)
  {
    assert(event);
    
    double EtaPhi_B[4]={-2.25, -1.93, 2.2, 2.5};
    double EtaPhi_C[4]={-3.489, -3.139, 2.237, 2.475};
    double EtaPhi_D[4]={-3.60, -3.139, 2.237, 2.475};
    
    const int njets = event->jets->size();

    //Needed if only some jets should be checked 
    /*
    int njet_cut = 3;
    if(njets<4) njet_cut = njets;
    */
    //

    for(int i=0; i < njets; i++){
      
      Jet* jet = &event->jets->at(i);// loop over all/some jets in event
      
      double jet_eta = jet->eta();
      double jet_phi = jet->phi();
      
      if(event->run < s_runB){
        if(jet_eta > EtaPhi_B[0] && jet_eta < EtaPhi_B[1] && jet_phi > EtaPhi_B[2] && jet_phi < EtaPhi_B[3]){
	  return false;
	}
      }
      else if(event->run > s_runB && event->run < s_runC){
	if(jet_eta > EtaPhi_C[0] && jet_eta < EtaPhi_C[1] && jet_phi > EtaPhi_C[2] && jet_phi < EtaPhi_C[3]){
	  return false;
	}
      }
      else if(event->run > s_runC && event->run < s_runD+1){
	if(jet_eta > EtaPhi_D[0] && jet_eta < EtaPhi_D[1] && jet_phi > EtaPhi_D[2] && jet_phi < EtaPhi_D[3]){
	  return false;
	}
      }
    }
    return true;
  }
  
  

  bool Selection::EtaPhiCleaning(uhh2::Event& evt)
  {
    assert(event);

    int n_bins_x = h_map->GetNbinsX();
    int n_bins_y = h_map->GetNbinsY();


    double xMin = h_map->GetXaxis()->GetXmin();
    double xWidth = h_map->GetXaxis()->GetBinWidth(1);

 
    double yMin = h_map->GetYaxis()->GetXmin();
    double yWidth = h_map->GetYaxis()->GetBinWidth(1);
    double cutValue=0;
    
    const int njets = event->jets->size();
    
    for(int i=0; i < njets; i++){
      int idx_x = 0;
      int idx_y = 0;
      Jet* jet = &event->jets->at(i);// loop over all jets in event
      
      while(jet->eta() > xMin+xWidth + idx_x * xWidth) idx_x++;
      while(jet->phi() > yMin+yWidth + idx_y * yWidth) idx_y++;
      
      cutValue = h_map->GetBinContent(idx_x+1, idx_y+1);
      
      if(cutValue > 0) break;
      
      
    }
    
    if(cutValue > 0) return false;
    
    return true;
  }
  
  bool Selection::EnergyEtaCut()
  {
    // cut away events with jets containing energy more than sqrt(s)/2
    assert(event);
    double probejet_eta = event->get(tt_probejet_eta);
    double ptave = event->get(tt_pt_ave);

    if(ptave*cosh(probejet_eta)>3250) return false; //3250 GeV =sqrt(s)/2 with s=13 TeV    
    return true;
}




Selection::~Selection()
{
}

}
