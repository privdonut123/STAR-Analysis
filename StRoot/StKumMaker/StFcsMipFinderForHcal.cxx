#include "StFcsMipFinderForHcal.h"
#include "StFwdTrackMaker/Common.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

using namespace std;

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"

#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFwdTrack.h"
#include "StEvent/StFstHit.h"
#include "StEventTypes.h"
#include "StEventUtilities/StEventHelper.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StChain/StChainOpt.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StBFChain/StBFChain.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include "StThreeVectorD.hh"
// #include "StContainers.h"


#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h"                            // Used in PID selection
#include "TRandom3.h"                            // Used in PID selection
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include <StRunInfo.h>
#include <StEventInfo.h>
#include <StHelix.hh>
#include "StThreeVector.hh"
#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"

//#include "StFwdTrackMaker/Common.h"
//#include "StEvent/StFwdTrackCollection.h"

#define NumberOfTH1F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH2F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH3F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTProfile  100                    // Number of Histograms (must be less than max in .h file)

#define MinTracks         0                    // Discard events with fewer than MinTracks (Flow with 2 tracks makes no sense :-)
#define MinEvents         0                    // Discard DST data when we have fewer than MinEvents on a file sequence
#define pi                TMath::Pi()
#define NMAX         10000

ClassImp(StFcsMipFinderForHcal)

StFcsMipFinderForHcal::StFcsMipFinderForHcal(const char* name) : StMaker(name)
  { // Initialize and/or zero all public/private data members here.
    /*
    for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointerOBs
      {
        h1list_mass_by_Ntower[i] = NULL;      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
        h1list_mass_by_Stower[i] = NULL;      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
        h1list_NEtower[i] = NULL;             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
        h1list_SEtower[i] = NULL;
      }
    */
    mMuDstMaker  = NULL   ; //Without this, maker is not defined in scope. Why?
    mHistogramOutput  =  NULL  ;                  // Zero the pointer to histogram output file

    mEventsStarted    =  0     ;
    mEventsProcessed  =  0     ;                  // Zero the Number of Events processed by the maker

    mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the .C macro
    
  }

StFcsMipFinderForHcal::~StFcsMipFinderForHcal()
  {
    // Destroy and/or zero out all public/private data members here.
  } 

Int_t StFcsMipFinderForHcal::Init()
  { // Do once at the start of the analysis, create histograms, etc.
      cout << termcolor::green << endl << "Initializing, please wait..." << endl << endl << termcolor::reset;
    mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); 

    
    if (!mFcsDb) 
      {
        if (mDebug > 0) {
        cerr << termcolor::red << "StFcsEventDisplay::Init Failed to get StFcsDbMaker" << termcolor::reset << endm;
        }
        return kStFatal;
      }

    for (int i = 0; i < 260; i++) 
      {
        char name_hist[50];
        char title_hist[100];
        sprintf(name_hist, "NEtower_%i", i);
        sprintf(title_hist, "North %i tower energy spectrum", i + 1);
        h1list_NEtower[i] = new TH1F(name_hist, title_hist, 200, 0, 20);
        h1list_NEtower[i]->SetXTitle("Energy [GeV]");
        sprintf(name_hist, "SEtower_%i", i);
        sprintf(title_hist, "South %i tower energy spectrum", i + 1);
        h1list_SEtower[i] = new TH1F(name_hist, title_hist, 200, 0, 20);
        h1list_SEtower[i]->SetXTitle("Energy [GeV]");
      }
    // h1_two_cluster_energy_nocut = new TH1F("h1_two_cluster_energy_nocut", "2 clusters energy(no cut)", bins, 0, 30);

    // h1_hcal_clusters_per_event = new TH1F("h1_hcal_clusters_per_event", "number of clusters per event on hcal", 10, -.5, 9.5);
    // h1_hcal_clusters_per_event->SetXTitle("number of clusters");
    // h1_hcal_clusters_per_event->SetYTitle("counts");
    h1_hcal_clusters_per_event = new TH1F("h1_hcal_clusters_per_event", "number of clusters per event on hcal", 101, -1, 100);
    h1_hcal_clusters_per_event->SetXTitle("number of clusters");
    h1_hcal_clusters_per_event->SetYTitle("counts");
    h1_ecal_clusters_per_event = new TH1F("h1_ecal_clusters_per_event", "number of clusters per event on ecal", 101, -1, 100);
    h1_ecal_clusters_per_event->SetXTitle("number of clusters");
    h1_ecal_clusters_per_event->SetYTitle("counts");
    h1_num_towers_in_hcal_cluster = new TH1F("h1_num_towers_in_hcal_cluster", "number of towers in hcal cluster", 101, -0.5, 100.5);
    h1_num_towers_in_hcal_cluster->SetXTitle("number of towers");
    h1_num_towers_in_hcal_cluster->SetYTitle("counts");
    h1_num_towers_in_ecal_cluster = new TH1F("h1_num_towers_in_ecal_cluster", "number of towers in ecal cluster", 101, -0.5, 100.5);
    h1_num_towers_in_ecal_cluster->SetXTitle("number of towers");
    h1_num_towers_in_ecal_cluster->SetYTitle("counts");
    h1_each_hcal_cluster_energy = new TH1F("h1_each_hcal_cluster_energy", "each cluster energy for HCal", bins, 0, 10);
    h1_each_hcal_cluster_energy->SetXTitle("Energy [GeV]");
    h1_each_hcal_cluster_energy->SetYTitle("counts");
    h1_each_ecal_cluster_energy = new TH1F("h1_each_ecal_cluster_energy", "each cluster energy for ECal", bins, 0, 10);
    h1_each_ecal_cluster_energy->SetXTitle("Energy [GeV]");
    h1_each_ecal_cluster_energy->SetYTitle("counts");
    h2_hcal_cluster_position = new TH2F("h2_hcal_cluster_position", "cluster_position", 300, -150, 150, 300, -150, 150);
    h2_hcal_cluster_position->SetXTitle("x [cm]");
    h2_hcal_cluster_position->SetYTitle("y [cm]");
    h1_num_of_fwd_tracks = new TH1F("h1_num_of_fwd_tracks", "number of forward tracks per event", 101, -1, 100);
    h1_num_of_fwd_tracks->SetXTitle("number of tracks");
    h1_num_of_fwd_tracks->SetYTitle("counts");
    h1_fwd_track_is_primary = new TH1F("h1_fwd_track_is_primary", "forward track is primary", 5, -2.5, 2.5);
    h1_fwd_track_is_primary->SetXTitle("is primary");
    h1_fwd_track_is_primary->SetYTitle("counts");
    h1_fwd_track_num_of_fit_points = new TH1F("h1_fwd_track_num_of_fit_points", "forward track number of fit points", 16, -1, 15);
    h1_fwd_track_num_of_fit_points->SetXTitle("number of fit points");
    h1_fwd_track_num_of_fit_points->SetYTitle("counts");
    h1_fwd_track_charge = new TH1F("h1_fwd_track_charge", "forward track charge", 50, -5, 5);
    h1_fwd_track_charge->SetXTitle("charge [e]");
    h1_fwd_track_charge->SetYTitle("counts");
    h1_fwd_track_pt = new TH1F("h1_fwd_track_pt", "forward track transverse momentum", bins, 0, 15);
    h1_fwd_track_pt->SetXTitle("p_{T} [GeV/c]");
    h1_fwd_track_pt->SetYTitle("counts");
    h1_fwd_track_chi2 = new TH1F("h1_fwd_track_chi2", "forward track chi2", 200, 0, 200);
    h1_fwd_track_chi2->SetXTitle("#chi^{2}");
    h1_fwd_track_chi2->SetYTitle("counts");
    h1_fwd_track_chi2_per_ndf = new TH1F("h1_fwd_track_chi2_per_ndf", "forward track chi2 per ndf", 200, 0, 20);
    h1_fwd_track_chi2_per_ndf->SetXTitle("#chi^{2}/ndf");
    h1_fwd_track_chi2_per_ndf->SetYTitle("counts");
    h1_num_of_track_matched_hcal_clusters = new TH1F("h1_num_of_track_matched_hcal_clusters", "number of track matched hcal clusters per event", 6, -0.5, 5.5);
    h1_num_of_track_matched_hcal_clusters->SetXTitle("number of clusters");
    h1_num_of_track_matched_hcal_clusters->SetYTitle("counts");
    h1_num_of_track_matched_ecal_clusters = new TH1F("h1_num_of_track_matched_ecal_clusters", "number of track matched ecal clusters per event", 6, -0.5, 5.5);
    h1_num_of_track_matched_ecal_clusters->SetXTitle("number of clusters");
    h1_num_of_track_matched_ecal_clusters->SetYTitle("counts");
    h1_track_match_cluster_hcal_energy = new TH1F("h1_track_match_cluster_hcal_energy", "track matched cluster energy on hcal", bins, 0, 10);
    h1_track_match_cluster_hcal_energy->SetXTitle("Energy [GeV]");
    h1_track_match_cluster_hcal_energy->SetYTitle("counts");
    h1_track_match_cluster_ecal_energy = new TH1F("h1_track_match_cluster_ecal_energy", "track matched cluster energy on ecal", bins, 0, 10);
    h1_track_match_cluster_ecal_energy->SetXTitle("Energy [GeV]");
    h1_track_match_cluster_ecal_energy->SetYTitle("counts");
    h1_single_cluster_id = new TH1F("h1_single_cluster_id", "single cluster id on hcal", 101, -0.5, 100.5);
    h1_single_cluster_id->SetXTitle("cluster id");
    h1_single_cluster_id->SetYTitle("counts");
    h1_num_towers_in_hcal_cluster_cut = new TH1F("h1_num_towers_in_hcal_cluster_cut", "number of towers in hcal cluster with cuts", 101, -0.5, 100.5);
    h1_num_towers_in_hcal_cluster_cut->SetXTitle("number of towers");
    h1_num_towers_in_hcal_cluster_cut->SetYTitle("counts");
    h1_num_towers_in_ecal_cluster_cut = new TH1F("h1_num_towers_in_ecal_cluster_cut", "number of towers in ecal cluster with cuts", 101, -0.5, 100.5);
    h1_num_towers_in_ecal_cluster_cut->SetXTitle("number of towers");
    h1_num_towers_in_ecal_cluster_cut->SetYTitle("counts");
    h1_track_match_cluster_hcal_energy_tower_cut = new TH1F("h1_track_match_cluster_hcal_energy_tower_cut", "track matched cluster energy on hcal with tower cut", bins, 0, 10);
    h1_track_match_cluster_hcal_energy_tower_cut->SetXTitle("Energy [GeV]");
    h1_track_match_cluster_hcal_energy_tower_cut->SetYTitle("counts");
    h1_track_match_cluster_ecal_energy_tower_cut = new TH1F("h1_track_match_cluster_ecal_energy_tower_cut", "track matched cluster energy on ecal with tower cut", bins, 0, 10);
    h1_track_match_cluster_ecal_energy_tower_cut->SetXTitle("Energy [GeV]");
    h1_track_match_cluster_ecal_energy_tower_cut->SetYTitle("counts");
    h1_hcal_hit_energy = new TH1F("h1_hcal_hit_energy", "hcal hit energy", bins, 0, 10);
    h1_hcal_hit_energy->SetXTitle("Energy [GeV]");
    h1_hcal_hit_energy->SetYTitle("counts");
    // h1_hcal_max_cluster_energy_per_event = new TH1F("h1_hcal_max_cluster_energy_per_event", "hcal max cluster energy per event", 100, 0, 100);
    // h1_hcal_max_cluster_energy_per_event->SetXTitle("Energy [GeV]");
    // h1_hcal_max_cluster_energy_per_event->SetYTitle("counts");
    // h1_Zgg_nocut_cluster = new TH1F("h1_Zgg_nocut_cluster", "Zgg without cut", bins, 0, 1);
    // h1_inv_mass_cluster_nocut = new TH1F("h1_inv_mass_cluster_nocut", "invariant mass plot (no cut)", bins, m_low, m_up);
    // h1_num_of_points = new TH1F("h1_num_of_points", "number of points", 10, -.5, 9.5);
    // h1_num_of_points->SetXTitle("number of points per event");
    // h1_num_of_points->SetYTitle("counts");

    // h1_ecal_clusters_per_event = new TH1F("h1_ecal_clusters_per_event", "number of clusters per event on ecal", 10, -.5, 9.5);
    // h1_ecal_clusters_per_event->SetXTitle("number of clusters");
    // h1_ecal_clusters_per_event->SetYTitle("counts");
    
    // h1_hcal_neigbors_per_cluster = new TH1F("h1_hcal_neigbors_per_cluster", "number of neighbors per cluster on hcal", 10, -.5, 9.5);
    // h1_hcal_neigbors_per_cluster->SetXTitle("number of neighbors");
    // h1_hcal_neigbors_per_cluster->SetYTitle("counts");

    


    

    
    h1_fwd_pt_res = new TH1F("h1_fwd_pt_res", "forward track p_{T} resolution", bins, -1, 1);
    h1_fwd_pt_res->SetXTitle("(p_{T} reconstructed - p_{T} geant_prim) / p_{T} geant_prim");
    h1_fwd_pt_res->SetYTitle("Counts");
    // h1_fwd_cal_hit_resolution = new TH1F("h1_fwd_cal_hit_resolution", "forward calorimeter hit energy resolution", bins, -1, 1);
    // h1_fwd_cal_hit_resolution->SetXTitle("(E_{reconstructed} - E_{geant_prim}) / E_{geant_prim}");
    // h1_fwd_cal_hit_resolution->SetYTitle("counts");
    // h1_fwd_cal_cluster_resolution = new TH1F("h1_fwd_cal_cluster_resolution", "forward calorimeter cluster energy resolution", bins, -1, 1);
    // h1_fwd_cal_cluster_resolution->SetXTitle("(E_{reconstructed} - E_{geant_prim}) / E_{geant_prim}");
    // h1_fwd_cal_cluster_resolution->SetYTitle("counts");
    
    
    h2_ecal_cluster_position = new TH2F("h2_ecal_cluster_position", "cluster_position", 300, -150, 150, 300, -150, 150);
    h2_ecal_cluster_position->SetXTitle("x [cm]");
    h2_ecal_cluster_position->SetYTitle("y [cm]");
    
    h2_ecal_hcal_cluster_energy = new TH2F("h2_ecal_hcal_cluster_energy", "ecal_hcal_summed_cluster_energy", bins, 0, 150, bins, 0, 150);
    h2_ecal_hcal_cluster_energy->SetXTitle("Ecal summed cluster energy [GeV]");
    h2_ecal_hcal_cluster_energy->SetYTitle("Hcal summed cluster energy [GeV]");
    h2_ecal_hcal_hit_energy = new TH2F("h2_ecal_hcal_hit_energy", "ecal_hcal_summed_hit_energy", bins, 0, 150, bins, 0, 150);
    h2_ecal_hcal_hit_energy->SetXTitle("Ecal summed hit energy [GeV]");
    h2_ecal_hcal_hit_energy->SetYTitle("Hcal summed hit energy [GeV]");
    h2_fwd_cal_energy_vs_track_pt_hit = new TH2F("h2_fwd_cal_energy_vs_track_pt_hit", "forward summed calorimeter hit energy vs track p_{T} adjusted for rapidity", bins, 0, 150, bins, 0, 150);
    // h2_fwd_cal_energy_vs_track_pt_hit = new TH2F("h2_fwd_cal_energy_vs_track_pt", "forward calorimeter energy vs track pt * Cosh ", bins, 0, 150, bins, 0, 150);
    h2_fwd_cal_energy_vs_track_pt_hit->SetXTitle("Track p_{T} * Cosh(#eta) [GeV/c]"); 
    h2_fwd_cal_energy_vs_track_pt_hit->SetYTitle("Calorimeter energy [GeV]");
    h2_fwd_cal_energy_vs_track_pt_cluster = new TH2F("h2_fwd_cal_energy_vs_track_pt_cluster", "forward summed calorimeter cluster energy vs track p_{T} adjusted for rapidity", bins, 0, 150, bins, 0, 150);
    h2_fwd_cal_energy_vs_track_pt_cluster->SetXTitle("Track p_{T} * Cosh(#eta) [GeV/c]");
    h2_fwd_cal_energy_vs_track_pt_cluster->SetYTitle("Calorimeter energy [GeV]");

    h2_fwd_pt_res_vs_eta = new TH2F("h2_fwd_pt_res_vs_eta", "forward track p_{T} resolution vs eta", bins, -1, 1, bins, 2.4, 4.1);
    h2_fwd_pt_res_vs_eta->SetXTitle("p_{T} resolution");
    h2_fwd_pt_res_vs_eta->SetYTitle("#eta");

    h2_ecal_hcal_cluster_energy_res_vs_eta = new TH2F("h2_ecal_hcal_cluster_energy_res_vs_eta", "ecal_hcal_cluster_energy_res_vs_eta", bins, -1, 1, bins, 2.4, 4.1);
    h2_ecal_hcal_cluster_energy_res_vs_eta->SetXTitle("cluster energy resolution");
    h2_ecal_hcal_cluster_energy_res_vs_eta->SetYTitle("#eta");
    h2_ecal_hcal_hit_energy_res_vs_eta = new TH2F("h2_ecal_hcal_hit_energy_res_vs_eta", "ecal_hcal_hit_energy_res_vs_eta", bins, -1, 1, bins, 2.4, 4.1);
    h2_ecal_hcal_hit_energy_res_vs_eta->SetXTitle("hit energy resolution");
    h2_ecal_hcal_hit_energy_res_vs_eta->SetYTitle("#eta");
    


    h1_fwd_cal_energy_hit = new TH1F("h1_fwd_cal_energy_hit", "forward calorimeter energy (hit)", bins, 0, 150);
    h1_fwd_cal_energy_hit->SetXTitle("Energy [GeV]");
    h1_fwd_cal_energy_hit->SetYTitle("counts");
    h1_fwd_cal_energy_cluster = new TH1F("h1_fwd_cal_energy_cluster", "forward calorimeter energy (cluster)", bins, 0, 150);
    h1_fwd_cal_energy_cluster->SetXTitle("Energy [GeV]");
    h1_fwd_cal_energy_cluster->SetYTitle("counts");


    return kStOK;
  }

Int_t StFcsMipFinderForHcal::Finish( )
  { // Do once at the end of the analysis, close files, etc.
    cout << termcolor::green << termcolor::bold << "Finishing, please wait..." << termcolor::reset << endl;
    // Write all histograms to disk
    mHistogramOutput = new TFile( mHistogramOutputFileName.Data() , "recreate" ) ;
    // h1_two_cluster_energy_nocut->Write();
    h1_hcal_clusters_per_event->Write();
    h1_ecal_clusters_per_event->Write();
    h1_num_towers_in_hcal_cluster->Write();
    h1_num_towers_in_ecal_cluster->Write();
    h1_each_hcal_cluster_energy->Write();
    h1_each_ecal_cluster_energy->Write();
    h2_hcal_cluster_position->Write();
    h2_ecal_cluster_position->Write();
    h1_num_of_fwd_tracks->Write();
    h1_fwd_track_is_primary->Write();
    h1_fwd_track_num_of_fit_points->Write();
    h1_fwd_track_charge->Write();
    h1_fwd_track_pt->Write();
    h1_fwd_track_chi2->Write();
    // h1_fwd_track_chi2_per_ndf->Write();
    h1_num_of_track_matched_hcal_clusters->Write();
    h1_num_of_track_matched_ecal_clusters->Write();
    h1_track_match_cluster_hcal_energy->Write();
    h1_track_match_cluster_ecal_energy->Write();
    // h1_single_cluster_id->Write();
    h1_num_towers_in_hcal_cluster_cut->Write();
    h1_num_towers_in_ecal_cluster_cut->Write();
    h1_track_match_cluster_hcal_energy_tower_cut->Write();
    h1_track_match_cluster_ecal_energy_tower_cut->Write();
    h1_hcal_hit_energy->Write();
    // h1_hcal_max_cluster_energy_per_event->Write();
    // h1_Zgg_nocut_cluster->Write();
    // h1_num_of_points->Write();
    // h1_inv_mass_cluster_nocut->Write();
    // h1_ecal_clusters_per_event->Write();
    
    // h1_hcal_neigbors_per_cluster->Write();
    
   
    h1_fwd_pt_res->Write();
    
    for (int i = 0; i < 260; i++) {
      // h1list_mass_by_Ntower[i]->Write();
      // h1list_mass_by_Stower[i]->Write();
      h1list_NEtower[i]->Write();
      h1list_SEtower[i]->Write();
   }
    // h1_geant_shower_proj_z->Write();
    // h2_geant_shower_proj_xy->Write();
    // h1_geant_parent_eta->Write();
    // h1_geant_parent_pt->Write();
    // h1_geant_parent_pz->Write();
    // h1_geant_track_id->Write();
    // h1_geant_particle_id->Write();
    // h1_geant_start_vertex->Write();
    // h1_geant_stop_vertex->Write();
    // h1_geant_num_of_tracks->Write();
    // h1_geant_primary_eta->Write();
    // h1_geant_primary_pt->Write();
    // h1_geant_primary_pz->Write();
    // h1_geant_primary_energy->Write();
    // h1_geant_vtx_z->Write();
    // h1_fwd_cal_hit_resolution->Write();
    // h1_fwd_cal_cluster_resolution->Write();
    
    
    
    h2_ecal_hcal_cluster_energy->Write();
    h2_ecal_hcal_hit_energy->Write();
    h2_fwd_cal_energy_vs_track_pt_hit->Write();
    h2_fwd_cal_energy_vs_track_pt_cluster->Write();
    
    h1_fwd_cal_energy_hit->Write();
    h1_fwd_cal_energy_cluster->Write();
    h2_fwd_pt_res_vs_eta->Write();
    h2_ecal_hcal_cluster_energy_res_vs_eta->Write();
    h2_ecal_hcal_hit_energy_res_vs_eta->Write();
    mHistogramOutput->Write();
    mHistogramOutput->Close();
    
    return kStOK;
 }


Int_t StFcsMipFinderForHcal::Make()
  { // Do every event
    if (mDebug > 0) 
      {
        cout << termcolor::green << "StFcsMipFinderForHcal::Make() called" << endl << termcolor::reset;
      }
    
    loadCollectionsAndDb();

    int total_nc = 0;
    int total_np = 0;


    for (int det = 0; det <= 3; det++) 
      {
        StSPtrVecFcsCluster& clusters = fcsColl->clusters(det);
        int nc = fcsColl->numberOfClusters(det);
        if (mDebug > 0) 
          {
            cout << termcolor::yellow << "Number of clusters: " << nc << termcolor::reset << endl;
          }
        // if (mDebug > 0) 
        //   LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;
        if(det == 0 || det == 1)
          {
            h1_ecal_clusters_per_event->Fill(nc);
          }
        else if(det == 2 || det == 3)
          {
            h1_hcal_clusters_per_event->Fill(nc);
          }
        // h1_hcal_clusters_per_event->Fill(nc);
        total_nc = nc + total_nc;
        // float maxEnergy = -1; // Initialize max energy to a very small value
        for (int i = 0; i < nc && nc > 0; i++) 
          {
            StFcsCluster* clu = clusters[i];
            // total_hcal_cluster_energy = 0;

            // float clu_energy = clu->energy();
            int nTowers_ecal = 0;
            int nTowers_hcal = 0;
            if(det == 0 || det == 1)
              {
                h1_num_towers_in_ecal_cluster->Fill(clu->nTowers());
                h1_each_ecal_cluster_energy->Fill(clu->energy());
                nTowers_ecal = clu->nTowers();
              }
            else if(det == 2 || det == 3)
              {
                h1_num_towers_in_hcal_cluster->Fill(clu->nTowers());
                h1_each_hcal_cluster_energy->Fill(clu->energy());
                nTowers_hcal = clu->nTowers();
              }
            // h1_each_hcal_cluster_energy->Fill(clu->energy());
            if (mDebug > 0)
              {          
                cout << termcolor::yellow << "Energy of cluster: " << clu->energy() << " GeV" << endl;
                cout << "In detector coordinates, cluster at x = " << clu->x() << ", y = " << clu->y() << termcolor::reset << endl;
                cout << termcolor::bright_yellow << "scale factor for x: " << mFcsDb->getXWidth(det) << ", y: " << mFcsDb->getYWidth(det) << termcolor::reset << endl;
              }
            if (nTowers_ecal <= 3 && nTowers_hcal <= 3) 
              {
                
              
                StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu->x(), clu->y());
                StLorentzVectorD p = mFcsDb->getLorentzVector(cluPos, clu->energy(), 0);
                if (mDebug > 0) 
                  {
                    cout << termcolor::underline << termcolor::yellow << "In physical coordinates, cluster at x = " << cluPos.x() << ", y = " << cluPos.y() << termcolor::reset << endl;
                  }
                if (det == 0 || det == 1)
                  {
                    h2_ecal_cluster_position->Fill(cluPos.x(), cluPos.y());
                  }
                else if (det == 2 || det == 3)
                  {
                    h2_hcal_cluster_position->Fill(cluPos.x(), cluPos.y());
                  }

                // if (det == 2 || det == 3)
                //   {
                //     h2_hcal_cluster_position->Fill(cluPos.x(), cluPos.y());
                //   }
                
                total_hcal_cluster_energy = total_hcal_cluster_energy + clu->energy();
                
                // Implement cut on number of towers in a cluster
                
                // if (clu->numberOfTowers() <= 2) 
                //   {
                //     h1_track_match_cluster_hcal_energy->Fill(clu->energy());
                //   }
                
                  // }
                // if (i == nc - 1) 
                //   continue;
                // for (int j = i + 1; j < nc; j++) 
                //   {
                //     StFcsCluster* cluj = clusters[j];
                //     float cluj_energy = cluj->energy();
                //     float cluj_x = cluj->x();
                //     float cluj_y = cluj->y();
                //     StThreeVectorD clujPos = mFcsDb->getStarXYZfromColumnRow(det, cluj_x, cluj_y);

                //     h1_two_cluster_energy_nocut->Fill(clu_energy + cluj_energy);
                //     float zgg = (abs(clu_energy - cluj_energy)) / (clu_energy + cluj_energy);
                //     h1_Zgg_nocut_cluster->Fill(zgg);
                //     StThreeVectorD xyzj = mFcsDb->getStarXYZfromColumnRow(det, cluj->x(), cluj->y());
                //     StLorentzVectorD pj = mFcsDb->getLorentzVector((xyzj), cluj->energy(), 0);
                //     h1_inv_mass_cluster_nocut->Fill((p + pj).m());
                //   }

              }
          }
      }
    h1_num_of_fwd_tracks -> Fill(ftc->numberOfTracks());  
    if (ftc->numberOfTracks() > 0) // CURRENT HACK OF FWDTRACKMAKER
      {
        for(int i = 0; i < ftc->numberOfTracks(); i++)
          {
            StFwdTrack* fwdTrack = ftc->tracks()[i];

            // if (mDebug > 0) 
            //   {
            //     cout << termcolor::blue << "Vertex Index: " << fwdTrack->vertexIndex() << termcolor::reset << endl;
            //   }
            if (fwdTrack->isPrimaryTrack() == true)
              {
                h1_fwd_track_is_primary->Fill(1);
              }
            else
              {
                h1_fwd_track_is_primary->Fill(-1);
              }

            h1_fwd_track_num_of_fit_points -> Fill(fwdTrack->numberOfFitPoints());
            // h1_fwd_track_charge -> Fill((int)fwdTrack->charge());
            // h1_fwd_track_pt -> Fill(fwdTrack->momentum().perp());
            // h1_fwd_track_chi2_per_ndf -> Fill(fwdTrack->chi2() / fwdTrack->ndf());
            h1_fwd_track_chi2 -> Fill(fwdTrack->chi2());
            if (mDebug > 0)
              {
                // cout << termcolor::blue << "chi2/ndf: " << fwdTrack->chi2() / fwdTrack->ndf() << endl;
                cout << termcolor::blue << "chi2: " << fwdTrack->chi2() << endl;
                cout << termcolor::blue << "ndf: " << fwdTrack->ndf() << endl;
                cout << termcolor::blue << "Number of fit points: " << fwdTrack->numberOfFitPoints() << endl;
                cout << termcolor::blue << "Charge: " << (int)fwdTrack->charge() << endl;
                cout << termcolor::blue << "Momentum: " << fwdTrack->momentum() << endl;
              }
            vector<StFwdTrack> fwdTrackClusterIdArray;
            // Issue with current fwdtrack maker, chi2/ndf is infinite
            // if (fwdTrack->isPrimaryTrack() == true && fwdTrack->numberOfFitPoints() >= 8 && fwdTrack->chi2() / fwdTrack->ndf() < 2.0) 
            if (fwdTrack->isPrimaryTrack() == true && fwdTrack->numberOfFitPoints() >= 4 && fwdTrack->chi2() < 20)
              {
                h1_fwd_track_charge -> Fill((int)fwdTrack->charge());
                h1_fwd_track_pt -> Fill(fwdTrack->momentum().perp());
                // h1_fwd_track_chi2_per_ndf -> Fill(fwdTrack->chi2() / fwdTrack->ndf());

                cout << termcolor::bold << termcolor::blue << "chi2: " << fwdTrack->chi2() << endl;
                cout << termcolor::bold << termcolor::blue << "ndf: " << fwdTrack->ndf() << endl;
                cout << termcolor::bold << termcolor::blue << "Charge: " << (int)fwdTrack->charge() << endl;
                cout << termcolor::bold << termcolor::blue << "Momentum: " << fwdTrack->momentum() << endl;
                  

                //TOF mult cut                                                                                                     
                int tofMult = 0;
                const StTriggerData* trgdata = event->triggerData();
                if(!trgdata && StMuDst::event()) 
                  {
                    trgdata = StMuDst::event()->triggerData();
                  }
                if(trgdata)
                  {
                    tofMult = trgdata->tofMultiplicity();
                    LOG_DEBUG<<"TOF mult="<<tofMult<<endm;
                    if (tofMult > 100) 
                      {
                        return kStOK;
                      }
                  }
                else
                  {
                    LOG_WARN << "No TriggerData found in StEvent nor Mudst. No TOFMult cut"<<endm;
                  }

                //TPC ZVERTEX
                float zTPC=-999.0;
                StPrimaryVertex* tpcvtx = event->primaryVertex();
                if(tpcvtx) 
                  {
                    zTPC=tpcvtx->position().z();
                  }
                else
                  {
                    if (StMuDst::numberOfPrimaryVertices() > 0)
                      {
                        StMuPrimaryVertex* mutpcvtx = StMuDst::primaryVertex();
                        if(mutpcvtx) 
                          {
                            zTPC=mutpcvtx->position().z();
                          }
                      }
                  }

                //BBC ZVERTEX
                float zBBC=-999.0;
                if(trgdata) zBBC = (4096 - trgdata->bbcTimeDifference())*0.016*30.0/2.0;      
                if(zBBC<-200 || zBBC>200) zBBC=-999;

                //VPD ZVERTEX from MuDst(TOF data)
                float zVPD=-999.0;
                if(StMuDst::btofHeader()) zVPD=StMuDst::btofHeader()->vpdVz();

                LOG_INFO << Form("ZTPX = %6.2f ZBBC = %6.2f ZVPD = %6.2f",zTPC,zBBC,zVPD) << endm;

                //test getLorentzVector
                StThreeVectorD xyz(20,0,720);     
                StLorentzVectorD pbbc,ptpc,pvpd;
                StLorentzVectorD p0 = mFcsDb->getLorentzVector((xyz), 10,    0);  LOG_INFO << "Zero " << p0 << endm;  	   
                if(zBBC>-200) {pbbc  = mFcsDb->getLorentzVector((xyz), 10, zBBC);	LOG_INFO << "BBC  " << pbbc << endm;}
                if(zTPC>-200) {ptpc  = mFcsDb->getLorentzVector((xyz), 10, zTPC);	LOG_INFO << "TPC  " << ptpc << endm;}
                if(zVPD>-200) {pvpd  = mFcsDb->getLorentzVector((xyz), 10, zVPD); LOG_INFO << "VPD  " << pvpd << endm;}

                float prim_trk_energy;
                StPtrVecFcsCluster& trackMatchEcalClusters = fwdTrack->ecalClusters();
                StPtrVecFcsCluster& trackMatchHcalClusters = fwdTrack->hcalClusters();

                // if (find(fwdTrackClusterIdArray.begin(), fwdTrackClusterIdArray.end(), 
                
                if (mDebug > 0)
                  {
                    cout << termcolor::green << "Number of track matched HCal clusters: " << trackMatchHcalClusters.size() << termcolor::reset << endl;
                  }
                
                for (int det = 0; det <= 3; det++)
                  {
                    if (det == 0 || det == 1)
                      {
                        h1_num_of_track_matched_ecal_clusters -> Fill(trackMatchEcalClusters.size());
                      }
                    else if (det == 2 || det == 3)
                      {
                        h1_num_of_track_matched_hcal_clusters -> Fill(trackMatchHcalClusters.size());
                      }
                  }

                // fwdTrackClusterIdArray.push_back(hcalCluster->id());

                for( int i = 0; i < trackMatchHcalClusters.size(); i++ )
                  {
                    StFcsCluster* hcalCluster = trackMatchHcalClusters[i];
                    StPtrVecFwdTrack& hcalClusterMatchedTracks = hcalCluster->tracks();
                    h1_num_towers_in_ecal_cluster_cut -> Fill(hcalCluster->nTowers());
                    for(int j = 0; j < trackMatchEcalClusters.size(); j++)
                      {
                        StFcsCluster* ecalCluster = trackMatchEcalClusters[j];
                        StPtrVecFwdTrack& ecalClusterMatchedTracks = ecalCluster->tracks();
                        h1_num_towers_in_hcal_cluster_cut -> Fill(ecalCluster->nTowers());
                    
                        // h1_single_cluster_id -> Fill(hcalCluster->id());
                        cout << termcolor::blue << "number of tracks matched to HCal cluster: " << hcalClusterMatchedTracks.size() << termcolor::reset << endl;
                        cout << termcolor::blue << "number of tracks matched to ECal cluster: " << ecalClusterMatchedTracks.size() << termcolor::reset << endl;
                        cout << termcolor::blue << "number of HCal clusters matched to track: " << trackMatchHcalClusters.size() << termcolor::reset << endl;
                        cout << termcolor::blue << "number of ECal clusters matched to track: " << trackMatchEcalClusters.size() << termcolor::reset << endl;
                        // This is a hack since we need to verify that the cluster matched tracks are primaries, clusterMatchedTracks.size() should not be 3
                        if (trackMatchHcalClusters.size() == 1 && trackMatchEcalClusters.size() == 1 && hcalClusterMatchedTracks.size() == 3 && ecalClusterMatchedTracks.size() == 3)
                          {
                            // NOTE: Check if cluster id is unique for each track
                            h1_track_match_cluster_hcal_energy -> Fill(hcalCluster->energy());
                            h1_track_match_cluster_ecal_energy -> Fill(ecalCluster->energy());
                            if (hcalCluster->nTowers() <= 3 && ecalCluster->nTowers() <= 3) 
                              {
                                if (mDebug > 0)
                                  {
                                    cout << termcolor::blue << "HCal Cluster: " << hcalCluster->id() 
                                          << ", energy: " << hcalCluster->energy() 
                                          << ", number of towers: " << hcalCluster->nTowers()
                                          << termcolor::reset << endl;
                                  }
                                h1_track_match_cluster_hcal_energy_tower_cut -> Fill(hcalCluster->energy());
                                h1_track_match_cluster_ecal_energy_tower_cut -> Fill(ecalCluster->energy());
                          
                                StPtrVecFcsHit& hits = hcalCluster->hits();
                                for(int k = 0; k < hits.size(); k++)
                                  {
                                    StFcsHit* hit = hits[k];
                                    int towerId = hit->id();
                                    h1_hcal_hit_energy -> Fill(hit->energy());
                                    if (hit->detectorId() == 2)
                                      {
                                        h1list_NEtower[towerId]->Fill(hit->energy());
                                      }
                                    else if (hit->detectorId() == 3)
                                      {
                                        h1list_SEtower[towerId]->Fill(hit->energy());
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
            
            total_ecal_hit_energy = 0;
            total_hcal_hit_energy = 0;      
            
            // StThreeVectorD north_ecal_corner1 = mFcsDb->getStarXYZfromColumnRow(0, 1, 1);
            // StThreeVectorD north_ecal_corner2 = mFcsDb->getStarXYZfromColumnRow(0, kFcsEcalNCol, kFcsEcalNRow);
            // StThreeVectorD south_ecal_corner1 = mFcsDb->getStarXYZfromColumnRow(1, 1, 1);
            // StThreeVectorD south_ecal_corner2 = mFcsDb->getStarXYZfromColumnRow(1, kFcsEcalNCol, kFcsEcalNRow);
            // if(mDebug > 0)
            //   {
            //     cout << termcolor::yellow << "North ECalorimeter dimensions: " << north_ecal_corner1.x() << ", " << north_ecal_corner1.y() << ", " << north_ecal_corner1.z() << termcolor::reset << endl;
            //     cout << termcolor::yellow << "North ECalorimeter dimensions: " << north_ecal_corner2.x() << ", " << north_ecal_corner2.y() << ", " << north_ecal_corner2.z() << termcolor::reset << endl;
            //     cout << termcolor::yellow << "South ECalorimeter dimensions: " << south_ecal_corner1.x() << ", " << south_ecal_corner1.y() << ", " << south_ecal_corner1.z() << termcolor::reset << endl;
            //     cout << termcolor::yellow << "South ECalorimeter dimensions: " << south_ecal_corner2.x() << ", " << south_ecal_corner2.y() << ", " << south_ecal_corner2.z() << termcolor::reset << endl;
            //   }
            // bool dead_zone_flag = false;
            St_g2t_track* trackTable = static_cast<St_g2t_track*>(GetDataSet("g2t_track"));
            St_g2t_vertex* vertexTable = static_cast<St_g2t_vertex*>(GetDataSet("g2t_vertex"));
            g2t_track_st* g2ttrk = 0;
            g2t_vertex_st* g2tvert = 0;

            if( !trackTable )
              { 
                cout << termcolor::red << "g2t_track Table not found" << termcolor::reset << std::endl; 
                // continue; 
              }
            else
              {
                const int nTrk = trackTable->GetNRows();
                // h1_geant_num_of_tracks -> Fill(nTrk);
                if( mDebug>0 )
                  { 
                    cout << termcolor::green << "g2t_track table has "<< nTrk << " tracks" << termcolor::reset << std::endl; 
                  }
                if( nTrk>0 )
                  {
                    g2ttrk = trackTable->GetTable();
                    g2tvert = vertexTable->GetTable();
                    if(mDebug > 0)
                      {
                        if( !g2ttrk ) 
                          { 
                            cout << termcolor::red << " g2t_track GetTable failed" << termcolor::reset << endl;
                            // continue; 
                          }
                        if( !g2tvert ) 
                          { 
                            cout << termcolor::red << " g2t_vertex GetTable failed" << termcolor::reset << endl;
                            // continue; 
                          }
                      }
                  }
              }
            g2t_track_st* primtrk = g2ttrk;
            // h1_geant_primary_eta->Fill(primtrk->eta);
            // h1_geant_primary_pt->Fill(primtrk->pt);
            // h1_geant_primary_pz->Fill(primtrk->p[2]);
            // h1_geant_primary_energy->Fill(primtrk->e);

            g2t_vertex_st* primvtx = g2tvert;
            // h1_geant_vtx_z->Fill(primvtx->ge_x[2]);

            h1_fwd_pt_res->Fill((fwdTrack->momentum().perp() - primtrk->pt) / primtrk->pt);
            h2_fwd_pt_res_vs_eta->Fill((fwdTrack->momentum().perp() - primtrk->pt) / primtrk->pt, fwdTrack->momentum().pseudoRapidity());

            // Calorimeter energy
            // for (int det = 0; det <= 3; det++) 
            //   {
            //     if (mDebug > 0)
            //       {
            //         cout << termcolor::yellow << termcolor::bold << "Detector: " << det << termcolor::reset << endl;
            //       }

            //     StSPtrVecFcsHit& hits = fcsColl->hits(det);
            //     int nh = fcsColl->numberOfHits(det);
            //     if (mDebug > 0)
            //       {
            //         cout << termcolor::yellow << "Number of hits: " << nh << termcolor::reset << endl;
            //       }
            //     for (int i = 0; i < nh; i++) 
            //       {
            //         //implement cut on energy
            //         StFcsHit* hit = hits[i];
            //         unsigned short hit_id = hit->id();
            //         if (hit->energy() > 0)
            //           {
            //             if (mDebug > 0) 
            //               {
            //                 cout << termcolor::yellow << "Hit ID: " << hit_id << termcolor::reset << endl;
            //               }
            //             float hit_energy = hit->energy();
            //             if (mDebug > 0)
            //               {
            //               cout << termcolor::yellow << "Energy of hit: " << hit->energy() << " GeV" << termcolor::reset << endl;
            //               }
            //             if(det == 0 || det == 1)
            //               {
            //                 total_ecal_hit_energy = total_ecal_hit_energy + hit_energy;
            //               }
            //             else if(det == 2 || det == 3)
            //               {
            //                 total_hcal_hit_energy = total_hcal_hit_energy + hit_energy;
            //               }
            //           }
            //       }

            //       }
                // if (mDebug > 0) 
                //   {
                //     cout << termcolor::green << "Finishing detector " << det << termcolor::reset << endl;
                //   }
                // if (maxEnergy > 0) 
                //   {
                //     if(det == 2 || det == 3)
                //       {
                //         h1_hcal_max_cluster_energy_per_event->Fill(maxEnergy);
                //       }
                //   }
              } 

            
            // //cut on energy, E_min
            // if (total_ecal_hit_energy > E_min && total_hcal_hit_energy > E_min) 
            //   {
            //     h2_ecal_hcal_hit_energy->Fill(total_ecal_hit_energy, total_hcal_hit_energy);
            //     h2_fwd_cal_energy_vs_track_pt_hit->Fill(fwdTrack->momentum().perp()*TMath::CosH(fwdTrack->momentum().pseudoRapidity()), total_ecal_hit_energy + total_hcal_hit_energy);
            //     h1_fwd_cal_energy_hit->Fill(total_ecal_hit_energy + total_hcal_hit_energy);
            //     // h1_fwd_cal_hit_resolution->Fill(((total_ecal_hit_energy + total_hcal_hit_energy) - primtrk->e) / primtrk->e);
            //     h2_ecal_hcal_hit_energy_res_vs_eta->Fill(((total_ecal_hit_energy + total_hcal_hit_energy) - primtrk->e) / primtrk->e, fwdTrack->momentum().pseudoRapidity());
            //   }
            // if (total_ecal_cluster_energy > E_min && total_hcal_cluster_energy > E_min) 
            //   {
            //     h2_ecal_hcal_cluster_energy->Fill(total_ecal_cluster_energy, total_hcal_cluster_energy);
            //     h2_fwd_cal_energy_vs_track_pt_cluster->Fill(fwdTrack->momentum().perp()*TMath::CosH(fwdTrack->momentum().pseudoRapidity()), total_ecal_cluster_energy + total_hcal_cluster_energy);
            //     h1_fwd_cal_energy_cluster->Fill(total_ecal_cluster_energy + total_hcal_cluster_energy);
            //     // h1_fwd_cal_cluster_resolution->Fill(((total_ecal_cluster_energy + total_hcal_cluster_energy) - primtrk->e) / primtrk->e);
            //     h2_ecal_hcal_cluster_energy_res_vs_eta->Fill(((total_ecal_cluster_energy + total_hcal_cluster_energy) - primtrk->e) / primtrk->e, fwdTrack->momentum().pseudoRapidity());
            //   }
      }
  
            //dead_zone_flag = false;
        
      return kStOK ;
  }
  

  
void StFcsMipFinderForHcal::loadCollectionsAndDb()
  {
    event = (StEvent*)GetInputDS("StEvent");
    if (!event) 
      {
        if (mDebug > 0)
          {
            cerr << termcolor::red << "StKumMaker::Make did not find StEvent" << endl << termcolor::reset;
          }
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got StEvent from input" << endl << termcolor::reset;
          }
      }
      
	  mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst");
    if (!mMuDstMaker)
      {
        if (mDebug > 0) 
          {
            cout << termcolor::bold << termcolor::red << "Could not retrieve MuDstMaker from the chain" << termcolor::reset << endl;
          }
      } 
    else
      {
        if (mDebug > 0) 
          {
            cout << termcolor::green << "Got MuDstMaker from the chain" << endl << termcolor::reset;
          }
      }
    
    StMuDst* muDst = mMuDstMaker->muDst();
    if ( !muDst )
      {
        if (mDebug > 0) 
          {
            cerr << termcolor::red << "Failed to get MuDst from input" << endl << termcolor::reset;
          }
      }
    else
      {
        if (mDebug > 0) 
          {
            cout  << termcolor::green << "Got MuDst from input" << endl << termcolor::reset;
          }
      }
      
    StMuEvent* muEvent =  muDst->event();
    if ( !muEvent )
      {
        if (mDebug > 0) 
          {
            cerr << termcolor::red << "Failed to get MuEvent from MuDst" << endl << termcolor::reset;
          }
      }
    else
      {
        if (mDebug > 0) 
          {
            cout << termcolor::green << "Got MuEvent from MuDst" << endl << termcolor::reset;
          }
      }

    StMuTriggerIdCollection* TrigMuColl = &(muEvent->triggerIdCollection());
    if(!TrigMuColl)
      {
        if (mDebug > 0) 
          {
            cerr << termcolor::red << "StFmsMaker::PopulateEventInfo - !TrigMuColl" << termcolor::reset << endl;
          }
      }
    else 
      {
        if (mDebug > 0) 
          {
            cout << termcolor::green << "Got TriggerIdCollection from MuEvent" << endl << termcolor::reset;
          }
      }

	  fcsColl = event->fcsCollection();
    if (!fcsColl)
      {
        if (mDebug > 0)
          {
            cerr << termcolor::red << "No FCS collection" << termcolor::reset << endl;
          }
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got FCS collection" << endl << termcolor::reset;
            // cout << termcolor::green << "Number of FwdTracks: " << ftc->numberOfTracks() << termcolor::reset << endl;
          }
      } 
      
    ftc = event->fwdTrackCollection();
    if (!ftc)
      {
        if (mDebug > 0)
          {
            cerr << termcolor::bold << termcolor::red << "No FwdTrack collection" << termcolor::reset << endl;
          }
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got FwdTrack collection" << endl << termcolor::reset;
          }
      }
    return;
  }

