#ifndef StFcsMipFinderForHcal_def
#define StFcsMipFinderForHcal_def

#include "StChain/StMaker.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <TTree.h>
// ROOT includes

#include "TNtuple.h"
#include "TTree.h"
// STL includes
#include <vector>
#include <memory>
#include <cmath>

#include "StEvent/StTriggerId.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "/direct/star+u/bmagh001/MuDst/include/termcolor.hpp"

// #include "Stiostream.h"
// #include "StObject.h"
// #include <vector>
// #include "StThreeVectorD.hh"
// #include "StContainers.h"
// #include <climits>

class StMuDstMaker ;
class StMuEvent    ;
class StMuTrack    ;
//class StSpinDbMaker;
class StEvent;
class StPicoEvent;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TProfile2D;
class TRandom;
class TRandom3;
class TCanvas;
class StFcsDb;
class StFcsDbMaker;
class StFcsCollection;
class StMuFwdTrackCollection;
class StFwdTrackCollection;
class StFcsEventDisplay;
//class StFcsTracker;
class StFwdTrack;
class StFstHitCollection;


class StFcsMipFinderForHcal : public StMaker
  {  
    private:
        StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions
        StEvent* event = 0;                              // StEvent pointer to access event data
        StPicoEvent* picoEvent = 0;                      // StPicoEvent pointer to access event data

        TH1F* h1list_mass_by_Ntower[748];      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
        TH1F* h1list_mass_by_Stower[748];      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
        TH1F* h1list_NEtower[260];             //h1list_NEtower: energy spectrum for north Ecal tower 
        TH1F* h1list_SEtower[260];             //h1list_SEtower: energy spectrum for south Ecal tower 
        
        // TH1F* h1_two_cluster_energy_nocut = 0;   //h1_two_point_energy_nocut:2 point energy(no cut)
        TH1F* h1_hcal_clusters_per_event = 0;
        TH1F* h1_ecal_clusters_per_event = 0;
        TH1F* h1_num_towers_in_hcal_cluster = 0; //h1_hcal_num_towers_in_cluster: number of towers in hcal cluster
        TH1F* h1_num_towers_in_ecal_cluster = 0; //h1_hcal_num_towers_in_ecal_cluster: number of towers in ecal cluster
        TH1F* h1_each_hcal_cluster_energy = 0;        //h1_each_cluster_energy:each cluster energy(no cut)
        TH1F* h1_each_ecal_cluster_energy = 0;        //h1_each_ecal_cluster_energy:each cluster energy on ecal(no cut)
        TH2F* h2_hcal_cluster_position = 0;  //h2_cluster_position_nocut
        TH1F* h1_num_of_fwd_tracks = 0;
        TH1F* h1_fwd_track_is_primary = 0;
        TH1F* h1_fwd_track_num_of_fit_points = 0;
        TH1F* h1_num_of_track_matched_hcal_clusters = 0; //h1_num_of_track_matched_hcal_clusters: number of track matched hcal clusters per event
        TH1F* h1_num_of_track_matched_ecal_clusters = 0;
        TH1F* h1_track_match_cluster_hcal_energy = 0; //h1_track_match_cluster_hcal_energy: track matched cluster energy on hcal
        TH1F* h1_track_match_cluster_ecal_energy = 0; //h1_track_match_cluster_ecal_energy: track matched cluster energy on ecal
        TH1F* h1_single_cluster_id = 0; //h1_single_cluster_id: single cluster id on hcal
        TH1F* h1_num_towers_in_hcal_cluster_cut = 0; //h1_num_towers_in_hcal_cluster_tower_cut: number of towers in hcal cluster with tower cut
        TH1F* h1_num_towers_in_ecal_cluster_cut = 0; //h1_num_towers_in_ecal_cluster_tower_cut: number of towers in ecal cluster with tower cut
        TH1F* h1_track_match_cluster_hcal_energy_tower_cut = 0; //h1_track_match_cluster_hcal_energy_tower_cut: track matched cluster energy on hcal with tower cut
        TH1F* h1_track_match_cluster_ecal_energy_tower_cut = 0; //h1_track_match_cluster_ecal_energy_tower_cut: track matched cluster energy on ecal with tower cut
        TH1F* h1_hcal_hit_energy = 0; //h1_hcal_hit_energy: hcal hit energy

        // TH1F* h1_hcal_max_cluster_energy_per_event = 0;    //h1_hcal_max_cluster_energy:hcal max cluster energy(no cut)
        // TH1F* h1_Zgg_nocut_cluster = 0;          //h1_Zgg_nocut_cluster:Zgg without cut
        // TH1F* h1_num_of_points = 0;
        // TH1F* h1_inv_mass_cluster_nocut = 0;
        // TH1F* h1_ecal_clusters_per_event = 0;
        
        // TH1F* h1_hcal_neigbors_per_cluster = 0;
        TH1F* h1_fwd_cal_energy_hit = 0;
        TH1F* h1_fwd_cal_energy_cluster = 0;
        // TH1F* h1_single_cluster_ecal_energy = 0; //h1_single_cluster_ecal_energy: single cluster energy on ecal
        TH1F* h1_fwd_track_pt = 0;
        TH1F* h1_fwd_pt_res = 0;
        
        TH1F* h1_fwd_track_charge = 0;
        TH1F* h1_fwd_track_chi2 = 0;
        TH1F* h1_fwd_track_chi2_per_ndf = 0;

        TH2F* h2_geant_shower_proj_xy = 0;
        TH2F* h2_ecal_cluster_position = 0;        
        
        TH2F* h2_ecal_hcal_cluster_energy = 0;        //h2_ecal_hcal_cluster_energy
        TH2F* h2_ecal_hcal_hit_energy = 0; 
        TH2F* h2_fwd_cal_energy_vs_track_pt_hit = 0;
        TH2F* h2_fwd_cal_energy_vs_track_pt_cluster = 0;
        TH2F* h2_fwd_pt_res_vs_eta = 0;
        TH2F* h2_ecal_hcal_cluster_energy_res_vs_eta = 0;
        TH2F* h2_ecal_hcal_hit_energy_res_vs_eta = 0;

        int bins = 150;
        float m_low = 0;
        float m_up = 0.4;
        float total_event_energy = 0;
        float total_ecal_cluster_energy = 0;
        float total_hcal_cluster_energy = 0;
        float total_ecal_hit_energy = 0;
        float total_hcal_hit_energy = 0;

        TFile* mHistogramOutput;     //  Histograms outputfile pointer
        TFile* spins;

        TTree* tree;
        
        ULong_t mEventsStarted;     //  Number of Events read
        ULong_t mEventsProcessed;     //  Number of Events processed and analyzed
        TString mHistogramOutputFileName;     //  Name of the histogram output file

        StFcsDb* mFcsDb = 0;
        StFcsCollection* fcsColl = 0;
        StMuFwdTrackCollection* mftc = 0;
        StFwdTrackCollection* ftc = 0;

        float E_min = 0;
        int mDebug = 1;
        
        

    protected: 
      
    public:
        //bool makerDebug();
        StFcsMipFinderForHcal(const char* name = "hello")   ;          //  Constructor
        virtual          ~StFcsMipFinderForHcal( )   ;          //  Destructor

        Int_t Init    ( ) ;                               //  Initialize the analysis tools ... done once
        Int_t Make    ( ) ;                               //  The main analysis that is done on each event
        Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

        void loadCollectionsAndDb();
        void setDebug(int v = 0) 
            {
                mDebug = v;
            }

        void set_outputfile(const char* name)
            {
                mHistogramOutputFileName = name;
            }
            
        ClassDef(StFcsMipFinderForHcal,1)                     //  Macro for CINT compatability
  };

#endif




