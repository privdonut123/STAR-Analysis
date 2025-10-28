#ifndef StHadronAnalysisMaker_def
#define StHadronAnalysisMaker_def

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


class StHadronAnalysisMaker : public StMaker
  {  
    private:
        StEvent* event = 0;                             // StEvent pointer to access event data
        StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions

        TH1F* h1list_mass_by_Ntower[748];      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
        TH1F* h1list_mass_by_Stower[748];      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
        TH1F* h1list_NEtower[748];             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
        TH1F* h1list_SEtower[748];             //h1list_SEtower: energy spectrum for south Ecal tower (no cut)
        
        TH1F* h1_two_cluster_energy_nocut = 0;   //h1_two_point_energy_nocut:2 point energy(no cut)
        TH1F* h1_each_cluster_energy = 0;        //h1_each_cluster_energy:each cluster energy(no cut)
        TH1F* h1_hcal_max_cluster_energy_per_event = 0;    //h1_hcal_max_cluster_energy:hcal max cluster energy(no cut)
        TH1F* h1_Zgg_nocut_cluster = 0;          //h1_Zgg_nocut_cluster:Zgg without cut
        TH1F* h1_num_of_points = 0;
        TH1F* h1_inv_mass_cluster_nocut = 0;
        TH1F* h1_ecal_clusters_per_event = 0;
        TH1F* h1_hcal_clusters_per_event = 0;
        TH1F* h1_hcal_neigbors_per_cluster = 0;
        TH1F* h1_fwd_cal_energy_hit = 0;
        TH1F* h1_fwd_cal_energy_cluster = 0;
        TH1F* h1_single_cluster_ecal_energy = 0; //h1_single_cluster_ecal_energy: single cluster energy on ecal
        TH1F* h1_single_cluster_hcal_energy = 0; //h1_single_cluster_hcal_energy: single cluster energy on hcal
        TH1F* h1_fwd_track_pt = 0;
        TH1F* h1_fwd_pt_res = 0;
        TH1F* h1_fwd_track_num_of_fit_points = 0;
        TH1F* h1_num_of_fwd_tracks = 0;
        TH1F* h1_fwd_track_is_primary = 0;
        TH1F* h1_fwd_track_charge = 0;
        TH1F* h1_fwd_track_chi2_per_ndf = 0;

        TH1F* h1_geant_shower_proj_z = 0;
        TH1F* h1_geant_track_id = 0;
        TH1F* h1_geant_vtx_z = 0;
        TH1F* h1_geant_start_vertex = 0;
        TH1F* h1_geant_stop_vertex = 0;
        TH1F* h1_geant_particle_id = 0;
        TH1F* h1_geant_num_of_tracks = 0;
        TH1F* h1_geant_primary_eta = 0;
        TH1F* h1_geant_primary_pt = 0;
        TH1F* h1_geant_primary_pz = 0;
        TH1F* h1_geant_primary_energy = 0;
        TH1F* h1_geant_parent_eta = 0;
        TH1F* h1_geant_parent_pt = 0;
        TH1F* h1_geant_parent_pz = 0;
        TH1F* h1_fwd_cal_hit_resolution = 0;
        TH1F* h1_fwd_cal_cluster_resolution = 0;
        TH1F* h1_matched_ecal_cluster_energy = 0;  // ECAL cluster energy for tracks matched to both ECAL and HCAL
        TH1F* h1_matched_hcal_cluster_energy = 0;  // HCAL cluster energy for tracks matched to both ECAL and HCAL
        TH1F* h1_matched_total_cluster_energy = 0; // Total cluster energy for tracks matched to both ECAL and HCAL
        TH1F* h1_track_matching_status = 0;        // Track matching status (0=no match, 1=ECAL only, 2=HCAL only, 3=both)
        TH1F* h1_all_track_fit_points = 0;         // Fit points distribution for all tracks
        TH1F* h1_track_selection_status = 0;       // Track selection status (0=rejected, 1=accepted)

        TH2F* h2_geant_shower_proj_xy = 0;
        TH2F* h2_ecal_cluster_position = 0;        
        TH2F* h2_hcal_cluster_position = 0;  //h2_cluster_position_nocut
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

        bool ecal_dead_zone = false;
        bool hcal_dead_zone = false;
        float E_min = 0;
        int mDebug = 1;
        
        // Statistics for track-cluster matching
        int mTotalTracks = 0;
        int mTracksWithBothMatches = 0;
        int mTracksWithEcalOnly = 0;
        int mTracksWithHcalOnly = 0;
        
        // Statistics for track selection criteria
        int mTotalTracksFound = 0;
        int mTracksPassingSelection = 0;
        int mTracksRejectedNotPrimary = 0;
        int mTracksRejectedFitPoints = 0;

    protected: 
      
    public:
        //bool makerDebug();
        StHadronAnalysisMaker(const char* name = "hello")   ;          //  Constructor
        virtual          ~StHadronAnalysisMaker( )   ;          //  Destructor

        Int_t Init    ( ) ;                               //  Initialize the analysis tools ... done once
        Int_t Make    ( ) ;                               //  The main analysis that is done on each event
        Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

        void setDebug(int v = 0) 
            {
                mDebug = v;
            }

        void set_outputfile(const char* name)
        {
            mHistogramOutputFileName = name;
        }
        ClassDef(StHadronAnalysisMaker,1)                     //  Macro for CINT compatability
  };

#endif




