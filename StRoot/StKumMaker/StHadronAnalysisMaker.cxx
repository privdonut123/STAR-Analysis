#include "StHadronAnalysisMaker.h"
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

ClassImp(StHadronAnalysisMaker)

StHadronAnalysisMaker::StHadronAnalysisMaker(const char* name) : StMaker(name)
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

StHadronAnalysisMaker::~StHadronAnalysisMaker()
  {
    // Destroy and/or zero out all public/private data members here.
  } 

Int_t StHadronAnalysisMaker::Init()
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

    for (int i = 0; i < 748; i++) 
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
    h1_two_cluster_energy_nocut = new TH1F("h1_two_cluster_energy_nocut", "2 clusters energy(no cut)", bins, 0, 30);
    h1_each_cluster_energy = new TH1F("h1_each_cluster_energy", "each cluster energy for FCS ECal", bins, 0, 25);
    h1_hcal_max_cluster_energy_per_event = new TH1F("h1_hcal_max_cluster_energy_per_event", "hcal max cluster energy per event", 100, 0, 100);
    h1_hcal_max_cluster_energy_per_event->SetXTitle("Energy [GeV]");
    h1_hcal_max_cluster_energy_per_event->SetYTitle("counts");
    h1_Zgg_nocut_cluster = new TH1F("h1_Zgg_nocut_cluster", "Zgg without cut", bins, 0, 1);
    h1_inv_mass_cluster_nocut = new TH1F("h1_inv_mass_cluster_nocut", "invariant mass plot (no cut)", bins, m_low, m_up);

    h1_ecal_clusters_per_event = new TH1F("h1_ecal_clusters_per_event", "number of clusters per event on ecal", 10, -.5, 9.5);
    h1_ecal_clusters_per_event->SetXTitle("number of clusters");
    h1_ecal_clusters_per_event->SetYTitle("counts");
    h1_hcal_clusters_per_event = new TH1F("h1_hcal_clusters_per_event", "number of clusters per event on hcal", 10, -.5, 9.5);
    h1_hcal_clusters_per_event->SetXTitle("number of clusters");
    h1_hcal_clusters_per_event->SetYTitle("counts");
    h1_hcal_neigbors_per_cluster = new TH1F("h1_hcal_neigbors_per_cluster", "number of neighbors per cluster on hcal", 10, -.5, 9.5);
    h1_hcal_neigbors_per_cluster->SetXTitle("number of neighbors");
    h1_hcal_neigbors_per_cluster->SetYTitle("counts");

    h1_geant_shower_proj_z = new TH1F("h1_geant_shower_proj_z", "z projection of geant shower", bins, -100, 100);
    h1_geant_shower_proj_z->SetXTitle("z [cm]");
    h1_geant_shower_proj_z->SetYTitle("counts");
    h2_geant_shower_proj_xy = new TH2F("h2_geant_shower_proj_xy", "xy projection of geant shower", 300, -150, 150, 300, -150, 150);
    h2_geant_shower_proj_xy->SetXTitle("x [cm]");
    h2_geant_shower_proj_xy->SetYTitle("y [cm]");
    h1_geant_parent_eta = new TH1F("h1_geant_parent_eta", "eta distribution of geant parent tracks", bins, -5, 5);
    h1_geant_parent_eta->SetXTitle("eta");
    h1_geant_parent_eta->SetYTitle("counts");
    h1_geant_parent_pt = new TH1F("h1_geant_parent_pt", "transverse momentum distribution of geant parent tracks", bins, 0, 5);
    h1_geant_parent_pt->SetXTitle("p_{T} [GeV/c]");
    h1_geant_parent_pt->SetYTitle("counts");
    h1_geant_parent_pz = new TH1F("h1_geant_parent_pz", "longitudinal momentum distribution of geant parent tracks", bins, 0, 5);
    h1_geant_parent_pz->SetXTitle("p_{z} [GeV/c]");
    h1_geant_parent_pz->SetYTitle("counts");

    h1_geant_primary_eta = new TH1F("h1_geant_primary_eta", "eta distribution of primary geant primary tracks", bins, -5, 5);
    h1_geant_primary_eta->SetXTitle("eta");
    h1_geant_primary_eta->SetYTitle("counts");
    h1_geant_primary_pt = new TH1F("h1_geant_primary_pt", "transverse momentum distribution of primary geant primary tracks", bins, 0, 5);
    h1_geant_primary_pt->SetXTitle("p_{T} [GeV/c]");
    h1_geant_primary_pt->SetYTitle("counts");
    h1_geant_primary_pz = new TH1F("h1_geant_primary_pz", "longitudinal momentum distribution of primary geant primary tracks", bins, 0, 100);
    h1_geant_primary_pz->SetXTitle("p_{z} [GeV/c]");
    h1_geant_primary_pz->SetYTitle("counts");

    h1_fwd_track_pt = new TH1F("h1_fwd_track_pt", "forward track transverse momentum", bins, 0, 5);
    h1_fwd_track_pt->SetXTitle("p_{T} [GeV/c]");
    h1_fwd_track_pt->SetYTitle("counts");
    h1_track_charge = new TH1F("h1_track_charge", "forward track charge", 50, -5, 5);
    h1_track_charge->SetXTitle("charge [e]");
    h1_track_charge->SetYTitle("counts");
    h2_ecal_cluster_position = new TH2F("h2_ecal_cluster_position", "cluster_position", 300, -150, 150, 300, -150, 150);
    h2_ecal_cluster_position->SetXTitle("x [cm]");
    h2_ecal_cluster_position->SetYTitle("y [cm]");
    h2_hcal_cluster_position = new TH2F("h2_hcal_cluster_position", "cluster_position", 300, -150, 150, 300, -150, 150);
    h2_hcal_cluster_position->SetXTitle("x [cm]");
    h2_hcal_cluster_position->SetYTitle("y [cm]");
    h2_ecal_hcal_cluster_energy = new TH2F("h2_ecal_hcal_cluster_energy", "ecal_hcal_summed_cluster_energy", bins, 0, 150, bins, 0, 150);
    h2_ecal_hcal_cluster_energy->SetXTitle("Ecal summed cluster energy [GeV]");
    h2_ecal_hcal_cluster_energy->SetYTitle("Hcal summed cluster energy [GeV]");
    h2_ecal_hcal_hit_energy = new TH2F("h2_ecal_hcal_hit_energy", "ecal_hcal_summed_hit_energy", bins, 0, 150, bins, 0, 150);
    h2_ecal_hcal_hit_energy->SetXTitle("Ecal summed hit energy [GeV]");
    h2_ecal_hcal_hit_energy->SetYTitle("Hcal summed hit energy [GeV]");

    return kStOK;
  }

Int_t StHadronAnalysisMaker::Make()
  { // Do every event
	  mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst");
    if (!mMuDstMaker)
      {
        if (mDebug > 0) 
          {
            cout << termcolor::bold << termcolor::red << "Could not retrieve MuDstMaker from the chain" << termcolor::reset << endl;
          }
          return kStWarn;
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
        return kStErr;
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
        return kStErr;
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
        return kStErr;
      }
    else 
      {
        if (mDebug > 0) 
          {
            cout << termcolor::green << "Got TriggerIdCollection from MuEvent" << endl << termcolor::reset;
          }
      }

	  StEvent* event = (StEvent*)GetInputDS("StEvent");
    if (!event) 
      {
        if (mDebug > 0)
          {
            cerr << termcolor::red << "StKumMaker::Make did not find StEvent" << endl << termcolor::reset;
          }
        return kStErr;
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got StEvent from input" << endl << termcolor::reset;
          }
      }
    
    // StFstHitCollection* fsthc = event->fstHitCollection();
    // if (!fsthc)
    //   {
    //     if (mDebug > 0)
    //       {
    //         cerr << termcolor::red << "No FST collection" << termcolor::reset << endl;
    //       }
    //     return kStErr;
    //   }
    // else
    //   {
    //     if (mDebug > 0)
    //       {
    //         cout << termcolor::green << "Got FST collection" << endl << termcolor::reset;
    //       }
    //   }

	  fcsColl = event->fcsCollection();
    if (!fcsColl)
      {
        if (mDebug > 0)
          {
            cerr << termcolor::red << "No FCS collection" << termcolor::reset << endl;
          }
        return kStErr;
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got FCS collection" << endl << termcolor::reset;
            // cout << termcolor::green << "Number of FwdTracks: " << ftc->numberOfTracks() << termcolor::reset << endl;
          }
      } 
      
    mftc = muDst->muFwdTrackCollection();
    if (!mftc)
      {
        if (mDebug > 0)
          {
            cerr << termcolor::bold << termcolor::red <<  "No muFwd collection" << termcolor::reset << endl;
          }
          return kStErr;
      }
    else
      {
        if (mDebug > 0)
          {
            cout << termcolor::green << "Got muFwd collection" << endl << termcolor::reset;
          }
      }

    int total_nc = 0;
    int total_np = 0;
    
    //cout << termcolor::magenta << "Number of FST hits: " << fsthc->getNumRawHits() << termcolor::reset << endl;
  
    if(mDebug > 0)
      {   
        cout << termcolor::red << "Number of FwdTracks: " << mftc->numberOfFwdTracks() << termcolor::reset << endl;
      }
    for (size_t iTrack = 0; iTrack < mftc->numberOfFwdTracks(); iTrack++)
          {
            StMuFwdTrack * muFwdTrack = mftc->getFwdTrack( iTrack );
            cout << termcolor::yellow << "StMuFwdTrack[ nProjections=" << muFwdTrack->mProjections.size() 
                                      << ", nFTTSeeds=" << muFwdTrack->mFTTPoints.size() 
                                      << ", nFSTSeeds=" << muFwdTrack->mFSTPoints.size() 
                                      << ", mPt=" << muFwdTrack->momentum().Pt() 
                                      << ", charge=" << (int)muFwdTrack->charge() << " ]" << termcolor::reset << endl;
            h1_track_charge -> Fill((int)muFwdTrack->charge());
            h1_fwd_track_pt -> Fill(muFwdTrack->momentum().Pt());
          }
    
    total_ecal_hit_energy = 0;
    total_hcal_hit_energy = 0;      
    total_ecal_cluster_energy = 0;
    total_hcal_cluster_energy = 0;
    StThreeVectorD north_ecal_corner1 = mFcsDb->getStarXYZfromColumnRow(0, 1, 1);
    StThreeVectorD north_ecal_corner2 = mFcsDb->getStarXYZfromColumnRow(0, kFcsEcalNCol, kFcsEcalNRow);
    StThreeVectorD south_ecal_corner1 = mFcsDb->getStarXYZfromColumnRow(1, 1, 1);
    StThreeVectorD south_ecal_corner2 = mFcsDb->getStarXYZfromColumnRow(1, kFcsEcalNCol, kFcsEcalNRow);
    if(mDebug > 0)
      {
        cout << termcolor::yellow << "North ECalorimeter dimensions: " << north_ecal_corner1.x() << ", " << north_ecal_corner1.y() << ", " << north_ecal_corner1.z() << termcolor::reset << endl;
        cout << termcolor::yellow << "North ECalorimeter dimensions: " << north_ecal_corner2.x() << ", " << north_ecal_corner2.y() << ", " << north_ecal_corner2.z() << termcolor::reset << endl;
        cout << termcolor::yellow << "South ECalorimeter dimensions: " << south_ecal_corner1.x() << ", " << south_ecal_corner1.y() << ", " << south_ecal_corner1.z() << termcolor::reset << endl;
        cout << termcolor::yellow << "South ECalorimeter dimensions: " << south_ecal_corner2.x() << ", " << south_ecal_corner2.y() << ", " << south_ecal_corner2.z() << termcolor::reset << endl;
      }
    // bool dead_zone_flag = false;
    for (int det = 0; det <= 3; det++) 
      {
        if (mDebug > 0)
          {
            cout << termcolor::yellow << termcolor::bold << "Detector: " << det << termcolor::reset << endl;
          }
        //int check_Cal = mFcsDb->ecalHcalPres(det); //Ecal North=0, Ecal South=1, Hcal North=2, Hcal South=3, Pres=4/5
        /*
        need to fix this check
        if (check_Cal != 0 || check_Cal != 1) 
          { 
            cout << termcolor::red << termcolor::bold << "No forward calorimeter found. Moving to the next" << termcolor::reset << endl;
            continue;
          }
        */
        StThreeVectorD cal_corner1 = mFcsDb->getStarXYZfromColumnRow(det, 1, 1);
        StThreeVectorD cal_corner2 = mFcsDb->getStarXYZfromColumnRow(det, kFcsEcalNCol, kFcsEcalNRow);
        vector<StThreeVectorD> corners = {cal_corner1, cal_corner2};

        // Create an array of the x and y components of the corners vector
        double corner_x[corners.size()];
        double corner_y[corners.size()];
        for (unsigned int i = 0; i < corners.size(); i++) 
          {
            corner_x[i] = corners[i].x();
            corner_y[i] = corners[i].y();
          }
        std::sort(corner_x, corner_x + corners.size());
        std::sort(corner_y, corner_y + corners.size());
        
        // no cut (points)
        /*
            if (mDebug>0 && (det == 0||det == 1))
              {
                cout << termcolor::yellow << "ECalorimeter dimensions: " << ecal_corner1.x() << ", " << ecal_corner1.y() << ", " << ecal_corner1.z() << termcolor::reset << endl;
                cout << termcolor::yellow << "ECalorimeter dimensions: " << ecal_corner2.x() << ", " << ecal_corner2.y() << ", " << ecal_corner2.z() << termcolor::reset << endl;
              }
              */
        StSPtrVecFcsPoint& points = fcsColl->points(det);
        int np = fcsColl->numberOfPoints(det);
        if (mDebug > 0) 
          {
            cout << termcolor::yellow << termcolor::bold << "Number of points: " << np << termcolor::reset << endl;
          }
          total_np = np + total_np;
        for( int ipoint=0; ipoint<np; ++ipoint )
          {
            StFcsPoint* point=points[ipoint];
            StFcsCluster* pointclus = point->cluster();
            St_g2t_track* trackTable = static_cast<St_g2t_track*>(GetDataSet("g2t_track"));
            //St_g2t_vertex* vertexTable = static_cast<St_g2t_vertex*>(GetDataSet("g2t_vertex"));
            g2t_track_st* g2ttrk = 0;
            //g2t_vertex_st* g2tvert = 0;
            float frac=0;
            int ntrk=0;
            if( !trackTable )
              { 
                cout << termcolor::red << "g2t_track Table not found" << termcolor::reset << std::endl; 
                continue; 
              }
            else
              {
                const int nTrk = trackTable->GetNRows();
                if( mDebug>0 )
                  { 
                    cout << termcolor::green << "g2t_track table has "<< nTrk << " tracks" << termcolor::reset << std::endl; 
                  }
                if( nTrk>0 )
                  {
                    g2ttrk = trackTable->GetTable();
                    if(mDebug > 0)
                      {
                        if( !g2ttrk ) 
                          { 
                            cout << termcolor::red << " g2t_track GetTable failed" << termcolor::reset << endl;
                            continue; 
                          }
                      }
                  }
              }
            // if( !vertexTable )
            //   { 
            //     std::cout<< "g2t_vertex Table not found" << std::endl; continue; 
            //   }
	          // else
            //   {
            //     const int nVertex = vertexTable->GetNRows();
            //     if( GetDebug()>0 )
            //       {
            //         std::cout << "g2t_vertex table has "<< nVertex << " vertices" << std::endl; 
            //       }
            //     if( nVertex>0 )
            //       {
            //         g2tvert = vertexTable->GetTable();
            //         if( !g2tvert) 
            //           { 
            //             std::cout << " g2t_vertex GetTable failed" << std::endl; continue; 
            //           }
	          //       }
	          //   }
            //std::cout << "|parenttrk|Id:"<<parenttrk->id << "|Pid:"<<parenttrk->ge_pid << "|E:"<<parenttrk->e << "|eta:"<<parenttrk->eta << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
            const g2t_track_st* primtrk = mFcsDb->getPrimaryG2tTrack(pointclus,g2ttrk,frac,ntrk);
            StThreeVectorD projshowerxyz1 = mFcsDb->projectTrackToEcal(primtrk);

            
            // bool ecal_dead_zone = ( projshowerxyz1.x() < corner_x[0] || projshowerxyz1.x() > corner_x[corners.size()-1] ||
            //             projshowerxyz1.y() < corner_y[0] || projshowerxyz1.y() > corner_y[corners.size()-1]  ); 
            if (mDebug > 0)
            {
              cout << "projshowerxyz1.x(): " << projshowerxyz1.x() << endl;
              cout << "projshowerxyz1.y(): " << projshowerxyz1.y() << endl;
              cout << "projshowerxyz1.z(): " << projshowerxyz1.z() << endl;
              cout << "corner_x[0]: " << corner_x[0] << endl;
              cout << "corner_x[1]: " << corner_x[1] << endl;
              cout << "corner_y[0]: " << corner_y[0] << endl;
              cout << "corner_y[1]: " << corner_y[1] << endl;
              cout << "projshowerxyz1.x() < corner_x[0]: " << (projshowerxyz1.x() < corner_x[0]) << endl;
              cout << "projshowerxyz1.x() > corner_x[1]: " << (projshowerxyz1.x() > corner_x[1]) << endl;
              cout << "projshowerxyz1.y() < corner_y[0]: " << (projshowerxyz1.y() < corner_y[0]) << endl;
              cout << "projshowerxyz1.y() > corner_y[1]: " << (projshowerxyz1.y() > corner_y[1]) << endl;
              cout << termcolor::blue << "Ecal geant shower projection: " << projshowerxyz1.x() << ", " << projshowerxyz1.y() << ", " << projshowerxyz1.z() << termcolor::reset << endl;  
            }

            // if (ecal_dead_zone && mDebug > 0)
            //   {
            //     cout << termcolor::red << "ECal dead zone tripped!" << termcolor::reset << endl;
            //     dead_zone_flag = true;
            //     continue;
            //   }
              
            h1_geant_shower_proj_z->Fill(projshowerxyz1.z());
            h2_geant_shower_proj_xy->Fill(projshowerxyz1.x(), projshowerxyz1.y());
            //StThreeVectorD projshowerxyz2 = mFcsDb->projectTrackToEcalSMax(primtrk,g2tvert);
            //std::cout << "|picotrk:"<<picotrk << "|primtrk:"<<primtrk << std::endl;
            h1_geant_primary_eta->Fill(primtrk->eta);
            h1_geant_primary_pt->Fill(sqrt(primtrk->p[0]*primtrk->p[0] + primtrk->p[1]*primtrk->p[1]));
            h1_geant_primary_pz->Fill(primtrk->p[2]);

            const g2t_track_st* parenttrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
            //StThreeVectorD projparentxyz = mFcsDb->projectTrackToEcalSMax(parenttrk,g2tvert);

            h1_geant_parent_eta->Fill(parenttrk->eta);
            h1_geant_parent_pt->Fill(sqrt(parenttrk->p[0]*parenttrk->p[0] + parenttrk->p[1]*parenttrk->p[1]));
            h1_geant_parent_pz->Fill(parenttrk->p[2]);
          }

        // if (dead_zone_flag && mDebug > 0)
        //   {
        //     cout << termcolor::red << "Dead zone tripped! Moving to next detector" << termcolor::reset << endl;
        //     continue;
        //   }



        StSPtrVecFcsHit& hits = fcsColl->hits(det);
        int nh = fcsColl->numberOfHits(det);
        if (mDebug > 0)
          {
            cout << termcolor::yellow << "Number of hits: " << nh << termcolor::reset << endl;
          }
        for (int i = 0; i < nh; i++) 
          {
            //implement cut on energy
            StFcsHit* hit = hits[i];
            unsigned short hit_id = hit->id();
            if (hit->energy() > 0)
              {
                if (mDebug > 0) 
                  {
                    cout << termcolor::yellow << "Hit ID: " << hit_id << termcolor::reset << endl;
                  }
                float hit_energy = hit->energy();
                if (mDebug > 0)
                  {
                  cout << termcolor::yellow << "Energy of hit: " << hit->energy() << " GeV" << termcolor::reset << endl;
                  }
                if(det == 0 || det == 1)
                  {
                    total_ecal_hit_energy = total_ecal_hit_energy + hit_energy;
                  }
                else if(det == 2 || det == 3)
                  {
                    total_hcal_hit_energy = total_hcal_hit_energy + hit_energy;
                  }
            /*
            if (det == 2) 
              {
                if (mDebug > 0) {
                cout << termcolor::green << "Filling North HCal histogram, tower " << hit_id << termcolor::reset << endl;
                }
                h1list_NEtower[hit_id]->Fill(hit_energy);
                if (mDebug > 0) {
                cout << termcolor::magenta << "Continuing" << termcolor::reset << endl;
                }
              } 
            else if (det == 3) 
              {
                if (mDebug > 0) {
                cout << termcolor::green << "Filling South HCal histogram, tower " << hit_id << termcolor::reset << endl;
                }
                h1list_SEtower[hit_id]->Fill(hit_energy);
                if (mDebug > 0) {
                cout << termcolor::magenta << "Continuing" << termcolor::reset << endl;
                }
              }  //fill in energy spectrum for tower
              */
              }
          }
        StSPtrVecFcsCluster& clusters = fcsColl->clusters(det);
        int nc = fcsColl->numberOfClusters(det);
        
        
        if (mDebug > 0) 
          {
            cout << termcolor::yellow << "Number of clusters: " << nc << termcolor::reset << endl;
          }
        if (mDebug > 0) 
          LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;
        if(det == 0 || det == 1)
          {
            h1_ecal_clusters_per_event->Fill(nc);
          }
        else if(det == 2 || det == 3)
          {
            h1_hcal_clusters_per_event->Fill(nc);
          }
        total_nc = nc + total_nc;
        //nc = 10;
        float maxEnergy = -1; // Initialize max energy to a very small value

        for (int i = 0; i < nc; i++) 
          {
            StFcsCluster* clu = clusters[i];
            if(det == 2 || det == 3)
              {
                h1_hcal_neigbors_per_cluster->Fill(clu->nNeighbor());
              }
            
            //implement cut on energy
            float clu_energy = clu->energy();
            if (mDebug > 0)
              {          
                cout << termcolor::yellow << "Energy of cluster: " << clu->energy() << " GeV" << endl;
                cout << "In detector coordinates, cluster at x = " << clu->x() << ", y = " << clu->y() << termcolor::reset << endl;
                cout << termcolor::bright_yellow << "scale factor for x: " << mFcsDb->getXWidth(det) << ", y: " << mFcsDb->getYWidth(det) << termcolor::reset << endl;
              }
            StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu->x(), clu->y());
            StLorentzVectorD p = mFcsDb->getLorentzVector(cluPos, clu_energy, 0);
            if (mDebug > 0) 
              {
                cout << termcolor::underline << termcolor::yellow << "In physical coordinates, cluster at x = " << cluPos.x() << ", y = " << cluPos.y() << termcolor::reset << endl;
              }
            if (clu->energy() > 0)
              {
                if (clu_energy > maxEnergy) 
                  {
                    maxEnergy = clu_energy; // Update max energy
                  }
                h1_each_cluster_energy->Fill(clu_energy);
                if(det == 0 || det == 1)
                  {
                    h2_ecal_cluster_position->Fill(cluPos.x(), cluPos.y());
                    total_ecal_cluster_energy = total_ecal_cluster_energy + clu_energy;
                  }
                else if(det == 2 || det == 3)
                  {
                    h2_hcal_cluster_position->Fill(cluPos.x(), cluPos.y());
                    total_hcal_cluster_energy = total_hcal_cluster_energy + clu_energy;
                  }
              }
            if (i == nc - 1) 
              continue;
            for (int j = i + 1; j < nc; j++) 
              {
                StFcsCluster* cluj = clusters[j];
                float cluj_energy = cluj->energy();
                float cluj_x = cluj->x();
                float cluj_y = cluj->y();
                StThreeVectorD clujPos = mFcsDb->getStarXYZfromColumnRow(det, cluj_x, cluj_y);

                h1_two_cluster_energy_nocut->Fill(clu_energy + cluj_energy);
                float zgg = (abs(clu_energy - cluj_energy)) / (clu_energy + cluj_energy);
                h1_Zgg_nocut_cluster->Fill(zgg);
                StThreeVectorD xyzj = mFcsDb->getStarXYZfromColumnRow(det, cluj->x(), cluj->y());
                StLorentzVectorD pj = mFcsDb->getLorentzVector((xyzj), cluj->energy(), 0);
                h1_inv_mass_cluster_nocut->Fill((p + pj).m());
              }
          }
        if (mDebug > 0) 
          {
            cout << termcolor::red << "Finishing detector " << det << termcolor::reset << endl;
          }
        if (maxEnergy > 0) 
          {
            if(det == 2 || det == 3)
              {
                h1_hcal_max_cluster_energy_per_event->Fill(maxEnergy);
              }
          }
      } 

    
    //cut on energy, E_min
    if (total_ecal_hit_energy > E_min && total_hcal_hit_energy > E_min) 
      {
        h2_ecal_hcal_hit_energy->Fill(total_ecal_hit_energy, total_hcal_hit_energy);
      }
    if (total_ecal_cluster_energy > E_min && total_hcal_cluster_energy > E_min) 
      {
        h2_ecal_hcal_cluster_energy->Fill(total_ecal_cluster_energy, total_hcal_cluster_energy);
      }
    //dead_zone_flag = false;
    return kStOK ;
  }
  

Int_t StHadronAnalysisMaker::Finish( )
  { // Do once at the end of the analysis, close files, etc.
    cout << termcolor::green << termcolor::bold << "Finishing, please wait..." << termcolor::reset << endl;
    // Write all histograms to disk
    mHistogramOutput = new TFile( mHistogramOutputFileName.Data() , "recreate" ) ;
    h1_two_cluster_energy_nocut->Write();
    h1_each_cluster_energy->Write();
    h1_hcal_max_cluster_energy_per_event->Write();
    h1_Zgg_nocut_cluster->Write();
    h1_inv_mass_cluster_nocut->Write();
    h1_ecal_clusters_per_event->Write();
    h1_hcal_clusters_per_event->Write();
    h1_hcal_neigbors_per_cluster->Write();
    h1_fwd_track_pt->Write();
    h1_track_charge->Write();
    h1_geant_shower_proj_z->Write();
    h2_geant_shower_proj_xy->Write();
    h1_geant_parent_eta->Write();
    h1_geant_parent_pt->Write();
    h1_geant_parent_pz->Write();
    h1_geant_primary_eta->Write();
    h1_geant_primary_pt->Write();
    h1_geant_primary_pz->Write();
    h2_ecal_cluster_position->Write();
    h2_hcal_cluster_position->Write();
    h2_ecal_hcal_cluster_energy->Write();
    h2_ecal_hcal_hit_energy->Write();
    mHistogramOutput->Write();
    mHistogramOutput->Close();
    
    return kStOK;
 }
