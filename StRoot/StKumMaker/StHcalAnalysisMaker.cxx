#include "StHcalAnalysisMaker.h"
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
#include "StEvent/StFcsHit.h"
#include "StEvent/StFwdTrack.h"
#include "StEventTypes.h"
#include "StEventUtilities/StEventHelper.h"
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

ClassImp(StHcalAnalysisMaker)

StHcalAnalysisMaker::StHcalAnalysisMaker(const char* name) : StMaker(name)
  { // Initialize and/or zero all public/private data members here.
    for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointerOBs
      {
        h1list_mass_by_Ntower[i] = NULL;      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
        h1list_mass_by_Stower[i] = NULL;      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
        h1list_NEtower[i] = NULL;             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
        h1list_SEtower[i] = NULL;
      }
    
    mMuDstMaker  = NULL   ; //Without this, maker is not defined in scope. Why?
    mHistogramOutput  =  NULL  ;                  // Zero the pointer to histogram output file

    mEventsStarted    =  0     ;
    mEventsProcessed  =  0     ;                  // Zero the Number of Events processed by the maker

    mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the .C macro
    
  }

StHcalAnalysisMaker::~StHcalAnalysisMaker()
  {
    // Destroy and/or zero out all public/private data members here.
  } 

Int_t StHcalAnalysisMaker::Init()
  { // Do once at the start of the analysis, create histograms, etc.
    cout << termcolor::green << endl << "Initializing, please wait..." << endl << endl << termcolor::reset;
    mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); 
    if (!mFcsDb) 
      {
        cerr << termcolor::red << "StFcsEventDisplay::Init Failed to get StFcsDbMaker" << termcolor::reset << endm;
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
    h1_Zgg_nocut_cluster = new TH1F("h1_Zgg_nocut_cluster", "Zgg without cut", bins, 0, 1);
    h1_inv_mass_cluster_nocut = new TH1F("h1_inv_mass_cluster_nocut", "invariant mass plot (no cut)", bins, m_low, m_up);
    h1_fwd_track_pt = new TH1F("h1_fwd_track_pt", "forward track transverse momentum", bins, 0, 5);
    h1_track_charge = new TH1F("h1_track_charge", "forward track charge", 50, -5, 5);
    h2_cluster_position = new TH2F("h2_cluster_position", "cluster_position", 300, -150, 150, 300, -150, 150);
    
    return kStOK;
  }

Int_t StHcalAnalysisMaker::Make()
  { // Do every event
	  mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst");
    if (!mMuDstMaker) 
      return kStWarn;
    
    StMuDst* muDst = mMuDstMaker->muDst();
    if ( !muDst )
      {
        cerr << termcolor::red << "Failed to get MuDst from input" << endl << termcolor::reset;
        return kStOK ;
      }
    else
      cout  << termcolor::green << "Got MuDst from input" << endl << termcolor::reset;

    StMuEvent* muEvent =  muDst->event();
    if ( !muEvent ) 
      return kStOK ;

    StMuTriggerIdCollection* TrigMuColl = &(muEvent->triggerIdCollection());
    if(!TrigMuColl)
      {
        cerr << termcolor::red << "StFmsMaker::PopulateEventInfo - !TrigMuColl" << termcolor::reset << endl;
        return false;
      }

    
	  StEvent* event = (StEvent*)GetInputDS("StEvent");
    if (!event) 
      {
        cerr << termcolor::red << "StKumMaker::Make did not find StEvent" << endl << termcolor::reset;
        return kStErr;
      }

	  mFcsColl = event->fcsCollection();
    if (!mFcsColl) 
      return kStOK;

    mFwdColl = muDst->muFwdTrackCollection();

    if (!mFwdColl)
      {
        cerr << termcolor::bold << termcolor::red <<  "No muFwd collection" << termcolor::reset << endl;
        return kStOK;
      }
    else
      {
        cout << termcolor::green << "Number of MuFwdTracks: " << mFwdColl->numberOfFwdTracks() << termcolor::reset << endl;
      }

    int total_nc = 0;
    int total_np = 0;
    int n_EcalMult = 0;
    int n_EcalClustMult = 0;
    int n_EcalClust_cut = 0;
    int n_Ecal_cut = 0;
    
    for (size_t iTrack = 0; iTrack < mFwdColl->numberOfFwdTracks(); iTrack++)
          {
            StMuFwdTrack * muFwdTrack = mFwdColl->getFwdTrack( iTrack );
            cout << termcolor::yellow << "StMuFwdTrack[ nProjections=" << muFwdTrack->mProjections.size() 
                                      << ", nFTTSeeds=" << muFwdTrack->mFTTPoints.size() 
                                      << ", nFSTSeeds=" << muFwdTrack->mFSTPoints.size() 
                                      << ", mPt=" << muFwdTrack->momentum().Pt() 
                                      << ", charge=" << (int)muFwdTrack->charge() << " ]" << termcolor::reset << endl;
            h1_track_charge -> Fill((int)muFwdTrack->charge());
            h1_fwd_track_pt -> Fill(muFwdTrack->momentum().Pt());
          }
    
    for (int det = 2; det <= 3; det++) 
      {
        int check_Hcal = mFcsDb->ecalHcalPres(det); //Ecal North=0, Ecal South=1, Hcal North=2, Hcal South=3, Pres=4/5
        if (check_Hcal != 1) 
          { 
            cout << termcolor::red << termcolor::bold << "No Hcal found. Moving to the next" << termcolor::reset << endl;
            continue;
          }

        // no cut (points)
        StSPtrVecFcsPoint& points = mFcsColl->points(det);
        int np = mFcsColl->numberOfPoints(det);
        cout << termcolor::yellow << termcolor::bold << "Number of points: " << np << termcolor::reset << endl;
        total_np = np + total_np;

        StSPtrVecFcsHit& hits = mFcsColl->hits(det);
        int nh = mFcsColl->numberOfHits(det);
        cout << termcolor::yellow << "Number of hits: " << nh << termcolor::reset << endl;
        for (int i = 0; i < nh; i++) 
          {
            //implement cut on energy
            StFcsHit* hit = hits[i];
            unsigned short hit_id = hit->id();
            cout << termcolor::yellow << "Hit ID: " << hit_id << termcolor::reset << endl;
            float hit_energy = hit->energy();
            cout << termcolor::yellow << "Energy of hit: " << hit->energy() << " GeV" << termcolor::reset << endl;
            if (det == 2) 
              {
                cout << termcolor::green << "Filling North HCal histogram, tower " << hit_id << termcolor::reset << endl;
                h1list_NEtower[hit_id]->Fill(hit_energy);
                cout << termcolor::magenta << "Continuing" << termcolor::reset << endl;
              } 
            else if (det == 3) 
              {
                cout << termcolor::green << "Filling South HCal histogram, tower " << hit_id << termcolor::reset << endl;
                h1list_SEtower[hit_id]->Fill(hit_energy);
                cout << termcolor::magenta << "Continuing" << termcolor::reset << endl;
              }  //fill in energy spectrum for tower
            if (hit_energy > E_min)
              {
                n_Ecal_cut++;
              }
            n_EcalMult++;
          }

        StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
        int nc = mFcsColl->numberOfClusters(det);
        cout << termcolor::yellow << "Number of clusters: " << nc << termcolor::reset << endl;
        if (mDebug > 0) 
          LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;
        total_nc = nc + total_nc;
        //nc = 10;
        for (int i = 0; i < nc; i++) 
          {
            StFcsCluster* clu = clusters[i];
            //implement cut on energy
            float clu_energy = clu->energy();          
            cout << termcolor::yellow << "Energy of cluster: " << clu->energy() << " GeV" << endl;
            cout << "In detector coordinates, cluster at x = " << clu->x() << ", y = " << clu->y() << termcolor::reset << endl;
            cout << termcolor::bright_yellow << "scale factor for x: " << mFcsDb->getXWidth(det) << ", y: " << mFcsDb->getYWidth(det) << termcolor::reset << endl;
            StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu->x(), clu->y());
            StLorentzVectorD p = mFcsDb->getLorentzVector(cluPos, clu_energy, 0);
            cout << termcolor::underline << termcolor::yellow << "In physical coordinates, cluster at x = " << cluPos.x() << ", y = " << cluPos.y() << termcolor::reset << endl;
            n_EcalClustMult++;
            if (clu->energy() > E_min) 
              {
                n_EcalClust_cut++;
              }
            h2_cluster_position->Fill(cluPos.x(), cluPos.y());
            h1_each_cluster_energy->Fill(clu_energy);
            
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
    }
    return kStOK ;
  }
  

Int_t StHcalAnalysisMaker::Finish( )
  { // Do once at the end of the analysis, close files, etc.
    cout << termcolor::green << termcolor::bold << "Finishing, please wait..." << termcolor::reset << endl;
    // Write all histograms to disk
    mHistogramOutput = new TFile( "./output/output_hadron.root" , "recreate" ) ;
    h1_two_cluster_energy_nocut->Write();
    h1_each_cluster_energy->Write();
    h1_Zgg_nocut_cluster->Write();
    h1_inv_mass_cluster_nocut->Write();
    h1_fwd_track_pt->Write();
    h1_track_charge->Write();
    h2_cluster_position->Write();
    mHistogramOutput->Write();
    mHistogramOutput->Close();
    
    return kStOK;
 }
