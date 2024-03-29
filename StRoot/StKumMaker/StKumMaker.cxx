#include "StKumMaker.h"
#include "StFwdTrackMaker/Common.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

//#include "StFcsPi0FinderForEcal.h"

#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFwdTrack.h"
#include "StEventTypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"

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
#include <StRunInfo.h>
#include <StEventInfo.h>
#include <StHelix.hh>
#include "StThreeVector.hh"
#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"

#include <string>




//#include "StFwdTrackMaker/Common.h"

#include "TMath.h"

#include <limits>
//#include <map>
#include <string>
#include <string>
#include <vector>

#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
//#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"


#include "TROOT.h"
#include "TLorentzVector.h"
//#include "StEvent/StFwdTrack.h"




















/*
//From fcspi0

#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEventTypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StSpinPool/StFcsQaMaker/StFcsQaMaker.h"
#include "StSpinPool/StFcsRawDaqReader/StFcsRawDaqReader.h"
#include "StThreeVectorF.hh"
#include "Stypes.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMarker.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"

*/




#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"


#define NumberOfTH1F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH2F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTH3F      50                    // Number of Histograms (must be less than max in .h file)
#define NumberOfTProfile  100                    // Number of Histograms (must be less than max in .h file)

#define MinTracks         0                    // Discard events with fewer than MinTracks (Flow with 2 tracks makes no sense :-)
#define MinEvents         0                    // Discard DST data when we have fewer than MinEvents on a file sequence
#define pi                TMath::Pi()

#define NMAX         10000
//float DataArray[NMAX][20];
//int ran_map[NMAX];
/*
const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
//const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");
*/
#include <stdio.h>

/*
void purple() {
  printf("\033[0;35m");
}

void reset_1() {
  printf("\033[0m");
}
*/
ClassImp(StKumMaker)

StKumMaker::StKumMaker(const char* name) : StMaker(name)
  { // Initialize and/or zero all public/private data members here.
    //histogram = NULL;
    //histogram[1] = NULL;
    
    for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointerOBs
      {
        h1list_mass_by_Ntower[i] = NULL;      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
        h1list_mass_by_Stower[i] = NULL;      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
        h1list_NEtower[i] = NULL;             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
        h1list_SEtower[i] = NULL;
      }
    
    mMuDstMaker  = NULL   ; //Without this, maker is not defined in scope. Why?
    mHistogramOutput  =  NULL  ;                  // Zero the pointer to histogram output file
    //evnt = NULL ;

    //mMuDstMaker       =  maker ;                  // Pass MuDst pointer to AnlysisMaker class member functions

    mEventsStarted    =  0     ;
    mEventsProcessed  =  0     ;                  // Zero the Number of Events processed by the maker

    mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the .C macro
    
  }
StKumMaker::~StKumMaker()
  {
    // Destroy and/or zero out all public/private data members here.
  } 

//bool StKumMaker::makerDebug()
  map <string,bool> makerDebug = 
    {
      {"fwdTrk", true},
      {"fcs", true}
    };
Int_t StKumMaker::Init()
 { 
  //makerDebug();
  // Do once at the start of the analysis, create histograms, etc.
    green(); bold(); printf("StKumMaker::Init() \n"); reset();
    LOG_INFO << endl << "Initializing, please wait..." << endl << endl;
	  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); 
    if (!mFcsDb) 
      {
        LOG_ERROR << "StFcsEventDisplay::Init Failed to get StFcsDbMaker" << endm;
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

  
  //mHistogramOutput = new TFile( "output.root" , "recreate" ) ;  // Name was set in "analysis".C macro
  //cout << "Creating ROOT File"  <<  endl;
  /*
  tree = new TTree("MyTree","Example Tree");
  tree->Branch("branch", "TH1F", &histogram, 128000, 0);
  cout << "tree created" << endl;
  tree->BuildIndex("Run","Event");
  TCanvas* nCanvas[2] ;
  */
    return kStOK;
  }
/*
 Int_t StKumMaker::InitRun(int runNo) {
   //gSystem->Exec("mkdir -p spins");
   //gSystem->cd("spins");
   //gSystem->Exec("pwd");
   //string fileName = to_string(runNo) + "_spins.root";
   //spins = new TFile( fileName.c_str() , "recreate" ) ;
   cout << " Updating Root file for " << runNo << endl;
   mH1_PolBlue   = new TH1F("pol_Blue",   "blue beam polarization",   120, -0.5, 119.5);
   mH1_PolYellow = new TH1F("pol_Yellow", "yellow beam polarization", 120, -0.5, 119.5);
   //histogram[2] = new TH1F("mTime", "Run Time ",10000 , 1491860037 ,1491861855)  ;

   cout << runNo << endl;
   // cout << "hey, hows it going" << endl;
  // Reading the polarization pattern
  mStSpinDbMaker = static_cast<StSpinDbMaker*>(GetMaker("spinDb"));
  mStSpinDbMaker->print(0);
  //cout << "testing, 1,2,3" << endl;
  
  for( int i=0; i<120; i++ ) {
    
    int bunch = mStSpinDbMaker->BXstarUsingBX48(i);
    int spin4 = mStSpinDbMaker->spin4usingBX48(i);

    int bpol = 0;
    int ypol = 0;
    if( spin4 & 0x1 ) ypol = +1;
    if( spin4 & 0x2 ) ypol = -1;
    if( spin4 & 0x4 ) bpol = +1;
    if( spin4 & 0x8 ) bpol = -1;

    //cout << ypol << endl;
    

    mH1_PolBlue->SetBinContent( bunch+1, bpol);
    //cout << "hello" << endl;
    mH1_PolYellow->SetBinContent( bunch+1, ypol);
  }
// 	
  mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst"); if (!mMuDstMaker) return kStWarn;

   StMuDst*   muDst   =  mMuDstMaker->muDst() ;  if ( !muDst )      return kStOK ;
   //StMuTrack* muTrack = (StMuTrack*)GetTracks.Next(); if (!mMuDstMaker) return kStWarn;
   StMuEvent* muEvent =  muDst->event() ;  if ( !muEvent )    return kStOK ;
   //histogram[2] -> Fill( muEvent -> eventInfo().time());
   start_time =  muEvent -> eventInfo().time();
   
   spins->Write() ;   // Write all histograms to disk
   spins->Close();
   //gSystem->Exec("cd ..");

   
   return kStOK ;
 }

Int_t StKumMaker::FinishRun(int oldRunNo)
{
  return kStOK;
}
*/
Int_t StKumMaker::Make()
  { // Do every event
    cout << endl << "Loading in makers:" << endl << endl;

    mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst");
    if (!mMuDstMaker) 
      return kStWarn;
    
    StMuDst* muDst = mMuDstMaker->muDst();
    if ( !muDst )
      {
        red(); printf("Failed to get MuDst from input\n"); reset();
        return kStOK ;
      }
    else
      red(); printf("Got MuDst from input\n"); reset();

    StMuEvent* muEvent =  muDst->event();
    if ( !muEvent ) 
      return kStOK ;

    StMuTriggerIdCollection* TrigMuColl = &(muEvent->triggerIdCollection());
    if(!TrigMuColl)
      {
        LOG_ERROR <<"StFmsMaker::PopulateEventInfo - !TrigMuColl" <<endl;
        return false;
      }

    
	  StEvent* event = (StEvent*)GetInputDS("StEvent");
    if (!event) 
      {
        LOG_ERROR << "\033[31mStKumMaker::Make did not find StEvent\033[m\n" << endm;
        return kStErr;
      }

	  mFcsColl = event->fcsCollection();
    if (!mFcsColl) 
      return kStOK;

//end block
//makerDebug().fwdTrk << endl;
    
    
    if(makerDebug["fwdTrk"])
      {
        //fwdColl = stEvent->fwdTrackCollection();
        mFwdColl = muDst->muFwdTrackCollection();

        if (!mFwdColl)
          {
            cout << "No muFwd collection" << endl;
            return kStOK;
          }
        else
          {
            green(); bold(); printf("Number of MuFwdTracks: %i\n", mFwdColl->numberOfFwdTracks()); reset();
          } 
      }
    /*  
    LOG_INFO << "Checking FcsCollection" << endm;
    StFcsCollection *fcs = stEvent->fcsCollection();
    if (!fcs) 
      cout << "No FcsCollection" << endl;
    */

    int total_nc = 0;
    int total_np = 0;
    int n_EcalMult = 0;
    int n_EcalClustMult = 0;
    int n_EcalClust_cut = 0;
    int n_Ecal_cut = 0;
    
    for (size_t iTrack = 0; iTrack < mFwdColl->numberOfFwdTracks(); iTrack++)
          {
            StMuFwdTrack * muFwdTrack = mFwdColl->getFwdTrack( iTrack );
            purple();
            printf("Track has a charge of: %x e\n", muFwdTrack->charge());
            reset();
            LOG_INFO << TString::Format("StMuFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f, charge=%d ]",
                                         muFwdTrack->mProjections.size(),
                                         muFwdTrack->mFTTPoints.size(),
                                         muFwdTrack->mFSTPoints.size(),
                                         muFwdTrack->momentum().Pt(),
                                         (int)muFwdTrack->charge()) << endm;
            //cout << "Charge is " << (float)muFwdTrack->charge() << endl;
            h1_track_charge -> Fill((int)muFwdTrack->charge());
            h1_fwd_track_pt -> Fill(muFwdTrack->momentum().Pt());
          }
    
    for (int det = 0; det < 2; det++) 
      {
        int check_Ecal = mFcsDb->ecalHcalPres(det); //Ecal North=0, Ecal South=1, Hcal North=2, Hcal South=3, Pres=4/5
        if (check_Ecal != 0) 
          { 
            red(); bold(); printf("No Ecal found. Moving to the next \n"); reset();
            continue;
          }

        
        /*
        int nt = mFwdColl->numberOfTracks();
        cout << green << "Number of tracks: " << nt << reset << endl;
        for(int i = 0; i < nt; i++)
          {
            StFwdTrack* track = tracks[i];
            cout << green << track->charge() << reset << endl;
            //StThreeVectorD 	mome
          }
        */
            //unsigned short track_id = hit->id();
            //cout << red << hit_id << reset << endl;
            //float hit_energy = track->;

        // Points are ignored for this section
        //if(makerDebug["fcs"])
          //{ 
            StSPtrVecFcsPoint& points = mFcsColl->points(det);
            int np = mFcsColl->numberOfPoints(det);
            //cout << green << "Number of points: " << np << reset << endl;
            green(); bold(); printf("Number of points: %i\n", np); reset();
            total_np = np + total_np;

            // no cut (hits)
            StSPtrVecFcsHit& hits = mFcsColl->hits(det);
            int nh = mFcsColl->numberOfHits(det);
            //cout << green << "Number of hits: " << nh << reset << endl;
            green(); printf("Number of hits: %i\n", nh); reset();
            for (int i = 0; i < 5/*nh*/; i++) 
              {
                //implement cut on energy
                StFcsHit* hit = hits[i];
                unsigned short hit_id = hit->id();
                //cout << red << hit_id << reset << endl;
                red(); printf("Hit ID: %i\n", hit_id); reset();
                float hit_energy = hit->energy();
                //cout << green << "Energy of hit: " << hit->energy() << " GeV" <<  endl;
                green(); printf("Energy of hit: %f GeV\n", hit->energy()); reset();
                if (det == 0) 
                  {
                    //cout << cyan << "Filling North ECal histogram, tower " << hit_id << reset << endl;
                    cyan(); printf("Filling North ECal histogram, tower %hu\n", hit_id); reset();
                    h1list_NEtower[hit_id]->Fill(hit_energy);
                    //cout << magenta << "Continuing" << reset << endl;
                    magenta(); printf("Continuing\n"); reset();
                  } 
                else if (det == 1) 
                  {
                    //cout << cyan << "Filling South ECal histogram, tower " << hit_id << reset << endl;
                    cyan(); printf("Filling South ECal histogram, tower %hu\n", hit_id); reset();
                    h1list_SEtower[hit_id]->Fill(hit_energy);
                    magenta(); printf("Continuing\n"); reset();
                    //cout << magenta << "Continuing" << reset << endl;
                  }  //fill in energy spectrum for tower
                if (hit_energy > E_min)
                  {
                    n_Ecal_cut++;
                  }
                n_EcalMult++;
              }

            //no cut (cluster)    
            StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
            int nc = mFcsColl->numberOfClusters(det);
            //cout << green << "Number of clusters: " << nc << reset << endl;
            green(); printf("Number of clusters: %i\n", nc); reset();
            if (mDebug > 0) 
              LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;
            total_nc = nc + total_nc;
            //nc = 10;
            for (int i = 0; i < nc; i++) 
              {
                StFcsCluster* clu = clusters[i];
                //cout << red << "ISSUE HERE?" << endl;
                //implement cut on energy
                float clu_energy = clu->energy();
                //cout << green << "Energy of cluster: " << clu->energy() << " GeV" << reset << endl;
                green(); printf("Energy of cluster: %f GeV\n", clu->energy()); reset();
                //cout << green << "In detector coordinates, cluster at x = " << clu->x() << ", y = " << clu->y() << reset << endl;
                green(); printf("In detector coordinates, cluster at x = %f, y = %f\n", clu->x(), clu->y()); reset(); 
                //cout << magenta << "scale factor for x: " << mFcsDb->getXWidth(det) << ", y: " << mFcsDb->getYWidth(det) << reset << endl;
                magenta(); printf("scale factor for x: %f, y: %f\n", mFcsDb->getXWidth(det), mFcsDb->getYWidth(det)); reset();
                StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu->x(), clu->y());
                StLorentzVectorD p = mFcsDb->getLorentzVector(cluPos, clu_energy, 0);
                float cluPos_x = cluPos.x();
                float cluPos_y = cluPos.y();
                //cout << cyan << "In physical coordinates, cluster at x = " << cluPos.x() << ", y = " << cluPos.y() << reset << endl;
                cyan(); printf("In physical coordinates, cluster at x = %f, y = %f\n", cluPos.x(), cluPos.y()); reset();
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
          //}
      //if (mDebug > 0)
        //LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;

    }

  //histogram[0] -> Fill( muEvent -> primaryVertexPosition().z() ) ;
  //histogram[1] -> Fill( muEvent -> primaryVertexPosition().x() );
  //histogram -> Fill( muEvent -> eventInfo().time());
  //histogram -> SetAxisRange(start_time - 10, muEvent -> eventInfo().time() + 10 ,"X");
  //histogram[1] -> Fill( muEvent -> refMult() );
  // cout << typeid( muEvent -> primaryVertexPosition().z()).name() << endl;

  //components for st spin db maker

   return kStOK ;
  }
  

Int_t StKumMaker::Finish( )
 { // Do once at the end of the analysis, close files, etc.
   //cout << endl << "Finishing, please wait..." << endl << endl;
   green(); bold(); printf("\nFinishing, please wait...\n\n"); reset();
   //tree->Print();
   //tree->Write();
  mHistogramOutput = new TFile( "./output/output_MuDst.root" , "recreate" ) ;
  h1_two_cluster_energy_nocut->Write();
  h1_each_cluster_energy->Write();
  h1_Zgg_nocut_cluster->Write();
  h1_inv_mass_cluster_nocut->Write();
  h1_fwd_track_pt->Write();
  h1_track_charge->Write();

  h2_cluster_position->Write();
  mHistogramOutput->Write() ;   // Write all histograms to disk
  mHistogramOutput->Close();
   //histogram[0] -> Draw() ; 
   return kStOK;
 }
