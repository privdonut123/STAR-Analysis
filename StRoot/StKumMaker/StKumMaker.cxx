#include "StKumMaker.h"


#include <iostream>
#include <fstream>
#include <iomanip>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

//#include "StFcsPi0FinderForEcal.h"

#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
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
float DataArray[NMAX][20];
int ran_map[NMAX];

const std::string red("\033[0;31m");
const std::string green("\033[1;32m");
//const std::string yellow("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string magenta("\033[0;35m");
const std::string reset("\033[0m");

ClassImp(StKumMaker)

StKumMaker::StKumMaker(const char* name) : StMaker(name)
{ // Initialize and/or zero all public/private data members here.
  //histogram = NULL;
  //histogram[1] = NULL;
  /*
  for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointerOBs
    {
      histogram[i] = NULL    ;
    }
  */
  mMuDstMaker  = 0   ; //Without this, maker is not defined in scope. Why?
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

Int_t StKumMaker::Init( )
 { // Do once at the start of the analysis, create histograms, etc.
   cout << endl << "Initializing, please wait..." << endl << endl;
	mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
if (!mFcsDb) {
      LOG_ERROR << "StFcsEventDisplay::Init Failed to get StFcsDbMaker" << endm;
      return kStFatal;}


      for (int i = 0; i < 748; i++) {
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
  
   //mHistogramOutput = new TFile( "output.root" , "recreate" ) ;  // Name was set in "analysis".C macro
   //cout << "Creating ROOT File"  <<  endl;
   //tree = new TTree("MyTree","Example Tree");
   //tree->Branch("branch", "TH1F", &histogram, 128000, 0);
   //cout << "tree created" << endl;
   //tree->BuildIndex("Run","Event");
   //TCanvas* nCanvas[2] ;
   
   //histogram[0] = new TH1F("pvpz", "Primary Vertex Positions, z", 100, -200 ,200)  ;
   //histogram[1] = new TH1F("pvpx", "Primary Vertex Positions, x", 100, -200 ,200)  ;
   //int start_of_run_17 = 1487656902;
   //int end_of_run_17 = 1498565074;
   //histogram = new TH1F("mTime", "Run Time ", end_of_run_17 - start_of_run_17 , start_of_run_17 , end_of_run_17)  ;
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
  /*
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
    cout << endl << "Hello World" << endl << endl;
    mMuDstMaker = (StMuDstMaker*)GetInputDS("MuDst");
    if (!mMuDstMaker) 
      return kStWarn;
    //cout << "error much?" <<endl; 

    StMuDst* muDst = mMuDstMaker->muDst();
    if ( !muDst )      
      return kStOK ;

    //StMuTrack* muTrack = (StMuTrack*)GetTracks.Next(); if (!mMuDstMaker) return kStWarn;
    StMuEvent* muEvent =  muDst->event();
    if ( !muEvent ) 
      return kStOK ;

	  StEvent* event = (StEvent*)GetInputDS("StEvent");
    if (!event) 
      {
        LOG_ERROR << "\033[31mStFcsPi0FinderForEcal::Make did not find StEvent\033[m\n" << endm;
        return kStErr;
      }
    
    
    cout << "\033[31m bolder red text\033[m\n" <<endl;

	  mFcsColl = event->fcsCollection();
    if (!mFcsColl) 
      return kStOK;

    StMuTriggerIdCollection* TrigMuColl = &(muEvent->triggerIdCollection());
    if(!TrigMuColl)
      {
        LOG_ERROR <<"StFmsMaker::PopulateEventInfo - !TrigMuColl" <<endl;
        return false;
      }

    //const StTriggerId& trgIDs = TrigMuColl->nominal();
    //Int_t ntrig = trgIDs.triggerIds().size();
    //cout << "\033[31m" << ntrig << "\033[m\n" << endl;
    cout << "\033[31mWho the hell is Steve Jobs?\033[m\n" << endl;
    //for(Int_t i = 0; i < ntrig ; ++i){cout << trgIDs.triggerIds().at(i) << endl;}
    //mFmsEventInfo->mNTrig = ntrig;
   //cout << runNo << endl;
   // Do 'event' analysis based on event pointers
   //StMuEvent* muEvent = mMuDstMaker -> muDst()-> event() ;
   //cout <<  "\033[32m" << muEvent -> eventInfo().time() << "\033[m\n" << endl;
    int total_nc = 0;
    int total_np = 0;
    int n_EcalMult = 0;
    int n_EcalClustMult = 0;
    int n_EcalClust_cut = 0;
    int n_Ecal_cut = 0;
    for (int det = 0; det < 2; det++) 
      {
        int check_Ecal = mFcsDb->ecalHcalPres(det);
        if (check_Ecal != 0) continue;

        StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
        int nc = mFcsColl->numberOfClusters(det);
        cout << "Number of clusters: " << nc << endl;
        total_nc = nc + total_nc;

        StSPtrVecFcsPoint& points = mFcsColl->points(det);
        int np = mFcsColl->numberOfPoints(det);

        total_np = np + total_np;

        StSPtrVecFcsHit& hits = mFcsColl->hits(det);
        int nh = mFcsColl->numberOfHits(det);
        cout << "Number of hits: " << nh << endl;
        if (mDebug > 0) 
        LOG_INFO << Form("StFcsEventDisplay Det=%1d nhit=%4d nclu=%3d", det, nh, nc) << endm;
        for (int i = 0; i < nh; i++) 
          {
            StFcsHit* hit = hits[i];
            unsigned short hit_id = hit->id();
            cout << red << hit_id << reset << endl;
            float hit_energy = hit->energy();
            cout << "\033[32mEnergy of hit: " + to_string(hit_energy) + " GeV\033[m\n" <<  endl;
            if (det == 0) 
              {
                cout << "\033[34mError after hit loop hit loop?\033[m\n" << endl;
                if (h1list_NEtower[hit_id])
                  {
                    cout << cyan << "Filling North ECal histogram, tower " << hit_id << reset << endl;
                  }
                h1list_NEtower[hit_id]->Fill(3);
                cout << magenta << "Continuing" << reset << endl;
              } 
            else if (det == 1) 
              {
                h1list_SEtower[hit_id]->Fill(3);
              }  //fill in energy spectrum for tower
            if (hit_energy > E_min)
              {
                n_Ecal_cut++;
              }
            n_EcalMult++;
          }
            

        			//no cut (cluster)
        for (int i = 0; i < nc; i++) 
          {
            StFcsCluster* clu = clusters[i];
            float clu_energy = clu->energy();
            float clu_x = clu->x();
            float clu_y = clu->y();
            StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu_x, clu_y);
            float cluPos_x = cluPos.x();
            float cluPos_y = cluPos.y();

            n_EcalClustMult++;
            if ((clu_energy > E_min)) 
              {
                n_EcalClust_cut++;
              }
          }

      }
  
   // histogram[0] -> Fill( muEvent -> primaryVertexPosition().z() ) ;
      //   histogram[1] -> Fill( muEvent -> primaryVertexPosition().x() );
      //histogram -> Fill( muEvent -> eventInfo().time());
      //histogram -> SetAxisRange(start_time - 10, muEvent -> eventInfo().time() + 10 ,"X");
   //histogram[1] -> Fill( muEvent -> refMult() );
   // cout << typeid( muEvent -> primaryVertexPosition().z()).name() << endl;

   //components for st spin db maker

   return kStOK ;
  }
  

Int_t StKumMaker::Finish( )
 { // Do once at the end of the analysis, close files, etc.
   cout << endl << "Finishing, please wait..." << endl << endl;
   //tree->Print();
   //tree->Write();
   //mHistogramOutput->Write() ;   // Write all histograms to disk
   //mHistogramOutput->Close();
   //histogram[0] -> Draw() ; 
   return kStOK;
 }
