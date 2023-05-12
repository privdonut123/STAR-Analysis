#ifndef StKumMaker_def
#define StKumMaker_def

#include "StMaker.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <TTree.h>
#include "StEvent/StTriggerId.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StFcsDbMaker/StFcsDb.h"
 

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
class StFwdTrackCollection;
class StFcsEventDisplay;
//class StFcsTracker;

#define MaxNumberOfTH1F      100
#define MaxNumberOfTH2F      100
#define MaxNumberOfTH3F      100
#define MaxNumberOfTProfile  200



class StKumMaker : public StMaker
  {  
    private:
      int start_time;
      //  StEvent* evnt;
      StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions
      //TH1F*         histogram   ;     //  1D Histograms
      //TH1F*         mH1_PolBlue;//[MaxNumberOfTH1F];
      //TH1F*         mH1_PolYellow;//[MaxNumberOfTH1F];
      TH1F* h1list_mass_by_Ntower[748];      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
      TH1F* h1list_mass_by_Stower[748];      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
      TH1F* h1list_NEtower[748];             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
      TH1F* h1list_SEtower[748];             //h1list_SEtower: energy spectrum for south Ecal tower (no cut)

      TH1F* h1_two_cluster_energy_nocut = 0;   //h1_two_point_energy_nocut:2 point energy(no cut)
      TH1F* h1_each_cluster_energy = 0;        //h1_each_cluster_energy:each cluster energy(no cut)
      TH1F* h1_Zgg_nocut_cluster = 0;          //h1_Zgg_nocut_cluster:Zgg without cut
      TH1F* h1_inv_mass_cluster_nocut = 0;

      TH2F* h2_cluster_position = 0;        //h2_cluster_position

      int bins = 150;
      float m_low = 0;
      float m_up = 0.4;
      //TH2F*         histogram2D[MaxNumberOfTH2F] ;     //  2D Histograms
      //TH3F*         histogram3D[MaxNumberOfTH3F] ;     //  3D Histograms
      
      //TProfile*     histogramTProfile[MaxNumberOfTProfile] ;     //  Profile Histograms
      // TProfile2D*   histogramTProfile2D[MaxNumberOfTProfile] ;     //  Profile Histograms
      
      TFile* mHistogramOutput;     //  Histograms outputfile pointer
      TFile* spins;

      TTree* tree;
      
      //int mRootS;
      
      ULong_t mEventsStarted;     //  Number of Events read
      ULong_t mEventsProcessed;     //  Number of Events processed and analyzed
      TString mHistogramOutputFileName;     //  Name of the histogram output file
      /*
      int mRootS;
      StRefMultCorr* mRefmultCorrUtil; //check out the *latest* most relevant version of StRefMultCorr with cvs co -r SLXXX StRoot/StRefMultCorr
      
      ULong_t       mEventsStarted               ;     //  Number of Events read
      ULong_t       mEventsProcessed             ;     //  Number of Events processed and analyzed
      TString       mHistogramOutputFileName     ;     //  Name of the histogram output file 
      */
      StFcsDb* mFcsDb = 0;
      StFcsCollection* mFcsColl = 0;
      StFwdTrackCollection* mFwdColl = 0;
    
    protected: 
      //StSpinDbMaker *mStSpinDbMaker;

    public:

      StKumMaker(const char* name = "hello")   ;          //  Constructor
      virtual          ~StKumMaker( )   ;          //  Destructor
      Int_t Init    ( ) ;                               //  Initialize the analysis tools ... done once
      //Int_t InitRun (int runNo);
      //Int_t FinishRun(int oldRunNo);
      Int_t Make    ( ) ;                               //  The main analysis that is done on each event
      Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

      //void SetOutputFileName(TString name) { mHistogramOutputFileName = name ; } // Make name available to member functions
      int mDebug = 0;
      float E_min = 1;
      ClassDef(StKumMaker,1)                     //  Macro for CINT compatability
  };

#endif




