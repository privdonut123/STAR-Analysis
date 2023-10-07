#ifndef StHcalAnalysisMaker_def
#define StHcalAnalysisMaker_def

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

#include "StEvent/StTriggerId.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "/direct/star+u/bmagh001/MuDst/include/termcolor.hpp"


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


class StHcalAnalysisMaker : public StMaker
  {  
    private:
      StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions

      TH1F* h1list_mass_by_Ntower[748];      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
      TH1F* h1list_mass_by_Stower[748];      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
      TH1F* h1list_NEtower[748];             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
      TH1F* h1list_SEtower[748];             //h1list_SEtower: energy spectrum for south Ecal tower (no cut)
      
      TH1F* h1_two_cluster_energy_nocut = 0;   //h1_two_point_energy_nocut:2 point energy(no cut)
      TH1F* h1_each_cluster_energy = 0;        //h1_each_cluster_energy:each cluster energy(no cut)
      TH1F* h1_Zgg_nocut_cluster = 0;          //h1_Zgg_nocut_cluster:Zgg without cut
      TH1F* h1_inv_mass_cluster_nocut = 0;
      TH1F* h1_fwd_track_pt = 0;
      TH1F* h1_track_charge = 0;

      TH2F* h2_cluster_position = 0;        //h2_cluster_position

      int bins = 150;
      float m_low = 0;
      float m_up = 0.4;

      TFile* mHistogramOutput;     //  Histograms outputfile pointer
      TFile* spins;

      TTree* tree;
      
      ULong_t mEventsStarted;     //  Number of Events read
      ULong_t mEventsProcessed;     //  Number of Events processed and analyzed
      TString mHistogramOutputFileName;     //  Name of the histogram output file

      StFcsDb* mFcsDb = 0;
      StFcsCollection* mFcsColl = 0;
      StMuFwdTrackCollection* mFwdColl = 0;
      StFwdTrackCollection* fwdColl = 0;

    protected: 
      
    public:
      //bool makerDebug();
      StHcalAnalysisMaker(const char* name = "hello")   ;          //  Constructor
      virtual          ~StHcalAnalysisMaker( )   ;          //  Destructor

      Int_t Init    ( ) ;                               //  Initialize the analysis tools ... done once
      Int_t Make    ( ) ;                               //  The main analysis that is done on each event
      Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

      int mDebug = 0;
      float E_min = 1;
      ClassDef(StHcalAnalysisMaker,1)                     //  Macro for CINT compatability
  };

#endif




