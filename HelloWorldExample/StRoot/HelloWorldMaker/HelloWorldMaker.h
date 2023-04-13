// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#ifndef HelloWorldMaker_def
#define HelloWorldMaker_def

#include "StMaker.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <TTree.h>

#include "StThreeVectorF.hh"

#include "StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h" 

class StMuDstMaker ;
class StMuEvent    ;
class StMuTrack    ;

class TFile        ;
class TH1F         ;
class TH2F         ;
class TH3F         ;
class TProfile     ;
class TProfile2D   ;
class TRandom      ;
class TRandom3      ;

#define MaxNumberOfTH1F      200
#define MaxNumberOfTH2F      100
#define MaxNumberOfTH3F      100
#define MaxNumberOfTProfile  200

class HelloWorldMaker : public StMaker
{
  
 private:

  //StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions
  /*
  TH1F*         histogram[MaxNumberOfTH1F]   ;     //  1D Histograms
  TH2F*         histogram2D[MaxNumberOfTH2F] ;     //  2D Histograms
  TH3F*         histogram3D[MaxNumberOfTH3F] ;     //  3D Histograms
  TProfile*     histogramTProfile[MaxNumberOfTProfile] ;     //  Profile Histograms
  TProfile2D*   histogramTProfile2D[MaxNumberOfTProfile] ;     //  Profile Histograms
  

  TFile*        mHistogramOutput             ;     //  Histograms outputfile pointer
 
  int mRootS;
  StRefMultCorr* mRefmultCorrUtil; //check out the *latest* most relevant version of StRefMultCorr with cvs co -r SLXXX StRoot/StRefMultCorr
  
  ULong_t       mEventsStarted               ;     //  Number of Events read
  ULong_t       mEventsProcessed             ;     //  Number of Events processed and analyzed
  TString       mHistogramOutputFileName     ;     //  Name of the histogram output file 

  bool          AcceptTrack   (StMuTrack*)   ;     //  Function to make cuts on track quality (ie. nhits or eta)


  bool          AcceptEvent   (StMuEvent*)   ;     //  Function to make cuts on event quanitites (ie. vertex position)
  bool          AcceptTrigger (StMuEvent*)   ;     //  Function to make cuts on the trigger words (Cusomize each year)

  bool          RejectRunNumbers(StMuEvent*);      // Function to reject bad runnumbers
  bool          CentralityCut(Double_t CentralityID);		//0-5% centrality

  TVector3      ConvertToTV3(StThreeVectorF inVector);

 protected:

  */
 public:

  HelloWorldMaker(const char* = "hello")   ;          //  Constructor
  virtual          ~HelloWorldMaker( )   ;          //  Destructor
  
  Int_t Init    ( ) ;                               //  Initiliaze the analysis tools ... done once
  Int_t Make    ( ) ;                               //  The main analysis that is done on each event
  Int_t Finish  ( ) ;                               //  Finish the analysis, close files, and clean up.

  //void SetOutputFileName(TString name) { mHistogramOutputFileName = name ; } // Make name available to member functions

  ClassDef(HelloWorldMaker,1)                     //  Macro for CINT compatability
    
};

#endif




