//cout << "ok" <<endl;
#include <iostream>

using namespace std;
const string red("\033[0;31m");
const string green("\033[1;32m");
const string yellow("\033[1;33m");
const string cyan("\033[0;36m");
const string magenta("\033[0;35m");
const string reset("\033[0m");
    

void Analysis(Int_t nFiles = 1,
                TString InputFileList = "test_prod.list",
                Int_t nevents=10, 
                Int_t pedLedPhy=2, 
                Int_t eventDisplay=1, 
                int readMuDst=1,
                Int_t debug=0)
{
  
    //TString opt = "in MakeEvent evout tpcDb trgd fcsDat fcsWFF fcsCluster fcsPoint";

// Load libraries

gROOT->Macro("Load.C");
gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
gSystem->Load("StEventMaker");
gSystem->Load("St_db_Maker");
gSystem->Load("StFcsDbMaker");
gSystem->Load("StFcsRawHitMaker");
gSystem->Load("StFcsWaveformFitMaker");
gSystem->Load("StFcsClusterMaker");
gSystem->Load("libMinuit");
gSystem->Load("StFcsPointMaker");
gSystem->Load("StEpdUtil");
gSystem->Load("StFcsEventDisplay");
gSystem->Load("StFcsTrackMaker")  ;  
gSystem -> Load("StKumMaker.so") ;

gMessMgr->SetLimit("I",0);
gMessMgr->SetLimit("Q",0);
gMessMgr->SetLimit("W",0);
 
TString OutputFileName;

cout << InputFileList.Data() << endl;
//OutputFileName = (InputFileList.Data()).erase((InputFileList.Data()).begin()+80, (InputFileList.Data()).end()-5);
cout << OutputFileName << endl;
//char edout[200];
StChain* chain = new StChain();
 
StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles);
//components for spin db maker
St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb");
if(dbMk){
	dbMk->SetAttr("blacklist", "tpc");
	dbMk->SetAttr("blacklist", "svt");
	dbMk->SetAttr("blacklist", "ssd");
	dbMk->SetAttr("blacklist", "ist");
	dbMk->SetAttr("blacklist", "pxl");
	dbMk->SetAttr("blacklist", "pp2pp");
	dbMk->SetAttr("blacklist", "ftpc");
	dbMk->SetAttr("blacklist", "emc");
	dbMk->SetAttr("blacklist", "eemc");
	dbMk->SetAttr("blacklist", "mtd");
	dbMk->SetAttr("blacklist", "pmd");
	dbMk->SetAttr("blacklist", "tof");
	dbMk->SetAttr("blacklist", "etof");
	dbMk->SetAttr("blacklist", "rhicf");
    }
//StSpinDbMaker* spindb = new StSpinDbMaker("spinDb");
StFcsDbMaker *fcsDbMkr= new StFcsDbMaker();
StFcsDb* fcsDb = (StFcsDb*) chain->GetDataSet("fcsDb");
    //fcsDb->setDbAccess(0);
    //fcsDb->setDebug(debug);

StEventMaker* eventMk = new StEventMaker(); 
StFcsRawHitMaker* hitmk = new StFcsRawHitMaker();  
    //hitmk->setDebug(debug);
    hitmk->setReadMuDst(readMuDst);
StFcsWaveformFitMaker *wff= new StFcsWaveformFitMaker();
    wff->setEnergySelect(13,13,1);
    //wff->SetDebug(debug);
    
StFcsClusterMaker *clu= new StFcsClusterMaker;
    //clu->setDebug(1);

//StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");
    //poi->setDebug(1);
    //poi->setShowerShape(3); 
//gSystem->Load("StVpdCalibMaker");
   //StVpdCalibMaker *vpdCalib = new StVpdCalibMaker();
// List of member links in the chain

//StFcsEventDisplay* fcsed = new StFcsEventDisplay();
/*
StFcsEventDisplay* fcsed;
    if(pedLedPhy>0 && eventDisplay>0){
        //gSystem->Load("StEpdUtil");
	    fcsed = new StFcsEventDisplay();
	    fcsed->setMaxEvents(10);
        string edout = "eventDisplay.png";
	    //sprintf(edout,"ligma.png");
	    fcsed->setFileName(edout.c_str());
	    fcsed->setFilter(1);
    }
    */
/*gSystem->Load("StFcsQaMaker");
    StFcsQaMaker *qaMkr=new StFcsQaMaker("FcsQa");  
    qaMkr->setRun(run);
    */    

StKumMaker* Hello = new StKumMaker();
Int_t nEvents = 5;
// Loop over the links in the chain
cout << "\033[1;31m" << "Chain initiating, please wait" << "\033[0m\n" << endl;
chain -> Init();
//chain-> InitRun();

cout <<  chain -> GetNTotal() << endl ;
chain -> EventLoop(0,nEvents);
cout << "\033[1;31m" << "Chain finishing, please wait" << "\033[0m\n" << endl;
chain -> Finish();
// Cleanup
cout << "Deleting chain" << endl;
delete chain ;
//exit();

}
