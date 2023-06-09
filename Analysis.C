#include <iostream>

using namespace std;
const string red("\033[0;31m");
const string green("\033[1;32m");
const string yellow("\033[1;33m");
const string cyan("\033[0;36m");
const string magenta("\033[0;35m");
const string reset("\033[0m");
    

void Analysis(	Int_t nFiles = 1,
                TString InputFileList = "/star/u/bmagh001/MuDst/input/input.list",
                Int_t nEvents=1, 
                Int_t pedLedPhy=2, 
                Int_t eventDisplay=1, 
                int readMuDst=1,
                Int_t debug=0)
	{
  

// Load libraries

//gROOT->Macro("Load.C");
gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");

gSystem->Load( "libStarRoot.so" );

if (gClassTable->GetID("TTable") < 0) {
		gSystem->Load("libStar");
		gSystem->Load("libPhysics");
	}  
	//gSystem->Load("libStarClassLibrary.so");
	gSystem->Load("libStarRoot.so");
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	
	gSystem->Load("StarMagField");
	gSystem->Load("StMagF");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StTpcDb");
	gSystem->Load("StDaqLib");
	gSystem->Load("StDbBroker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("St_db_Maker");

	gSystem->Load("StEvent");
	gSystem->Load("StEventMaker");
	//gSystem->Load("StarMagField");
 
	gSystem->Load("libGeom");
	gSystem->Load("St_g2t");
	
	// Added for Run16 And beyond
	gSystem->Load("libGeom.so");
	
	gSystem->Load("St_base.so");
	gSystem->Load("StUtilities.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("StarAgmlUtil.so");
	gSystem->Load("StarAgmlLib.so");
	gSystem->Load("libStarGeometry.so");
	gSystem->Load("libGeometry.so");
	
	gSystem->Load("xgeometry");
 
	gSystem->Load("St_geant_Maker");


	// needed since I use the StMuTrack
	//gSystem->Load("StarClassLibrary");
	gSystem->Load("StStrangeMuDstMaker");
	gSystem->Load("StMuDSTMaker");
	gSystem->Load("StBTofCalibMaker");
	gSystem->Load("StVpdCalibMaker");
	gSystem->Load("StBTofMatchMaker");
	gSystem->Load("StBTofUtil");	


gSystem->Load("StEventMaker");
gSystem->Load("St_db_Maker");
gSystem->Load("StFcsDbMaker");
gSystem->Load("StFcsRawHitMaker");
gSystem->Load("StFcsWaveformFitMaker");
gSystem->Load("StFcsClusterMaker");
gSystem->Load("libMinuit");
gSystem->Load("StFcsPointMaker");
gSystem->Load("StEpdUtil");
//gSystem->Load("StFcsEventDisplay");
//gSystem->Load("StFwdTrackMaker")  ;

 
gSystem -> Load("StKumMaker.so") ;


//gMessMgr->SetLimit("I",0); //Turn off log info messages
//gMessMgr->SetLimit("Q",0); //turn off log warn messages
//gMessMgr->SetLimit("W",0);

TString OutputFileName;

cout << InputFileList.Data() << endl;
//OutputFileName = (InputFileList.Data()).erase((InputFileList.Data()).begin()+80, (InputFileList.Data()).end()-5);
cout << OutputFileName << endl;
//char edout[200];
StChain* chain = new StChain();
 
StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",InputFileList.Data(),"MuDst.root",nFiles);
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
//StFwdTrackMaker* track = new StFwdTrackMaker();
    //wff->SetDebug(debug);
    
StFcsClusterMaker *clu= new StFcsClusterMaker;
    //clu->setDebug(1);

StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");
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
chain->AddMaker(Hello);
//Int_t nEvents = 5;
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
