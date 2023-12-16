//#include <iostream>
//#include "/direct/star+u/bmagh001/MuDst/include/termcolor.hpp"
using namespace std;
const string red("\033[0;31m");
const string green("\033[1;32m");
const string yellow("\033[1;33m");
const string cyan("\033[0;36m");
const string magenta("\033[0;35m");
const string reset("\033[0m");

    
void hadron_Analysis(	Int_t nFiles = 1,
                		TString InputFileList = "/direct/star+u/bmagh001/MuDst/input/st_fwd_23074018_raw_1000004.MuDst.root",
                		Int_t nEvents=1000, 
                		Int_t pedLedPhy=2, 
						Int_t eventDisplay=1, 
						int readMuDst=1,
						Int_t debug=0	)
	{
  
//gROOT->ProcessLine("#include \"include/termcolor.hpp\"" );
// Load libraries

gSystem->Load( "libStarRoot.so" );

if (gClassTable->GetID("TTable") < 0) {
		gSystem->Load("libStar");
		gSystem->Load("libPhysics");
	}  
gSystem->Load("libStarRoot.so");

gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
loadSharedLibraries();

gSystem->Load("StMagF");
gSystem->Load("StDetectorDbMaker");
gSystem->Load("StTpcDb");
gSystem->Load("StDaqLib");
gSystem->Load("StDbBroker");
gSystem->Load("StDbUtilities");
gSystem->Load("St_db_Maker");

gSystem->Load("StEvent");
gSystem->Load("StEventMaker");

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
gSystem->Load("StarClassLibrary");
gSystem->Load("StStrangeMuDstMaker");
gSystem->Load("StMuDSTMaker");
gSystem->Load("StBTofCalibMaker");
gSystem->Load("StVpdCalibMaker");
gSystem->Load("StBTofMatchMaker");
gSystem->Load("StBTofUtil");	

gSystem->Load("StFcsDbMaker");
gSystem->Load("StFcsRawHitMaker");
gSystem->Load("StFcsWaveformFitMaker");
gSystem->Load("StFcsClusterMaker");
gSystem->Load("libMinuit");
gSystem->Load("StFcsPointMaker");
gSystem->Load("StEpdUtil");
//gSystem->Load("StFcsEventDisplay");

// Getting libraries for forward tracking
gSystem->Load("libStarAgmlUtil.so");
gSystem->Load("libStarAgmlLib.so");
gSystem->Load("libGeometry.so");
gSystem->Load("libStarGeometry.so");
gSystem->Load("libVMC.so");
gSystem->Load("libxgeometry.so");
gSystem->Load("libStBTofUtil.so");
gSystem->Load("libStFstUtil.so");
gSystem->Load("libTree.so");
gSystem->Load("libStDbBroker.so");
gSystem->Load("libSt_db_Maker.so");
gSystem->Load("libStMagF.so");
gSystem->Load("libStDetectorDbMaker.so");
gSystem->Load("libStDbUtilities.so");
gSystem->Load("libStFcsDbMaker.so");
gSystem->Load("libStFttDbMaker.so");
gSystem->Load("libStFstDbMaker.so");
gSystem->Load("libStEventMaker.so");
gSystem->Load("libStFstRawHitMaker.so");
gSystem->Load("libStFstClusterMaker.so");
gSystem->Load("libStFstHitMaker.so");
gSystem->Load("libStFcsRawHitMaker.so");
gSystem->Load("libStFcsWaveformFitMaker.so");
gSystem->Load("libStFcsClusterMaker.so");
gSystem->Load("libStFcsPointMaker.so");
gSystem->Load("libMinuit.so");
gSystem->Load("libStFttRawHitMaker.so");
gSystem->Load("libStFttHitCalibMaker.so");
gSystem->Load("libStFttClusterMaker.so");
gSystem->Load("libStFttPointMaker.so");
gSystem->Load("libStFttQAMaker.so");

gSystem->Load("libXMLIO.so");
gSystem->Load("libgenfit2.so");
gSystem->Load("libKiTrack.so");
gSystem->Load("libStarGeneratorUtil.so");
gSystem->Load("libMathMore.so");
gSystem->Load("libStEpdUtil.so");
gSystem->Load("libStFwdTrackMaker.so");
// make five lines that have gSystem->Load("")

gSystem->Load("StFwdTrackMaker")  ;

 
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
    wff->setEnergySelect(0);
//StFwdTrackMaker* track = new StFwdTrackMaker();
    //wff->SetDebug(debug);
    
StFcsClusterMaker *clu= new StFcsClusterMaker;
    //clu->setDebug(1);

StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");

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

cout << cyan << "InputFileList.Data(): " << InputFileList << reset << endl;

StHcalAnalysisMaker* Hello = new StHcalAnalysisMaker();
TString out_file=Form("./output/output_hadron.root");
Hello->set_outputfile(out_file.Data());
chain->AddMaker(Hello);
cout << green << "StHcalAnalysisMaker added to the chain" << reset << endl;
cout << green << "Chain initiating, please wait" << reset << endl;
chain -> Init();
//chain-> InitRun();

cout <<  chain -> GetNTotal() << endl ;
chain -> EventLoop(nEvents);
cout << green << "Chain finishing, please wait" << reset << endl;
chain -> Finish();
// Cleanup
cout << green << "Deleting chain"  << reset << endl;
delete chain;

}
