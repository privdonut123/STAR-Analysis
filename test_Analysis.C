void test_Analysis(	Int_t nFiles = 10,
                TString InputFileList = "/star/u/bmagh001/temp/checkpoint-07-06-2023/st_fwd_23074018_raw_1000002.event.root",
                Int_t nEvents=1000, 
                Int_t pedLedPhy=2, 
                Int_t eventDisplay=1, 
                int readMuDst=1,
                Int_t debug=0,
                const char *geom = "y2023")
{
    TString _chain;
    //gSystem->Load( "libStarRoot.so" );

    cout << "LL0" << endl;
	//gSystem->Load("libStarClassLibrary.so");
	//gSystem->Load("libStarRoot.so");
	cout << "LL1" << endl;
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	cout << "LL2" << endl;

    bool rerun_tracking = false;

    string fwdTrackOpt = "fwdTrack";
    if (rerun_tracking == false)
        fwdTrackOpt = "";
        
    _chain = Form("in, %s, useXgeom, AgML, db, StEvent, MakeEvent, MuDST, trgd, fcs, fst, ftt, fttQA, fwdTrack, %s", geom, fwdTrackOpt.c_str());
    // "in, y2023, useXgeom, AgML, db, StEvent, MakeEvent"
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");


    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, InputFileList);

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
//StFcsDbMaker *fcsDbMkr= new StFcsDbMaker();
//StFcsDb* fcsDb = (StFcsDb*) chain->GetDataSet("fcsDb");
    //fcsDb->setDbAccess(0);
    //fcsDb->setDebug(debug);


//StEventMaker* eventMk = new StEventMaker(); 
//StFcsRawHitMaker* hitmk = new StFcsRawHitMaker();  
    //hitmk->setDebug(debug);
    //hitmk->setReadMuDst(readMuDst);
//StFcsWaveformFitMaker *wff= new StFcsWaveformFitMaker();
//    wff->setEnergySelect(13,13,1);
StFwdTrackMaker* track = new StFwdTrackMaker();
    //wff->SetDebug(debug);
    
//StFcsClusterMaker *clu= new StFcsClusterMaker;
    //clu->setDebug(1);

//StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");
//gSystem->Load("StFwdTrackMaker");
//StFwdTrackMaker *fwdTrack = (StFwdTrackMaker *)chain->GetMaker("StFwdTrackMaker");

    gSystem->Load( "libStFcsDbMaker.so" ) ;
    StFcsDbMaker * fcsDb = new StFcsDbMaker();
    chain->AddMaker(fcsDb);
	fcsDb->SetDebug();

    //gSystem->Load("StFwdUtils.so");
    //StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
    //chain->AddMaker(fwdAna);
    gSystem->Load("StKumMaker.so");
    StModMaker * kumar = new StModMaker();
    //StModMaker * mod = new StModMaker();
    // Initialize the chain
    chain->Init();

    for (int i = 0; i < nEvents; i++) {
        chain->Clear();
        if (kStOK != chain->Make())
            break;
    }
}