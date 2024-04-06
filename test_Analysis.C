void test_Analysis(	Int_t nFiles = 10,
                TString InputFileList = "/gpfs01/star/pwg_tasks/FwdCalib/DAQ/st_fwd_23074018_raw_1000002.daq",
                Int_t nEvents=10, 
                Int_t pedLedPhy=2, 
                Int_t eventDisplay=1, 
                int readMuDst=1,
                Int_t debug=1,
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
        
    _chain = Form("in, %s, useXgeom, AgML, db, StEvent, MakeEvent, MuDST, trgd, fcs, fcsDat, fst, ftt, fttQA, fwdTrack, %s", geom, fwdTrackOpt.c_str());
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
StFcsDbMaker *fcsDbMkr= new StFcsDbMaker();
StFcsDb* fcsDb = (StFcsDb*) chain->GetDataSet("fcsDb");
    //fcsDb->setDbAccess(0);
    //fcsDb->setDebug(debug);
StMuDstMaker muDstMaker(0, 0, "", InputFileList.Data(), "st:MuDst.root", 10); // set up maker in read mode
  //                      0, 0                        this means read mode
  //                           dir                    read all files in this directory
  //                               file               bla.lis read all file in this list, if (file!="") dir is ignored
  //                                    filter        apply filter to filenames, multiple filters are separated by ':'
  //                                          10      maximum number of file to read

//StEventMaker* eventMk = new StEventMaker(); 
StFcsRawHitMaker* hitmk = (StFcsRawHitMaker*) chain->GetMaker("fcsHit");
    hitmk->setReadMuDst(readMuDst);
cout << hitmk << endl;
StFcsWaveformFitMaker *wff= (StFcsWaveformFitMaker*) chain->GetMaker("StFcsWaveformFitMaker");
    //wff->setDebug(debug);
    wff->setAnaWaveform(true);
    //wff->setEnergySelect(13,13,1);
    
    // cout << wff << endl;
    StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
    fwdTrack->setConfigForData( );
    fwdTrack->setGeoCache( "fGeom.root" );
    fwdTrack->setSeedFindingWithFst();
    fwdTrack->SetGenerateTree( false );
    fwdTrack->SetGenerateHistograms( false );
    fwdTrack->SetDebug();
    //wff->SetDebug(debug);
    
//StFcsClusterMaker *clu= new StFcsClusterMaker;
    //clu->setDebug(1);

//StFcsPointMaker *poi=(StFcsPointMaker *)chain->GetMaker("StFcsPointMaker");
//gSystem->Load("StFwdTrackMaker");
//StFwdTrackMaker *fwdTrack = (StFwdTrackMaker *)chain->GetMaker("StFwdTrackMaker");

    // gSystem->Load( "libStFcsDbMaker.so" ) ;
    // StFcsDbMaker * fcsDb = new StFcsDbMaker();
    // chain->AddMaker(fcsDb);
	// fcsDb->SetDebug();

    //gSystem->Load("StFwdUtils.so");
    //StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
    //chain->AddMaker(fwdAna);
    gSystem->Load("StKumMaker.so");
    StHadronAnalysisMaker * kumar = new StHadronAnalysisMaker();
    kumar->setDebug(1);
    TString out_file=Form("./output/output_hadron.root");
    kumar->set_outputfile(out_file.Data());
    //StModMaker * mod = new StModMaker();
    // Initialize the chain
    chain->Init();

    for (int i = 0; i < nEvents; i++) {
        chain->Clear();
        if (kStOK != chain->Make())
            break;
    }
}