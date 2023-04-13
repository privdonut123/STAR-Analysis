void HelloWorld(Int_t nFiles, TString InputFileList)
{
// Load libraries
gROOT -> Macro("loadMuDst.C");
gSystem -> Load("HelloWorldMaker.so") ;
gMessMgr->SetLimit("I",0);
gMessMgr->SetLimit("Q",0);
gMessMgr->SetLimit("W",0);

// List of member links in the chain
StChain* chain = new StChain;
StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles);
HelloWorldMaker* Hello = new HelloWorldMaker( );
Int_t nEvents = 1;
// Loop over the links in the chain
chain -> Init();
chain -> EventLoop(1,nEvents);
chain -> Finish();
// Cleanup
delete chain ;
 exit(1);
}
