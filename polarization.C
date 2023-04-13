#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"


void polarization()
{
  
  fstream list_file;
   list_file.open("input.list",ios::in); //open a file to perform read operation using file object
   if (list_file.is_open()){ //checking whether the file is open
      string root_file_name;
      int nFiles = 0;
	while(getline(list_file, root_file_name) && nFiles < 1)      //read data from file object and put it into string.
	{
	  TFile *root_file = new TFile(root_file_name.c_str());
	  TTree *t1 = (TTree*)root_file->Get("MuDst");
	  Int_t mTime;
	  ((TLeaf*)t1->GetBranch("MuEvent")->GetLeaf("MuEvent.mEventInfo.mTime"))->SetAddress(&mTime);
   
	  TH1F* time = 0;
	  
	   //t1->SetBranchAddress(,&time1);
	  const int NumEntries = t1->GetEntries();
	  cout << NumEntries << endl;


	  for( Int_t iEvent = 0; iEvent<NumEntries; ++iEvent )
	    {
	      t1->GetEntry(iEvent);
      
      
	    }
	  //->At(0))->SetAddress(&DetId);
	  //t1->SetBranchAddress("MuEvent.mEventInfo.mTime",&mTime);

	  //t1->Print();
	  //cout << root_file_name << "\n"; //print the data of the string
	  nFiles++;
	}
      list_file.close(); //close the file object.
   }
   //TFile *f=new TFile("output.root");
  

   
}
