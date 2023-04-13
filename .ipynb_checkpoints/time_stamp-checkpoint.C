void time_stamp()

{

    // if you want to use root.exe instead of root4star, uncomment block below:

//    gSystem->AddDynamicPath("/usr/lib/mysql");
//    gSystem->AddDynamicPath("/usr/lib64/mysql");
//    gSystem->AddDynamicPath("$OPTSTAR/lib/mysql/");
//    gSystem->Load("libmysqlclient");

    // base libraries
    gSystem->Load("St_base");
    gSystem->Load("StChain");
    gSystem->Load("StUtilities");
    gSystem->Load("StIOMaker");
    gSystem->Load("StarClassLibrary");

    // db-related libraries
    gSystem->Load("St_Tables");
    gSystem->Load("StDbLib");
    gSystem->Load("StDbBroker");
    gSystem->Load("St_db_Maker");

    St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    dbMk->SetDebug(0);
    dbMk->SetDateTime(20161220,1); // event or run start time, set to your liking
    dbMk->SetFlavor("ofl");

    dbMk->Init();
    //dbMk->InitRun();
    dbMk->Make();

    TDataSet *DB = 0;
    DB = dbMk->GetDataBase("RunLog/onl/L0TriggerInfo");
    //    DB = dbMk->GetDataBase("RunLog/onl/starMagAvg");
    if (!DB) {
        std::cout << "ERROR: no table found in db, or malformed local db config" << std::endl;
	return;
    }
    /*
    St_fmsGain *dataset = 0;
    dataset = (St_fmsGain*) DB->Find("fmsGain");

    
    St_starMagAvg *dataset = 0;
    dataset = (St_starMagAvg*) DB->Find("starMagAvg");
    */
    St_L0TriggerInfo *dataset = 0;
    dataset = (St_L0TriggerInfo*) DB->Find("L0TriggerInfo");

if (dataset) {

        Int_t rows = dataset->GetNRows();
	if (rows > 1) {
    	    std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        }

	TDatime val[2];
        dbMk->GetValidity((TTable*)dataset,val);
        std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
	<< val[1].GetDate() << "." << val[1].GetTime() << " ] "
        << std::endl;

        L0TriggerInfo_st *table = dataset->GetTable();
        for (Int_t i = 0; i < rows; i++) {
    	    // sample output of first member variable
	  std::cout << i << "th row : " << table[i]->runNumber << std::endl;
        }
    } else {
        std::cout << "ERROR: dataset does not contain requested table" << std::endl;
    }
    
    /*
    if (dataset) {

        Int_t rows = dataset->GetNRows();
	if (rows > 1) {
    	    std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        }

	TDatime val[2];
        dbMk->GetValidity((TTable*)dataset,val);
        std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
	<< val[1].GetDate() << "." << val[1].GetTime() << " ] "
        << std::endl;

        starMagAvg_st *table = dataset->GetTable();
        for (Int_t i = 0; i < rows; i++) {
    	    // sample output of first member variable
            std::cout << i << "th row : " << table[i]->runNumber << std::endl;
        }
    } else {
        std::cout << "ERROR: dataset does not contain requested table" << std::endl;
	}*/

}
