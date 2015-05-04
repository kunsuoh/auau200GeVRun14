void runPicoNpeAnaMaker(TString Npelist = "test_npe.list",TString outFileName = "test_Ana.root")
{
    outFileName= "Out/Ana/out_"+Npelist+".root";
    //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
    string SL_version = "SL15c";
    string env_SL = getenv ("STAR");
    if(env_SL.find(SL_version)==string::npos)
    {
        cout<<"Environment Star Library does not match the requested library in runPicoNpeEventMaker.C. Exiting..."<<endl;
        exit(1);
    }
    
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StPicoPrescales");
    gSystem->Load("StPicoNpeEventMaker");
    gSystem->Load("StPicoNpeAnaMaker");
    
    chain = new StChain();
    
    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0,"script/" + Npelist + ".list","picoDstMaker");
    StPicoNpeAnaMaker*  picoNpeAnaMaker = new StPicoNpeAnaMaker("picoNpeAnaMaker","fileLists/Run14/AuAu/200GeV/picoNpeLists/runs/"+Npelist+".list",outFileName.Data(),picoDstMaker);
    
    // -------------- USER variables -------------------------
    // add your cuts here.
    
    chain->Init();
    int nEntries = picoNpeAnaMaker->getEntries();
    cout << " Total entries = " << nEntries << endl;

    for(int iEvent = 0; iEvent<nEntries; ++iEvent)
    {
        if(iEvent%1000==0)
            cout << "Working on eventNumber " << iEvent << endl;

        chain->Clear();
        int iret = chain->Make();
        if (iret) { cout << "Bad return code!" << iret << endl; break;}
    }
    
    chain->Finish();
    delete chain;
    
    // delete list of picos
    //  command = "rm -f correspondingPico.list";
    // gSystem->Exec(command.Data());
    
}
