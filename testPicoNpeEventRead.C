void testPicoD0EventRead(TString filename="test.picoD0.root")
{
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StPicoD0EventMaker");
    
    TFile* f = new TFile(filename.Data());
    TTree* T = (TTree*)f->Get("T");
    StPicoD0Event* event = new StPicoD0Event();
    T->SetBranchAddress("dEvent",&event);
    
    TFile ff("read_test.root","RECREATE");
    TNtuple* nt = new TNtuple("nt","","m:pt:eta:phi:theta:"
                              "decayL:kDca:pDca:dca12:cosThetaStar");
    
    StElectronPair* kp = 0;
    
    for(Int_t i=0;i<100000;++i)
    {
        T->GetEntry(i);
        
        TClonesArray* arrKPi = event->electronPairArray();
        
        for(int idx=0;idx<event->nElectronPair();++idx)
        {
            kp = (StElectronPair*)arrKPi->At(idx);
            
            nt->Fill(kp->m(),kp->pt(),kp->eta(),kp->phi(),kp->pointingAngle(),
                     kp->decayLength(),kp->e1Dca(),kp->e2Dca(),kp->dcaDaughters(),kp->cosThetaStar());
        }
    }
    
    nt->Write();
    ff.Close();
}
