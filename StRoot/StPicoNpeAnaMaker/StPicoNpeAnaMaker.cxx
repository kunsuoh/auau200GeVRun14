#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "THnSparse.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoNpeEventMaker/StPicoNpeEvent.h"
#include "StPicoNpeEventMaker/StElectronPair.h"
#include "StPicoNpeEventMaker/StElectronTrack.h"
#include "StPicoNpeAnaMaker.h"
#include "StCuts.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "StBTofUtil/tofPathLength.hh"

ClassImp(StPicoNpeAnaMaker)

StPicoNpeAnaMaker::StPicoNpeAnaMaker(char const * name,char const * inputFilesList,
                                     char const * outName,StPicoDstMaker* picoDstMaker):
StMaker(name),mPicoDstMaker(picoDstMaker),mPicoNpeEvent(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
mOutputFile(NULL), mChain(NULL), mEventCounter(0)
{}

Int_t StPicoNpeAnaMaker::Init()
{
    mPicoNpeEvent = new StPicoNpeEvent();
    
    mChain = new TChain("T");
    std::ifstream listOfFiles(mInputFileList.Data());
    if (listOfFiles.is_open())
    {
        std::string file;
        while (getline(listOfFiles, file))
        {
            LOG_INFO << "StPicoNpeAnaMaker - Adding :" << file << endm;
            mChain->Add(file.c_str());
        }
    }
    else
    {
        LOG_ERROR << "StPicoNpeAnaMaker - Could not open list of files. ABORT!" << endm;
        return kStErr;
    }
    
    mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
    mChain->SetBranchAddress("dEvent", &mPicoNpeEvent);
    
    mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
    mOutputFile->cd();
    
    
    // -------------- USER VARIABLES -------------------------
    hEvent = new TH1F("hEvent","hEvent",10,0,10);
    hRefMult = new TH1F("hRefMult","hRefMult",1000,0,1000);
    hZDCx = new TH1F("hZDCx","hZDCx",10000,0,1000000);
    hHFTInnerOuter = new TH2F("hHFTInnerOuter","hHFTInnerOuter",1000,0,30000,1000,0,30000);
    hHFTInner = new TH1F("hHFTInner","hHFTInner",1000,0,30000);
    hHFTOuter = new TH1F("hHFTOuter","hHFTOuter",1000,0,30000);
    hTrigger = new TH1I("hTrigger","hTrigger",40,0,40);
    
    int bins[4] =   {20, 1000,   289, 800};
    double mins[4] = {-1, 0,      -13, -0.2};
    double maxs[4] = {1, 10,     13,  0.6};
    hsTrackPid = new THnSparseF("hsTrackPid", "hsTrackPid", 4, bins, mins, maxs);
    
    tIncPion = new TTree("tIncPion","tree for Pion form PicoDst");
    tIncPion->Branch("dca",&dca,"dca/F");
    tIncPion->Branch("pt",&pt,"pt/F");
    tIncPion->Branch("eta",&eta,"eta/F");
    
    tIncPion->Branch("nsige",&nsige,"nsige/F");
    tIncPion->Branch("beta",&beta,"beta/F");
    tIncPion->Branch("e",&e,"e/F");
    tIncPion->Branch("e0",&e0,"e0/F");
    tIncPion->Branch("e1",&e1,"e1/F");
    tIncPion->Branch("e2",&e2,"e2/F");
    tIncPion->Branch("e3",&e3,"e3/F");
    tIncPion->Branch("neta",&neta,"neta/b");
    tIncPion->Branch("nphi",&nphi,"nphi/b");
    tIncPion->Branch("phiDist",&phiDist,"phiDist/F");
    tIncPion->Branch("zDist",&zDist,"zDist/F");
    tIncPion->Branch("etaTowDist",&etaTowDist,"etaTowDist/F");
    tIncPion->Branch("phiTowDist",&phiTowDist,"phiTowDist/F");
    
    tIncPion->Branch("mZDCx",&mZDCx,"mZDCx/s");
    tIncPion->Branch("mRefMult",&mRefMult,"mRefMult/s");
    
    
    tInc = new TTree("tInc","tree for Inclusive Elctron");
    tInc->Branch("dca",&dca,"dca/F");
    tInc->Branch("pt",&pt,"pt/F");
    tInc->Branch("eta",&eta,"eta/F");
    
    tInc->Branch("nsige",&nsige,"nsige/F");
    tInc->Branch("beta",&beta,"beta/F");
    tInc->Branch("e",&e,"e/F");
    tInc->Branch("e0",&e0,"e0/F");
    tInc->Branch("e1",&e1,"e1/F");
    tInc->Branch("e2",&e2,"e2/F");
    tInc->Branch("e3",&e3,"e3/F");
    tInc->Branch("neta",&neta,"neta/b");
    tInc->Branch("nphi",&nphi,"nphi/b");
    tInc->Branch("phiDist",&phiDist,"phiDist/F");
    tInc->Branch("zDist",&zDist,"zDist/F");
    tInc->Branch("etaTowDist",&etaTowDist,"etaTowDist/F");
    tInc->Branch("phiTowDist",&phiTowDist,"phiTowDist/F");
    
    tInc->Branch("mZDCx",&mZDCx,"mZDCx/s");
    tInc->Branch("mRefMult",&mRefMult,"mRefMult/s");
    
    
    tPhE = new TTree("tPhE","tree for Photonic Elctron");
    tPhE->Branch("dca",&dca,"dca/F");
    tPhE->Branch("pt",&pt,"pt/F");
    tPhE->Branch("parnter_pt",&partner_pt,"partner_pt/F");
    tPhE->Branch("eta",&eta,"eta/F");
    
    tPhE->Branch("nsige",&nsige,"nsige/F");
    tPhE->Branch("partner_nsige",&partner_nsige,"partner_nsige/F");
    tPhE->Branch("beta",&beta,"beta/F");
    tPhE->Branch("e",&e,"e/F");
    tPhE->Branch("e0",&e0,"e0/F");
    tPhE->Branch("e1",&e1,"e1/F");
    tPhE->Branch("e2",&e2,"e2/F");
    tPhE->Branch("e3",&e3,"e3/F");
    tPhE->Branch("neta",&neta,"neta/b");
    tPhE->Branch("nphi",&nphi,"nphi/b");
    tPhE->Branch("phiDist",&phiDist,"phiDist/F");
    tPhE->Branch("zDist",&zDist,"zDist/F");
    tPhE->Branch("etaTowDist",&etaTowDist,"etaTowDist/F");
    tPhE->Branch("phiTowDist",&phiTowDist,"phiTowDist/F");
    
    tPhE->Branch("pairAngle3d",&pairAngle3d,"pairAngle3d/F");
    tPhE->Branch("pairAnglePhi",&pairAnglePhi,"pairAnglePhi/F");
    tPhE->Branch("pairAngleTheta",&pairAngleTheta,"pairAngleTheta/F");
    tPhE->Branch("pairMass",&pairMass,"pairMass/F");
    tPhE->Branch("pairCharge",&pairCharge,"pairCharge/B");
    tPhE->Branch("pairDca",&pairDca,"pairDca/F");
    tPhE->Branch("pairPositionX",&pairPositionX,"pairPositionX/F");
    tPhE->Branch("pairPositionY",&pairPositionY,"pairPositionY/F");
    
    tPhE->Branch("mZDCx",&mZDCx,"mZDCx/s");
    tPhE->Branch("mRefMult",&mRefMult,"mRefMult/s");
    
    
    tPureE = new TTree("tPureE","tree for Pure Elctron");
    tPureE->Branch("dca",&dca,"dca/F");
    tPureE->Branch("pt",&pt,"pt/F");
    tPureE->Branch("parnter_pt",&partner_pt,"partner_pt/F");
    tPureE->Branch("eta",&eta,"eta/F");
    
    tPureE->Branch("nsige",&nsige,"nsige/F");
    tPureE->Branch("partner_nsige",&partner_nsige,"partner_nsige/F");
    tPureE->Branch("beta",&beta,"beta/F");
    tPureE->Branch("e",&e,"e/F");
    tPureE->Branch("e0",&e0,"e0/F");
    tPureE->Branch("e1",&e1,"e1/F");
    tPureE->Branch("e2",&e2,"e2/F");
    tPureE->Branch("e3",&e3,"e3/F");
    tPureE->Branch("neta",&neta,"neta/b");
    tPureE->Branch("nphi",&nphi,"nphi/b");
    tPureE->Branch("phiDist",&phiDist,"phiDist/F");
    tPureE->Branch("zDist",&zDist,"zDist/F");
    tPureE->Branch("etaTowDist",&etaTowDist,"etaTowDist/F");
    tPureE->Branch("phiTowDist",&phiTowDist,"phiTowDist/F");
    
    tPureE->Branch("pairAngle3d",&pairAngle3d,"pairAngle3d/F");
    tPureE->Branch("pairAnglePhi",&pairAnglePhi,"pairAnglePhi/F");
    tPureE->Branch("pairAngleTheta",&pairAngleTheta,"pairAngleTheta/F");
    tPureE->Branch("pairMass",&pairMass,"pairMass/F");
    tPureE->Branch("pairCharge",&pairCharge,"pairCharge/B");
    tPureE->Branch("pairDca",&pairDca,"pairDca/F");
    tPureE->Branch("pairPositionX",&pairPositionX,"pairPositionX/F");
    tPureE->Branch("pairPositionY",&pairPositionY,"pairPositionY/F");
    
    
    return kStOK;
}
//-----------------------------------------------------------------------------
StPicoNpeAnaMaker::~StPicoNpeAnaMaker()
{
    /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Finish()
{
    LOG_INFO << " StPicoNpeAnaMaker - writing data and closing output file " <<endm;
    mOutputFile->cd();
    
    hEvent->Write();
    hRefMult->Write();
    hZDCx->Write();
    hHFTInnerOuter->Write();
    hHFTInner->Write();
    hHFTOuter->Write();
    hTrigger->Write();
    
    tInc->Write();
    tIncPion->Write();
    tPhE->Write();
    tPureE->Write();
    
    mOutputFile->Close();
    
    return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoNpeAnaMaker::Make()
{
    hEvent->Fill(0);
    readNextEvent();
    
    if (!mPicoDstMaker)
    {
        LOG_WARN << " StPicoNpeAnaMaker - No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(1);
    
    StPicoDst const* picoDst = mPicoDstMaker->picoDst();
    
    if (!picoDst)
    {
        LOG_WARN << "StPicoNpeAnaMaker - No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    hEvent->Fill(2);
    
    if(mPicoNpeEvent->runId() != picoDst->event()->runId() ||
       mPicoNpeEvent->eventId() != picoDst->event()->eventId())
    {
        LOG_ERROR <<" runId() " << mPicoNpeEvent->runId() << " " << picoDst->event()->runId()  <<endm;
        LOG_ERROR <<" eventId() " << mPicoNpeEvent->eventId() << " " << picoDst->event()->eventId()  <<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
        LOG_ERROR <<" StPicoNpeAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoNpeEvent are not in sync."<<endm;
        return kStWarn;
    }
    hEvent->Fill(3);
    
    mRefMult = std::numeric_limits<Int_t>::quiet_NaN();
    mZDCx = std::numeric_limits<Int_t>::quiet_NaN();
    
    mRefMult = picoDst->event()->refMult();
    mZDCx = picoDst->event()->ZDCx();
    
    hRefMult->Fill(mRefMult);
    hZDCx->Fill(mZDCx);
    hHFTInnerOuter->Fill(picoDst->event()->numberOfPxlInnerHits(),picoDst->event()->numberOfPxlOuterHits());
    hHFTInner->Fill(picoDst->event()->numberOfPxlInnerHits());
    hHFTOuter->Fill(picoDst->event()->numberOfPxlOuterHits());
    
    for (int i=0;i<40;i++) if (picoDst->event()->triggerWord()>>i & 1)  hTrigger->Fill(i);
    
    
    // -------------- USER ANALYSIS -------------------------
    StThreeVectorF const pVtx = picoDst->event()->primaryVertex();
    
    // Inclusive hadrons with StPicoTrack
    UInt_t nTracks = picoDst->numberOfTracks();
    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
        StPicoTrack* track = picoDst->track(iTrack);
        if (!track) continue;
        if (isGoodPion(track)) {
            StPhysicalHelixD eHelix = track->dcaGeometry().helix();
            
            dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
            pt = track->gPt();
            eta = track->gMom().pseudoRapidity();
            
            nsige = track->nSigmaElectron();
            
            beta = std::numeric_limits<float>::quiet_NaN();
            e = std::numeric_limits<float>::quiet_NaN();
            e0 = std::numeric_limits<float>::quiet_NaN();
            e1 = std::numeric_limits<float>::quiet_NaN();
            e2 = std::numeric_limits<float>::quiet_NaN();
            e3 = std::numeric_limits<float>::quiet_NaN();
            nphi = std::numeric_limits<unsigned short>::quiet_NaN();
            neta = std::numeric_limits<unsigned short>::quiet_NaN();
            phiDist = std::numeric_limits<float>::quiet_NaN();
            zDist = std::numeric_limits<float>::quiet_NaN();
            etaTowDist = std::numeric_limits<float>::quiet_NaN();
            phiTowDist = std::numeric_limits<float>::quiet_NaN();
            if (isGoodTofTrack(track))  {
                StPicoBTofPidTraits * Tof = picoDst->btofPidTraits(track->bTofPidTraitsIndex());
                
                // start global beta calculation
                double newBeta = Tof->btofBeta();
                if(newBeta<1e-4 && fabs(dca)<3.) { 
                    StThreeVectorF btofHitPos = Tof->btofHitPos();
                    float L = tofPathLength(&pVtx, &btofHitPos, eHelix.curvature());
                    float tof = Tof->btof();
                	if (tof>0) newBeta = L/(tof*(C_C_LIGHT/1.e9));
                }
                beta = 1./newBeta;
                // end global beta calculation
            }
            if (isGoodEmcTrack(track)) {
                StPicoEmcPidTraits * Emc =  picoDst->emcPidTraits(track->emcPidTraitsIndex());
                e = Emc->e();
                e0 = Emc->e0();
                e1 = Emc->e1();
                e2 = Emc->e2();
                e3 = Emc->e3();
                nphi = Emc->nPhi();
                neta = Emc->nEta();
                phiDist = Emc->phiDist();
                zDist = Emc->zDist();
                etaTowDist = Emc->etaTowDist();
                phiTowDist = Emc->phiTowDist();
            }
            
            
            tIncPion->Fill();
        }
        
    }
    
    // Inclusive Electron with StPicoNpeEvents
    TClonesArray const * aElectron = mPicoNpeEvent->electronArray();
    for (int idx = 0; idx < aElectron->GetEntries(); ++idx)
    {
        StElectronTrack const* electron  = (StElectronTrack*)aElectron->At(idx);
        StPicoTrack const* etrack = picoDst->track(electron->electronIdx());
        if (!etrack || !isGoodElectron(etrack)) continue;
        if (!isGoodTofTrack(etrack) && !isGoodEmcTrack(etrack)) continue;
        
        StPhysicalHelixD eHelix = etrack->dcaGeometry().helix();
        
        dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
        //     dca = eHelix.geometricSignedDistance(pVtx);
        pt = etrack->gPt();
        eta = etrack->gMom().pseudoRapidity();
        nsige = etrack->nSigmaElectron();
        
        beta = std::numeric_limits<float>::quiet_NaN();
        e = std::numeric_limits<float>::quiet_NaN();
        e0 = std::numeric_limits<float>::quiet_NaN();
        e1 = std::numeric_limits<float>::quiet_NaN();
        e2 = std::numeric_limits<float>::quiet_NaN();
        e3 = std::numeric_limits<float>::quiet_NaN();
        nphi = std::numeric_limits<unsigned short>::quiet_NaN();
        neta = std::numeric_limits<unsigned short>::quiet_NaN();
        phiDist = std::numeric_limits<float>::quiet_NaN();
        zDist = std::numeric_limits<float>::quiet_NaN();
        etaTowDist = std::numeric_limits<float>::quiet_NaN();
        phiTowDist = std::numeric_limits<float>::quiet_NaN();
        
        if (isGoodTofTrack(etrack))  {
            StPicoBTofPidTraits * Tof = picoDst->btofPidTraits(etrack->bTofPidTraitsIndex());
                // start global beta calculation
                double newBeta = Tof->btofBeta();
                if(newBeta<1e-4 && fabs(dca)<3.) { 
                    StThreeVectorF btofHitPos = Tof->btofHitPos();
                    float L = tofPathLength(&pVtx, &btofHitPos, eHelix.curvature());
                    float tof = Tof->btof();
                	if (tof>0) newBeta = L/(tof*(C_C_LIGHT/1.e9));
                }
                beta = 1./newBeta;
                // end global beta calculation
        }
        if (isGoodEmcTrack(etrack)) {
            StPicoEmcPidTraits * Emc =  picoDst->emcPidTraits(etrack->emcPidTraitsIndex());
            e = Emc->e();
            e0 = Emc->e0();
            e1 = Emc->e1();
            e2 = Emc->e2();
            e3 = Emc->e3();
            nphi = Emc->nPhi();
            neta = Emc->nEta();
            phiDist = Emc->phiDist();
            zDist = Emc->zDist();
            etaTowDist = Emc->etaTowDist();
            phiTowDist = Emc->phiTowDist();
        }
        tInc->Fill();
    }
    
    float const bField = picoDst->event()->bField();
    
    // Photonic Electron
    TClonesArray const * aElectronPair = mPicoNpeEvent->electronPairArray();
    for (int idx = 0; idx < aElectronPair->GetEntries(); ++idx)
    {
        // this is an example of how to get the ElectronPair pairs and their corresponsing tracks
        StElectronPair const* epair = (StElectronPair*)aElectronPair->At(idx);
        if(!isGoodPair(epair)) continue;
        
        StPicoTrack const* electron = picoDst->track(epair->electronIdx());
        StPicoTrack const* partner = picoDst->track(epair->partnerIdx());
        
        // calculate Lorentz vector of electron-partner pair
        StPhysicalHelixD electronHelix = electron->dcaGeometry().helix();
        StPhysicalHelixD partnerHelix = partner->dcaGeometry().helix();
        pair<double,double> ss = electronHelix.pathLengths(partnerHelix);
        StThreeVectorF const electronMomAtDca = electronHelix.momentumAt(ss.first, bField * kilogauss);
        StThreeVectorF const partnerMomAtDca = partnerHelix.momentumAt(ss.second, bField * kilogauss);
        StLorentzVectorF const electronFourMom(electronMomAtDca, electronMomAtDca.massHypothesis(M_ELECTRON));
        StLorentzVectorF const partnerFourMom(partnerMomAtDca, partnerMomAtDca.massHypothesis(M_ELECTRON));
        StLorentzVectorF const epairFourMom = electronFourMom + partnerFourMom;
        StPhysicalHelixD eHelix = electron->dcaGeometry().helix();
        
        pairMass = epairFourMom.m();
        pairAngle3d = electronMomAtDca.angle(partnerMomAtDca);
        pairAnglePhi = fabs(electronMomAtDca.phi() - partnerMomAtDca.phi());
        pairAngleTheta = fabs(electronMomAtDca.theta() - partnerMomAtDca.theta());
        pairCharge = electron->charge()+partner->charge();
        pairDca = epair->pairDca();
        pairPositionX = epair->positionX();
        pairPositionY = epair->positionY();
        
        dca = eHelix.curvatureSignedDistance(pVtx.x(),pVtx.y());
        pt = electron->gPt();
        partner_pt = partner->gPt();
        eta = electron->gMom().pseudoRapidity();
        nsige = electron->nSigmaElectron();
        partner_nsige = partner->nSigmaElectron();
        
        beta = std::numeric_limits<float>::quiet_NaN();
        e = std::numeric_limits<float>::quiet_NaN();
        e0 = std::numeric_limits<float>::quiet_NaN();
        e1 = std::numeric_limits<float>::quiet_NaN();
        e2 = std::numeric_limits<float>::quiet_NaN();
        e3 = std::numeric_limits<float>::quiet_NaN();
        nphi = std::numeric_limits<unsigned short>::quiet_NaN();
        neta = std::numeric_limits<unsigned short>::quiet_NaN();
        phiDist = std::numeric_limits<float>::quiet_NaN();
        zDist = std::numeric_limits<float>::quiet_NaN();
        etaTowDist = std::numeric_limits<float>::quiet_NaN();
        phiTowDist = std::numeric_limits<float>::quiet_NaN();
        
        if (isGoodTofTrack(electron))  {
            StPicoBTofPidTraits * Tof = picoDst->btofPidTraits(electron->bTofPidTraitsIndex());
                // start global beta calculation
                double newBeta = Tof->btofBeta();
                if(newBeta<1e-4 && fabs(dca)<3.) { 
                    StThreeVectorF btofHitPos = Tof->btofHitPos();
                    float L = tofPathLength(&pVtx, &btofHitPos, eHelix.curvature());
                    float tof = Tof->btof();
                	if (tof>0) newBeta = L/(tof*(C_C_LIGHT/1.e9));
                }
                beta = 1./newBeta;
                // end global beta calculation
        }
        if (isGoodEmcTrack(electron)) {
            StPicoEmcPidTraits * Emc =  picoDst->emcPidTraits(electron->emcPidTraitsIndex());
            e = Emc->e();
            e0 = Emc->e0();
            e1 = Emc->e1();
            e2 = Emc->e2();
            e3 = Emc->e3();
            nphi = Emc->nPhi();
            neta = Emc->nEta();
            phiDist = Emc->phiDist();
            zDist = Emc->zDist();
            etaTowDist = Emc->etaTowDist();
            phiTowDist = Emc->phiTowDist();
        }
        if(isGoodPureElectron(epair)) tPureE->Fill();
        if(!isGoodTofTrack(electron) && !isGoodEmcTrack(electron)) continue;
        tPhE->Fill();
    }
    return kStOK;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPureElectron(StElectronPair const* const epair) const
{
    if(!epair) return false;
    
    StPicoTrack const* electron = mPicoDstMaker->picoDst()->track(epair->electronIdx());
    StPicoTrack const* partner = mPicoDstMaker->picoDst()->track(epair->partnerIdx());
    
    return
    isGoodElectron(electron) &&
    isGoodElectron(partner) &&
    epair->pairMass() < cutsAna::pureElectronMass &&
    epair->pairDca() < cutsAna::pureElectronDca
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPair(StElectronPair const* const epair) const
{
    if(!epair) return false;
    
    StPicoTrack const* electron = mPicoDstMaker->picoDst()->track(epair->electronIdx());
    StPicoTrack const* partner = mPicoDstMaker->picoDst()->track(epair->partnerIdx());
    
    return
    isGoodElectron(electron) &&
    isGoodPartner(partner) &&
    epair->pairMass() < cutsAna::pairMass &&
    epair->pairDca() < cutsAna::pairDca
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodElectron(StPicoTrack const * const trk) const
{
    return
    (!cutsAna::electronRequireHFT || trk->isHFTTrack()) &&
    trk->nHitsFit() >= cutsAna::electronNHitsFit &&
    trk->nHitsDedx() >= cutsAna::electronNhitsDedx &&
    //  trk->hitRatio() >= cutsAna::electronHitRatio &&
    fabs(trk->gMom().pseudoRapidity()) <= cutsAna::electronEta &&
    trk->gPt() >= cutsAna::electronPt &&
    trk->nSigmaElectron() <= cutsAna::nSigElectron
    ;
}//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPion(StPicoTrack const * const trk) const
{
    return
    (!cutsAna::pionRequireHFT || trk->isHFTTrack()) &&
    trk->nHitsFit() >= cutsAna::pionNHitsFit &&
    trk->nHitsDedx() >= cutsAna::pionNhitsDedx &&
    //  trk->hitRatio() >= cutsAna::pionHitRatio &&
    fabs(trk->gMom().pseudoRapidity()) <= cutsAna::pionEta &&
    trk->gPt() >= cutsAna::pionPt
    //   trk->nSigmaPion() <= cutsAna::nSigPion
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodPartner(StPicoTrack const * const trk) const
{
    return
    trk->nHitsFit() >= cutsAna::partnerNHitsFit &&
    fabs(trk->gMom().pseudoRapidity()) <= cutsAna::partnerEta &&
    trk->gPt() >= cutsAna::partnerPt &&
    abs(trk->nSigmaElectron()) <= cutsAna::partnerNSigElectron
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodTofTrack(StPicoTrack const * const trk) const
{
    if (trk->bTofPidTraitsIndex() < 0) return false;
    StPicoBTofPidTraits * Tof = mPicoDstMaker->picoDst()->btofPidTraits(trk->bTofPidTraitsIndex());
    
    return
    // TMath::Abs(1-1/Tof->btofBeta()) < cutsAna::tofBeta
    true
    ;
}

//-----------------------------------------------------------------------------
bool StPicoNpeAnaMaker::isGoodEmcTrack(StPicoTrack const * const trk) const
{
    if (trk->emcPidTraitsIndex() < 0) return false;
    StPicoEmcPidTraits * Emc =  mPicoDstMaker->picoDst()->emcPidTraits(trk->emcPidTraitsIndex());
    
    return
    /*   Emc->nPhi() > cutsAna::emcNPhi &&
     Emc->nEta() > cutsAna::emcNEta &&
     TMath::Abs(Emc->phiDist()) < cutsAna::emcPhiDist &&
     TMath::Abs(Emc->zDist()) < cutsAna::emcZDist &&
     TMath::Sqrt(Emc->phiTowDist()*Emc->phiTowDist() + Emc->etaTowDist()*Emc->etaTowDist()) < cutsAna::emcAssDist &&
     Emc->e0()/trk->dcaGeometry().momentum().mag() < cutsAna::emcEoverPHigh &&
     Emc->e0()/trk->dcaGeometry().momentum().mag() > cutsAna::emcEoverPLow
     */
    true
    ;
}
