#include <vector>
#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"
#include "../StPicoDstMaker/StPicoDst.h"
#include "../StPicoDstMaker/StPicoDstMaker.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoDstMaker/StPicoTrack.h"
#include "../StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoNpeEvent.h"
#include "StPicoNpeEventMaker.h"
#include "StPicoNpeHists.h"
#include "StCuts.h"

ClassImp(StPicoNpeEventMaker)

//-----------------------------------------------------------------------------
StPicoNpeEventMaker::StPicoNpeEventMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName)
: StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mPicoNpeHists(NULL)
{
    mPicoNpeEvent = new StPicoNpeEvent();
    
    TString baseName(fileBaseName);
    mOutputFile = new TFile(Form("%s.picoNpe.root",fileBaseName), "RECREATE");
    mOutputFile->SetCompressionLevel(1);
    int BufSize = (int)pow(2., 16.);
    int Split = 1;
    mTree = new TTree("T", "T", BufSize);
    mTree->SetAutoSave(1000000); // autosave every 1 Mbytes
    mTree->Branch("dEvent", "StPicoNpeEvent", &mPicoNpeEvent, BufSize, Split);
    
    mPicoNpeHists = new StPicoNpeHists(fileBaseName);
}

//-----------------------------------------------------------------------------
StPicoNpeEventMaker::~StPicoNpeEventMaker()
{
    /* mTree is owned by mOutputFile directory, it will be destructed once
     * the file is closed in ::Finish() */
    delete mPicoNpeHists;
}

//-----------------------------------------------------------------------------
Int_t StPicoNpeEventMaker::Init()
{
    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StPicoNpeEventMaker::Finish()
{
    mOutputFile->cd();
    mOutputFile->Write();
    mOutputFile->Close();
    mPicoNpeHists->closeFile();
    return kStOK;
}
//-----------------------------------------------------------------------------
void StPicoNpeEventMaker::Clear(Option_t *opt)
{
    mPicoNpeEvent->clear("C");
}

//-----------------------------------------------------------------------------
Int_t StPicoNpeEventMaker::Make()
{
    if (!mPicoDstMaker)
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    
    StPicoDst const * picoDst = mPicoDstMaker->picoDst();
    if (!picoDst)
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    
    mPicoEvent = picoDst->event();
    mPicoNpeEvent->addPicoEvent(*mPicoEvent);
    
    if (isGoodEvent())
    {
        UInt_t nTracks = picoDst->numberOfTracks();
        
        std::vector<unsigned short> idxPicoTaggedEs;
        std::vector<unsigned short> idxPicoPartnerEs;
        
        unsigned int nHftTracks = 0;
        
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
        {
            StPicoTrack* trk = picoDst->track(iTrack);
            
            if (!trk || !isGoodTrack(trk)) continue;
            ++nHftTracks;
            
            if (isElectron(trk)) {
                idxPicoTaggedEs.push_back(iTrack);
                StElectronTrack electronTrack((StPicoTrack const *)trk, iTrack);
                mPicoNpeEvent->addElectron(&electronTrack);
            }
            if (isPartnerElectron(trk)) idxPicoPartnerEs.push_back(iTrack);
            
        } // .. end tracks loop

        float const bField = mPicoEvent->bField();
        StThreeVectorF const pVtx = mPicoEvent->primaryVertex();
        
        mPicoNpeEvent->nElectrons(idxPicoTaggedEs.size());
        mPicoNpeEvent->nPartners(idxPicoPartnerEs.size());
        
        for (unsigned short ik = 0; ik < idxPicoTaggedEs.size(); ++ik)
        {
            
            StPicoTrack const * electron = picoDst->track(idxPicoTaggedEs[ik]);
            
            // make electron pairs
            for (unsigned short ip = 0; ip < idxPicoPartnerEs.size(); ++ip)
            {
                
                if (idxPicoTaggedEs[ik] == idxPicoPartnerEs[ip]) continue;
                
                StPicoTrack const * partner = picoDst->track(idxPicoPartnerEs[ip]);
                
                StElectronPair electronPair(electron, partner, idxPicoTaggedEs[ik], idxPicoPartnerEs[ip], pVtx, bField);
                
                
                if (!isGoodElectronPair(electronPair, electron->gPt())) continue;
                
                mPicoNpeEvent->addElectronPair(&electronPair);
                
                if(electron->charge() * partner->charge() <0) // fill histograms for unlike sign pairs only
                {
                    bool fillMass = isGoodQaElectronPair(&electronPair, *electron, *partner);
                    mPicoNpeHists->addElectronPair(&electronPair, electron->gPt(), fillMass);
                }
                
            } // .. end make electron pairs
        } // .. end of tagged e loop
        
        mPicoNpeHists->addEvent(*mPicoEvent,*mPicoNpeEvent,nHftTracks);
        idxPicoTaggedEs.clear();
        idxPicoPartnerEs.clear();
     //   printf("%d\n",mPicoEvent->triggerWord());
    } //.. end of good event fill

    // This should never be inside the good event block
    // because we want to save header information about all events, good or bad
    mTree->Fill();
    mPicoNpeEvent->clear("C");
    
    return kStOK;
}

//-----------------------------------------------------------------------------
bool StPicoNpeEventMaker::isGoodEvent()
{
    return
    fabs(mPicoEvent->primaryVertex().z()) < cuts::vz &&
    fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeEventMaker::isGoodTrack(StPicoTrack const * const trk) const
{
    // Require at least one hit on every layer of PXL and IST.
    // It is done here for tests on the preview II data.
    // The new StPicoTrack which is used in official production has a method to check this
    return
    trk->nHitsFit() >= cuts::nHitsFit;
}

//-----------------------------------------------------------------------------
bool StPicoNpeEventMaker::isElectron(StPicoTrack const * const trk) const
{
    return
    (!cuts::requireHFT || trk->isHFTTrack()) &&
    fabs(trk->nSigmaElectron()) < cuts::nSigmaElectron &&
    trk->gPt() > cuts::pt
    ;
}
//-----------------------------------------------------------------------------
bool StPicoNpeEventMaker::isPartnerElectron(StPicoTrack const * const trk) const
{
    return
    fabs(trk->nSigmaElectron()) < cuts::nSigmaPartnerElectron;
}
//-----------------------------------------------------------------------------
bool StPicoNpeEventMaker::isGoodElectronPair(StElectronPair const & epair, float pt) const
{
    return
    ((epair.pairMass() < cuts::pairMass && pt < 1) || (epair.pairMass() < cuts::pairMassHigh && pt > 1)) &&
    epair.pairDca() < cuts::pairDca &&
    fabs(epair.positionX()) < cuts::positionX &&
    fabs(epair.positionY()) < cuts::positionY &&
    fabs(epair.positionZ()) < cuts::positionZ
    ;
}
//-----------------------------------------------------------------------------
bool  StPicoNpeEventMaker::isGoodQaElectronPair(StElectronPair  const& epair, StPicoTrack const& electron, StPicoTrack const& partner)
{
    return
    electron.nHitsFit() >= cuts::qaNHitsFit &&
    partner.nHitsFit() >= cuts::qaNHitsFit &&

    electron.nSigmaElectron() < cuts::qaNSigmaElectronMax && electron.nSigmaElectron() > cuts::qaNSigmaElectronMin &&
    partner.nSigmaElectron() < cuts::qaNSigmaElectronMax && partner.nSigmaElectron() > cuts::qaNSigmaElectronMin &&
    
    epair.pairDca() < cuts::qaPairDca &&
    epair.pairMass() < cuts::qaPairMass
    ;
}






