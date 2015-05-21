#ifndef StPicoNpeAnaMaker_h
#define StPicoNpeAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoNpeEvent
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */
 

#include "TChain.h"
#include "StMaker.h"
#include "TH2F.h"

class TString;
class TFile;
class StPicoNpeEvent;
class StElectronTrack;
class StElectronPair;
class StPicoDstMaker;
class StPicoTrack;
class StPicoDst;

class StPicoNpeAnaMaker : public StMaker
{
  public:
    StPicoNpeAnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker);
    virtual ~StPicoNpeAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;
    
    
  private:
    StPicoNpeAnaMaker() {}
    void readNextEvent();

    bool isGoodElectron(StPicoTrack const*) const;
    bool isGoodPion(StPicoTrack const*) const;
    bool isGoodPartner(StPicoTrack const*) const;
    bool isGoodTofTrack(StPicoTrack const*) const;
    bool isGoodEmcTrack(StPicoTrack const*) const;
    bool isGoodPair(StElectronPair const*) const;
    bool isGoodPureElectron(StElectronPair const*) const;
    bool isGoodEvent() const;
    void setTree(TTree *, TString);
    void initVariables();
    void setVariables(StPicoTrack *);
    void setVariables(StElectronPair *);

    StPicoDstMaker* mPicoDstMaker;
    StPicoNpeEvent* mPicoNpeEvent;
    StPicoDst * picoDst;
    
    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;
    
    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    TH1F * hEvent;
    TH1F * hRefMult;
    TH1F * hZDCx;
    TH1I * hTrigger;
    
    TH2F * hHFTInnerOuter;
    TH1F * hHFTInner;
    TH1F * hHFTOuter;
    
    TTree * tIncPion;
    TTree * tInc;
    TTree * tPhE;
    TTree * tPureE;

    Float_t pairAngle3d;
    Float_t pairAnglePhi;
    Float_t pairAngleTheta;
    Float_t pairMass;
    char pairCharge;
    Float_t pairDca;
    Float_t pairPositionX;
    Float_t pairPositionY;
    
    
    Float_t dca;
    Float_t pt;
    Float_t partner_pt;
    Float_t eta;
    Float_t nsige;
    Float_t partner_nsige;
    Float_t beta;
    Float_t e;
    Float_t e0;
    Float_t e1;
    Float_t e2;
    Float_t e3;
    unsigned char neta;
    unsigned char nphi;
    Float_t phiDist;
    Float_t zDist;
    Float_t etaTowDist;
    Float_t phiTowDist;
    
    UShort_t mZDCx;
    UShort_t mRefMult;

    ClassDef(StPicoNpeAnaMaker, 1)
};

inline int StPicoNpeAnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoNpeAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

#endif
