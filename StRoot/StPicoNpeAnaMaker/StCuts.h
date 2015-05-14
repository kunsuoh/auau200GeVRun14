#ifndef CUTS_H
#define CUTS_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */


namespace cutsAna
{
    // event
    extern float const vz;
    extern float const vzVpdVz;
    extern unsigned char const trigger;
    
    // hadron selectron
    extern float const nSigPion;
    extern bool const pionRequireHFT;
    extern int const pionNHitsFit;
    extern int const pionNhitsDedx;
    extern float const pionHitRatio;
    extern float const pionEta;
    extern float const pionPt;

    // electron
    extern bool const electronRequireHFT;
    extern int const electronNHitsFit;
    extern int const electronNhitsDedx;
    extern float const electronHitRatio;
    extern float const electronEta;
    extern float const electronPt;
    
    
    // partner
    extern int const partnerNHitsFit;
    extern float const partnerEta;
    extern float const partnerPt;
    // partner pid
    extern float const partnerNSigElectron;
    
    
    // electron + partner pair cuts
    extern float const pairMass;
    extern float const pairDca;
    
    
    // pure electron
    extern int const pureElectronNHitsFit;
    extern int const pureElectronNhitsDedx;
    extern float const pureElectronHitRatio;
    extern float const pureElectronEta;
    extern float const pureElectronPt;
    // pure electron pair cut
    extern float const pureElectronMass;
    extern float const pureElectronDca;
    

    //  pid
    extern float const nSigElectron;
    extern int const emcNEta;
    extern int const emcNPhi;
    extern float const emcEoverPLow;
    extern float const emcEoverPHigh;
    extern float const emcPhiDist;
    extern float const emcZDist;
    extern float const emcAssDist;
    extern float const tofBeta;

}
#endif
