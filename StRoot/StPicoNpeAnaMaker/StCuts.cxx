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

#include "StCuts.h"

namespace cutsAna
{
    // event
    float const vz = 6.0;// cm.
    float const vzVpdVz = 3.0; // 3 cm.
    unsigned char const trigger = 19; // 19: BHT1, 21: BHT2, 23: BHT3
    unsigned char const triggerLength = 0x3;
    
    // hadron selectron
    float const nSigPion = 2;
    bool const pionRequireHFT = true;
    int const pionNHitsFit = 20;
    int const pionNhitsDedx = 15;
    float const pionHitRatio = 0.52;
    float const pionEta = 0.7;
    float const pionPt = 1.5;

    // electron
    bool const electronRequireHFT = true;
    int const electronNHitsFit = 20;
    int const electronNhitsDedx = 15;
    float const electronHitRatio = 0.52;
    float const electronEta = 0.7;
    float const electronPt = 1.5;

    
    // partner
    int const partnerNHitsFit = 15;
    float const partnerEta = 1.;
    float const partnerPt = 0.2;
    
    
    // electron + partner pair cuts
    float const pairMass = 0.4;
    float const pairDca = 3.;
    

    
    // pure electron pair cut
    float const pureElectronMass = 0.01;
    float const pureElectronDca = 0.01;
    
    
    // pid
    float const nSigElectron = 3;
    float const partnerNSigElectron = 13;
    int const emcNEta = 0;//1;
    int const emcNPhi = 0;//1;
    float const emcEoverPLow = 0;//0.5;
    float const emcEoverPHigh = 10;//1.7;
    float const emcPhiDist = 10;//0.013;
    float const emcZDist = 100;//2;
    float const emcAssDist = 10;//0.05;
    float const tofBeta = 1;//0.1;
    
    
}
