#include "Langevin.h"

void Paths::initLangevinThermostat()
{
  double LTTau0 = 1.0; // <------ May need to change

  // LT Omegas
  // See eq 20, Ref 3
  LTOmega.set_size(nBead);
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    LTOmega(iBead) = 2.0 * wp * sqrt(nBead) * sin(iBead*pi/(nBead*1.0));
  }

  // LT Gammas
  // See eq 36, Ref 3
  LTGamma.set_size(nBead);
  LTGamma(0) = 1.0/LTTau0;
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    LTGamma(iBead) = 1.0;//2.0 * LTOmega(iBead);
  }

  // LT Constants
  // See eq 28, Ref 3
  // Not sure about whether mass should be m or M(iBead).
  LTC1.set_size(nBead);
  LTC2.set_size(nBead);
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    LTC1(iBead) = exp((-dt/2.0)*LTGamma(iBead));
    LTC2(iBead) = sqrt(M(iBead)/beta) * sqrt(1.0 - LTC1(iBead)*LTC1(iBead));
  }

  // LT C-Matrix
  // See eq 18, Ref 3

}

void Paths::LangevinThermostat()
{ 
		
  // See eqs 28, Ref 3
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) = LTC1(iBead) * P(iPart,iBead) + LTC2(iBead) * normRand(0.0,1.0);
    }
  }
  
}
