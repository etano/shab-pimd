#include "NoseHoover.h"

void Paths::initNoseHoover()
{
  // Masses (Nose-Hoover)
  // See eq 12.6.14, Ref 2 
  NHM.set_size(nNH);
  NHM(0) = kT/wp2;  
  for (unsigned int iNH = 1; iNH < nNH; iNH += 1)
    NHM(iNH) = kT/wp2;
  oneOverNHM = 1.0/NHM;
  oneOverM = 1.0/M;

  // Positions (Nose-Hoover)
  NHR.zeros(nPart,nBead,nNH);

  // Momenta (Nose-Hoover)
  NHP.zeros(nPart,nBead,nNH);

  // Forces (Nose-Hoover)
  NHF.zeros(nPart,nBead,nNH);
  
  // Suzuki-Yoshida Scheme Weights (Nose-Hoover) NHd
  if (SYOrder == 4) {
    // See eq 4.11.11, Ref 2
    nSY = 3;
    NHd.set_size(nSY);
    NHd(0) = 1.0/(2.0 - pow(2.0,1.0/3.0));
    NHd(2) = NHd(0);
    NHd(1) = 1.0 - NHd(0) - NHd(2);
  } else if (SYOrder == 6) {
    // See eq 4.11.12, Ref 2
    nSY = 7;
    NHd.set_size(nSY);
    NHd(0) = 0.784513610477560;
    NHd(6) = NHd(0);
    NHd(1) = 0.235573213359357;
    NHd(5) = NHd(1);
    NHd(2) = -1.17767998417887;
    NHd(4) = NHd(2);   
    NHd(3) = 1.0 - NHd(0) - NHd(1) - NHd(2)  - NHd(4)  - NHd(5)  - NHd(6);
  } else {
    std::cout << "\nSuzuki-Yoshida Order must be 4 or 6!\n";
    exit(0);
  }
  for (unsigned int iSY = 0; iSY < nSY; iSY += 1) {
    NHd(iSY) *= dt/(nNHsteps*1.0);
  }

}

void Paths::NHThermostat()
{ 
		
  // See eqs 12.6.14 & 4.11.17, Ref 2
  for (unsigned int iNHsteps = 0; iNHsteps < nNHsteps; iNHsteps += 1) {          
    for (unsigned int iSY = 0; iSY < nSY; iSY += 1) {
    
      NHF.slice(nNH-1) = NHP.slice(nNH-2) % NHP.slice(nNH-2) * oneOverNHM(nNH-2) - kT;
      NHP.slice(nNH-1) += 0.25 * NHd(iSY) * NHF.slice(nNH-1);
      
      for (unsigned int iNH = nNH-2; iNH > 0; iNH -= 1) {
        NHP.slice(iNH) %= exp(-0.125 * NHd(iSY) * NHP.slice(iNH+1) * oneOverNHM(iNH+1));
        NHF.slice(iNH) = NHP.slice(iNH-1) % NHP.slice(iNH-1) * oneOverNHM(iNH-1) - kT;
        NHP.slice(iNH) += 0.25 * NHd(iSY) * NHF.slice(iNH);
        NHP.slice(iNH) %= exp(-0.125 * NHd(iSY) * NHP.slice(iNH+1) * oneOverNHM(iNH+1));
      }
      
      NHP.slice(0) %= exp(-0.125 * NHd(iSY) * NHP.slice(1) * oneOverNHM(1));
      
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
          NHF(iPart,iBead,0) = dot( P(iPart,iBead) , P(iPart,iBead) ) * oneOverM(iBead);
        }
      }
      NHF.slice(0) -= nD * kT;

      NHP.slice(0) += 0.25 * NHd(iSY) * NHF.slice(0);
      NHP.slice(0) %= exp(-0.125 * NHd(iSY) * NHP.slice(1) * oneOverNHM(1));
      
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
          for (unsigned int iD = 0; iD < nD; iD += 1) {
            P(iPart,iBead)(iD) *= exp(-0.5 * NHd(iSY) * NHP(iPart,iBead,0) * oneOverNHM(0));  
          }
        }
      }
          
      NHP.slice(0) %= exp(-0.125 * NHd(iSY) * NHP.slice(1) * oneOverNHM(1));
      
      for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
        for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
          NHF(iPart,iBead,0) = dot( P(iPart,iBead) , P(iPart,iBead) ) * oneOverM(iBead);
        }
      }
      NHF.slice(0) -= nD * kT;

      NHP.slice(0) += 0.25 * NHd(iSY) * NHF.slice(0);
      NHP.slice(0) %= exp(-0.125 * NHd(iSY) * NHP.slice(1) * oneOverNHM(1));
      
      for (unsigned int iNH = 1; iNH < nNH-1; iNH += 1) {
        NHP.slice(iNH) %= exp(-0.125 * NHd(iSY) * NHP.slice(iNH+1) * oneOverNHM(iNH+1));
        NHF.slice(iNH) = NHP.slice(iNH-1) % NHP.slice(iNH-1) * oneOverNHM(iNH-1) - kT;
        NHP.slice(iNH) += 0.25 * NHd(iSY) * NHF.slice(iNH);
        NHP.slice(iNH) %= exp(-0.125 * NHd(iSY) * NHP.slice(iNH+1) * oneOverNHM(iNH+1));
      }
      
      NHF.slice(nNH-1) = NHP.slice(nNH-2) % NHP.slice(nNH-2) * oneOverNHM(nNH-2) - kT;
      NHP.slice(nNH-1) += 0.25 * NHd(iSY) * NHF.slice(nNH-1);
      
    }   
  }   
  
}
