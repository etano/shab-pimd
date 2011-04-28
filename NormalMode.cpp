#include "Staging.h"

// Initiate Vector Fields
void Paths::initNormalMode()
{
  // Positions
  W.set_size(nPart,nBead);   
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    W(iPart,0).zeros(nD);
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      W(iPart,iBead).zeros(nD);
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        W(iPart,iBead)(iD) = normRand(0.0,1.0/(beta*M(iBead)*wp2)); //Quantum Free-Particle Distribution
      }
    }
  }
  WtoRNormal();

  // Forces
  UpdateFNormal(F);
}

// The Verlet time-stepping algorithm, 'dt' is the time step.
// See eq 116, Ref 1
void Paths::takeStepNormal()
{    
  if (useNH) NHThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      W(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
    }
  }
  
  WtoRNormal();
  UpdateFNormal(F);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  if (useNH) NHThermostat();  
}

// Assign actual positions, going from wi's to xi's
// See eq 12.6.6, Ref 2
void Paths::WtoRNormal()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    R(iPart, 0) = W(iPart, 0);
    PutInBox( R(iPart, 0) );
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      R(iPart, iBead) = W(iPart, iBead) + ((1.0*iBead)/(iBead + 1.0))*R(iPart, bL(iBead+1)) + (1.0/(iBead + 1.0))*W(iPart, 0);
      PutInBox( R(iPart,iBead) );
    }  
  }
}

// Update the Force for every bead of every particle
// See eqs 12.6.12 & 12.6.13, Ref 2
void Paths::UpdateFNormal( field<rowvec>& FX )
{
  rowvec gradVNormalA(nD), gradVNormalB(nD);

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // 0th bead (note that the mass is zero, so there is no first term)
    gradVNormalA.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)
      gradVNormalA += getgradV(iPart, iBead);
    FX(iPart,0) = -1.0 * oneOvernBead * gradVNormalA;
    gradVNormalB = gradVNormalA;
    
    // Other beads
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      gradVNormalA = getgradV(iPart, iBead) + ((iBead - 1.0)/(1.0*iBead))*gradVNormalB;
      FX(iPart,iBead) = -M(iBead) * wp2 * W(iPart,iBead) - oneOvernBead * gradVNormalA;
      gradVNormalB = gradVNormalA;
    }    
    
  }
}
