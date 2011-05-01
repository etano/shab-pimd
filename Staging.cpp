#include "Staging.h"

// Initiate Vector Fields
void Paths::initStaging()
{
  // Masses
  for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {  
    M(iBead) = ((iBead + 1.0)/(1.0*iBead)) * m;
  }

  // Positions
  U.set_size(nPart,nBead);   
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    U(iPart,0).zeros(nD);
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      U(iPart,iBead).zeros(nD);
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        U(iPart,iBead)(iD) = normRand(0.0,1.0/(beta*M(iBead)*wp2)); //Quantum Free-Particle Distribution
      }
    }
  }
  UtoRStage();

  // Forces
  UpdateFStage(F);
}

// The Verlet time-stepping algorithm, 'dt' is the time step.
// See eq 116, Ref 1
void Paths::takeStepStage()
{    
  if(thermostat) ApplyThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      U(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
    }
  }
  
  UtoRStage();
  UpdateFStage(F);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  if(thermostat) ApplyThermostat();
}

// Assign actual positions, going from ui's to xi's
// See eq 12.6.6, Ref 2
void Paths::UtoRStage()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    R(iPart, 0) = U(iPart, 0);
    PutInBox( R(iPart, 0) );
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      R(iPart, iBead) = U(iPart, iBead) + ((1.0*iBead)/(iBead + 1.0))*R(iPart, bL(iBead+1)) + (1.0/(iBead + 1.0))*U(iPart, 0);
      PutInBox( R(iPart,iBead) );
    }  
  }
}

// Update the Force for every bead of every particle
// See eqs 12.6.12 & 12.6.13, Ref 2
void Paths::UpdateFStage( field<rowvec>& FX )
{
  rowvec gradVStageA(nD), gradVStageB(nD);

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // 0th bead (note that the mass is zero, so there is no first term)
    gradVStageA.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)
      gradVStageA += getgradV(iPart, iBead);
    FX(iPart,0) = -1.0 * oneOvernBead * gradVStageA;
    gradVStageB = gradVStageA;
    
    // Other beads
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      gradVStageA = getgradV(iPart, iBead) + ((iBead - 1.0)/(1.0*iBead))*gradVStageB;
      FX(iPart,iBead) = -M(iBead) * wp2 * U(iPart,iBead) - oneOvernBead * gradVStageA;
      gradVStageB = gradVStageA;
    }    
    
  }
}
