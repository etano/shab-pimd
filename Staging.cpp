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
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      U(iPart,iBead).zeros(nD);
    }
  }
  RtoUStage();

  // Momenta
  field<rowvec> Q = P;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    P(iPart, 0) = Q(iPart, 0);
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      P(iPart, iBead) = Q(iPart, iBead) - ((1.0*iBead)/(iBead + 1.0))*Q(iPart, bL(iBead+1)) - (1.0/(iBead + 1.0))*P(iPart, 0);
    }  
  }

  // Forces
  UpdateFStage();
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
  UpdateFStage();
  
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

// Assign staging positions, going from xi's to ui's
// See eq 12.6.6, Ref 2
void Paths::RtoUStage()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    U(iPart, 0) = R(iPart, 0);
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      U(iPart, iBead) = R(iPart, iBead) - ((1.0*iBead)/(iBead + 1.0))*R(iPart, bL(iBead+1)) - (1.0/(iBead + 1.0))*U(iPart, 0);
    }  
  }
}

// Update the Force for every bead of every particle
// See eqs 12.6.12 & 12.6.13, Ref 2
void Paths::UpdateFStage()
{
  if(interaction) UpdateVint();

  rowvec gradVStageA(nD), gradVStageB(nD);

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // 0th bead (note that the mass is zero, so there is no first term)
    gradVStageA.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)
      gradVStageA += getgradV(iPart, iBead);
    F(iPart,0) = -1.0 * oneOvernBead * gradVStageA;
    gradVStageB = gradVStageA;
    
    // Other beads
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      gradVStageA = getgradV(iPart, iBead) + ((iBead - 1.0)/(1.0*iBead))*gradVStageB;
      F(iPart,iBead) = -M(iBead) * wp2 * U(iPart,iBead) - oneOvernBead * gradVStageA;
      gradVStageB = gradVStageA;
    }    
    
  }
}
