#include "NormalMode.h"

// Initiate Vector Fields
void Paths::initNormalMode()
{
  // Create Orthogonal Matrix
  mat NormA(nBead,nBead);
  NormA.zeros();
  for (unsigned int i = 0; i < nBead; i += 1) {
    for (unsigned int j = 0; j < nBead; j += 1) {
      if (i==j) {
        NormA(i,j) = 2;
      } else if (i==j-1) {
        NormA(i,j) = -1;
      } else if (i==j+1) {
        NormA(i,j) = -1;
      }
    }
  }
  NormA(0,nBead-1) = -1;
  NormA(nBead-1,0) = -1;

  vec lambda;
  eig_sym(lambda, NormO, NormA);
  NormO *= sqrt(nBead);
  NormO.print();
  NormA.print();
  lambda.print();

  // Masses
  for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {  
    M(iBead) = m * nBead * lambda(iBead);
  }

  // Positions
  W.set_size(nPart,nBead);   
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      W(iPart,iBead).zeros(nD);
    }
  }
  RtoWNormal();
  
  // Momenta
  field<rowvec> Q = P;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
        P(iPart,iBead) += NormO(jBead,iBead) * Q(iPart,iBead);
      }
      P(iPart,iBead) *= oneOvernBead;
    }
  }

  // Forces
  UpdateFNormal();
}

// The Verlet time-stepping algorithm, 'dt' is the time step.
// See eq 116, Ref 1
void Paths::takeStepNormal()
{    
  if(thermostat) ApplyThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      W(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
    }
  }
  
  WtoRNormal();
  UpdateFNormal();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  if(thermostat) ApplyThermostat();
}

// Assign actual positions, going from wi's to xi's
// See eq 12.6.17, Ref 2
void Paths::WtoRNormal()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      R(iPart,iBead).zeros();
      for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
        R(iPart,iBead) += NormO(iBead,jBead) * W(iPart,jBead);
      }
      PutInBox( R(iPart,iBead) );
    }
  }
}

// Assign actual positions, going from wi's to xi's
// See eq 12.6.17, Ref 2
void Paths::RtoWNormal()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      W(iPart,iBead).zeros();
      for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
        W(iPart,iBead) += NormO(jBead,iBead) * R(iPart,jBead);
      }
      W(iPart,iBead) *= oneOvernBead;
    }
  }
}

// Update the Force for every bead of every particle
// See eqs 12.6.20, Ref 2
void Paths::UpdateFNormal()
{
  if(interaction) UpdateVint();

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // 0th bead (note that the mass is zero, so there is no first term)
    F(iPart,0).zeros();
    for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
      F(iPart,0) -= getgradV(iPart, jBead);
    }
    F(iPart,0) *= oneOvernBead;
    
    // Other beads
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      F(iPart,iBead).zeros();
      for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
        F(iPart,iBead) -= getgradV(iPart, jBead) * NormO(jBead,iBead);
      }
      F(iPart,iBead) *= oneOvernBead;
      F(iPart,iBead) -= M(iBead) * wp2 * W(iPart,iBead);
    }    
    
  }
}
