#include "Staging.h"

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

  /*double pi = math::pi();
  vec lambda(nBead);
  if ((nBead % 2)==0) {
    lambda(0) = 0.0;
    for (unsigned int i = 1; i < nBead/2; i += 1) {
      lambda(2*i - 1) = 2.0*(1.0 - cos(2.0*pi*i/(nBead*1.0)));
      lambda(2*i) = lambda(2*i - 1);
    }
    lambda(nBead - 1) = 2.0*(1.0 - cos(pi*nBead/(nBead*1.0)));
  } else {
    lambda(0) = 0.0;
    for (unsigned int i = 1; i < nBead; i += 2) {
      lambda(i) = 2.0*(1.0 - cos(2.0*pi*i/(nBead*1.0)));
      lambda(i+1) = lambda(i);
    }
  }*/

  // Masses
  for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {  
    M(iBead) = m * nBead * lambda(iBead);
  }

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

// Update the Force for every bead of every particle
// See eqs 12.6.20, Ref 2
void Paths::UpdateFNormal( field<rowvec>& FX )
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // 0th bead (note that the mass is zero, so there is no first term)
    FX(iPart,0).zeros();
    for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
      FX(iPart,0) -= getgradV(iPart, jBead);
    }
    FX(iPart,0) *= oneOvernBead;
    
    // Other beads
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      FX(iPart,iBead).zeros();
      for (unsigned int jBead = 0; jBead < nBead; jBead += 1) {
        FX(iPart,iBead) -= getgradV(iPart, jBead) * NormO(jBead,iBead);
      }
      FX(iPart,iBead) *= oneOvernBead;
      FX(iPart,iBead) -= M(iBead) * wp2 * W(iPart,iBead);
    }    
    
  }
}
