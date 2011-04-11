
#include "Paths.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const bool stageIn )
{
  // Set constants
  nPart = nPartIn;
  nD = nDIn;
  nBead = nBeadIn;
  beta = betaIn;
  dt = dtIn;
  L = LIn;    
  stage = stageIn;
  
  m = 1.0; // Default (harmonic oscillator units)
  hbar = 1.0; // Default (harmonic oscillator units)
  w = 1.0; // Default (harmonic oscillator units)
  
  wp = sqrt(nBead*1.0)/(beta*hbar);
  wp2 = wp*wp;
  mwp2 = m*wp*wp;
  mw2 = m*w*w;
  
  oneOvernBead = 1.0/(1.0*nBead);
  oneOvernPartnBead = 1.0/(1.0*nPart*nBead);
  nDnPartOver2Beta = 1.0*nD*nPart/(2.0*beta);
  nDnPartnBeadOver2Beta = 1.0*nD*nPart*nBead/(2.0*beta);
  
  kT = 1.0/beta;
  
  // Set bead loop
  bL.set_size(2*nBead);
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    bL(iBead) = iBead;
    bL(iBead + nBead) = iBead;
  }
  
  // Set masses
  M.set_size(nBead);
  if (useNH) Q.set_size(nBead); // N-H masses
  if (stage) {
    M(0) = m; // See eq 12.6.10, Ref 2
    if (useNH) Q(0) = kT*beta*beta/(nBead*nBead); // See eq 12.6.14, Ref 2
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      M(iBead) = ((iBead + 1.0)/(1.0*iBead)) * m; // See eq 12.6.10, Ref 2
      if (useNH) Q(iBead) = kT/wp2; // See eq 12.6.14, Ref 2
    }
  } else {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      M(iBead) = m;
    }
  }
  
  // Size vector fields  
  
  // Positions R
  R.set_size(nPart,nBead);
  if (useNH) NHR.set_size(nPart,nBead);
  
  // Positions U
  if (stage) U.set_size(nPart,nBead);
  
  // Momenta
  P.set_size(nPart,nBead);
  if (useNH) NHP.set_size(nPart,nBead);
  
  // Forces
  F.set_size(nPart,nBead);
  nF.set_size(nPart,nBead);
  
  // Size vectors
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {  

      // Positions R
      R(iPart,iBead).zeros(nD);
      
      // Positions U
      if (stage) U(iPart,iBead).zeros(nD);

      // Momenta
      P(iPart,iBead).zeros(nD); 

      // Forces 
      F(iPart,iBead).zeros(nD);
      nF(iPart,iBead).zeros(nD); 

    }
  }
  
  // Initialize system
  if (stage) {
    InitPosition(U);
    UtoRStage();
  } else {
    InitPosition(R);
  }
  InitMomentum(1.0/beta);

  // Compute First Force
  if (stage) {
    UpdateFStage(F); 
  } else {
    UpdateF(F);
  }

}

// Paths destructor
Paths::~Paths()
{
}

//////////////////////
/* Intialize System */
//////////////////////

// Initialize Positions
void Paths::InitPosition( field<rowvec>& X )
{  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      X(iPart,iBead).zeros();
    }
  }
}

// Initialize Momentum
void Paths::InitMomentum( double T )
{
  rowvec newP(nD), netP(nD), avgP(nD);
  double netE, pScale;
    
  netP.zeros();
  netE = 0.0;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {  
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        newP(iD) = unifRand() - 0.5;
      }
      netP += newP;
      netE += dot(newP,newP);
      P(iPart,iBead) = newP;
    }
  }

  avgP = netP/(1.0*nPart*nBead);
  if (nPart < 2) avgP.zeros();
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      pScale = sqrt(M(iBead)*nD*nPart*nBead*T/netE);
      P(iPart,iBead) = (P(iPart,iBead) - avgP)*pScale;      
    }
  }
}

////////////////////////////////////////////////////
/* Routines to ensure periodic boundary coditions */
////////////////////////////////////////////////////

void Paths::PutInBox(rowvec& Ri)
{  
  /*for(unsigned int i = 0; i < nD; i++) {
    while(Ri(i) > L/2.0) Ri(i) -= L;
    while(Ri(i) < -L/2.0) Ri(i) += L;
  }*/
}
  
double Paths::Distance(rowvec& Ri, rowvec& Rj)
{  
  double d = norm( Displacement(Ri, Rj) , 2 );  
  return d;
}
  
rowvec Paths::Displacement(rowvec& Ri, rowvec& Rj)
{
  rowvec dR = Ri - Rj;  
  PutInBox(dR);  
  return dR;
}

/////////////////
/* OBSERVABLES */
/////////////////

// Get Primitive Energy Estimator
// See eq. 12.6.40, Ref. 2
double Paths::getPE()
{
  rowvec dR(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dR = Displacement( R(iPart,iBead) , R(iPart,bL(iBead+1)) );
      total += 0.5 * mwp2 * dot(dR,dR) - oneOvernBead * getV(iPart, iBead);
    }
  }
  
  return nDnPartnBeadOver2Beta - total;
}

// Get Virial Energy Estimator
// See eq. 12.6.40, Ref. 2
double Paths::getVE()
{
  rowvec rc(nD), dR(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    rc.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) rc += R(iPart,iBead);
    rc = rc * oneOvernBead;
    PutInBox(rc);
    
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dR = Displacement( R(iPart,iBead) , rc );
      total += 0.5 * dot(dR,dR) * getdV(iPart, iBead) + getV(iPart, iBead);
    }
  }
  
  return nDnPartOver2Beta + oneOvernBead*total;
}

// Position Estimator
double Paths::getR()
{
  rowvec Rsum(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {  
    Rsum.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)
      Rsum += R(iPart,iBead);

    total += Rsum(0);
  }
  
  return oneOvernPartnBead * total;
}

// Absolute Position Estimator
double Paths::getR2()
{
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      total += norm( R(iPart,iBead) , 2 );
    }
  }
  
  return oneOvernPartnBead * total;
}

/////////////////////////
/* Potential Functions */
/////////////////////////

// Get Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// V = (1/2) m w^2 r^2
double Paths::getV( const int iPart, const int iBead )
{
  return 0.5 * mw2 * dot( R(iPart,iBead) , R(iPart,iBead) );
}

// Get Derivative of Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// dV/dr = m w^2 r
double Paths::getdV( const int iPart, const int iBead )
{
  return mw2 * norm( R(iPart,iBead) , 2 );
}

// Get Derivative of Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// dV/dr = m w^2 r
rowvec Paths::getgradV( const int iPart, const int iBead )
{
  return mw2 * R(iPart,iBead);
}

//////////////////////////////////
/* Molecular Dynamics Functions */
//////////////////////////////////

// The Verlet time-stepping algorithm, 'dt' is the time step.
void Paths::takeStep()
{  
  /*
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      R(iPart,iBead) += P(iPart,iBead)*dt/M(iBead) + 0.5*F(iPart,iBead)*dt*dt/M(iBead);
      PutInBox( R(iPart,iBead) );
    }
  }
    
  UpdateF(nF);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*(F(iPart,iBead) + nF(iPart,iBead))*dt;
    }
  }
  
  F = nF;
  */
    
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      R(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
      PutInBox( R(iPart,iBead) );
    }
  }
    
  UpdateF(F);  
}

// Update the Force for every bead of every particle
void Paths::UpdateF( field<rowvec>& FX )
{
  rowvec dR1, dR2;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // Handle iBead = 0 case
    dR1 = R(iPart,bL(1)) - R(iPart,0);
    dR2 = R(iPart,nBead-1) - R(iPart,0);
    PutInBox(dR1);
    PutInBox(dR2);
    FX(iPart,0) = M(0)*wp2*(dR1 + dR2) - oneOvernBead*getgradV(iPart,0);
    
    // Do the rest of the time slices
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      dR1 = R(iPart,bL(iBead+1)) - R(iPart,iBead);
      dR2 = R(iPart,iBead-1) - R(iPart,iBead);
      PutInBox(dR1);
      PutInBox(dR2);
      FX(iPart,iBead) = M(iBead)*wp2*(dR1 + dR2) - oneOvernBead*getgradV(iPart,iBead);
    }
    
  }
}

///////////////////////
/* Staging Functions */
///////////////////////

// The Verlet time-stepping algorithm, 'dt' is the time step.
// See eq 116, Ref 1
void Paths::takeStepStage()
{  
  /*
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      U(iPart,iBead) += P(iPart,iBead)*dt/M(iBead) + 0.5*F(iPart,iBead)*dt*dt/M(iBead);
      PutInBox( U(iPart,iBead) );
    }
  }
  
  UtoRStage();
  UpdateFStage(nF);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*(F(iPart,iBead) + nF(iPart,iBead))*dt;
    }
  }
  
  F = nF;
  */
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      U(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
      PutInBox( U(iPart,iBead) );
    }
  }
  
  UtoRStage();
  UpdateFStage(F);
}

// Assign actual positions, going from ui's to xi's
// See eq 12.6.6, Ref 2
void Paths::UtoRStage()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    R(iPart, 0) = U(iPart, 0);
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      R(iPart, iBead) = U(iPart, iBead) + ((1.0*iBead)/(iBead + 1.0))*R(iPart, bL(iBead+1)) + (1.0/(iBead + 1.0))*U(iPart, 0);
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

