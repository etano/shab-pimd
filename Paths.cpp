
#include "Paths.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const bool useStageIn , const bool useNHIn , const int nNHIn )
{
  // Set constants
  nPart = nPartIn; // # of Particles
  nD = nDIn; // # of Dimensions
  nBead = nBeadIn; // # of Beads
  beta = betaIn; // Inverse Temperature
  dt = dtIn; // MD Time Step
  L = LIn; // Size of Simulation Box
  useStage = useStageIn; // Use Staging
  useNH = useNHIn; // Use Nose-Hoover Thermostat
  nNH = nNHIn; // Length of Nose-Hoover Thermostat
  
  m = 1.0; // Particle Mass, Default 1
  hbar = 1.0; // Plank's Constant/2pi, Default 1
  k = 1.0; // Boltzmann Constant, Default 1
   
  // Harmonic Oscillator Constants
  w = 1.0;
  mw2 = m*w*w;
  
  // System Constants
  kT = 1.0/beta;
  wp = sqrt(nBead*1.0)/(beta*hbar);
  wp2 = wp*wp;
  mwp2 = m*wp*wp;  
  oneOvernBead = 1.0/(1.0*nBead);
  oneOvernPartnBead = 1.0/(1.0*nPart*nBead);
  nDnPartOver2Beta = 1.0*nD*nPart/(2.0*beta);
  nDnPartnBeadOver2Beta = 1.0*nD*nPart*nBead/(2.0*beta);  
  
  // Set bead loop
  bL.set_size(2*nBead);
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    bL(iBead) = iBead;
    bL(iBead + nBead) = iBead;
  }
  
  ///////////////////////
  /* Initialize masses */
  ///////////////////////
  
  // Masses M
  // See eq 12.6.10, Ref 2
  M.set_size(nBead);  
  M(0) = m;
  for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {  
    if (useStage) M(iBead) = ((iBead + 1.0)/(1.0*iBead)) * m;
    else M(iBead) = m;  
  }
  
  // Masses (Nose-Hoover) Q
  // See eq 12.6.14, Ref 2 
  if (useNH) {
    Q.set_size(nBead);
    Q(0) = kT/(w*w);  
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1)
      Q(iBead) = kT/wp2;
  }
  
  //////////////////////////////
  /* Initialize Vector Fields */
  //////////////////////////////
  
  // Positions R
  R.set_size(nPart,nBead);
  InitPosition(R);
  
  // Positions (Staging) U
  if (useStage) {  
    U.set_size(nPart,nBead);
    InitPosition(U);
    UtoRStage();
  }
  
  // Positions (Nose-Hoover) NHR
  if (useNH) {
    NHR = new field<rowvec>[nBead];
    for (unsigned int iNH = 0; iNH < nNH; iNH += 1) {
      NHR[iNH].set_size(nPart,nBead);
      InitPosition(NHR[iNH]);
    }
  }
  
  // Momenta P
  P.set_size(nPart,nBead);
  InitMomentum(P, M);
  
  // Momenta (Nose-Hoover) NHP
  if (useNH) {
    NHP = new field<rowvec>[nBead];
    for (unsigned int iNH = 0; iNH < nNH; iNH += 1) {
      NHP[iNH].set_size(nPart,nBead);      
      InitMomentum(NHP[iNH], Q);  
    }
  }
  
  // Forces F
  F.set_size(nPart,nBead);
  if (useStage) UpdateFStage(F); 
  else UpdateF(F);

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
      X(iPart,iBead).zeros(nD);
    }
  }
}

// Initialize Momentum
void Paths::InitMomentum( field<rowvec>& Mom , rowvec& Mass )
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
      Mom(iPart,iBead) = newP;
    }
  }

  avgP = netP/(1.0*nPart*nBead);
  if (nPart < 2) avgP.zeros();
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      pScale = sqrt(Mass(iBead)*nD*nPart*nBead/(beta*netE));
      Mom(iPart,iBead) = (Mom(iPart,iBead) - avgP)*pScale;      
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

// Position Squared Estimator
double Paths::getR2()
{
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      total += dot( R(iPart,iBead) , R(iPart,iBead) );
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

