
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
  
  mnBeadOver2Beta2hbar2 = 1.0*nBead/(2.0*beta*beta*hbar*hbar);
  oneOvernBead = 1.0/(1.0*nBead);
  oneOver2Beta = 1.0/(2.0*beta);
  mOmega2 = m*w*w;
  nBeadOver2Beta = 1.0*nBead/(2.0*beta);
  wp = sqrt(nBead*1.0)/(beta*hbar);
  
  // Set bead loop
  try {  
    bL = new int[2*nBead];
  } catch (std::bad_alloc xa) { 
    std::cout << "Allocation Failure\n"; 
    delete[] bL;
    throw;
  }   
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    bL[iBead] = iBead;
    bL[iBead + nBead] = iBead;
  }
  
  // Set arbitrary masses
  M.set_size(nBead);
  if (stage) {
    M(0) = m;
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      M(iBead) = ((iBead + 1.0)/(1.0*iBead)) * m;
    }
  }
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    M(iBead) = m;
  }
  
  // Size vector fields  
  
  // Positions R
  R.set_size(nPart,nBead);
  nR.set_size(nPart,nBead);
  
  // Positions U
  if (stage) {
    U.set_size(nPart,nBead);
    nU.set_size(nPart,nBead);    
  }
  
  // Momenta
  P.set_size(nPart,nBead);
  nP.set_size(nPart,nBead);
  
  // Forces
  F.set_size(nPart,nBead);
  nF.set_size(nPart,nBead);
  
  // Size vectors
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {  

      // Positions R
      R(iPart,iBead).zeros(nD);
      nR(iPart,iBead).zeros(nD);
      
      // Positions U
      if (stage) {
        U(iPart,iBead).zeros(nD);
        nU(iPart,iBead).zeros(nD);
      }

      // Momenta
      P(iPart,iBead).zeros(nD); 
      nP(iPart,iBead).zeros(nD);

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
    UpdateFStage(F, U); 
  } else {
    UpdateF(F, R);
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

  avgP = netP/(1.0*nPart);
  if (nPart < 2) avgP.zeros();
  pScale = sqrt(m*nD*nPart*nBead*T/netE);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
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
// (See eq. 74 from Ref. 1)
double Paths::getPE()
{
  rowvec dR(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dR = Displacement( R(iPart,iBead) , R(iPart,bL[iBead+1]) );
      total += mnBeadOver2Beta2hbar2 * dot(dR,dR) - oneOvernBead * getV(iPart, iBead);
    }
  }
  
  return nBeadOver2Beta - total;
}

// Get Virial Energy Estimator
// (See eqs. 75, 76 from Ref. 1)
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
  
  return oneOver2Beta + oneOvernBead*total;
}

/////////////////////////
/* Potential Functions */
/////////////////////////

// Get Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// V = (1/2) m w^2 r^2
double Paths::getV( const int iPart, const int iBead )
{
  return 0.5 * mOmega2 * dot( R(iPart,iBead) , R(iPart,iBead) );
}

// Get Derivative of Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// dV/dr = m w^2 r
double Paths::getdV( const int iPart, const int iBead )
{
  return mOmega2 * norm( R(iPart,iBead) , 2 );
}

// Get Derivative of Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// dV/dr = m w^2 r
rowvec Paths::getgradV( const int iPart, const int iBead )
{
  return mOmega2 * R(iPart,iBead);
}

//////////////////////////////////
/* Molecular Dynamics Functions */
//////////////////////////////////

// The Verlet time-stepping algorithm, 'dt' is the time step.
void Paths::takeStep()
{  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      nR(iPart,iBead) = R(iPart,iBead) + P(iPart,iBead)*dt/M(iBead) + 0.5*F(iPart,iBead)*dt*dt/M(iBead);
      PutInBox( nR(iPart,iBead) );
    }
  }
  
  UpdateF(nF, nR);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      nP(iPart,iBead) = P(iPart,iBead) + 0.5*(F(iPart,iBead) + nF(iPart,iBead))*dt;
    }
  }
  
  R = nR;  
  P = nP;
  F = nF;
}

// Update the Force for every bead of every particle
void Paths::UpdateF( field<rowvec>& FX , field<rowvec>& RX )
{
  rowvec dR1, dR2;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // Handle iBead = 0 case
    dR1 = RX(iPart,bL[1]) - RX(iPart,0);
    dR2 = RX(iPart,nBead-1) - RX(iPart,0);
    PutInBox(dR1);
    PutInBox(dR2);
    FX(iPart,0) = M(0)*wp*wp*(dR1 + dR2) - oneOvernBead*getgradV(iPart,0);
    
    // Do the rest of the time slices
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      dR1 = RX(iPart,bL[iBead+1]) - RX(iPart,iBead);
      dR2 = RX(iPart,iBead-1) - RX(iPart,iBead);
      PutInBox(dR1);
      PutInBox(dR2);
      FX(iPart,iBead) = M(iBead)*wp*wp*(dR1 + dR2) - oneOvernBead*getgradV(iPart,iBead);
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
      nU(iPart,iBead) = U(iPart,iBead) + P(iPart,iBead)*dt/M(iBead) + 0.5*F(iPart,iBead)*dt*dt/M(iBead);
      PutInBox( nU(iPart,iBead) );
    }
  }
  
  UpdateFStage(nF, nU);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      nP(iPart,iBead) = P(iPart,iBead) + 0.5*(F(iPart,iBead) + nF(iPart,iBead))*dt;
    }
  }
  
  U = nU;  
  P = nP;
  F = nF;

  UtoRStage();
}

// Assign actual positions, going from ui's to xi's
// See eq 87, Ref 1
void Paths::UtoRStage()
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    R(iPart, 0) = U(iPart, 0);
    for (unsigned int iBead = nBead-1; iBead > 0; iBead -= 1) {
      R(iPart, iBead) = U(iPart, iBead) + ((1.0*iBead)/(iBead + 1.0))*R(iPart, bL[iBead+1]) + (1.0/(iBead + 1.0))*U(iPart, 0);
    }  
  }
}

// Update the Force for every bead of every particle
// See eq 93, Ref 1
void Paths::UpdateFStage( field<rowvec>& FX , field<rowvec>& UX )
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    FX(iPart,0) = -getgradVStage(iPart,0);
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      FX(iPart,iBead) = -M(iBead)*wp*wp*UX(iPart,iBead) - getgradVStage(iPart,iBead);
    }    
  }
}

// Get Derivative of Potential for iPart, iBead
// See eq 94, Ref 1
rowvec Paths::getgradVStage( const int iPart, const int iBead )
{
  rowvec total(nD);
  if(iBead==0) {
    total.zeros();
    for (unsigned int jBead = 0; jBead < nBead; jBead += 1)
      total += getgradV(iPart, jBead);
    //total = oneOvernBead * total;
  } else {
    total = getgradV(iPart, iBead) + ((iBead - 1.0)/(1.0*iBead))*getgradV(iPart, iBead-1);
  }
  return total;
}
