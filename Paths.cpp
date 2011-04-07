
#include "Paths.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn )
{
  // Set constants
  nPart = nPartIn;
  nD = nDIn;
  nBead = nBeadIn;
  beta = betaIn;
  dt = dtIn;
  L = LIn;    
  
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
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    M(iBead) = m;
  }
  
  // Size vector fields  
  
  // Positions
  R.set_size(nPart,nBead);
  nR.set_size(nPart,nBead);
  
  // Velocities
  V.set_size(nPart,nBead);
  nV.set_size(nPart,nBead);
  
  // Accelerations
  A.set_size(nPart,nBead);
  nA.set_size(nPart,nBead);
  
  // Forces
  F.set_size(nPart,nBead);
  nF.set_size(nPart,nBead);
  
  // Size vectors
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {  
      R(iPart,iBead).zeros(nD); // Positions
      nR(iPart,iBead).zeros(nD); // Positions
      V(iPart,iBead).zeros(nD); // Velocities
      nV(iPart,iBead).zeros(nD); // Velocities
      A(iPart,iBead).zeros(nD); // Accelerations
      nA(iPart,iBead).zeros(nD); // Accelerations
      F(iPart,iBead).zeros(nD); // Forces  
      nF(iPart,iBead).zeros(nD); // Forces  
    }
  }
  
  // Initialize system
  InitPosition(R);
  InitVelocity(V, 1.0/beta);
  
}

// Paths destructor
Paths::~Paths()
{
}

//////////////////////
/* Intialize System */
//////////////////////

// Initialize Positions
void Paths::InitPosition( field<rowvec>& R )
{  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      R(iPart,iBead).zeros();
    }
  }
}

// Initialize Velocities
void Paths::InitVelocity( field<rowvec>& V , double T )
{
  rowvec newP(nD), netP(nD), avgP(nD);
  double netE, vScale;
  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
  
    netE = 0.0;
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        newP(iD) = unifRand() - 0.5;
      }
      netP += newP;
      netE += dot(newP,newP);
      V(iPart,iBead) = newP;
    }

    avgP = netP/(1.0*nPart);
    if (nPart < 2) avgP.zeros();
    vScale = sqrt(nD*nPart*T/(m*netE*nBead));
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      V(iPart,iBead) = (V(iPart,iBead) - avgP)*vScale;
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
  UpdateF(F, R);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      A(iPart,iBead) = F(iPart,iBead)/m;
      nR(iPart,iBead) = R(iPart,iBead) + V(iPart,iBead)*dt + 0.5*A(iPart,iBead)*dt*dt;
      PutInBox( nR(iPart,iBead) );
    }
  }
  
  UpdateF(nF, nR);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      nA(iPart,iBead) = nF(iPart,iBead)/m;
      nV(iPart,iBead) = V(iPart,iBead) + 0.5*(A(iPart,iBead) + nA(iPart,iBead))*dt;
    }
  }
  
  R = nR;  
  V = nV;
}

// Update the Force for every bead of every particle
// !! RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// !! F_i = -m_i * w_p * w_p ( r_i+1 + r_i-1 - 2 * r_i ) - (1/P) * (dV/dr_i)
void Paths::UpdateF( field<rowvec>& F , field<rowvec>& R )
{
  rowvec dR1, dR2;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // Handle iBead = 0 case
    dR1 = R(iPart,bL[1]) - R(iPart,0);
    dR2 = R(iPart,nBead-1) - R(iPart,0);
    PutInBox(dR1);
    PutInBox(dR2);
    F(iPart,0) = -M(0)*wp*wp*(dR1 + dR2) - oneOvernBead*getgradV(iPart,0);
    
    // Do the rest of the time slices
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      dR1 = R(iPart,bL[iBead+1]) - R(iPart,iBead);
      dR2 = R(iPart,iBead-1) - R(iPart,iBead);
      PutInBox(dR1);
      PutInBox(dR2);
      F(iPart,iBead) = -M(iBead)*wp*wp*(dR1 + dR2) - oneOvernBead*getgradV(iPart,iBead);
    }
    
  }
}
