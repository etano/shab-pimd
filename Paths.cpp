#include "Paths.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn )
{
  // Set constants
  nPart = nPartIn;
  nD = nDIn;
  nBead = nBeadIn;
  beta = betaIn;    
  dt = dtIn;
  
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
  
  // Set up matrices  
  R.zeros(nBead,nD,nPart); // Positions
  nR.zeros(nBead,nD,nPart); // Positions
  V.zeros(nBead,nD,nPart); // Velocities
  nV.zeros(nBead,nD,nPart); // Velocities
  A.zeros(nBead,nD,nPart); // Accelerations
  nA.zeros(nBead,nD,nPart); // Accelerations
  F.zeros(nBead,nD,nPart); // Forces  
  nF.zeros(nBead,nD,nPart); // Forces  
  
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
void Paths::InitPosition( cube& R )
{  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        R(iBead,iD,iPart) = 0.0;    
      }
    }
  }
}

// Initialize Velocities
void Paths::InitVelocity( cube& V , double T )
{
  rowvec netP(nD);
  double netE, newP, vscale;
  
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
  
    netP.zeros();
    netE = 0.0;
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        newP = unifRand() - 0.5;
        netP(iD) += newP;
        netE += newP*newP;
        V(iBead,iD,iPart) = newP;        
      }
    }

    netP = netP*(1.0/(1.0*nPart));
    vscale = sqrt(nD*nPart*T/(m*netE));
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        V(iBead,iD,iPart) = (V(iBead,iD,iPart) - netP(iD))*vscale;
      }
    }
    
  }
}

/////////////////
/* OBSERVABLES */
/////////////////

// Get Primitive Energy Estimator
// (See eq. 74 from Ref. 1)
double Paths::getPE()
{
  rowvec dr(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dr = R.slice(iPart).row(iBead) - R.slice(iPart).row(bL[iBead+1]);
      total += mnBeadOver2Beta2hbar2 * norm(dr,2) - oneOvernBead * getV(iPart, iBead);
    }
  }
  
  return nBeadOver2Beta - total;
}

// Get Virial Energy Estimator
// (See eqs. 75, 76 from Ref. 1)
double Paths::getVE()
{
  rowvec rc(nD), dr(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    rc.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) rc += R.slice(iPart).row(iBead);
    rc = rc * oneOvernBead;
    
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dr = R.slice(iPart).row(iBead) - rc;
      total += 0.5 * norm(dr,2) * getdV(iPart, iBead) + getV(iPart, iBead);
    }
  }
  
  return oneOver2Beta + oneOvernBead*total;
}

/////////////////////////
/* Potential Functions */
/////////////////////////

// Get Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL, i.e. (1/2) m w^2 r^2 <---------------------
double Paths::getV( const int iPart, const int iBead )
{
  return 0.5 * mOmega2 * norm(R.slice(iPart).row(iBead),2);
}

// Get Derivative of Potential for iPart, iBead
// RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL, i.e. m w^2 |r| <---------------------
double Paths::getdV( const int iPart, const int iBead )
{
  return mOmega2 * norm(R.slice(iPart).row(iBead),1);
}

//////////////////////////////////
/* Molecular Dynamics Functions */
//////////////////////////////////

// The Verlet time-stepping algorithm, 'dt' is the time step.
void Paths::takeStep()
{
  UpdateF(F, R);
  A = F/m;
  nR = R + V*dt + 0.5*A*dt*dt;
  UpdateF(nF, nR);
  nA = nF/m;
  nV = V + 0.5*(A + nA)*dt;
  R = nR;
  V = nV;
}

// Update the Force for every bead of every particle
// !! RIGHT NOW THIS IS ONLY FOR A HARMONIC POTENTIAL
// !! F_i = -m_i * w_p * w_p ( r_i+1 + r_i-1 - 2 * r_i ) - (1/P) * (dV/dr_i)
void Paths::UpdateF( cube& F , cube& R )
{
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    // Handle iBead = 0 case
    F.slice(iPart).row(0) = -M(0) * wp * wp *
                            (R.slice(iPart).row(1) + R.slice(iPart).row(nBead-1) - 2.0*R.slice(iPart).row(0)) - 
                            oneOvernBead * getdV(iPart,0);
    // Do the rest of the time slices
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      F.slice(iPart).row(iBead) = -M(iBead) * wp * wp * 
                                  (R.slice(iPart).row(bL[iBead+1]) + R.slice(iPart).row(iBead-1) - 2.0*R.slice(iPart).row(iBead)) - 
                                  oneOvernBead*getdV(iPart,iBead);
    }
  }
}
