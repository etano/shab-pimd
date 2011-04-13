
#include "Paths.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const bool useStageIn , const bool useNHIn , const int nNHIn , const int nSYIn )
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
  nSY = nSYIn; // Number of Suzuki-Yoshida Weights
  
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
  
  //////////////////////////////
  /* Initialize Vector Fields */
  //////////////////////////////
  
  // Positions R
  R.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      R(iPart,iBead).zeros(nD);
    }
  }
  
  // Positions (Staging)
  if (useStage) {  
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
  }
  
  // Momenta P
  P.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    P(iPart,0).zeros(nD);
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      P(iPart,iBead).zeros(nD);
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        P(iPart,iBead)(iD) = normRand(0.0,M(iBead)*kT); // Maxwell-Boltzmann Distribution
      }
    }    
  }
  
  // Forces F
  F.set_size(nPart,nBead);
  if (useStage) UpdateFStage(F);
  else UpdateF(F);
  
  //////////////////////////////////
  /* Initialize Nose-Hoover Chain */
  //////////////////////////////////  
    
  if (useNH) {
  
    // Masses (Nose-Hoover)
    // See eq 12.6.14, Ref 2 
    NHM.set_size(nNH);
    NHM(0) = kT/wp2;  
    for (unsigned int iNH = 1; iNH < nNH; iNH += 1)
      NHM(iNH) = kT/wp2;
  
    // Positions (Nose-Hoover)
    NHR.set_size(nPart,nBead);
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
        NHR(iPart,iBead).zeros(nNH);
      }
    }
  
    // Momenta (Nose-Hoover)
    NHP.set_size(nPart,nBead);
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
        NHP(iPart,iBead).zeros(nNH);
        for (unsigned int iNH = 0; iNH < nNH; iNH += 1) {
          NHP(iPart,iBead)(iNH) = normRand(0.0,NHM(iNH)*kT); // Maxwell-Boltzmann Distribution
        }
      }
    }
  
    // Forces (Nose-Hoover)
    NHF.set_size(nPart,nBead);
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
        NHF(iPart,iBead).zeros(nNH);
      }
    }
    
    // Suzuki-Yoshida Scheme Weights (Nose-Hoover) NHd
    NHd.set_size(nSY);
    for (unsigned int iSY = 0; iSY < nSY; iSY += 1) {
      if (iSY==2) NHd(iSY) = 1.0 - 4.0/(4.0 - pow(4.0,1.0/3.0));
      else NHd(iSY) = 1.0/(4.0 - pow(4.0,1.0/3.0));
      NHd(iSY) *= dt;
    }
    
  }

}

// Paths destructor
Paths::~Paths()
{
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
		
  if (useNH) NHThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      R(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
      PutInBox( R(iPart,iBead) );
    }
  }
    
  UpdateF(F); 
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }  
  
  if (useNH) NHThermostat();
  
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
  if (useNH) NHThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      U(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
      PutInBox( U(iPart,iBead) );
    }
  }
  
  UtoRStage();
  UpdateFStage(F);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }
  
  if (useNH) NHThermostat();
  
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

////////////////////////////
/* Nose-Hoover Thermostat */
////////////////////////////

void Paths::NHThermostat()
{ 
		
  // See eqs 12.6.14 & 4.11.17, Ref 2
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {      
      for (unsigned int iSY = 0; iSY < nSY; iSY += 1) {
      
        NHF(iPart,iBead)(nNH-1) = NHP(iPart,iBead)(nNH-2) * NHP(iPart,iBead)(nNH-2) / NHM(nNH-2) - kT;
        NHP(iPart,iBead)(nNH-1) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(nNH-1);
        
        for (unsigned int iNH = nNH-2; iNH > 0; iNH -= 1) {
          NHP(iPart,iBead)(iNH) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(iNH+1) / NHM(iNH+1));
          NHF(iPart,iBead)(iNH) = NHP(iPart,iBead)(iNH-1) * NHP(iPart,iBead)(iNH-1) / NHM(iNH-1) - kT;
          NHP(iPart,iBead)(iNH) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(iNH);
          NHP(iPart,iBead)(iNH) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(iNH+1) / NHM(iNH+1));
        }
        
        NHP(iPart,iBead)(0) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(1) / NHM(1));
        NHF(iPart,iBead)(0) = dot( P(iPart,iBead) , P(iPart,iBead) ) / M(iBead) - nD * kT;
        NHP(iPart,iBead)(0) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(0);
        NHP(iPart,iBead)(0) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(1) / NHM(1));
        
        for (unsigned int iD = 0; iD < nD; iD += 1) {
          P(iPart,iBead)(iD) *= exp(-0.5 * NHd(iSY) * NHP(iPart,iBead)(0) / NHM(0));  
        }
        
        // Don't really need this! Would only want it if checking conserved quantity.
        //for (unsigned int iNH = 0; iNH < nNH; iNH += 1) {
        //  NHR(iPart,iBead)(iNH) += -0.5 * NHd(iSY) * NHP(iPart,iBead)(iNH) / NHM(iNH);
        //}
        
        NHP(iPart,iBead)(0) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(1) / NHM(1));
        NHF(iPart,iBead)(0) = dot( P(iPart,iBead) , P(iPart,iBead) ) / M(iBead) - nD * kT;
        NHP(iPart,iBead)(0) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(0);
        NHP(iPart,iBead)(0) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(1) / NHM(1));
        
        for (unsigned int iNH = 1; iNH < nNH-1; iNH += 1) {
          NHP(iPart,iBead)(iNH) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(iNH+1) / NHM(iNH+1));
          NHF(iPart,iBead)(iNH) = NHP(iPart,iBead)(iNH-1) * NHP(iPart,iBead)(iNH-1) / NHM(iNH-1) - kT;
          NHP(iPart,iBead)(iNH) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(iNH);
          NHP(iPart,iBead)(iNH) *= exp(-0.125 * NHd(iSY) * NHP(iPart,iBead)(iNH+1) / NHM(iNH+1));
        }
        
        NHF(iPart,iBead)(nNH-1) = NHP(iPart,iBead)(nNH-2) * NHP(iPart,iBead)(nNH-2) / NHM(nNH-2) - kT;
        NHP(iPart,iBead)(nNH-1) += 0.25 * NHd(iSY) * NHF(iPart,iBead)(nNH-1);
        
      }      		    
    }  
  }
  
}

