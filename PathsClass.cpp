#include "PathsClass.h"

// Paths constructor
Paths::Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const int transformationIn , const int thermostatIn , const int interactionIn , const int nNHIn , const int SYOrderIn , const int nNHstepsIn )
{
  // Set constants
  nPart = nPartIn; // # of Particles
  nD = nDIn; // # of Dimensions
  nBead = nBeadIn; // # of Beads
  beta = betaIn; // Inverse Temperature
  dt = dtIn; // MD Time Step
  L = LIn; // Size of Simulation Box

  transformation = transformationIn; // 0 - No Transformation, 1 - Staging, 2 - Normal Mode
  thermostat = thermostatIn; // 0 - No Thermostat, 1 - Nose-Hoover, 2 - Langevin
  interaction = interactionIn; // 0 - No Interaction, 1 - Lennard Jones, 2 - Coulomb

  nNH = nNHIn; // Length of Nose-Hoover Thermostat
  SYOrder = SYOrderIn; // Order of Suzuki-Yoshida Factorization
  nNHsteps = nNHstepsIn; // Number of Nose-Hoover Steps per sweep
  
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
  oneOvernPart = 1.0/(1.0*nPart);
  oneOvernBead = 1.0/(1.0*nBead);
  oneOvernPartnBead = 1.0/(1.0*nPart*nBead);
  nDnPartOver2Beta = 1.0*nD*nPart/(2.0*beta);
  nDnPartnBeadOver2Beta = 1.0*nD*nPart*nBead/(2.0*beta); 

  // Lennard-Jones-Constants
	rcut = L; // Set to system Size Needs to be adjusted for Periodic Boundary Calculation
	r0 = 1; // Length scale of interaction
	e0 = 1; // Strength of interaction
	ecut = 0; // Adjust here  
  
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
    M(iBead) = m;  
  }
  
  //////////////////////////////
  /* Initialize Vector Fields */
  //////////////////////////////
  
  // Positions R
  InitR();
  
  // Momenta P
  P.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead).zeros(nD);
      for (unsigned int iD = 0; iD < nD; iD += 1) {
        P(iPart,iBead)(iD) = normRand(0.0,M(iBead)*kT); // Maxwell-Boltzmann Distribution
      }
    }    
  }

  // Interaction 
  Vint.set_size(nPart,nBead);
  Vint.zeros();

  GradVint.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      GradVint(iPart,iBead).zeros(nD);
    }
  }
  
  // Forces F
  F.set_size(nPart,nBead);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      F(iPart,iBead).zeros(nD);
    }
  }
  UpdateF();
    
  ///////////////////////////////
  /* Initialize Transformation */
  /////////////////////////////// 

  useStage = 0;
  useNormal = 0;

  if (transformation == 1) useStage = 1;
  if (transformation == 2) useNormal = 1;

  if (useStage) initStaging();
  if (useNormal) initNormalMode();
  
  ///////////////////////////
  /* Initialize Thermostat */
  ///////////////////////////  

  useNH = 0;
  useLT = 0;
    
  if (thermostat == 1) useNH = 1;
  if (thermostat == 2) useLT = 1;

  if (useNH) initNoseHoover();
  if (useLT) initLangevinThermostat();

}

// Paths destructor
Paths::~Paths()
{
}

///////////////////////////
/* System Initialization */
///////////////////////////

void Paths::InitR()
{
  R.set_size(nPart,nBead);

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      R(iPart,iBead).zeros(nD);
    }
  }

  int nCube = 1;
  while (pow(nCube,nD) < nPart) nCube += 1;
  double rs = nPart*1.0/(1.0*nCube);
  double rOffset = (1.0 - rs)/2.0;

  int iPart;
  for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
    iPart = 0;
    if (nD == 1) {
      for (unsigned int x = 0; x < nCube; x += 1) {
        if (iPart < nPart) {
          R(iPart,iBead)(0) = rs*x - rOffset;
          iPart++;
        }
      }
    }
    if (nD == 2) {
      for (unsigned int x = 0; x < nCube; x += 1) {
        for (unsigned int y = 0; y < nCube; y += 1) {
          if (iPart < nPart) {
            R(iPart,iBead)(0) = rs*x - rOffset;
            R(iPart,iBead)(1) = rs*y - rOffset;
            iPart++;
          }
        }
      }
    }
    if (nD == 3) {
      for (unsigned int x = 0; x < nCube; x += 1) {
        for (unsigned int y = 0; y < nCube; y += 1) {
          for (unsigned int z = 0; z < nCube; z += 1) {
            if (iPart < nPart) {
              R(iPart,iBead)(0) = rs*x - rOffset;
              R(iPart,iBead)(1) = rs*y - rOffset;
              R(iPart,iBead)(2) = rs*z - rOffset;
              iPart++;
            }
          }
        }
      }
    }
  }

}

/////////////////////////
/* Potential Functions */
/////////////////////////

// Update Gradient of Interacting Potential
void Paths::UpdateVint()
{
  rowvec dRij, GradVintTemp;
  double VintTemp;
  for (unsigned int iBead = 0; iBead < nBead; iBead +=1) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      GradVint(iPart,iBead).zeros();
      Vint(iPart,iBead) = 0.0;
    }
    
    for (unsigned int iPart = 0; iPart < nPart-1; iPart += 1) {
      for (unsigned int jPart = iPart+1; jPart < nPart; jPart += 1) {
        dRij = Displacement( R(iPart,iBead) , R(jPart,iBead) );
        // Grad Vint
        GradVintTemp = -1.0 * InteractionForce(dRij);
        GradVint(iPart,iBead) += GradVintTemp;
        GradVint(jPart,iBead) -= GradVintTemp;
        // Vint
        VintTemp = InteractionPotential( norm(dRij,2) );
        Vint(iPart,iBead) += 0.5 * VintTemp;
        Vint(jPart,iBead) += 0.5 * VintTemp;
      }
    }
  }
}

//Calculate Pair-Potential-Force; 1-LJ, 2-Coulomb, 
rowvec Paths::InteractionForce( rowvec& dRij )
{
  if (interaction==0) return 0.0 * dRij;

	double r = norm(dRij,2); 
	double r2 = r*r; 
	if (interaction==1) {
    double r02i = 1.0/(r0*r0);
		r2 = r2*r02i;  //Distance in units of r0
		double r2i = 1.0/r2; //Inverse Distance
		double r6i = r2i*r2i*r2i;
    
		double ff = 48.0*e0*r2i*r6i*(r6i-0.5)*r02i;
		return ff * dRij;
	}
	if (interaction == 2) {
		return -(e0/r2/r) * dRij;
	}
}

//Calculate Pair-Potential; 1-LJ, 2-Coulomb, 
double Paths::InteractionPotential( const double dRij )
{
  if (interaction==0) return 0.0;

  double r2 = dRij * dRij;
	if (interaction==1) {
		r2 = r2/r0;  //Distance in units of r0
		double r2i = 1.0/r2; //Inverse Distance
		double r6i = r2i*r2i*r2i;
    
		double en = 4.0*e0*r6i*(r6i-1.0);
		return en * dRij;
	}
	if (interaction == 2) {
		return e0/dRij;
	}
}

//////////////////////////////////
/* Molecular Dynamics Functions */
//////////////////////////////////

// The Verlet time-stepping algorithm, 'dt' is the time step.
void Paths::takeStep()
{      
		
  if(thermostat) ApplyThermostat();
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
      R(iPart,iBead) += P(iPart,iBead)*dt/M(iBead);
      PutInBox( R(iPart,iBead) );
    }
  }
    
  UpdateF(); 
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      P(iPart,iBead) += 0.5*F(iPart,iBead)*dt;
    }
  }  
  
  if(thermostat) ApplyThermostat();
  
}

// Update the Force for every bead of every particle
void Paths::UpdateF()
{
  if(interaction) UpdateVint();

  rowvec dR1, dR2;
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
  
    // Handle iBead = 0 case
    dR1 = R(iPart,bL(1)) - R(iPart,0);
    dR2 = R(iPart,nBead-1) - R(iPart,0);
    PutInBox(dR1);
    PutInBox(dR2);    
    F(iPart,0) = M(0)*wp2*(dR1 + dR2) - oneOvernBead*getgradV(iPart,0);;
    
    // Do the rest of the time slices
    for (unsigned int iBead = 1; iBead < nBead; iBead += 1) {
      dR1 = R(iPart,bL(iBead+1)) - R(iPart,iBead);
      dR2 = R(iPart,iBead-1) - R(iPart,iBead);
      PutInBox(dR1);
      PutInBox(dR2);
      F(iPart,iBead) = M(iBead)*wp2*(dR1 + dR2) - oneOvernBead*getgradV(iPart,iBead);
    }
    
  }
}

//////////////////////////////////
/* Molecular Dynamics Functions */
//////////////////////////////////

void Paths::ApplyThermostat()
{
  if (useNH) NHThermostat();
  else if (useLT) LangevinThermostat();
}

