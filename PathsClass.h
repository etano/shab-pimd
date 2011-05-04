#ifndef Paths_H
#define Paths_H

#include "StandardLibs.h"       // Standard libraries
#include "armadillo"
#include "RNG.h"

using namespace arma;

class Paths
{
public:
  Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const int transformationIn , const int thermostatIn , const int nNHIn , const int SYOrderIn , const int nNHstepsIn );  // Constructor
  ~Paths(); // Destructor

  // Observables Functions
  double getPE(); // Get "Primitive" Energy estimator
  double getVE(); // Get "Virial" Energy estimator
  double getR(); // Get Position estimator
  double getR2(); // Get Position Squared estimator
  double getBS(); // Get Bead Spread estimator
  rowvec getR(int iPart); // Get Position estimator for single particle
  void UpdateDensity(field<rowvec>& RDen); // Update Density
  void UpdateGrr(vec& Grr); // Update Pair Correlation
  
  // Molecular Dynamics Functions
  void takeStep(); // Take a step
  void takeStepStage(); // Take a step (using staging)
  void takeStepNormal(); // Take a step (using normal mode)
  
protected:
  //protected things
private:   

  // Constants
  int nPart; // Number of particles
  int nD; // Number of dimensions
  int nBead; // Number of beads
  double beta; // Inverse temperature
  double dt; // Time step  
  double L; // Simulation Box Size
  double wp; // Defined frequency
  double oneOvernPart, oneOvernBead, oneOvernPartnBead, nDnPartOver2Beta, nDnPartnBeadOver2Beta, wp2, mwp2, kT; // These are context clear (check Paths constructor)
  
  double m; // Particle Mass, Default 1
  double hbar; // Plank's Constant/2pi, Default 1
  double k; // Boltzmann Constant, Default 1
  
  double w; // Harmonic Oscillator Frequency
  double mw2; // m * w^2
  
  int transformation; // 0 - No Transformation, 1 - Staging, 2 - Normal Mode
  int thermostat; // 0 - No Thermostat, 1 - Nose-Hoover, 2 - Langevin

  // 3D Matrices
  field<rowvec> R; // Positions R
  field<rowvec> P; // Momenta
  field<rowvec> F, nF; // Forces
    
  // Periodic Boundary Conditions
  void PutInBox( rowvec& Ri );
  double Distance( rowvec& Ri , rowvec& Rj );
  rowvec Displacement( rowvec& Ri , rowvec& Rj );  
  
  // Potential Functions
  double getV( const int iPart, const int iBsead ); // Get Potential for iPart, iBead
  double getdV( const int iPart, const int iBead ); // Get Derivative of Potential for iPart, iBead
  rowvec getgradV( const int iPart, const int iBead ); // Get Gradient of Potential for iPart, iBead
  
  // Molecular Dynamics Functions
  void UpdateF( field<rowvec>& FX ); // Update Force
  
  // Thermostat
  void ApplyThermostat();
  
  // Masses
  rowvec M;
  
  // Bead Loop
  ivec bL; 

  // Staging
  bool useStage; // 1 - Use staging, 0 - Don't use staging
  field<rowvec> U, nU; // Positions (Staging)
  void initStaging(); // Initiate Staging
  void UtoRStage(); // Switch from U to R
  void UpdateFStage( field<rowvec>& FX ); // Update Force with Staging
  
  // Normal Mode
  bool useNormal; // 1 - Use normal mode, 0 - Don't use normal mode
  mat NormO; // Orthgonal Transformation Matrix
  field<rowvec> W, nW; // Positions (Normal Mode)
  void initNormalMode(); // Initiate Normal Mode
  void WtoRNormal(); // Switch from W to R
  void UpdateFNormal( field<rowvec>& FX ); // Update Force with Normal Mode
  
  // Nose-Hoover Thermostat
  bool useNH; // 1 - Use Nose-Hoover, 0 - Don't use Nose-Hoover
  int nNH; // Length of Nose-Hoover chain
  field<rowvec> NHR; // Positions (Nose-Hoover)
  field<rowvec> NHP; // Momenta (Nose-Hoover)
  field<rowvec> NHF; // Forces (Nose-Hoover)
  rowvec NHM; // Masses (Nose-Hoover)
  int SYOrder; // Order of Suzuki-Yoshida factorization
  int nSY; // Number of Suzuki-Yoshida weights
  int nNHsteps; // Number of Nose-Hoover Steps
  rowvec NHd; // Suzuki-Yoshia Weight Factors 
  void initNoseHoover(); // Initiate Nose-Hoover
  void NHThermostat(); // Thermostatting function

  // Langevin Thermostat
  bool useLT; // 1 - Use Langevin Thermostat, 0 - Don't
  rowvec LTOmega; // Normal Mode Frequencies
  rowvec LTGamma; // Langevin Friction Constants
  rowvec LTC1; // Langevin Momentum Multiplier
  rowvec LTC2; // Langevin Force Multiplier
  void initLangevinThermostat(); // Thermostatting function
  void LangevinThermostat(); // Thermostatting function

};

#endif
