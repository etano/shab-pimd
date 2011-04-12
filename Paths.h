#ifndef Paths_H
#define Paths_H

#include "StandardLibs.h"       // Standard libraries
#include "armadillo"
#include "RNG.h"

using namespace arma;

class Paths
{
public:
  Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const bool useStageIn , const bool useNHIn , const int nNHIn );  // Constructor
  ~Paths(); // Destructor

  // Observables Functions
  double getPE(); // Get "Primitive" Energy estimator
  double getVE(); // Get "Virial" Energy estimator
  double getR(); // Get Position estimator
  double getR2(); // Get Position Squared estimator
  
  // Molecular Dynamics Functions
  void takeStep(); // Take a step
  void takeStepStage(); // Take a step (using staging)
  
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
  double oneOvernBead, oneOvernPartnBead, nDnPartOver2Beta, nDnPartnBeadOver2Beta, wp2, mwp2, kT; // These are context clear (check Paths constructor)
  
  double m; // Particle Mass, Default 1
  double hbar; // Plank's Constant/2pi, Default 1
  double k; // Boltzmann Constant, Default 1
  
  double w; // Harmonic Oscillator Frequency
  double mw2; // m * w^2
  
  // 3D Matrices
  field<rowvec> R; // Positions R
  field<rowvec> P; // Momenta
  field<rowvec> F, nF; // Forces
  
  // System Initialization
  void InitPosition( field<rowvec>& X );
  void InitMomentum( field<rowvec>& Mom , rowvec& Mass );
  
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
  
  // Masses
  rowvec M;
  
  // Bead Loop
  ivec bL; 

  // Staging
  bool useStage; // 1 - Use staging, 0 - Don't use staging
  field<rowvec> U, nU; // Positions U (staging positions)
  void UtoRStage(); // Switch from U to R
  void UpdateFStage( field<rowvec>& FX ); // Update Force with Staging
  
  // Nose-Hoover Thermostat
  bool useNH; // 1 - Use Nose-Hoover, 0 - Don't use Nose-Hoover
  int nNH; // Length of Nose-Hoover chain
  field<rowvec> *NHR; // Nose-Hoover Positions R
  field<rowvec> *NHP; // Nose-Hoover Momenta P  
  rowvec Q; // Nose-Hoover Masses Q

};

#endif
