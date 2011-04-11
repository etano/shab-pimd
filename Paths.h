#ifndef Paths_H
#define Paths_H

#include "StandardLibs.h"       // Standard libraries
#include "armadillo"
#include "RNG.h"

using namespace arma;

class Paths
{
public:
  Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn , const bool stageIn );  // Constructor
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
  double oneOvernBead, oneOvernPartnBead, nDnPartOver2Beta, nDnPartnBeadOver2Beta, mw2, wp2, mwp2, kT; // These are context clear (check Paths constructor)
  
  double m; // Default (harmonic oscillator units)
  double hbar; // Default (harmonic oscillator units)
  double w; // Default (harmonic oscillator units)
  
  // 3D Matrices
  field<rowvec> R; // Positions R
  field<rowvec> P; // Momenta
  field<rowvec> F, nF; // Forces
  
  // System Initialization
  void InitPosition( field<rowvec>& RX );
  void InitMomentum( double T );
  
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
  bool stage; // 1 - Use staging, 0 - Don't use staging
  field<rowvec> U, nU; // Positions U (staging positions)
  void UtoRStage(); // Switch from U to R
  void UpdateFStage( field<rowvec>& FX ); // Update Force with Staging
  
  // Nose-Hoover Thermostat
  bool useNH; // 1 - Use Nose-Hoover, 0 - Don't use Nose-Hoover
  field<rowvec> NHR; // Nose-Hoover Positions R
  field<rowvec> NHP; // Nose-Hoover Momenta P  
  rowvec Q; // Nose-Hoover Masses Q

};

#endif
