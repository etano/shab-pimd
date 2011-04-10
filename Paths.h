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
  
  // Molecular Dynamics Functions
  void takeStep(); // Take a step
  void takeStepStage(); // Take a step (using staging)
  
  // 3D Matrices
  field<rowvec> R, nR; // Positions R
  field<rowvec> U, nU; // Positions U (staging positions)
  field<rowvec> P, nP; // Velocities
  field<rowvec> F, nF; // Forces
  
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
  double mnBeadOver2Beta2hbar2, oneOvernBead, oneOver2Beta, mOmega2, nBeadOver2Beta; // These are context clear (check Paths constructor)
  double wp; // Defined frequency
  
  double m; // Default (harmonic oscillator units)
  double hbar; // Default (harmonic oscillator units)
  double w; // Default (harmonic oscillator units)
  
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
  void UpdateF( field<rowvec>& FX , field<rowvec>& RX );
  
  // Masses
  rowvec M;
  
  int *bL; // Bead Loop

  // Staging
  bool stage; // 1 - Use staging, 0 - Don't use staging
  void UtoRStage(); // Switch from U to R
  rowvec getgradVStage( const int iPart , const int iBead ); // Get Gradient of Potential for iPart, iBead (using staging)
  void UpdateFStage( field<rowvec>& FX , field<rowvec>& UX );

};

#endif
