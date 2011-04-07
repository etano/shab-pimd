#ifndef Paths_H
#define Paths_H

#include "StandardLibs.h"       // Standard libraries
#include "armadillo"
#include "RNG.h"

using namespace arma;

class Paths
{
public:
  Paths( const int nPartIn , const int nDIn , const int nBeadIn , const double betaIn , const double dtIn , const double LIn );  // Constructor
  ~Paths(); // Destructor

  // Observables Functions
  double getPE(); // Get "Primitive" Energy estimator
  double getVE(); // Get "Virial" Energy estimator
  
  // Molecular Dynamics Functions
  void takeStep(); // Take a step
  
  // 3D Matrices
  cube R, nR; // Positions
  cube V, nV; // Velocities
  cube A, nA; // Accelerations
  cube F, nF; // Forces
  
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
  void InitPosition( cube& R );
  void InitVelocity( cube& V , double T );
  
  // Periodic Boundary Conditions
  void PutInBox( rowvec& Ri );
  double Distance( rowvec Ri , rowvec Rj );
  rowvec Displacement( rowvec Ri , rowvec Rj );  
  
  // Potential Functions
  double getV( const int iPart, const int iBead ); // Get Potential for iPart, iBead
  double getdV( const int iPart, const int iBead ); // Get Derivative of Potential for iPart, iBead
  
  // Molecular Dynamics Functions
  void UpdateF( cube& F , cube& R );
  
  // Masses
  rowvec M;
  
  int *bL; // Bead Loop
  
};

#endif
