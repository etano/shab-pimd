#include "shab-pimd.h"

using namespace std;

int main (int argc, char* argv[])
{
  // Inputs
  int nPart = atoi(argv[1]); // Number of particles
  int nD = atoi(argv[2]); // Dimension
  int nBead = atoi(argv[3]); // Number of time slices
  double beta = atof(argv[4]); // Inverse temperature (kb = 1)
  double dt = atof(argv[5]); // Time step of simulation
  double L = atof(argv[6]); // Simulation box size
  bool stage = atoi(argv[7]); // Simulation box size

  cout << "\nSimulation Settings:";
  cout << "\nN: " << nPart << "\nD: " << nD << "\nM: " << nBead << "\nBeta: " << beta << "\nTime Step (s): " << dt << "\nBox Size: " << L << "\nStaging?: " << stage << "\n";

  // Random Seed
  srand ( time(NULL) );
  
  // Measurements
  bool measureScalars = 1;

  // Output files
  fstream scalarTrace;   
  char scalarFormat[] = "data/traces/scalarTrace-%d-%d-%d-%g-%g.dat";
  char scalarFile[sizeof scalarFormat+100];
  sprintf(scalarFile,scalarFormat,nPart,nD,nBead,beta,dt);  
  scalarTrace.open (scalarFile, ios::out | ios::trunc);
  
  // COUT FORMATTING
  cout << scientific << setprecision(4); 
    
  // Intialise paths
  Paths path(nPart,nD,nBead,beta,dt,L,stage);   
  
  // Initialize observables;
  double PE, VE; 
  PE = path.getPE();
  VE = path.getVE();
  cout << "\nInitial Measurements:\n";
  cout << "PE: " << PE << " VE: " << VE << " " << "\n"; 
  scalarTrace << "t PE VE R00\n";
  
  // Main Simulation Loop
  int eSteps = 0;
  int rSteps = 100000;
  int totSteps = eSteps + rSteps;
  int perSkip = totSteps/10;
  cout << "\nRunning Simulation:\n";
  for(unsigned int t = 0; t < totSteps; t++) {
    
    // Verlet Step
    if (stage) path.takeStepStage();
    else path.takeStep();
  
    if(measureScalars && t>eSteps) {
      // Compute Scalars
      PE = path.getPE();
      VE = path.getVE();
    
      // Output Scalars
      scalarTrace << t << " " << PE << " " << VE << " " << path.R(0,0)(0) << "\n";
    }
    
    // Percentage Complete
    if((t%perSkip)==0) cout << 100.0*t/rSteps << "%\n";
  }
  
  // Close files
  scalarTrace.close();
  
  // Output results
  PE = path.getPE();
  VE = path.getVE();
  cout << "\nFinal Measurements:\n";
  cout << "PE: " << PE << " VE: " << VE << " " << "\n"; 
  path.R(0,0).print();

  cout << "\nDone.\n\n"; 
  
  return 0;
}
