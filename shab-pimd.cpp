#include "shab-pimd.h" // Import headers, set global constants, and build system

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

  cout << "\nN: " << nPart << "\nD: " << nD << "\nM: " << nBead << "\nBeta: " << beta << "\nTime Step (s): " << dt << "\nBox Size: " << L << "\n";

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
  Paths path(nPart,nD,nBead,beta,dt,L);   
  
  // Equilibration
  int eSteps = 0;
  int perSkip = eSteps/10;
  cout << "Equilibrating...\n";
  for(unsigned int t = 0; t < eSteps; t++) {
  
    // Verlet Step
    path.takeStep();

    // Percentage Complete    
    if((t%perSkip)==0) cout << 100.0*t/eSteps << "%\n";
  }
  
  // Initialize observables;
  double PE, VE;
  
  // Start Recording Data
  int rSteps = 100;
  perSkip = rSteps/10;
  cout << "Recording Data...\n";
  for(unsigned int t = 0; t < rSteps; t++) {
    
    // Verlet Step
    path.takeStep();
  
    if(measureScalars) {
      // Compute Scalars
      PE = path.getPE();
      VE = path.getVE();
    
      // Output Scalars
      scalarTrace << t << " " << PE << " " << VE << "\n";
      cout << t << " " << PE << " " << VE << "\n";
    }
    
    // Percentage Complete
    if((t%perSkip)==0) cout << 100.0*t/rSteps << "%\n";
  }
  
  // Close files
  scalarTrace.close();
  
  // Output results
  PE = path.getPE();
  VE = path.getVE();
  cout << PE << " " << VE << " " << "\n"; 

  cout << "\n"; 
  
  return 0;
}
