#include "shab-pimd.h"

using namespace std;

int main (int argc, char* argv[])
{
  // Inputs
  cout << "\nSimulation Settings:\n";
  int nPart = atoi(argv[1]); // Number of particles
  cout << "Number of Particles: " << nPart << "\n";
  int nD = atoi(argv[2]); // Dimension
  cout << "Dimension: " << nD << "\n";
  int nBead = atoi(argv[3]); // Number of time slices
  cout << "Number of Time Slices: " << nBead << "\n";
  double beta = atof(argv[4]); // Inverse temperature (kb = 1)
  cout << "Inverse temperature (beta): " << beta << "\n";
  double dt = atof(argv[5]); // Time step of simulation
  cout << "Simulation Time Step: " << dt << "\n";
  double L = atof(argv[6]); // Simulation box size
  cout << "Simulation Box Size: " << L << "\n";
  double eSteps = atoi(argv[7]); // Number of Equilibration Sweeps
  cout << "Number of Equilibration Sweeps: " << eSteps << "\n";
  double rSteps = atoi(argv[8]); // Number of Recording Sweeps
  cout << "Number of Recording Sweeps: " << rSteps << "\n";

  bool useStage = atoi(argv[9]); // Use Staging
  cout << "Staging?: " << useStage << "\n";

  bool useNH = atoi(argv[10]); // Use Nose-Hoover Thermostat
  cout << "Nose-Hoover?: " << useNH << "\n";
  int nNH = atoi(argv[11]); // Length of Nose-Hoover Thermostat
  cout << "Nose-Hoover Length: " << nNH << "\n";
  int SYOrder = atoi(argv[12]); // Order of Suzuki-Yoshida Factorization
  cout << "Order of Suzuki-Yoshida Factorization: " << SYOrder << "\n";
  int nNHsteps = atoi(argv[13]); // Number of Nose-Hoover Steps
  cout << "Number of Nose-Hoover Steps: " << nNHsteps << "\n";

  // Random Seed
  srand ( time(NULL) );
  
  // Measurements
  bool measureScalars = 1;

  // Output files
  fstream scalarTrace;   
  char scalarFormat[] = "data/traces/scalarTrace-%d-%d-%d-%g-%g-%g-%d-%d-%d-%d.dat";
  char scalarFile[sizeof scalarFormat+100];
  sprintf(scalarFile,scalarFormat,nPart,nD,nBead,beta,dt,L,eSteps,rSteps,useStage,useNH,nNH,SYOrder,nNHsteps);  
  scalarTrace.open (scalarFile, ios::out | ios::trunc);
  
  // COUT FORMATTING
  cout << scientific << setprecision(4);
    
  // Intialise paths
  Paths path(nPart,nD,nBead,beta,dt,L,useStage,useNH,nNH,SYOrder,nNHsteps);   
  
  // Initialize observables;
  double PE, VE, R, R2; 
  PE = path.getPE();
  VE = path.getVE();
  R = path.getR();
  R2 = path.getR2();
  cout << "\nInitial Measurements:\n";
  cout << "PE: " << PE << " VE: " << VE << " R: " << R << " R2: " << R2 << "\n"; 
  scalarTrace << "t PE VE R R2\n";
           
  // Reset Scalars
  PE = 0.0;
  VE = 0.0;
  R = 0.0;
  R2 = 0.0;
  
  // Main Simulation Loop
  int totSteps = eSteps + rSteps;
  int block = 1;
  int nBlock = 0;
  int perSkip = totSteps/10;
  cout << "\nRunning Simulation:\n";
  for(unsigned int t = 0; t < totSteps; t++) {
    
    // Verlet Step
    if (useStage) path.takeStepStage();
    else path.takeStep();
  
    if(measureScalars && t>=eSteps) {
    
      // Compute Scalars
      PE += path.getPE();
      VE += path.getVE();
      R += path.getR();
      R2 += path.getR2();  
            
      // Blocking      
      if ((t%block)==0) {  
         
        // Output Scalars
        scalarTrace << t/block << " " << PE/block << " " << VE/block << " " << R/block << " " << R2/block << "\n"; 
           
        // Reset Scalars
        PE = 0.0;
        VE = 0.0;
        R = 0.0;
        R2 = 0.0;
           
        nBlock += 1;
      }
    }
    
    // Percentage Complete
    if((t%perSkip)==0) cout << 100.0*t/totSteps << "%\n";
  }
  
  // Close files
  scalarTrace.close();
  
  // Compute and Output stats
  if(measureScalars) statsScalars(scalarFile,nBlock);

  cout << "\nDone.\n\n"; 
  
  return 0;
}
