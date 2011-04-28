#include "shab-pimd.h"

using namespace std;

int main (int argc, char* argv[])
{
  // Input File
  char* inputFile = argv[1];
  cout << "\nReading input from " << inputFile << ".\n";
  ifstream inputStream;
  inputStream.open(inputFile);
  if (!inputStream) {
    cout << "There was a problem opening the input file " << inputFile << " for reading.\n";
  }

  // Inputs
  int nPart; // Number of particles
  int nD; // Dimension
  int nBead; // Number of time slices
  double beta; // Inverse temperature (kb = 1)
  double dt; // Time step of simulation
  double L; // Simulation box size
  double eSteps; // Number of Equilibration Sweeps
  double rSteps; // Number of Recording Sweeps
  bool useStage; // Use Staging
  bool useNormal; // Use Normal Mode
  bool useNH; // Use Nose-Hoover Thermostat
  int nNH; // Length of Nose-Hoover Thermostat
  int SYOrder; // Order of Suzuki-Yoshida Factorization
  int nNHsteps; // Number of Nose-Hoover Steps

  int nLineSkip;
  if(argv[2] == NULL) nLineSkip = 1;
  else nLineSkip = atoi(argv[2]);
  for (unsigned int iLine = 0; iLine < nLineSkip; iLine += 1) {
    inputStream >> nPart;
    inputStream >> nD;
    inputStream >> nBead;
    inputStream >> beta;
    inputStream >> dt;
    inputStream >> L;
    inputStream >> eSteps;
    inputStream >> rSteps;
    inputStream >> useStage;
    inputStream >> useNormal;
    inputStream >> useNH;
    inputStream >> nNH;
    inputStream >> SYOrder;
    inputStream >> nNHsteps;    
  }
  inputStream.close();

  cout << "Number of Particles: " << nPart << "\n";
  cout << "Dimension: " << nD << "\n";
  cout << "Number of Time Slices: " << nBead << "\n";
  cout << "Inverse temperature (beta): " << beta << "\n";
  cout << "Simulation Time Step: " << dt << "\n";
  cout << "Simulation Box Size: " << L << "\n";
  cout << "Number of Equilibration Sweeps: " << eSteps << "\n";
  cout << "Number of Recording Sweeps: " << rSteps << "\n";
  cout << "Staging?: " << useStage << "\n";
  cout << "Normal Mode?: " << useNormal << "\n";
  cout << "Nose-Hoover?: " << useNH << "\n";
  cout << "Nose-Hoover Length: " << nNH << "\n";
  cout << "Order of Suzuki-Yoshida Factorization: " << SYOrder << "\n";
  cout << "Number of Nose-Hoover Steps: " << nNHsteps << "\n";

  // Random Seed
  srand ( time(NULL) );
  
  // Measurements
  bool measureScalars = 1;

  // Output files
  fstream scalarTrace;   
  char scalarFormat[] = "data/traces/scalarTrace-%d-%d-%d-%g-%g-%g-%g-%d-%d-%d-%d.dat";
  char scalarFile[sizeof scalarFormat+100];
  sprintf(scalarFile,scalarFormat,nPart,nD,nBead,beta,dt,L,eSteps,rSteps,useStage,useNormal,useNH,nNH,SYOrder,nNHsteps);  
  scalarTrace.open (scalarFile, ios::out | ios::trunc);
  
  // COUT FORMATTING
  cout << scientific << setprecision(4);
    
  // Intialise paths
  Paths path(nPart,nD,nBead,beta,dt,L,useStage,useNormal,useNH,nNH,SYOrder,nNHsteps);   
  
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
    else if (useNormal) path.takeStepNormal();
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
