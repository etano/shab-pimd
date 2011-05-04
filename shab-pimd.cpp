#include "shab-pimd.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main (int argc, char* argv[])
{
  // Input File
  char* inputFile = argv[1];
  ifstream inputStream;
  inputStream.open(inputFile);
  if (!inputStream) cout << "There was a problem opening the input file " << inputFile << " for reading.\n";
  string inputFileLabel;
  getline(inputStream, inputFileLabel); // Skip first line, get label
  int nLineSkip;
  if(argv[2] == NULL) nLineSkip = 1;
  else nLineSkip = atoi(argv[2]);

  cout << "\nRunning simulation " << nLineSkip << " from " << inputFile << " ( " << inputFileLabel << " ) ...\n\n";

  // Inputs
  int nPart; // Number of particles
  int nD; // Dimension
  int nBead; // Number of time slices
  double beta; // Inverse temperature (kb = 1)
  double dt; // Time step of simulation
  double L; // Simulation box size
  int eSteps; // Number of Equilibration Sweeps
  int rSteps; // Number of Recording Sweeps
  int transformation; // Variable Transformation: 0 - None, 1 - Staging, 2 - Normal Mode
  int thermostat; // 0 - No Thermostat, 1 - Nose-Hoover, 2 - Langevin
  int nNH; // Length of Nose-Hoover Thermostat
  int SYOrder; // Order of Suzuki-Yoshida Factorization
  int nNHsteps; // Number of Nose-Hoover Steps

  for (unsigned int iLine = 0; iLine < nLineSkip; iLine += 1) {
    inputStream >> nPart;
    inputStream >> nD;
    inputStream >> nBead;
    inputStream >> beta;
    inputStream >> dt;
    inputStream >> L;
    inputStream >> eSteps;
    inputStream >> rSteps;
    inputStream >> transformation;
    inputStream >> thermostat;
    inputStream >> nNH;
    inputStream >> SYOrder;
    inputStream >> nNHsteps;    
  }
  inputStream.close();

  bool useStage = 0; // Use Staging
  if (transformation==1) useStage = 1;
  bool useNormal = 0; // Use Normal Mode
  if (transformation==2) useNormal = 1;

  cout << "Number of Particles: " << nPart << "\n";
  cout << "Dimension: " << nD << "\n";
  cout << "Number of Time Slices: " << nBead << "\n";
  cout << "Inverse temperature (beta): " << beta << "\n";
  cout << "Simulation Time Step: " << dt << "\n";
  cout << "Simulation Box Size: " << L << "\n";
  cout << "Number of Equilibration Sweeps: " << eSteps << "\n";
  cout << "Number of Recording Sweeps: " << rSteps << "\n";
  cout << "Variable Transformation: " << transformation << "\n";
  cout << "Thermostat: " << thermostat << "\n";
  cout << "Nose-Hoover Length: " << nNH << "\n";
  cout << "Order of Suzuki-Yoshida Factorization: " << SYOrder << "\n";
  cout << "Number of Nose-Hoover Steps: " << nNHsteps << "\n";
  
  // COUT FORMATTING
  cout << scientific << setprecision(4);
    
  // Intialise paths
  Paths path(nPart,nD,nBead,beta,dt,L,transformation,thermostat,nNH,SYOrder,nNHsteps);     

  // Form Output String
  char outputFormat[] = "-%d-%d-%d-%3.1f-%3.1f-%3.1f-%d-%d-%d-%d-%d-%d-%d.dat";
  char outputFile[sizeof outputFormat];
  sprintf(outputFile,outputFormat,nPart,nD,nBead,beta,dt,L,eSteps,rSteps,transformation,thermostat,nNH,SYOrder,nNHsteps); 
  string outputString(outputFile);

  // Scalar Observables Output File
  string scalarFile = "data/traces/scalarTrace" + outputString;
  cout << "\nOutputting scalar data to " << scalarFile << ".\n";
  fstream scalarTrace;
  scalarTrace.open (scalarFile.c_str(), ios::out | ios::trunc);
  scalarTrace << "t PE VE R R2 BS\n";
  
  // Initialize Scalar Observables;
  bool measureScalars = 1;
  double PE = 0.0;
  double VE = 0.0;
  double R = 0.0;
  double R2 = 0.0;
  double BS = 0.0;

  // Pair Correlation Output File
  string GrrFile = "data/traces/GrrTrace" + outputString;
  cout << "\nOutputting Grr data to " << GrrFile << "\n";
  fstream GrrTrace;
  GrrTrace.open (GrrFile.c_str(), ios::out | ios::trunc);
  GrrTrace << "Grr\n";

  // Initialize Pair Correlation
  bool measureGrr = 1;
  const int nPair = choose(nPart, 2); // nPart choose 2
  vec Grr(nPair);
  Grr.zeros();

  // Density Output File
  string RDenFile = "data/traces/RDenTrace" + outputString;
  cout << "\nOutputting Density data to " << RDenFile << "\n";
  fstream RDenTrace;
  RDenTrace.open (RDenFile.c_str(), ios::out | ios::trunc);
  for (unsigned int iD = 0; iD < nD; iD += 1) {
    RDenTrace << "r ";
  }
  RDenTrace << "\n";

  // Initialize Density
  bool measureDensity = 1;
  field<rowvec> RDen(nPart);
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    RDen(iPart).zeros(nD);
  }
  
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
  
    // Compute and Update Values
    if (t > eSteps-1) {
    
      // Scalars
      if (measureScalars) {
        PE += path.getPE();
        VE += path.getVE();
        R += path.getR();
        R2 += path.getR2();  
        BS += path.getBS();
      }

      // Density
      if (measureDensity) path.UpdateDensity(RDen);
      
      // Pair Correlation
      if(measureGrr) path.UpdateGrr(Grr);
      
      // Output and Reset Values      
      if ((t%block)==0) {
         
        // Scalars
        if (measureScalars) {
          scalarTrace << t/block << " " << PE/block << " " << VE/block << " " << R/block << " " << R2/block << " " << BS/block << "\n"; 
          PE = 0.0;
          VE = 0.0;
          R = 0.0;
          R2 = 0.0;
          BS = 0.0;
        }
      
        // Density
        if (measureDensity) {
          for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
            for (unsigned int iD = 0; iD < nD; iD += 1) {
              RDenTrace << RDen(iPart)(iD)/block << " ";
            }
            RDenTrace << "\n";
            RDen(iPart).zeros();
          }
        }

        // Pair Correlation
        if (measureGrr) {
          for (unsigned int iPair = 0; iPair < nPair; iPair += 1) {
            GrrTrace << Grr(iPair)/block << "\n";
          }
          Grr.zeros();
        }         
      
        nBlock += 1;
      }
    }
    
    // Percentage Complete
    if((t%perSkip)==0) cout << 100.0*t/totSteps << "%\n";
  }
  
  // Close files
  scalarTrace.close();
  
  // Compute and Output Scalar Stats
  if(measureScalars) statsScalars(scalarFile.c_str(),nBlock);
  
  cout << "\nDone.\n\n"; 

  return 0;
}
