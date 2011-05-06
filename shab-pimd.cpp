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
  int nBin; // Number of bins per D
  int eSteps; // Number of Equilibration Sweeps
  int rSteps; // Number of Recording Sweeps
  int transformation; // Variable Transformation: 0 - None, 1 - Staging, 2 - Normal Mode
  int thermostat; // 0 - No Thermostat, 1 - Nose-Hoover, 2 - Langevin
  int interaction; // 0 - No Interaction, 1 - Lennard Jones, 2 - Coulomb
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
    inputStream >> nBin;
    inputStream >> eSteps;
    inputStream >> rSteps;
    inputStream >> transformation;
    inputStream >> thermostat;
    inputStream >> interaction;
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
  cout << "# of Bins per D: " << nBin << "\n";
  cout << "Number of Equilibration Sweeps: " << eSteps << "\n";
  cout << "Number of Recording Sweeps: " << rSteps << "\n";
  cout << "Variable Transformation: " << transformation << "\n";
  cout << "Thermostat: " << thermostat << "\n";
  cout << "Interaction: " << interaction << "\n";
  cout << "Nose-Hoover Length: " << nNH << "\n";
  cout << "Order of Suzuki-Yoshida Factorization: " << SYOrder << "\n";
  cout << "Number of Nose-Hoover Steps: " << nNHsteps << "\n";
  
  // COUT FORMATTING
  cout << scientific << setprecision(4);
    
  // Intialise paths
  Paths path(nPart,nD,nBead,beta,dt,L,transformation,thermostat,interaction,nNH,SYOrder,nNHsteps);     

  // Form Output String
  char outputFormat[] = "-%d-%d-%d-%3.1f-%.1g-%3.1f-%d-%d-%d-%d-%d-%d-%d-%d-%d.dat";
  char outputFile[sizeof outputFormat];
  sprintf(outputFile,outputFormat,nPart,nD,nBead,beta,dt,L,nBin,eSteps,rSteps,transformation,thermostat,interaction,nNH,SYOrder,nNHsteps); 
  string outputString(outputFile);

  // Scalar Observables Output File
  string scalarFile = "data/traces/scalarTrace" + outputString;
  cout << "\nOutputting scalar data to " << scalarFile << ".\n";
  fstream scalarTrace;
  scalarTrace.open (scalarFile.c_str(), ios::out | ios::trunc);
  scalarTrace << "n t PE VE R R2 BS\n";
  
  // Initialize Scalar Observables;
  bool measureScalars = 1;
  double PE = 0.0;
  double VE = 0.0;
  double R = 0.0;
  double R2 = 0.0;
  double BS = 0.0;

  // Binning
  double dBin = L*1.0/nBin;

  // Pair Correlation Output File
  string GrrFile = "data/traces/GrrTrace" + outputString;
  cout << "\nOutputting Grr data to " << GrrFile << "\n";
  fstream GrrTrace;
  GrrTrace.open (GrrFile.c_str(), ios::out | ios::trunc);
  for (unsigned int iBin = 0; iBin < nBin; iBin += 1) {
    GrrTrace << iBin*dBin + dBin/2.0 << " ";
  }
  GrrTrace << "\n";

  // Initialize Pair Correlation
  bool measureGrr = 1;
  vec Grr(nBin);
  Grr.zeros();

  // Density Output File
  string RDenFile = "data/traces/RDenTrace" + outputString;
  cout << "\nOutputting Density data to " << RDenFile << "\n";
  fstream RDenTrace;
  RDenTrace.open (RDenFile.c_str(), ios::out | ios::trunc);
  for (unsigned int iBin = 0; iBin < nBin; iBin += 1) {
    RDenTrace << iBin*dBin + dBin/2.0 << " ";
  }
  RDenTrace << "\n";

  // Initialize Density
  bool measureDensity = 1;
  vec RDen(nBin);
  RDen.zeros();

  // Main Simulation Loop
  int totSteps = eSteps + rSteps;
  int block = rSteps/100;
  int nBlock = 0;
  int perSkip = totSteps/10;
  cout << "\nRunning Simulation:\n";

  // Timer
  time_t start, end;
  double timeDif;

  for(unsigned int n = 0; n < totSteps; n++) {

    // Verlet Step
    if (useStage) path.takeStepStage();
    else if (useNormal) path.takeStepNormal();
    else path.takeStep();

    // Compute and Update Values
    if (n > eSteps-1) {

      // Scalars
      if (measureScalars) {
        PE += path.getPE();
        VE += path.getVE();
        R += path.getR();
        R2 += path.getR2();  
        BS += path.getBS();
      }

      // Density
      if (measureDensity) path.UpdateDensity(RDen,dBin);
      
      // Pair Correlation
      if (measureGrr) path.UpdateGrr(Grr,dBin);
      
      // Output and Reset Values      
      if ((n%block)==0) {
         
        // Timer
        if (!nBlock) {
          timeDif = 0.0;
          start = clock();
        }
        else {
          end = clock();
          timeDif = end - start;
          start = end;
        }

        // Scalars
        if (measureScalars) {
          scalarTrace << n-eSteps << " " << timeDif << " "  << PE/block << " " << VE/block << " " << R/block << " " << R2/block << " " << BS/block << "\n"; 
          PE = 0.0;
          VE = 0.0;
          R = 0.0;
          R2 = 0.0;
          BS = 0.0;
        }
      
        // Density
        if (measureDensity) {
          for (unsigned int iBin = 0; iBin < nBin; iBin += 1) {
            RDenTrace << RDen(iBin) << " ";
            RDen(iBin) = 0.0;
          }
          RDenTrace << "\n";
        }

        // Pair Correlation
        if (measureGrr) {
          for (unsigned int iBin = 0; iBin < nBin; iBin += 1) {
            GrrTrace << Grr(iBin) << " ";
            Grr(iBin) = 0.0;
          }
          GrrTrace << "\n";
        }         
      
        nBlock += 1;
      }
    }
    
    // Percentage Complete
    if((n%perSkip)==0) cout << 100.0*n/totSteps << "%\n";
  }
  
  // Close files
  scalarTrace.close();
  
  // Compute and Output Scalar Stats
  if(measureScalars) statsScalars(scalarFile.c_str(),nBlock);
  
  cout << "\nDone.\n\n"; 

  return 0;
}
