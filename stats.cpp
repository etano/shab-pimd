#include "stats.h"

double getMean ( const std::vector<double>& data ) {
  int N = data.size();
  double mean = 0.0;
  for (unsigned int i = 0; i < N; i += 1) mean += data[i]/(1.0*N);
  return mean;
}

double getVar ( const std::vector<double>& data ) {
  int N = data.size();
  std::vector<double> data2;
  for (unsigned int i = 0; i < N; i += 1) data2.push_back(data[i]*data[i]);
  double mean = getMean(data);
  double mean2 = getMean(data2);
  double var = (1.0*N/(N-1.0))*(mean2 - mean*mean);
  return var;
}

double getC( const std::vector<double>& data , int k , int N , double mean , double var ) {
  double delta[N];
  double tot = 0.0;
  for (unsigned int i = 0; i < N; i += 1) delta[i] = data[i] - mean;
  for (unsigned int i = 0; i < N-k; i += 1) tot += delta[i]*delta[i+k];
  return tot/((N-k)*var);
}

double getKappa( const std::vector<double>& data ) {
  int N = data.size();
  bool done = 0;
  int k = 1;
  double mean = getMean(data);
  double var = getVar(data);
  double c, csum = 1.0;
  while ((!done) && (k < N)) {
    c = getC(data,k,N,mean,var);
    if (c < 0.0) done = 1;
    else csum += 2.0*c;
    k += 1;
  }
  return csum;
}

double getError( const std::vector<double>& data ) {
  int N = data.size();
  double var = getVar(data);
  double kappa = getKappa(data);
  return sqrt(kappa*var/N);
}

void statsScalars ( char* scalarFile , int nblock ) {

  std::ifstream scalarTrace;

  int count; 
  double t, PE, VE, R00, R002;
  std::vector<double> tvec, PEvec, VEvec, R00vec, R002vec;
    
  scalarTrace.open (efile);
  
  count = 0;
  while ( count < nblock ) { //!energyTrace.eof() ) { // keep reading until end-of-file !indata.eof()
    scalarTrace >> t; 
    scalarTrace >> PE;
    scalarTrace >> VE;
    scalarTrace >> R00;
    scalarTrace >> R002;
    
    tvec.push_back(t);
    PEvec.push_back(PE);
    VEvec.push_back(PE);
    R00vec.push_back(R00);
    R002vec.push_back(R002);
    
    count += 1;
  }    
  
  std::cout << "\nScalar Estimates:\n";
  std::cout << "PE: " << getMean(PEvec) << " (" << getError(PEvec) << ")\n";
  std::cout << "VE: " << getMean(VEvec) << " (" << getError(VEvec) << ")\n";
  std::cout << "R00: " << getMean(R00vec) << " (" << getError(R00vec) << ")\n";
  std::cout << "R002: " << getMean(R002vec) << " (" << getError(R002vec) << ")\n";
  
  // Close files
  scalarTrace.close();
 
}
