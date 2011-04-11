#include "Stats.h"

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
  double t, PE, VE, R, R2;
  std::string st, sPE, sVE, sR, sR2;
  std::vector<double> tvec, PEvec, VEvec, Rvec, R2vec;
    
  scalarTrace.open(scalarFile);
  
  count = 0;
  while ( count < nblock ) { //!energyTrace.eof() ) { // keep reading until end-of-file !indata.eof()
    
    if (count > 0) {
      scalarTrace >> t; 
      scalarTrace >> PE;
      scalarTrace >> VE;
      scalarTrace >> R;
      scalarTrace >> R2;
      tvec.push_back(t);
      PEvec.push_back(PE);
      VEvec.push_back(VE);
      Rvec.push_back(R);
      R2vec.push_back(R2);      
    } else {
      scalarTrace >> st; 
      scalarTrace >> sPE;
      scalarTrace >> sVE;
      scalarTrace >> sR;
      scalarTrace >> sR2;      
    }
    
    count += 1;
  }    
  
  std::cout << "\nScalar Estimates:\n";
  std::cout << "PE: " << getMean(PEvec) << " (" << getError(PEvec) << ")\n";
  std::cout << "VE: " << getMean(VEvec) << " (" << getError(VEvec) << ")\n";
  std::cout << "R: " << getMean(Rvec) << " (" << getError(Rvec) << ")\n";
  std::cout << "R2: " << getMean(R2vec) << " (" << getError(R2vec) << ")\n";
  
  // Close files
  scalarTrace.close();
 
}
