#ifndef STATS_H
#define STATS_H

#include "StandardLibs.h"
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>


int factorial(int n);

double choose ( const int N , const int M );

double getMean ( const std::vector<double>& data );

double getVar ( const std::vector<double>& data );

double getC( const std::vector<double>& data , int k , int N , double mean , double var );

double getKappa( const std::vector<double>& data );

double getError( const std::vector<double>& data );

void statsScalars ( const char* efile , int nblock );

#endif
