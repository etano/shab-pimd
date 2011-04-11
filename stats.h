#ifndef STATS_H
#define STATS_H

#include "StandardLibs.h"
#include <iomanip>

double getMean ( const std::vector<double>& data );

double getVar ( const std::vector<double>& data );

double getC( const std::vector<double>& data , int k , int N , double mean , double var );

double getKappa( const std::vector<double>& data );

double getError( const std::vector<double>& data );

void statsScalars ( char* efile , int ntype , int nblock );

#endif
