#include "RNG.h"

int seed = (int)time(0); // random seed
CRandomSFMT1 RanGen(seed); // make instance of random number generator
StochasticLib1 sto(seed); // make instance of random library

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
  return RanGen.Random();
}

// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
  return (b-a)*RanGen.Random() + a;
}

// Generate a random integer between 1 and a given value.
// param n the largest value 
// return a uniform random value in [1,...,n]
long unifRand(long n)
{
  return RanGen.IRandom(1,n);
}

// Generate a box muller normal distribution random number
// with mean m and standard deviation s
double normRand(double m, double s)
{				        
  return sto.Normal(m,s);
}
