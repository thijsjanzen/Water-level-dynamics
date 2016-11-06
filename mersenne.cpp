/**************************   mersenne.cpp   **********************************
* Author:        Agner Fog
* Date created:  2001
* Last modified: 2008-11-16
* Project:       randomc.h
* Platform:      Any C++
* Description:
* Random Number generator of type 'Mersenne Twister'
*
* This random number generator is described in the article by
* M. Matsumoto & T. Nishimura, in:
* ACM Transactions on Modeling and Computer Simulation,
* vol. 8, no. 1, 1998, pp. 3-30.
* Details on the initialization scheme can be found at
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*
* Further documentation:
* The file ran-instructions.pdf contains further documentation and 
* instructions.
*
* Copyright 2001-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

#include "randomc.h"
#include <cmath>

CRandomMersenne rndgen(5); //the one random number generator

void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}

uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((double)(uint32_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}

double CRandomMersenne::normal(double m, double s)
{
	double normal_x1;                   // first random coordinate (normal_x2 is member of class)
	double w;                           // radius
		
	if (normal_x2_valid) {              // we have a valid result from last call
		normal_x2_valid = 0;
		return normal_x2 * s + m;
	}
		
	// make two normally distributed variates by Box-Muller transformation
	do {
		normal_x1 = 2. * Random() - 1.;
		normal_x2 = 2. * Random() - 1.;
		w = normal_x1*normal_x1 + normal_x2*normal_x2;
	} while (w >= 1. || w < 1E-30);
		
	w = std::sqrt(std::log(w)*(-2./w));
	normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
	normal_x2_valid = 1;                // save normal_x2 for next call
	return normal_x1 * s + m;
}


void set_seed(int seed)
{
	rndgen.RandomInit(seed);
}

double uniform()
{
	return rndgen.Random();
}

int random_number(int n)
{
	return rndgen.IRandom(0,n-1);
}

double normal(double m, double s)
{
	return rndgen.normal(m,s);
}


double Expon(double lambda)
{
	return log(1 - rndgen.Random()) / (-1.0 * lambda);
}








