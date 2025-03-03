/*
    
	ran2_mod.cpp - This module is taken from NUMERICAL RECIPES IN C by Press,
            Teukolsky, Vetterling, and Flannery (PTVF); 1992; page 282.

       Long period (> 2 * 10^18) random number generator of L'Ecuyer with Bays-Durham
       shuffle and added safeguards. Returns a uniform random deviate between 0.0 and
       1.0 (exclusive of the endpoint values). Call with idum a negative integer to
       initialize; thereafter, do not alter idum between successive deviates in a
       sequence. RNMX should approximate the largest floating value that is less than 1.

     Notes on modifications:

       - All "double" variables were "float" in PTVF.
       - PTVF uses function prototyping (ie. argument types are specified
        in the parameter list) while
        this module specifies variable formats after the function
        declaration.

     The difference between ran2_mod.c and ran2.c is that *idum is the code is replaced with local_idum

*/

#include "functions.h"

#include <iostream>
#include <iomanip>
#include <ios>
#include <math.h>
#include <stdio.h>
#include <string>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2_mod(long &idum)
{
  int j;
  long k;
  long local_idum;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  //local_idum=*idum;
  local_idum=idum;

  if (local_idum <= 0) {
    if (-(local_idum) < 1) 
        local_idum=1;
    else 
        local_idum = -(local_idum);

    idum2=(local_idum);

    for (j=NTAB+7;j>=0;j--) {
      k=(local_idum)/IQ1;
      local_idum = IA1*(local_idum-k*IQ1)-k*IR1;
      if (local_idum < 0) 
          local_idum += IM1;
      
      if (j < NTAB) 
          iv[j] = local_idum;
    }
    iy=iv[0];
  }
  k = (local_idum)/IQ1;
  local_idum = IA1*(local_idum-k*IQ1)-k*IR1;
  if (local_idum < 0)
      local_idum += IM1;
  
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) 
      idum2 += IM2;
  
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = local_idum;
  if (iy < 1) 
      iy += IMM1;

  //*idum = local_idum;
  idum = local_idum;

  if ((temp=AM*iy) > RNMX) 
      return RNMX;
  else 
      return temp;
}
