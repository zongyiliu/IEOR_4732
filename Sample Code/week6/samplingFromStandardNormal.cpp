#include "functions.h"

#include <iostream>
#include <iomanip>
#include <ios>
#include <math.h>
#include <stdio.h>
#include <sys/time.h> 
#include <cstdlib>

// sampling from a standard normal distribution
// not an efficient code

using namespace std;

int main(int argc, char*argv[])
{ 

	const double PI = 4.0*atan(1.0);

	long idum = -126876761;
	int i, iI, jJ;
	double   u_tmp;
	double u1, u2, u3;
	double z, x1, x2, v1, v2, r;
	int flag = 0;

	int nN = 200000000;

	double *z1 = new double[nN];
	double *z2 = new double[nN];
	double *z3 = new double[nN];
	double *z4 = new double[nN];

	struct timeval startTime, endTime; 
	long seconds, useconds; 
	double mtime;


	int flag2 = 1;
	if (flag2 == 1){
		// Box-Muller
		gettimeofday(&startTime, NULL); 
		for (i=0;i<nN/2;i++){
			
			u1 = ran2_mod(&idum);
			u2 = ran2_mod(&idum);

			z1[i]      = sqrt(-2*log(u1))*cos(2*PI*u2);
			z1[i+nN/2] = sqrt(-2*log(u1))*sin(2*PI*u2);
		}
		gettimeofday(&endTime, NULL);
		seconds = endTime.tv_sec - startTime.tv_sec; 
		useconds = endTime.tv_usec - startTime.tv_usec; 
		mtime = ((seconds) * 1000 + useconds/1000.0); 
		cout << " Box-Muller is done. Time elapsed was (milliseconds): " << mtime << endl;

		// Marsaglia Polar
		gettimeofday(&startTime, NULL); 
		for (i=0;i<nN/2;i++){
			do
			{
				v1 = 2*ran2_mod(&idum) - 1;
				v2 = 2*ran2_mod(&idum) - 1;
				r = v1*v1 + v2*v2;
			}
			while (r > 1.0);

			r = sqrt(-2*log(r)/r);
			z2[i] = r*v1;
			z2[i+nN/2] = r*v2;
		}
		gettimeofday(&endTime, NULL); 
		seconds = endTime.tv_sec - startTime.tv_sec; 
		useconds = endTime.tv_usec - startTime.tv_usec; 
		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5; 
		cout << " Marsaglia Polar is done. Time elapsed was (milliseconds): " << mtime << endl;

		// rational approximation
		gettimeofday(&startTime, NULL); 
		for (i=0;i<nN;i++){
			z4[i] = invertingStandardNormalDistribution(ran2_mod(&idum));
		}
		gettimeofday(&endTime, NULL); 
		seconds = endTime.tv_sec - startTime.tv_sec; 
		useconds = endTime.tv_usec - startTime.tv_usec; 
		mtime = ((seconds) * 1000 + useconds/1000.0); 
		cout << " Rational Approximation is done. Time elapsed was (milliseconds): " << mtime << endl;

		// Acceptance-Rejection (not efficient)
		gettimeofday(&startTime, NULL); 
		for (i=0;i<nN;i++){
			flag = 0;
			while (flag==0){
				u1 = ran2_mod(&idum);
				u2 = ran2_mod(&idum);
				x1 = -log(u1);
				x2 = -log(u2);
				if (x2 >(x1-1.0)*(x1-1.0)/2.0)
					flag = 1;
			}
			u3 = ran2_mod(&idum);
			if (u3<0.5){
				z3[i] = x1;
			}
			else{
				z3[i] = -x1;
			}
		}
		gettimeofday(&endTime, NULL); 
		seconds = endTime.tv_sec - startTime.tv_sec; 
		useconds = endTime.tv_usec - startTime.tv_usec; 
		mtime = ((seconds) * 1000 + useconds/1000.0); 
		cout << " Acceptance-Rejection is done. Time elapsed was (milliseconds): " << mtime << endl;
	}
	delete [] z1;
	delete [] z2;
	delete [] z3;
	delete [] z4;

}
