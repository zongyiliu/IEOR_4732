#include "functions.h"

#include <iostream>
#include <iomanip>
#include <ios>
#include <math.h>
#include <stdio.h>
#include <sys/time.h> 
#include <cstdlib>

// sampling from a standard uniform
// test code

using namespace std;

int main(int argc, char*argv[])
{ 

	struct timeval startTime, endTime; 
	long seconds, useconds; 
	double mtime;

	long idum = -26876761;
	
	/* initialize random seed: */
    srand (time(NULL));
  
	double u1, u2;
	
	int iI, jJ;


	int nN = 20000000;
	int numBins = 20;

	int *freq1 = new int[numBins];
	int *freq2 = new int[numBins];
	for (int i=0;i<numBins;i++){
		freq1[i] = 0;
		freq2[i] = 0;
	}


	gettimeofday(&startTime, NULL);
	// binning the samples
	// counting the frequencies in each bin
	for(int i=0;i<nN;i++){
		u1 = ran2_mod(&idum);
		//u2 = (double)rand() / ((double)RAND_MAX+1.0);
		//
		iI = (int) floor(u1*numBins);
		//jJ = (int) floor(u2*numBins);
		//
		++freq1[iI];
		//++freq2[jJ];
	}
	gettimeofday(&endTime, NULL);
	seconds = endTime.tv_sec - startTime.tv_sec; 
	useconds = endTime.tv_usec - startTime.tv_usec; 
	mtime = ((seconds) * 1000 + useconds/1000.0); 
	cout << " ran2_mod is done. Time elapsed was (milliseconds): " << mtime << endl;
	
	gettimeofday(&startTime, NULL);
	// binning the samples
	// counting the frequencies in each bin
	for(int i=0;i<nN;i++){
		//u1 = ran2_mod(&idum);
		u2 = (double)rand() / ((double)RAND_MAX+1.0);
		//
		//iI = (int) floor(u1*numBins);
		jJ = (int) floor(u2*numBins);
		//
		//++freq1[iI];
		++freq2[jJ];
	}
	gettimeofday(&endTime, NULL);
	seconds = endTime.tv_sec - startTime.tv_sec; 
	useconds = endTime.tv_usec - startTime.tv_usec; 
	mtime = ((seconds) * 1000 + useconds/1000.0); 
	cout << " rand() is done. Time elapsed was (milliseconds): " << mtime << endl;
	
	
	for(int i=0;i<numBins;i++){
		cout << freq1[i]<< " " << freq2[i] << endl;
	}

	delete [] freq1;
	delete [] freq2;

}
