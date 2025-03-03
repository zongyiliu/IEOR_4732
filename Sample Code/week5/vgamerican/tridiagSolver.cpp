
#include <iostream>
#include <iomanip>
#include <ios>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

//convention 1, ..., nSize
//zero left empty

void triDiagSolver(double *lowerD ,double *diag, double *upperD, double *B, int nSize){

	double Xmult;
	int i;
	double *L, *D, *U;
	int index;

	L = new double[nSize+1];
	D = new double[nSize+1];
	U = new double[nSize+1];

	for (i = 1; i <= nSize; i++) {
		L[i] = lowerD[i];
		D[i] = diag[i];
		U[i] = upperD[i];
	}

	for (i = nSize-1; i >= 1; i--) {
		Xmult = U[i] / D[i+1];
		D[i]  = D[i] - Xmult * L[i+1];
		B[i]  = B[i] - Xmult * B[i+1];
	}
	B[1] = B[1] / D[1];

	for (i = 2; i <= nSize; i++)
		B[i] = (B[i] - L[i] * B[i-1]) / D[i];

	delete [] L;
	delete [] D;
	delete [] U;

}