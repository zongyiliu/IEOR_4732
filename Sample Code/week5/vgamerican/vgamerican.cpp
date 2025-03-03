#include "stdafx.h"

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

#define MAXIT 100
#define EULER .577215664901532860606512
#define FPMIN 1.0e-30
#define EPS 1.0e-7

int  tridiagPut(double *, double *, double *, double *, double *, double, int);
double max(double, double);
double expint(double);

int main(){

double   *L, *D, *U, *B, *W, *x;
double   *ei_n, *ei_p, *exp_n, *exp_p;
double   A;
double   *s_tau;
double   cs0;
double   r = .0541, q = .012;
double   *t;
double   Dx, Dt;
int      i, j, k, ir;
int      index, counter = 0;
double   cn, cp, c;
double   Bp, Bn;

double   sig=0.20722, nu=0.50215, theta=-0.22898;
double   lambda_n, lambda_p;
double   omega;

double   T=0.56164;
double   S0=1369.41, K=1300.0; 
double   SMIN = 10, SMAX = 2000;

int  N=2000, M=100;


FILE * fout1 ;
fout1=fopen("prices2.dat","w"); 

FILE * fout2 ;
fout2=fopen("freeboundary2.dat","w");    

     Dx = (log(SMAX)-log(SMIN))/N;
     Dt = T/M;

	 L = new double[N];
	 D = new double[N];
	 U = new double[N];
	 B = new double[N];
	 W = new double[N+1];
	 x = new double[N+1];

	 s_tau = new double[M+1];

	 exp_n = new double[N];
     exp_p = new double[N];

	 ei_n = new double[N];
	 ei_p = new double[N];

	 t = new double[M+1];

	 omega    = log(1-theta*nu-sig*sig*nu/2)/nu;
	 lambda_n = sqrt(theta*theta/pow(sig,4)+2/(nu*sig*sig)) + theta/(sig*sig);
     lambda_p = sqrt(theta*theta/pow(sig,4)+2/(nu*sig*sig)) - theta/(sig*sig);

	 cn = Dt/(lambda_n*nu*Dx);
	 cp = Dt/(lambda_p*nu*Dx);

	 c = Dt/nu;

	 Bn = cn*(1-exp(-lambda_n*Dx));
     Bp = cp*(1-exp(-lambda_p*Dx));

	 A = (r-q+omega)*Dt/(2*Dx);

	 // Pre-calculated Vectors
	 for ( k = 1; k <= N-1; k++ ){
		 ei_n[ k] = expint(k*Dx*lambda_n);
		 ei_p[ k] = expint(k*Dx*lambda_p);
		 exp_n[k] = exp(-k*Dx*lambda_n);
         exp_p[k] = exp(-k*Dx*lambda_p);
	 }
     // 
     for ( i = 0; i <= N; i++){
		 x[i] = log(SMIN)+i*Dx;
		 if(exp(x[i])<K)
			 W[i] = K - exp(x[i]);
		 else{
			 if (counter == 0){
				 index=i;
				 counter=1;
			 }
			 W[i] = 0.0;
		 }
	 }
     s_tau[0] = K;
	 t[0] = 0.0;
	 for (j = 1; j <= M; j++) {
		 cout << M-j << "\n" << endl;
		 t[j]=Dt*j;
		 // Boundary Conditions 
		 W[0] = K - exp(x[0]);
		 W[N] = 0;
         // 
		 for (i=1;i<=N-1;i++){
			 B[i]=W[i];
		 }
		 for (i = 1; i <=N-1; i++){
			 if (i==1){
				 D[i] = 1 + r*Dt + Bp + Bn + c*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx));
		         U[i] = -A - Bp;
				 for (k=1;k<=N-i-1;k++)
					 B[i] += c*((W[i+k]-W[i])-k*(W[i+k+1]-W[i+k]))*(ei_p[k]-ei_p[k+1]) + cp*(W[i+k+1]-W[i+k])*(exp_p[k]-exp_p[k+1]);
				 B[i] += c*( K*expint(i*Dx*lambda_n) - exp(x[i])*expint(i*Dx*(lambda_n+1)) );
				 B[i] += -(A - Bn)*( K - exp(x[0]) );  // Boundary Condition Effect 

			 }
		     else if (i==N-1){
				 L[i] = A - Bn;
                 D[i] = 1 + r*Dt + Bp + Bn + c*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx));
				 for (k=1;k<=i-1;k++)
					 B[i] += c*((W[i-k]-W[i])-k*(W[i-k-1]-W[i-k]))*(ei_n[k]-ei_n[k+1]) + cn*(W[i-k-1]-W[i-k])*(exp_n[k]-exp_n[k+1]);
				 B[i] += c*(K*expint(i*Dx*lambda_n) - exp(x[i])*expint(i*Dx*(lambda_n+1)));
			 }
		     else{
				 L[i] =  A - Bn ;
                 D[i] =  1.0 + r*Dt + Bp + Bn + c*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx));
                 U[i] = -A - Bp;
				 for (k=1;k<=N-i-1;k++)
					 B[i] += c*((W[i+k]-W[i])-k*(W[i+k+1]-W[i+k]))*(ei_p[k]-ei_p[k+1]) + cp*(W[i+k+1]-W[i+k])*(exp_p[k]-exp_p[k+1]);
                 for (k=1;k<=i-1;k++)
					 B[i] += c*((W[i-k]-W[i])-k*(W[i-k-1]-W[i-k]))*(ei_n[k]-ei_n[k+1]) + cn*(W[i-k-1]-W[i-k])*(exp_n[k]-exp_n[k+1]);
				 B[i] += c*(K*expint(i*Dx*lambda_n) - exp(x[i])*expint(i*Dx*(lambda_n+1)));
			 }
			 if ( x[i] < x[index]  ){
				 B[i] += Dt*(r*K - q*exp(x[i]));
				 for (k=index-i;k<=N-i-1;k++)
					 B[i] -= c*((W[i+k]-W[i])-k*(W[i+k+1]-W[i+k]))*(ei_p[k]-ei_p[k+1]) + cp*(W[i+k+1]-W[i+k])*(exp_p[k]-exp_p[k+1]);
				 B[i] += c * ( K*expint((index-i)*Dx*lambda_p) - exp(x[i])*expint((index-i)*Dx*(lambda_p-1)) );
			 }
		 }
		 for (i=1;i<=N-1;i++){
			 W[i]=B[i];
		 }
	     index = tridiagPut(L, D, U, W, x, K, N-1);
		 s_tau[j] = exp(x[index]);
	 }
     for(i = 0; i <= N-1; i++){
		 if(x[i]>log(S0)){
			 ir=i;
		 break;
		 }
	 }
     for(i = 0; i <= N; i++)
		 fprintf(fout1, "%g %g \n",exp(x[i]),W[i]);

	 cs0 = (W[ir]-W[ir-1])*(log(S0)-x[ir-1])/Dx + W[ir-1];


     cout << "\n Stock Price  "    << S0          << endl;
	 cout << "\n American Put Value  " << cs0 << "\n" << endl;

	 for(j = 0; j <=M; j++)
		 fprintf(fout2, "%g %g \n",t[j],s_tau[j]);

	 fclose(fout1);
     fclose(fout2);

	 delete [] L;
	 delete [] D;
	 delete [] U;
	 delete [] W;
	 delete [] B;
     delete [] x;

	 delete [] exp_n;
     delete [] exp_p;
	 delete [] ei_n;
	 delete [] ei_p;

	 delete [] s_tau;

	 delete [] t;

     return(0);
     }

int tridiagPut(double *LL ,double *DD, double *UU, double *B, double *x, double K, int NN){

   double Xmult;
   int i,Ic;
   double *L,*D,*U;
   int index;

   L = new double[NN+1];
   D = new double[NN+1];
   U = new double[NN+1];

   for (i = 1; i <= NN; i++) {
     L[i] = LL[i];
     D[i] = DD[i];
     U[i] = UU[i];
   }

     for (i = NN-1; i >= 1; i--) {
       Xmult = U[i] / D[i+1];
       D[i]  = D[i] - Xmult * L[i+1];
       B[i]  = B[i] - Xmult * B[i+1];
	 }
	 Ic=1;
	 B[Ic] = B[Ic] / D[Ic];
	 //
	 while ( B[Ic] <= max(K - exp(x[Ic]),0) ){
		 B[Ic] = max(K - exp(x[Ic]), 0);
		 B[Ic+1]=(B[Ic+1] - L[Ic+1] * B[Ic]) / D[Ic+1];
		 Ic++;
	 }
	 index = Ic;

     for (i = Ic + 1; i <= NN; i++)
		 B[i] = (B[i] - L[i] * B[i-1]) / D[i];

	 delete [] L;
	 delete [] D;
	 delete [] U;

     return index;
}

double max(double a, double b)
{
	double temp;
	if (a>=b)
		temp = a;
	else
		temp = b;
	return temp;
}

double expint(double x)
{
	int n = 1;
	int i,ii,nm1;
	double a, b, c, d, del, fact, h, psi, ans;

	nm1 = n-1;
	if(n<0 || x <0.0 || (x==0 && (n==0 || n==1)))
		cout << "bad arguments in expint" << endl;
	else{
		if (n==0){ 
			ans=exp(-x)/x;
		}
		else {
			if (x==0.0){
				ans=1.0/nm1;
			}
			else{
				if(x>1.0){
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++){
						a = -i*(nm1+i);
						b += 2.0;
						d = 1.0/(a*d+b);
						c = b+a/c;
						del = c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans = h*exp(-x);
							return ans;
						}
					}
					cout << "continued fraction failed in expint" << endl;
				} else {
					ans = (nm1 != 0 ? 1.0/nm1 : -log(x)-EULER);
					fact = 1.0;
					for (i=1;i<MAXIT;i++) {
						fact *= -x/i;
						if ( i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;
							for (ii=1; ii<=nm1; ii++) 
								psi += 1.0/ii;

							del=fact*(-log(x)+psi);
						}
						ans +=del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					cout << "series failed in expint" << endl;
				}
			}
		}
	}
	return ans;
}
