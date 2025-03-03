#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXIT 100
#define EULER .577215664901532860606512
#define FPMIN 1.0e-30
#define EPS 1.0e-7

void   tridiag(double *, double *, double *, double *, int);
double expint(double);

main(){

double   *L, *D, *U, *B, *W, *x;
double   *ei_n, *ei_p, *exp_n, *exp_p;
double   Barrier;
double   Strike;
double   *A;
double   UOC;
double   *r, *q;
double   Dx, Dt;
int      i, j, k, ir;
int      ii, jj;
double   cn, cp, c;
double   Bp, Bn;

double   T = 0.92273; 
double   S0 = 139.0906; 
double   sig=0.1952, nu=0.6140, theta=-0.1994;
double   rfr = .0654, div = 0.0116;
double   lambda_n, lambda_p;
double   omega;
 
double   SMIN = 1.0;
double   Rebate = 0.0;

int      N=600, M=370;


FILE * fout ;
fout=fopen("VG5.dat","w"); 

for  (ii=0; ii<15; ii++){
	for (jj = 0; jj <24; jj++){

		cout << ii << " " << jj << endl;

		Barrier = 180-2.5*ii;
		Strike = 120+2.5*jj;

     Dx = (log(Barrier)-log(SMIN))/N;
     Dt = T/M;

	 L = new double[N];
	 D = new double[N];
	 U = new double[N];
	 B = new double[N];
	 W = new double[N+1];
	 x = new double[N+1];

	 exp_n = new double[N];
     exp_p = new double[N];

	 ei_n = new double[N];
	 ei_p = new double[N];


     A = new double[M];
	 r = new double[M];
	 q = new double[M];


	 omega    = log(1-theta*nu-sig*sig*nu/2)/nu;
	 lambda_n = sqrt(theta*theta/pow(sig,4)+2/(nu*sig*sig)) + theta/(sig*sig);
     lambda_p = sqrt(theta*theta/pow(sig,4)+2/(nu*sig*sig)) - theta/(sig*sig);

	 cn = Dt/(lambda_n*nu*Dx);
	 cp = Dt/(lambda_p*nu*Dx);

	 c = Dt/nu;

	 Bn = Dt*(1-exp(-lambda_n*Dx))/(nu*lambda_n*Dx);
     Bp = Dt*(1-exp(-lambda_p*Dx))/(nu*lambda_p*Dx);

	 // Precaculated Vectors
	 for ( k = 1; k <= N-1; k++ ){
		 ei_n[ k] = expint(k*Dx*lambda_n);
		 ei_p[ k] = expint(k*Dx*lambda_p);
		 exp_n[k] = exp(-k*Dx*lambda_n);
         exp_p[k] = exp(-k*Dx*lambda_p);
	 }
     // 


     for ( i = 0; i <= N; i++){
		 x[i] = log(SMIN)+i*Dx;
		 if(i==N)
			 W[i] = Rebate;
		 else{
			 if(exp(x[i])>Strike)
				 W[i] = exp(x[i])-Strike;
		     else
			     W[i]=0.0;
		 }
	 }


	 for (j = M-1; j >= 0; j--) {
		 //cout << j << endl;
		 r[j]=rfr;
		 q[j]=div;
		 A[j]=(r[j]-q[j]+omega)*Dt/(2*Dx);
		 for (i=1;i<=N-1;i++){
			 B[i]=W[i];
		 }

		 for (i = 1; i <=N-1; i++){
			 if (i==1){
				 D[i] = 1 + r[j]*Dt + Bp + Bn + Dt*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx))/nu;
		         U[i] = -A[j] - Bp;

				 for (k=1;k<=N-i-1;k++)
					 B[i] += c*((W[i+k]-W[i])-k*(W[i+k+1]-W[i+k]))*(ei_p[k]-ei_p[k+1]) + cp*(W[i+k+1]-W[i+k])*(exp_p[k]-exp_p[k+1]);
			 }
		     else if (i==N-1){
				 L[i] = A[j] - Bn;
                 D[i] = 1 + r[j]*Dt + Bp + Bn + Dt*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx))/nu;

				 for (k=1;k<=i-1;k++)
					 B[i] += c*((W[i-k]-W[i])-k*(W[i-k-1]-W[i-k]))*(ei_n[k]-ei_n[k+1]) + cn*(W[i-k-1]-W[i-k])*(exp_n[k]-exp_n[k+1]);
			 }
		     else{
				 L[i] =  A[j] - Bn ;
                 D[i] =  1 + r[j]*Dt + Bp + Bn + Dt*(expint((N-i)*lambda_p*Dx)+expint(i*lambda_n*Dx))/nu;
                 U[i] = -A[j] - Bp;

				 for (k=1;k<=N-i-1;k++)
					 B[i] += c*((W[i+k]-W[i])-k*(W[i+k+1]-W[i+k]))*(ei_p[k]-ei_p[k+1]) + cp*(W[i+k+1]-W[i+k])*(exp_p[k]-exp_p[k+1]);
                 for (k=1;k<=i-1;k++)
					 B[i] += c*((W[i-k]-W[i])-k*(W[i-k-1]-W[i-k]))*(ei_n[k]-ei_n[k+1]) + cn*(W[i-k-1]-W[i-k])*(exp_n[k]-exp_n[k+1]);
			 }
		 }
		 B[N-1] = B[N-1] + Rebate * (A[j]+Bp);
		 for (i=1;i<=N-1;i++){
			 W[i] = B[i];
		 }
	     tridiag(L,D,U,W,N-1);
	 }

	 ir=0;
     for(i = 0; i <= N-1; i++){
		 if(x[i]>log(S0)){
			 ir=i;
		 break;
		 }
	 }
     //for(i = 0; i <= N; i++)
	 //	 fprintf(fout, "%g %g \n",exp(x[i]),W[i]);
	 UOC = (W[ir]-W[ir-1])*(log(S0)-x[ir-1])/Dx + W[ir-1];

	 if (jj == 0 && ii != 0){
		 fprintf(fout, "\n");
	 }

     fprintf(fout, "%g ", UOC);


     //cout << "\n Stock Price  "    << S0          << endl;
	 //cout << "\n UOC Call Value  " << cs0 << "\n" << endl;

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

	 delete [] r;
	 delete [] q;
     delete [] A;
	}
}

     fclose(fout);

     return(0);
     }

void tridiag(double *LL ,double *DD, double *UU, double *B, int NN){

   double Xmult;
   int i,Ic;
   double *L,*D,*U;

   L = new double[NN+1];
   D = new double[NN+1];
   U = new double[NN+1];

   for (i=1;i<=NN;i++) {
     L[i] = LL[i];
     D[i] = DD[i];
     U[i] = UU[i];
   }

     for (i = NN-1; i >= 1; i--) {
       Xmult = U[i] / D[i+1];
       D[i]  = D[i] - Xmult * L[i+1];
       B[i]  = B[i] - Xmult * B[i+1];
	 }
	 B[1] = B[1]/D[1];
	 Ic = 1;
     for (i = Ic + 1; i <= NN; i++)
     B[i] = (B[i] - L[i] * B[i-1]) / D[i];

	 delete [] L;
	 delete [] D;
	 delete [] U;

     return;
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
		if (n==0) ans=exp(-x)/x;
		else {
			if (x==0.0) ans=1.0/nm1;

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
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
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


