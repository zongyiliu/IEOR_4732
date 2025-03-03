import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#import time

a = 0
b = 100

T = 1

N = 100
M = 100

dX = (b-a)/N
dT = T/M

kappa = 2
rho = kappa*dT/(dX*dX)



# Explicit Scheme
def setUpDiagonals(N):
    
    # lower, diag, upper
    l = np.zeros((N,1))
    d = np.zeros((N,1))
    u = np.zeros((N,1))
    
    for i in range(N):
        l[i] = -rho
        d[i] = 1 + 2*rho
        u[i] = -rho
        
        
    return l, d, u

# Implicit Scheme
# Tri-diagonal solver
def diagSolver(l, d, u, b):
    
    '''
    tri-diagonal solver
    '''
    N = len(d)
    lL, dD, uU, bB = map(np.array, (l, d, u, b)) # copy arrays
    for i in range(1, N):
        tmp = lL[i-1]/dD[i-1]
        dD[i] = dD[i] - tmp*uU[i-1] 
        bB[i] = bB[i] - tmp*bB[i-1]
        	    
    x = dD
    x[-1] = bB[-1]/dD[-1]

    for il in range(N-2, -1, -1):
        x[il] = (bB[il]-uU[il]*x[il+1])/dD[il]

    return x

def formTheMatrix(N):
    
    A = np.zeros((N,N))
    
    for i in range(N):
        A[i,i] = 1-2*rho
        if (i<N-1):
            A[i,i+1] = rho
        
        if (i>0):
            A[i,i-1] = rho
            
    return A
    
    
#initial condition    
x = np.zeros((N,1))
u = np.zeros((N,1))
for i in range(N):
    x[i] = i
    if (i<N/2) & (i>N/4):
        u[i] = i-N/4.0
    elif (i>=N/2) & (i<3*N/4):
        u[i] = 3.0*N/4.0-i
        
 
lL, dD, uU = setUpDiagonals(N) 

A = formTheMatrix(N)
    
#plt.figure(figsize=(20,10))  
plt.ion()
plt.figure()

for j in range(M):
    #print('iter = %s' % (str(j)))
    u = np.matmul(A, u)
    #u = diagSolver(lL, dD, uU, u)
    plt.show(block=False)
    plt.clf() # clear figure
    plt.plot(x, u)
    #plt.axis((0,100, -.5, 25))
    plt.grid(True)
    plt.pause(0.01)
    #plt.xlabel('$S_T$')
    #plt.ylabel('lognormal density $f(S_T|S_0)$')

    

    

    
    


