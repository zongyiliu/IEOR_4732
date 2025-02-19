# 這個代碼如果用intergral 的數的話, 能跑出正確的 FFT 結果 (符合課本給的數值, 這個我都測了)

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import cmath
import math
import time

# Fixed Parameters
S0 = 1900
K = 2100
k = math.log(K)
r = 0.02
q = 0.0187

# Parameters for FFT and FrFFT 

n_FFT = 9
N_FFT = 2**n_FFT

n_FrFFT = 9
N_FrFFT = 2**n_FrFFT

N = 512

#step-size
eta = 0.25
# damping factor
alpha = 0.4

# step-size in log strike space
lda_FFT = (2*math.pi/N_FFT)/eta # lda is fixed under FFT
lda_FrFFT = 0.001 # lda is an adjustable parameter under FrFFT, 


#Choice of beta
beta = np.log(S0)-N*lda_FFT/2
#beta = np.log(S0)-N*lda_FrFFT/2
#beta = np.log(K)

#model-specific Parameters
model = 'BlackScholes'

params = []     
if (model == 'GBM'):
    
    sig = 0.30
    params.append(sig);

elif model == 'BlackScholes':  
    sigma = 0.36  # Volatility  
    params.append(sigma)  
    
elif (model == 'VG'):
    
    sig = 0.3
    nu = 0.5
    theta = -0.4
    #
    params.append(sig);
    params.append(nu);
    params.append(theta);


    
elif (model == 'Heston'):
    
    kappa = 2.0
    theta = 0.05
    sig = 0.30
    rho = -0.70
    v0 = 0.04
    #
    params.append(kappa)
    params.append(theta)
    params.append(sig)
    params.append(rho)
    params.append(v0)

def generic_CF(u, params, S0, r, q, T, model):
    
    if (model == 'GBM'):
        
        sig = params[0]
        mu = np.log(S0) + (r-q-sig**2/2)*T
        a = sig*np.sqrt(T)
        phi = np.exp(1j*mu*u-(a*u)**2/2)

    elif (model == 'BlackScholes'):
        sigma = params[0]  # Volatility

        mu = np.log(S0) + (r - q - 0.5 * sigma**2) * T  # Drift
        a = sigma * np.sqrt(T)  # Standard deviation over the maturity period
        phi = np.exp(1j * mu * u - 0.5 * a**2 * u**2)  # Characteristic function for Black-Scholes

        
    elif(model == 'Heston'):
        
        kappa  = params[0]
        theta  = params[1]
        sigma  = params[2]
        rho    = params[3]
        v0     = params[4]
        
        tmp = (kappa-1j*rho*sigma*u)
        g = np.sqrt((sigma**2)*(u**2+1j*u)+tmp**2)
        
        pow1 = 2*kappa*theta/(sigma**2)
        
        numer1 = (kappa*theta*T*tmp)/(sigma**2) + 1j*u*T*r + 1j*u*math.log(S0)
        log_denum1 = pow1 * np.log(np.cosh(g*T/2)+(tmp/g)*np.sinh(g*T/2))
        tmp2 = ((u*u+1j*u)*v0)/(g/np.tanh(g*T/2)+tmp)
        log_phi = numer1 - log_denum1 - tmp2
        phi = np.exp(log_phi)
        
        #g = np.sqrt((kappa-1j*rho*sigma*u)**2+(u*u+1j*u)*sigma*sigma)
        #beta = kappa-rho*sigma*1j*u
        #tmp = g*T/2
        
        #temp1 = 1j*(np.log(S0)+(r-q)*T)*u + kappa*theta*T*beta/(sigma*sigma)
        #temp2 = -(u*u+1j*u)*v0/(g/np.tanh(tmp)+beta)
        #temp3 = (2*kappa*theta/(sigma*sigma))*np.log(np.cosh(tmp)+(beta/g)*np.sinh(tmp))
        
        #phi = np.exp(temp1+temp2-temp3);
        

    elif (model == 'VG'):
        
        sigma  = params[0];
        nu     = params[1];
        theta  = params[2];

        if (nu == 0):
            mu = np.log(S0) + (r-q - theta -0.5*sigma**2)*T
            phi  = np.exp(1j*u*mu) * np.exp((1j*theta*u-0.5*sigma**2*u**2)*T)
        else:
            mu  = np.log(S0) + (r-q + np.log(1-theta*nu-0.5*sigma**2*nu)/nu)*T
            phi = np.exp(1j*u*mu)*((1-1j*nu*theta*u+0.5*nu*sigma**2*u**2)**(-T/nu))

    return phi
def evaluateIntegral(params, S0, K, r, q, T, alpha, eta, N, model):
    
    # Just one strike at a time
    # no need for Fast Fourier Transform
    
    # discount factor
    df = math.exp(-r*T)
    
    sum1 = 0
    for j in range(N):
        nuJ = j*eta
        psi_nuJ = df*generic_CF(nuJ-(alpha+1)*1j, params, S0, r, q, T, model)/((alpha + 1j*nuJ)*(alpha+1+1j*nuJ))
        if j == 0:
            wJ = (eta/2)
        else:
            wJ = eta
        sum1 += np.exp(-1j*nuJ*k)*psi_nuJ*wJ
        
    cT_k = (np.exp(-alpha*k)/math.pi)*sum1
    
    return np.real(cT_k) 

def genericFFT(params, S0, K, r, q, T, alpha, eta, n, model):
    
    N = 2**n
    
    # step-size in log strike space
    lda = (2*np.pi/N)/eta
    
    #Choice of beta
    #beta = np.log(S0)-N*lda/2
    #beta = np.log(K)
    
    # forming vector x and strikes km for m=1,...,N
    km = np.zeros((N))
    xX = np.zeros((N))
    
    # discount factor
    df = math.exp(-r*T)
    
    nuJ = np.arange(N)*eta
    psi_nuJ = generic_CF(nuJ-(alpha+1)*1j, params, S0, r, q, T, model)/((alpha + 1j*nuJ)*(alpha+1+1j*nuJ))
    
    for j in range(N):  
        km[j] = beta+j*lda
        if j == 0:
            wJ = (eta/2)
        else:
            wJ = eta
            
        xX[j] = np.exp(-1j*beta*nuJ[j])*df*psi_nuJ[j]*wJ
     
    yY = np.fft.fft(xX)
    cT_km = np.zeros((N))  
    for i in range(N):
        multiplier = np.exp(-alpha*km[i])/math.pi
        cT_km[i] = multiplier*np.real(yY[i])
    
    return km, cT_km

def genericFrFFT(params, S0, K, r, q, T, alpha, eta, n, lda, model):
    
    N = 2**n
    gamma = eta*lda/(2*math.pi)

    #Choice of beta
    #beta = np.log(S0)-N*lda/2
    beta = np.log(K)

    # initialize x, y, z, and cT_km
    km = np.zeros((N))
    x = np.zeros((N))
    y = np.zeros((2*N), dtype=np.complex128)
    z = np.zeros((2*N), dtype=np.complex128)
    cT_km = np.zeros((N)) 

    # discount factor
    df = math.exp(-r*T)

    # compute x
    nuJ = np.arange(N)*eta
    psi_nuJ = generic_CF(nuJ-(alpha+1)*1j, params, S0, r, q, T, model)/((alpha + 1j*nuJ)*(alpha+1+1j*nuJ))

    for j in range(N):  
        km[j] = beta+j*lda
        if j == 0:
            wJ = (eta/2)
        else:
            wJ = eta
        x[j] = np.exp(-1j*beta*nuJ[j])*df*psi_nuJ[j]*wJ

    # set up y
    for i in range(N):
        y[i] = np.exp(-1j*math.pi*gamma*i**2)*x[i]
    y[N:] = 0

    # set up z
    for i in range(N):
        z[i] = np.exp(1j*math.pi*gamma*i**2)
    z[N:] = z[:N][::-1]

    # compute xi_hat
    xi_hat = np.fft.ifft(np.fft.fft(y) * np.fft.fft(z))

    # compute call prices
    for i in range(N):
        cT_km[i] = np.exp(-alpha*(beta + i*lda))/math.pi * (np.exp(-1j*math.pi*gamma*i**2)*xi_hat[i]).real

    return km, cT_km

print(' ')
print('===================')
print('Model is %s' % model)
print('-------------------')
    
T = 0.25

# FFT
print(' ')
start_time = time.time()
km, cT_km = genericFFT(params, S0, K, r, q, T, alpha, eta, n_FFT, model)
#cT_k = cT_km[0]
cT_k = np.interp(k, km, cT_km)

elapsed_time = time.time() - start_time
    
#cT_k = np.interp(np.log(k), km, cT_km)
print("Option via FFT: for strike %s the option premium is %6.4f" % (np.exp(k), cT_k))
#print("Option via FFT: for strike %s the option premium is %6.4f" % (np.exp(k), cT_km[0]))
print('FFT execution time was %0.7f' % elapsed_time)

# FrFFT
print(' ')
start_time = time.time()
km, cT_km = genericFrFFT(params, S0, K, r, q, T, alpha, eta, n_FrFFT, lda_FrFFT, model)
#cT_k = cT_km[0]
cT_k = np.interp(k, km, cT_km)

elapsed_time = time.time() - start_time
    
#cT_k = np.interp(np.log(), km, cT_km)
print("Option via FrFFT: for strike %s the option premium is %6.4f" % (np.exp(k), cT_k))
#print("Option via FFT: for strike %s the option premium is %6.4f" % (np.exp(k), cT_km[0]))
print('FrFFT execution time was %0.7f' % elapsed_time)


# Integral
print(' ')
start_time = time.time()
cT_k = evaluateIntegral(params, S0, K, r, q, T, alpha, eta, N, model)
elapsed_time = time.time() - start_time
print("Option via Integration: for strike %s the option premium is %6.4f" % (np.exp(k), cT_k))
print('Evaluation of integral time was %0.7f' % elapsed_time)
