"""
Calaix de Sastre
"""

from math import tanh, tan, cos, pi, log10
import numpy as np
from numpy.linalg import solve
from scipy.special import jv, yv
import subprocess
import matplotlib.pyplot as plt
import cmath

ro = 1000
g = 9.81

def Dispersion_T(L, const):
    # Linear dispersion equation, f=0, for travelling modes
    T,h= const
    A= g*T**2/2/pi
    B= 2*pi*h
    f= L-A*tanh(B/L)
    df= 1-(A*B*(tanh(B/L)**2-1))/L**2
    return f,df

def Dispersion_E(L, const):
    # Linear dispersion equation, f=0, for evanescent modes
    T,h= const
    A= g*T**2/2/pi
    B= 2*pi*h
    f= L+A*tan(B/L)
    df= 1-A*B/cos(B/L)**2/L**2
    return f,df

def Newton(func, const, x = 1, feat= (100,1e-6,1e-6)):
    """
    Newton method for solving zeros of non-linear equations
    Where:
    - func=f(x)=0 is the non-linear equation
    - const is a tuple consisting of the f constants
    - x is the first iterate
    - feat is a tuple consisting of the max. number of iterates, the tolerance in x and f
    """
    kmax,Tol_x,Tol_f= feat
    found= 0
    k= 0
    while found == 0 and k < kmax:
        k+= 1 # Update k
        # Evaluating f and df
        f,df= func(x,const)
        # Avoiding dividing by zero
        if abs(df) < Tol_f*1e-3:
            found= 1
            bl= (0,'Local minimum. log|f| = {:.2f}'.format(log10(abs(f)+1e-99)))
        else:
            # Update x
            x0= x
            x+= -f/df
            # Check on convergence
            if abs(x-x0) < Tol_x and abs(f) < Tol_f:
                found= 1
                bl= (1,'Converged in {} iterations and with log|f| = {:.2f}'.format(k,log10(abs(f)+1e-99)))
    # Did not converge
    if found == 0:
        bl= (0,'Did not converge. log|f| = {:.2f}'.format(log10(abs(f)+1e-99)))
    return x,bl

def Bisection(func, const, a, b, feat= (1000,1e-6)) :
    """
    """
    kmax,Tol = feat
    k = 0
    fa = func(a, const)[0]
    fb = func(b, const)[0]
    fc = 1e6
    while abs(fc) > Tol and k < kmax :
        k += 1
        c = (float(a)+float(b))/2.
        fc = func(c, const)[0]
        if (fc < 0. and fb < 0.) or (fc > 0. and fb > 0.) :
            b = c
            fb = fc
        elif (fc < 0. and fa < 0.) or (fc > 0. and fa > 0.) :
            a = c
            fa = fc
        else :
            bl = (0,'fa*fb = {:}'.format(fa*fb))
            break
    if abs(fc) < Tol :
        bl = (1, 'converged')
    elif k >= kmax :
        bl = (0, 'max number of iterations')
    return c, bl

def WNumber(period, depth):
    """
    Calculates wave-number, k0, from wave-period and water depth (T,h).
    Dispersion relationship for travelling modes
    Inputs:
    - period: a float or a list of float numbers. Wave periods
    - depth: float(). It is the water depth.
    Outputs:
    """
    try :
        L = np.array([Newton(Dispersion_T, (p, depth))[0] for p in period], dtype=float)
    except TypeError:
        L = Newton(Dispersion_T, (period, depth))[0]
    return 2*pi/L

def WNumber_E(period, depth, L, step = .001, kmax = 1000):
    """
    Calculates wave-numbers, ke, from wave-period and water depth (T,h).
    Dispersion relationship for evanescent modes
    Inputs:
    - period: a float or a list of float numbers. Wave periods
    - depth: float(). It is the water depth.
    - L (integer):  number of evanescent modes (L).
    Outputs:
    """
    Le = np.zeros((len2(period), L), dtype = float)
    for ip in range(len2(period)) :
        if len2(period) > 1 :
            p = period[ip]
        else :
            p = period
        for l in range(L) :
            k0 = np.pi/2./depth*(2*l+1.)*(1.+1e-6)
            kf = np.pi/2./depth*(2*l+2.)
            a = 2*np.pi/k0
            b = 2*np.pi/kf
            sol = Bisection(Dispersion_E, (p, depth), a, b)
            if sol[1][0] == 1 : # Solution found in (a,b)
                Le[ip,l] = sol[0]
            else :
                print('Fail evanescent')
                break
    return 2*pi/Le

def len2(x):
    try:
        if type(x) == int or type(x) == float or \
        type(x) == np.int32 or type(x) == np.float64: return 1
        else : return len(x)
    except TypeError:
        raise TypeError('Attempted len2(x) but type(x) = {} is not supported in this version'.format(type(x)))

def CylWaveField(X, Y, Z, a, fr, k, h, coord, disregard, convention='W', plane=0):
    Nb = len2(coord)
    Nm = (len2(a)/Nb-1)/2
    mode = range(-Nm,Nm+1)
    XX , YY = np.meshgrid(X, Y, indexing='ij')
    Phi = np.zeros((Nb, XX.shape[0], YY.shape[1]), dtype=complex)
    for i in range(Nb):
        Xi, Yi= coord[i]
        ri = np.sqrt((XX-Xi)**2+(YY-Yi)**2)
        thi = np.arctan2(YY-Yi,XX-Xi)
        Phi[i][ri<=disregard] = np.nan
        ai = a[(2*Nm+1)*i:(2*Nm+1)*(i+1)]
        for m in mode[Nm:]:
            H = jv(m,k*ri)-1j*yv(m,k*ri)
            if convention == 'N' :
                H = np.conj(H)
            if plane == 1:
                Phi[i] += ai[Nm+m]*H.real*np.exp(1j*m*thi)
                if m > 0:
                    Phi[i] += ai[Nm-m]*(-1)**m*H.real*np.exp(1j*-m*thi)
            else :
                Phi[i] += ai[Nm+m]*H*np.exp(1j*m*thi)
                if m > 0:
                    Phi[i] += ai[Nm-m]*(-1)**m*H*np.exp(1j*-m*thi)
    if convention == 'N' :
        Phi *= -1j*g/fr
    else :
        Phi *= 1j*g/fr
    if len2(Z) == 1:
        return np.cosh(k*(Z+h))/np.cosh(k*h)*Phi
    else :
        # Velocity potential
        return np.array([np.cosh(k*(z+h))/np.cosh(k*h)*Phi for z in Z], dtype=complex)

def execute(command):
    """
    execute(['ping', 'localhost'])
    """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''

    # Poll process for new output until finished
    # for line in iter(process.stdout.readline, ""):
    #     if process.poll()!=None:
    #         break
    #     print(line)
    #     output += line.decode()
    while process.poll() is None:
        line = process.stdout.readline()
        print(line)
        output += line.decode()

    print(process.stdout.read())


    process.wait()
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise Exception(command, exitCode, output)

def MotionFreq(Mass, Damping, Stiffness, Force, freq, convention = 'W'):
    """
    """
    if convention == 'N' :
        H = -freq**2*Mass-1j*freq*Damping+Stiffness
    else :
        H = -freq**2*Mass+1j*freq*Damping+Stiffness
    return solve(H,Force)

def PowerFreq(CPTO, Mass, Damping, Stiffness, Force, freq, convention = 'W') :
    """
    """
    if len(Force.shape) > 1 : # more than 1 column (more than 1 wave direction)
        Nd = Force.shape[-1]
    X = MotionFreq(Mass, CPTO + Damping, Stiffness, Force, freq, convention)
    if convention == 'N' :
        V = -1j*freq*X
    else :
        V = 1j*freq*X
    pdir = 0.5*np.dot(V.real.T,np.dot(CPTO,V.real)) + 0.5*np.dot(V.imag.T,np.dot(CPTO,V.imag))
    if len(Force.shape) > 1 :
        pdir = pdir[range(Nd),range(Nd)]
    return pdir

def ZeroDownCrossing(t, y, discard=0., intgr=False) :
    """
    Perfrom a zero-down-crossing analysis of the time signal y.

    t (1D numpy.array): time
    y (1D): time signal of length Number of records (i.e. len(t))
    discard (float): takes waves bigger than maximum wave height times discard
    """
    ##
    cross = np.zeros(y.shape, dtype = bool)
    cross[1:] = (y[:-1]*y[1:] < 0)*(y[:-1] > 0) + (y[1:] == 0.)*(y[:-1] > 0)
    # single wave heights and periods
    ind = np.array(range(cross.shape[0]), dtype = int)[cross]
    T = np.array([t[nf]-t[n0] for n0, nf in zip(ind[:-1],ind[1:])], dtype = float)
    if not intgr:
        H = np.array([y[n0:nf].max()-y[n0:nf].min() for n0, nf in zip(ind[:-1],ind[1:])], dtype = float)
    else:
        Intgr = np.array([(y[n0+1:nf]**2+y[n0:nf-1]**2).sum()*.5*(t[1]-t[0]) for n0, nf in zip(ind[:-1],ind[1:])], dtype = float)
        H = 2.*np.sqrt(2.*Intgr/T)
    ## discard waves which are smaller than the maximum wave height times discard
    cond = H >= H.max()*discard
    ind0, indf = ind[:-1][cond], ind[1:][cond]
    crossind = np.array([ind0,indf], dtype = int)
    auxt = [t[ind0[i]:indf[i]] for i in range(len(ind0))]
    auxy = [y[ind0[i]:indf[i]] for i in range(len(ind0))]
    tmod = np.array([itemt for subauxt in auxt for itemt in subauxt], dtype = float)
    ymod = np.array([itemy for subauxy in auxy for itemy in subauxy], dtype = float)
    return H[cond], T[cond], crossind, tmod, ymod

def Time2Freq(y, inct, tap=0.05, Nsub=1) :
    """
    np.fft.fft() does F[m] = sum_{n=0}^{Nt-1}(f[n]*np.exp(-1j*n*2*pi/Nt*m)) for m in range(len(t))

    I = int_{-T/2}^{T/2}(f*exp(-1j*2*pi/T*t*m)*dt)
    F[m] = 1/T*I
    np.fft.fft() outputs I/dt, so to get F[m] multiply by dt/T or,
    what is the same, multiply by 1./Nt, since T = Nt*dt

    Compute amplitudes for wave harmonics given the time signal y.

    y (1D): time signal of length Number of records

    frec (int) : frequency of recording, frec = Number of records / second

    tap (float) : if tap > 0, a taper filter will be used on
                 the time signal. tap corresponds to the fraction
                 of the total number of records which will be fully subjected to the
                 tapper filter. e.g number of records equal to 100 and
                 tap of 0.1 will take 100*tap = 10 records from the beggining
                 and from the ending of the signal where a taper filter
                 will be applied. Of course 0.0<tap<0.5. If Nsub > 1 the taper is the
                 fraction of the total number of records within the subserie
                 that will be fully subjected to the tapper filter!!!
    Nsub (int): Number of subseries, which will be as much overlapped as
                tap
    """
    #############
    
    #############
    
    ##
    Nrec = len(y) # Number of records
    Nrecsub = int(Nrec/(Nsub+tap*(Nsub+1)))
    Ntap = int(Nrecsub*tap)
    while Nrecsub*Nsub+Ntap*(Nsub+1) <= Nrec:
        Nrecsub += 1
        Ntap = int(Nrecsub*tap)
    Nrecsub -= 1
    Ntap = int(Nrecsub*tap)
    ypp = np.zeros((Nsub, Ntap+Nrecsub+Ntap))
    ypp[:, Ntap:] = y[Ntap: Nrecsub*Nsub+Ntap*(Nsub+1)].reshape((Nsub, -1))
    ypp[0, :Ntap] = y[:Ntap]
    ypp[1:, :Ntap] = ypp[:-1, Ntap+Nrecsub:]
    ## taper filter
    Np = Ntap+Nrecsub+Ntap
    filt = np.ones(Np, dtype = float)
    if tap > 0. :
        law = np.cos(np.linspace(0., pi*.5, Ntap))
        filt[:Ntap] = law[::-1]
        filt[-Ntap:] = law
    ##
    F = np.zeros((Nsub, Np/2), complex)
    for ind, ysub in enumerate(ypp):
        F[ind] = np.fft.fft(ysub.T*filt).T[:Np/2] # Nrec/2 due to aliasing
    F *= 1./Np
    ## associated frequencies in Hz
    fr = np.arange(Np/2)/float(Np)/inct
    return F, fr

def AverageFilter(y, Np, weights) :
    """
    y (1D): time signal of length Number of records
    Np (int) : 2*Np+1 number of points window for the average filter
    weights (1D numpy array): weights for half window the length should be Np+1
    """
    if Np > 0 :
        weights = np.array(list(weights[1:][::-1])+list(weights), float)
        aux = np.arange(Np, len(y)-Np)
        ind = np.array([aux+i for i in range(-Np,Np+1)], dtype = int).T
        y[Np:-Np] = (y[ind]*weights).sum(axis = 1)/weights.sum()
    return y

def convolute(h, f, order=0):
    """
    x = int_{-inf}^t(h(t-tau)*f(tau)dtau)
    x /= dtau
    """
    n = len( f )
    x = np.zeros( n, dtype=type( f[0] ) )
    for tau in range( n ):
        if order == 0:
            x[tau] = np.dot(h[:tau+1][::-1], f[:tau+1])
        elif order == 1:
            intgrnd = h[:tau+1][::-1]*f[:tau+1]
            x[tau] = .5*(intgrnd[1:]+intgrnd[:-1]).sum()
    return x


if __name__ == '__main__' :
    # L = 32
    # freq = np.linspace(.5, 1.5, 4)
    # period = 2*pi/freq
    # depth = 20.
    # evanescent = [0.03, 5., L]
    # Le = 2*pi/WNumber_E(period, depth, evanescent, step = .001, kmax = 1000)
    # ke = 2*pi/Le

    L = 1
    # period = 2*pi/0.3
    freq = np.linspace(0.3, 2, 20)
    period = 2*pi/freq
    depth = 20.
    k = WNumber(period, depth)
    Le = 2*pi/WNumber_E(period, depth, L, step = .001, kmax = 1000)
    k_l = 2*pi/Le
    # print(k_l)

    # Création du graphique
    # plt.figure(figsize=(8, 5))
    # plt.plot(freq, 2*pi/k, marker="o", label='Progressive', color='blue')
    # plt.plot(freq, Le, marker="o", label='Evanescent', color='orange', linestyle='--')

    # # Personnalisation
    # plt.title('')
    # plt.xlabel('Frequency')
    # plt.ylabel('Wave Length')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()

    # # Sauvegarde du graphique
    # plt.savefig("Dispersion_lambda.png", dpi=300)
    # plt.close()

    with open(f"Data_progressive_waves.dat", "w") as f:
            f.write(f"frequence - nombre d'onde - longueur d'onde  \n  ")
            for i, w in enumerate(freq) : 
                f.write(f"{w:8.4f}     {k[i]:8.4f}    {2*pi/k[i]:8.4f}  \n  ")