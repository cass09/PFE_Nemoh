"""Computing the potential and velocity of cylindrical wave."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import numpy as np
from capytaine.tools.lists_of_points import _normalize_points, _normalize_free_surface_points
from scipy.special import jv, iv
def cylindrical_waves_potential(points, pb):
    """Compute the potential for cylindrical waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest
        !!!!!!!!!! can be progressive -> m mode, k_0 
                        or evanescent --> m mode, l mode, k_l
    Returns
    -------
    array of shape (1) or (N x 1)
        The potential
    """
    points, output_shape = _normalize_points(points)

    x, y, z = points.T
    k = pb.wavenumber
    h = pb.water_depth
    r=np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        # sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        # sih = np.exp(k*z)
    if pb.modeL==0 :
        phi = cih *jv(pb.modeM, k*r)  *np.exp(1j * pb.modeM * theta)
    else :
        kl=WNumber_E(pb.period, h, pb.modeL, One=True) 
        phi = np.cos(kl*(z+h))* iv(pb.modeM, kl*r) * np.exp(1j * pb.modeM * theta)
    return phi.reshape(output_shape)


def cylindrical_waves_velocity(points, pb):
    """Compute the fluid velocity for cylindrical waves at a given point (or array of points).

    Parameters
    ----------
    points: array of shape (3) or (N x 3)
        coordinates of the points in which to evaluate the potential.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    array of shape (3) or (N x 3)
        the velocity vectors
    """

    points, output_shape = _normalize_points(points)

    x, y, z = points.T
    k = pb.wavenumber
    h = pb.water_depth
    r=np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
   
    if 0 <= k*h < 20:
        cih = np.cosh(k*(z+h))/np.cosh(k*h)
        sih = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        cih = np.exp(k*z)
        sih = np.exp(k*z)
    # print(pb.modeM, pb.modeL)

    if pb.modeL==0 :
        v = np.exp(1j * pb.modeM * theta) * \
        np.array([(k*x/r*jv(pb.modeM-1, k*r)-pb.modeM/r**2*(x+1j*y)*jv(pb.modeM, k*r)) * cih, 
                  (k*y/r*jv(pb.modeM-1, k*r)-pb.modeM/r**2*(y-1j*x)*jv(pb.modeM, k*r)) * cih, 
                  k*jv(pb.modeM, k*r)* sih])
    else :
        kl=WNumber_E(pb.period, h, pb.modeL, One=True) 
        v = np.exp(1j * pb.modeM * theta) * \
        np.array([(kl*x/r*iv(pb.modeM-1, kl*r)-pb.modeM/r**2*(x+1j*y)*iv(pb.modeM, kl*r)) * np.cos(kl*(z+h)), 
                  (kl*y/r*iv(pb.modeM-1, kl*r)-pb.modeM/r**2*(y-1j*x)*iv(pb.modeM, kl*r)) * np.cos(kl*(z+h)), 
                  -kl*iv(pb.modeM, kl*r)* np.sin(kl*(z+h))])

    return v.T.reshape((*output_shape, 3))


def cylindrical_waves_pressure(points, pb):
    return 1j * pb.omega * pb.rho * cylindrical_waves_potential(points, pb)


def froude_krylov_force_cyl(pb):
    return pb.body.integrate_pressure(cylindrical_waves_pressure(pb.body.mesh.faces_centers, pb))


def cylindrical_waves_free_surface_elevation(points, pb):
    """Compute the free surface elevation at points of the undisturbed cylindrical waves

    Parameters
    ----------
    points: array of shape (3) or (N × 3) or (2) or (N × 2)
        coordinates of the points in which to evaluate the potential.
        If only two coordinates are passed, the last one is filled with zeros.
    pb: DiffractionProblem
        problem with the environmental conditions (g, rho, ...) of interest

    Returns
    -------
    complex-valued array of shape (1,) or (N,)
        the free surface elevations
    """
    points, output_shape = _normalize_free_surface_points(points)
    return 1j * pb.omega / pb.g * cylindrical_waves_potential(points, pb).reshape(output_shape)


def Dispersion_E(L, const):
    # Linear dispersion equation, f=0, for evanescent modes
    T,h= const
    g=9.81
    A= g*T**2/2/np.pi
    B= 2*np.pi*h
    f= L+A*np.tan(B/L)
    df= 1-A*B/np.cos(B/L)**2/L**2
    return f,df

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

def WNumber_E(period, depth, L, One=False):
    """
    Calculates wave-numbers, ke, from wave-period and water depth (T,h).
    Dispersion relationship for evanescent modes
    Inputs:
    - period: a float. Wave periods
    - depth: float(). It is the water depth.
    - L (integer):  evanescent mode (L).
    Outputs:
    """
    if One :
        p = period
        l=L
        k0 = np.pi/2./depth*(2*l+1.)*(1.+1e-6)
        kf = np.pi/2./depth*(2*l+2.)
        a = 2*np.pi/k0
        b = 2*np.pi/kf
        sol = Bisection(Dispersion_E, (p, depth), a, b)
        if sol[1][0] == 1 : # Solution found in (a,b)
            Le = sol[0]
        else :
            print('Fail evanescent')
                
    else :
        Le = np.zeros((len(period), L), dtype = float)
        for ip in range(len(period)) :
            # if len(period) > 1 :
            p = period[ip]
            # else :
            #     p = period
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
    return 2*np.pi/Le