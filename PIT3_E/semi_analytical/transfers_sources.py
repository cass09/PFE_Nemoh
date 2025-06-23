#!/usr/bin/python2.7
# encoding: utf-8
"""
This module contains the methods to derive the diffraction and force
transfer matrices, which are required input of the direct matrix
method (Kagemoto and Yue, 1986). These are computed from the method
presented by McNatt using standard BEM outputs.

.. module:: transfers
   :platform: Windows
   :synopsis: Numerical model of WEC builder

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
"""

from math import log
import numpy as np
from scipy.special import jv, yv, iv, kv
from toolbox.CalaixSastre import WNumber, WNumber_E

def transfers_sources(water_depth,
              directions,
              periods,
              BodyMesh,
              source_scat,
              source_rad,
              fex,
              tol,
              convention,
              Evanescent,
              Nmodes_E):
    """ Computes cylindrical amplitude coefficients from the sources terms
      values on the isolated body.

    :param water_depth: water depth (m)
    :type water_depth: float
    :param directions: wave heading angles (rad)
    :type directions: list or 1D numpy array
    :param periods: wave periods (s)
    :type periods: list or 1D numpy array
    :param BodyMesh: (centers(r,theta,z), areas)
    :type discrete_cyl: tuple
    :param source_scat: complex amplitude of the term source for each discrete panel
                  of the body and for the inputted water depth, wave periods and wave,
                  directions for the diffraction problems. shape (number wave frequencies,
                  number wave directions, number of panels)
    :type source_scat: 3D numpy array
    :param source_rad: complex amplitude of the term source for each discrete panel
                  of the body and for the inputted water depth, wave periods and degrees
                  of freedom, for the diffraction problems. shape (number wave frequencies,
                  number of degrees of freedom, number of panels)
    :type source_rad: 3D numpy array
    :param fex: complex amplitude of the diffraction excitation force () for the inputted
                water depth, wave periods, wave directions and degrees of freedom.
                shape (number wave frequencies, number of wave directions, number of
                 degrees of freedom)
    :type fex: 3D numpy array
    :param tol: 1e-(number of significant decimals) to be used in max_trunc_order()
    :type tol: float
    """
    targ_order = int((len(directions)-1)/2)
    act_order = np.zeros((len(periods), 2), dtype=int)
    decimals = np.zeros((len(periods), 2), dtype=int)
    
    diffmat = np.zeros((len(periods), 2*targ_order+1, 2*targ_order+1), dtype=complex)
    frcmat = np.zeros((len(periods), 2*targ_order+1, fex.shape[-1]), dtype=complex)
    a_s_rad = np.zeros((len(periods), fex.shape[-1], 2*targ_order+1), dtype=complex)
    a_s_scat = np.zeros((len(periods), source_scat.shape[1], 2*targ_order+1), dtype=complex)
    dirs, modes = np.meshgrid(directions, range(-targ_order, targ_order+1),
                              indexing='ij', sparse=True)
    print("tranfers")
    for ind, per in enumerate(periods):
        Kw=(2.*np.pi/per)**2/9.81
        k0= WNumber(per, water_depth)
        C0=(Kw**2-k0**2)/((k0**2-Kw**2)*water_depth+Kw)
        wave_cond = (water_depth, 2.*np.pi/per, k0)
        int_scat = integral_sources(wave_cond, BodyMesh, source_scat[ind], targ_order)
        int_rad = integral_sources(wave_cond, BodyMesh, source_rad[ind], targ_order)
        a_s_rad[ind]=1j/2*C0*np.cosh(k0*water_depth)*int_rad
        # print(a_s_rad)
        diffmat[ind] = 1j/2*C0*np.cosh(k0*water_depth)*int_scat
        # frcmat[ind] = np.linalg.lstsq(psi_i, fex[ind], rcond=None)[0]
        frcmat[ind]=fex[ind]
        # act_order[ind, 0], decimals[ind, 0] = max_trunc_order(a_s_scat[ind], targ_order, tol)
        act_order[ind, 0] = targ_order
        act_order[ind, 1], decimals[ind, 1] = max_trunc_order(a_s_rad[ind], targ_order, tol)
    
    # print(diffmat.shape, a_i_plane_E.shape, coef_scat.shape)
    # print(frcmat.shape, coef_rad.shape)
    # print("G", frcmat[0], "aR", coef_rad[0])
    # with open("D.dat", 'w') as f:
    #     for row in diffmat[0]:
    #         line = "\t".join(f"{val.real:.6e}+{val.imag:.6e}j" for val in row)
    #         f.write(line + "\n")
    # with open("G.dat", 'w') as f:
    #     for row in frcmat[0]:
    #         line = "\t".join(f"{val.real:.6e}+{val.imag:.6e}j" for val in row)
    #         f.write(line + "\n")
    # with open("aR.dat", 'w') as f:
    #     for row in coef_rad[0]:
    #         line = "\t".join(f"{val.real:.6e}+{val.imag:.6e}j" for val in row)
    #         f.write(line + "\n")


    # act_order[:]=targ_order
    # Shrink G, D and AR according to the truncation order Nm
    ini = targ_order-act_order.max()
    fin = ini+2*act_order.max()+1
    # print(diffmat.shape, ini, fin)
    # print(act_order.shape, act_order)
   
    return (diffmat[:, ini:fin, ini:fin],
                frcmat[:, ini:fin, :],
                a_s_rad[:, :, ini:fin],
                a_s_scat[:, :, ini:fin],
                act_order.max(axis=0),
                act_order)
   
    # act_order[:]=targ_order
    # return (diffmat.round(decimals.max()),
    #         frcmat.round(decimals.max()),
    #         a_s_rad.round(decimals.max()),
    #         act_order.max(axis=0),
    #         act_order)

def max_trunc_order(a_prob,
                    targ_order,
                    tol):
    """ Selects the biggest wave mode fulfiling a_prob_mode > abs(a_prob).max()*tol

    :param a_prob: cylindrical amplitude coefficients. First dimension is taken for
                   wave directions, if diffraction problem, or for degrees of freedom,
                   if radiation problem.
    :type a_prob: 2D numpy array
    :param targ_order: maximum possible truncation order
    :type targ_order: int
    :param tol: 1e-(number of significant decimals)
    :type tol: float
    """
    tol = max([tol*abs(a_prob).max(), 1e-99])
    decimals = int(abs(log(tol, 10)))
    act_order = 0
    for n_mode in (range(targ_order)):
        if any(abs(a_prob[:, n_mode]) > tol):
            act_order = targ_order-n_mode
            break
    return (act_order, decimals)

def bem2cyl(wave_cond,
            discrete_cyl,
            vpot_cyl,
            trunc_ord,
            convention):
    """ Computes cylindrical amplitude coefficients from the velocity
    potential values on a cylinder.

    :param wave_cond: (water_depth, cfreq, wnum)
    :type wave_cond: tuple
    :param water_depth: water depth (m)
    :type water_depth: float
    :param cfreq: cyclic wave frequency (rad/s)
    :type cfreq: float
    :param wnum: wave number (rad/m) associated to that wave frequency and water depth
    :type wnum: float
    :param discrete_cyl: (radius_cyl, azimuth_cyl, axial_cyl)
    :type discrete_cyl: tuple
    :param radius_cyl: radius (m) of the cylinder
    :type radius_cyl: float
    :param azimuth_cyl: azimuthal discretization (rad) of the cylinder. An equispaced
                        discretization is assumed. Each element of the array is
                        an azimuth. Same axial discretization is assumed for each
                        azimuth
    :type azimuth_cyl: 1D numpy array
    :param axial_cyl: axial discretization (m) of the cylinder. Each element of the array is
                      the z-axial coordinate of a discrete point of the cylinder. Same
                      azimuthal discretization is assumed for each z-coordinate
    :type axial_cyl: 1D numpy array
    :param vpot_cyl: complex amplitude of the velocity potential (m**2/s) for each discrete point
                     of the cylinder and for the inputted water depth and wave frequency. shape
                     (extra axis, number of axial discretization, number of azimuthal
                     discretization)
    :type vpot_cyl: 3D numpy array
    :param trunc_ord: truncation order for the number of wave modes included.
                      Total number of wave modes is 2*trunc_ord+1
    :type trunc_ord: int
    """
    (water_depth, cfreq, wnum) = wave_cond
    (radius_cyl, azimuth_cyl, axial_cyl) = discrete_cyl
    dz = axial_cyl[1:]-axial_cyl[:-1]
    dth = azimuth_cyl[1]-azimuth_cyl[0] # equispaced is assumed
    rightz = all(dz > 0)
    rightth = dth > 0
    rightth2 = azimuth_cyl[-1] == 2.*np.pi + azimuth_cyl[0] or azimuth_cyl[-1] == azimuth_cyl[0]
    # Build integration domain
    (z_cyl, th_cyl) = np.meshgrid(axial_cyl, azimuth_cyl, indexing='ij')
    # Initialize
    a_s = np.zeros((vpot_cyl.shape[0], 2*trunc_ord+1), dtype=complex)
    for n_mode, mode in enumerate(range(-trunc_ord, trunc_ord+1)):
        integrand = vpot_cyl*np.cosh(wnum*(z_cyl+water_depth))*np.exp(-1j*mode*th_cyl)
        # Integrate along th
        int_th = (integrand[:, :, 1:]+integrand[:, :, :-1]).sum(axis=2)*.5*dth
        if not rightth2: # add last paralepipede
            int_th += (integrand[:, :, 0]+integrand[:, :, -1])*.5*dth
        if not rightth:
            int_th *= -1
        # Integrate I_th along z
        int_th_z = ((int_th[:, 1:]+int_th[:, :-1])*dz).sum(axis=1)*.5
        if not rightz:
            int_th_z *= -1
        # Cm
        cntm = -1j*cfreq/(2*np.pi*9.809)
        if convention == 'N':
            cntm *= -1
        cntm *= 2*np.cosh(wnum*water_depth)
        cntm /= water_depth*(1+np.sinh(2*wnum*water_depth)/(2*wnum*water_depth))
        if convention == 'N':
            cntm /= jv(mode, wnum*radius_cyl)+1j*yv(mode, wnum*radius_cyl)
        else:
            cntm /= jv(mode, wnum*radius_cyl)-1j*yv(mode, wnum*radius_cyl)
        # amplitude coefficients
        a_s[:, n_mode] = cntm*int_th_z
    return a_s


def integral_sources(wave_cond, BodyMesh, source, trunc_ord):
    # sigma_all, centers, areas, k, k0, d, M, p_values=None):
    """
    Calcule une matrice de projection :
        A_{mp} = ∫_S σ_p(r,θ,z) * J_m(k0 r) * cosh(k0(z + d)) * exp(-imθ) dS

    pour tous les m ∈ [-M, M] et p donné (mode q ou dof k).

    :param wave_cond: (water_depth, cfreq, wnum)
    :type wave_cond: tuple
    :param water_depth: water depth (m)
    :type water_depth: float
    :param cfreq: cyclic wave frequency (rad/s)
    :type cfreq: float
    :param wnum: wave number (rad/m) associated to that wave frequency and water depth
    :type wnum: float
    :param BodyMesh: (centers(r,theta,z), areas)
    :param source: complex amplitude of the term source for each discrete panel
                  of the body and for the inputted water depth, wave periods and wave,
                  directions for the diffraction or radiation problems. shape
                   (number wave directions or dof, number of panels)
    :type source: 2D numpy array
    :param trunc_ord: truncation order for the number of wave modes included.
                      Total number of wave modes is 2*trunc_ord+1
    :type trunc_ord: int
    Retour
    ------
    A : np.ndarray
        Matrice complexe de forme (2M+1, len(sigma_all))
    """
    centers=BodyMesh[0]
    areas=BodyMesh[1]
    r = centers[:, 0]
    theta = centers[:, 1]
    z = centers[:, 2]
    (water_depth, cfreq, k0) = wave_cond

    Np = len(source[:,0])
    A = np.zeros((Np, 2 * trunc_ord + 1), dtype=complex)
    coshz = np.cosh(k0 * (z + water_depth))

    # Pour chaque m, calcule les termes communs
    for mi, mode in enumerate(range(-trunc_ord, trunc_ord + 1)):
        Jm_kr = jv(mode, k0 * r)
        exp_m_theta = np.exp(-1j * mode * theta)

        for pi in range(Np):
            sigma_p = source[pi]
            integrand = sigma_p * Jm_kr * coshz * exp_m_theta
            A[pi, mi] = np.sum(integrand * areas)

    return A