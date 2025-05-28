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

def transfers(water_depth,
              directions,
              periods,
              discrete_cyl,
              vpot_scat,
              vpot_rad,
              fex,
              tol,
              convention,
              Evanescent,
              Nmodes_E):
    """ Computes cylindrical amplitude coefficients from the velocity
    potential values on a cylinder.

    :param water_depth: water depth (m)
    :type water_depth: float
    :param directions: wave heading angles (rad)
    :type directions: list or 1D numpy array
    :param periods: wave periods (s)
    :type periods: list or 1D numpy array
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
    :param vpot_scat: complex amplitude of the velocity potential (m**2/s) for each discrete point
                  of the cylinder and for the inputted water depth, wave periods and wave,
                  directions for the diffraction problems. shape (number wave frequencies,
                  number wave directions, number of axial discretization, number of azimuthal
                  discretization)
    :type vpot_scat: 4D numpy array
    :param vpot_rad: complex amplitude of the velocity potential (m**2/s) for each discrete point
                  of the cylinder and for the inputted water depth, wave periods and degrees
                  of freedom, for the diffraction problems. shape (number wave frequencies,
                  number of degrees of freedom, number of axial discretization, number of
                  azimuthal discretization)
    :type vpot_rad: 4D numpy array
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
    diffmat_E = np.zeros((len(periods), (2*targ_order+1)*Nmodes_E, (2*targ_order+1)*Nmodes_E), dtype=complex)
    frcmat = np.zeros((len(periods), 2*targ_order+1, fex.shape[-1]), dtype=complex)
    frcmat_E = np.zeros((len(periods), (2*targ_order+1)*Nmodes_E, fex.shape[-1]), dtype=complex)
    a_s_rad = np.zeros((len(periods), fex.shape[-1], 2*targ_order+1), dtype=complex)
    a_s_scat = np.zeros((len(periods), vpot_scat.shape[1], 2*targ_order+1), dtype=complex)
    b_s_rad = np.zeros((len(periods), fex.shape[-1], (2*targ_order+1)*Nmodes_E), dtype=complex)
    b_s_scat = np.zeros((len(periods), vpot_scat.shape[1], (2*targ_order+1)*Nmodes_E), dtype=complex)
    dirs, modes = np.meshgrid(directions, range(-targ_order, targ_order+1),
                              indexing='ij', sparse=True)
    
    wave_number_e=WNumber_E(periods, water_depth, Nmodes_E)
    for ind, per in enumerate(periods):
        wave_cond = (water_depth, 2.*np.pi/per, WNumber(per, water_depth))
        a_s_scat[ind] = bem2cyl(wave_cond, discrete_cyl, vpot_scat[ind], targ_order, convention)
        a_s_rad[ind] = bem2cyl(wave_cond, discrete_cyl, vpot_rad[ind], targ_order, convention)
        
        # for evanescent waves
        if Evanescent : 
            wave_cond_E = (water_depth, 2.*np.pi/per, wave_number_e[ind])
            b_s_scat[ind] = bem2cyl_ev(wave_cond_E, discrete_cyl, vpot_scat[ind], targ_order, convention)
            # b_s_scat = np.zeros((11, (2*targ_order+1)*Nmodes_E), dtype=complex)
            b_s_rad[ind] = bem2cyl_ev(wave_cond_E, discrete_cyl, vpot_rad[ind], targ_order, convention)
        
        if convention == 'N':
            a_i_plane = np.exp(1j*modes*(np.pi/2.-dirs))
            if Evanescent : 
                a_i_plane_E = np.zeros((a_i_plane.shape[0], a_i_plane.shape[1] * Nmodes_E), dtype=complex)
                for l in range(Nmodes_E):
                    start = l * a_i_plane.shape[1]
                    end = start + a_i_plane.shape[1]
                    a_i_plane_E[:, start:end] = a_i_plane
        else:
            a_i_plane = np.exp(-1j*modes*(np.pi/2.+dirs))
        diffmat[ind] = np.linalg.lstsq(a_i_plane, a_s_scat[ind], rcond=None)[0]
        frcmat[ind] = np.linalg.lstsq(a_i_plane, fex[ind], rcond=None)[0]
        if Evanescent : 
            diffmat_E[ind] = np.linalg.lstsq(a_i_plane_E, b_s_scat[ind], rcond=None)[0]
            frcmat_E[ind] = np.linalg.lstsq(a_i_plane_E, fex[ind], rcond=None)[0]
        # find maximum truncation order
        act_order[ind, 0], decimals[ind, 0] = max_trunc_order(a_s_scat[ind], targ_order, tol)
        act_order[ind, 1], decimals[ind, 1] = max_trunc_order(a_s_rad[ind], targ_order, tol)

    # Shrink G, D and AR according to the truncation order Nm
    ini = targ_order-act_order.max()
    fin = ini+2*act_order.max()+1

    # print(diffmat.shape, diffmat[0])
    if Evanescent : 
        # Matrice réduite de D_e
        diffmat_E_reduced = np.zeros((len(periods), (2*act_order.max()+1)*Nmodes_E, (2*act_order.max()+1)*Nmodes_E), dtype=complex)
        frcmat_E_reduced = np.zeros((len(periods), (2*act_order.max()+1)*Nmodes_E, fex.shape[-1]), dtype=complex)
        # On coupe bloc par bloc
        for t in range(len(periods)):
            for i in range(Nmodes_E):
                block_f= frcmat_E[t, i*(2*targ_order+1):(i+1)*(2*targ_order+ 1), :]
                frcmat_E_reduced[t, i*(2*act_order.max()+1):(i+1)*(2*act_order.max()+1), :] = block_f[ini:fin, :]
                for j in range(Nmodes_E):
                    block = diffmat_E[t, i*(2*targ_order+1):(i+1)*(2*targ_order+ 1), j*(2*targ_order+ 1):(j+1)*(2*targ_order+1)]
                    diffmat_E_reduced[t, i*(2*act_order.max()+1):(i+1)*(2*act_order.max()+1), j*(2*act_order.max()+1):(j+1)*(2*act_order.max()+1)] = block[ini:fin, ini:fin]

        # Matrice D global = progressive + evanescent 
        Nm_red=(2*act_order.max()+1)
        D_global = np.zeros((len(periods), Nm_red*(1+Nmodes_E), Nm_red*(1+Nmodes_E)), dtype=complex)
        G_global = np.zeros((len(periods), Nm_red*(1+Nmodes_E), fex.shape[-1]), dtype=complex)
        coef_rad = np.zeros((len(periods), fex.shape[-1], (2*targ_order+1)*(1+Nmodes_E)), dtype=complex)
        coef_scat = np.zeros((len(periods), vpot_scat.shape[1], (2*targ_order+1)*(1+Nmodes_E)), dtype=complex)
        for t in range(len(periods)):
            # Partie centrale : D
            D_global[t, :Nm_red, :Nm_red] = diffmat[t]
            G_global[t, :Nm_red, :] = frcmat[t]
            coef_rad[t, :, :2*targ_order+1]=a_s_rad[t, :, :]
            coef_scat[t, :, :2*targ_order+1]=a_s_scat[t, :, :]
            # Partie évanescente : D_e
            start = Nm_red
            end = Nm_red * (1 + Nmodes_E)
            D_global[t, start:end, start:end] = diffmat_E_reduced[t]
            G_global[t, start:end, :] = frcmat_E_reduced[t]
            coef_rad[t, :, 2*targ_order+1:(2*targ_order+1)* (1 + Nmodes_E)]=b_s_rad[t, :, :]
            coef_scat[t, :, 2*targ_order+1:(2*targ_order+1)* (1 + Nmodes_E)]=b_s_scat[t, :, :]
        
        # with open("D.dat", 'w') as f:
        #     for row in D_global[0]:
        #         line = "\t".join(f"{val.real:.6e}+{val.imag:.6e}j" for val in row)
        #         f.write(line + "\n")
        return (D_global.round(decimals.max()),
                G_global.round(decimals.max()),
                coef_rad.round(decimals.max()),
                coef_scat.round(decimals.max()),
                act_order.max(axis=0),
                act_order)
    else : 
        return (diffmat[:, ini:fin, ini:fin].round(decimals.max()),
                frcmat[:, ini:fin, :].round(decimals.max()),
                a_s_rad[:, :, ini:fin].round(decimals.max()),
                a_s_scat[:, :, ini:fin].round(decimals.max()),
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

    
def bem2cyl_ev(wave_cond,
            discrete_cyl,
            vpot_cyl,
            trunc_ord,
            convention):
    """ Computes cylindrical amplitude coefficients from the velocity
    potential values on a cylinder.

    :param wave_cond: (water_depth, cfreq, wnum for evanescent waves)
    :type wave_cond: tuple
    :param water_depth: water depth (m)
    :type water_depth: float
    :param cfreq: cyclic wave frequency (rad/s)
    :type cfreq: float
    :param wnum: wave number (rad/m) associated to that wave frequency and water depth
                for evanescent waves
    :type wnum: 1D numpy array for Nmodes_Evanescent
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
    Ne=len(wnum)
    (radius_cyl, azimuth_cyl, axial_cyl) = discrete_cyl
    dz = axial_cyl[1:]-axial_cyl[:-1]
    dth = azimuth_cyl[1]-azimuth_cyl[0] # equispaced is assumed
    rightz = all(dz > 0)
    rightth = dth > 0
    rightth2 = azimuth_cyl[-1] == 2.*np.pi + azimuth_cyl[0] or azimuth_cyl[-1] == azimuth_cyl[0]
    # Build integration domain
    (z_cyl, th_cyl) = np.meshgrid(axial_cyl, azimuth_cyl, indexing='ij')
    # Initialize
    a_s = np.zeros((vpot_cyl.shape[0], (2*trunc_ord+1)* Ne), dtype=complex)
    for l_mode in range(0,  Ne):    
        for n_mode, mode in enumerate(range(-trunc_ord, trunc_ord+1)):
            integrand = vpot_cyl*np.cos(wnum[l_mode]*(z_cyl+water_depth))*np.exp(-1j*mode*th_cyl)
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
            cntm *= 2
            cntm /= water_depth*(1+np.sin(2*wnum[l_mode]*water_depth)/(2*wnum[l_mode]*water_depth))
            cntm /= kv(mode, wnum[l_mode]*radius_cyl)
            # amplitude coefficients
            a_s[:, l_mode*(2*trunc_ord+1)+n_mode] = cntm*int_th_z
    return a_s

    
