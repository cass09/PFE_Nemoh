"""Experimental function to compute the Response Amplitude Operator."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import xarray as xr
from capytaine.post_pro.impedance import rao_transfer_function

LOG = logging.getLogger(__name__)


def rao(dataset, wave_direction=None, dissipation=None, stiffness=None):
    """Response Amplitude Operator.

    Parameters
    ----------
    dataset: xarray Dataset
        The hydrodynamical dataset.
        This function supposes that variables named 'inertia_matrix' and 'hydrostatic_stiffness' are in the dataset.
        Other variables can be computed by Capytaine, by those two should be manually added to the dataset.
    wave_direction: float, optional
        Select a wave directions for the computation. (Not recommended, kept for legacy.)
        Default: all wave directions in the dataset.
    dissipation: array, optional
        An optional dissipation matrix (e.g. Power Take Off) to be included in the RAO.
        Default: none.
    stiffness: array, optional
        An optional stiffness matrix (e.g. mooring stiffness) to be included in the RAO.
        Default: none.

    Returns
    -------
    xarray DataArray
        The RAO as an array depending of omega and the degree of freedom.
    """

    # ASSEMBLE MATRICES
    H = rao_transfer_function(dataset, dissipation, stiffness)
    fex = dataset.excitation_force

    LOG.info("Compute RAO.")

    # SOLVE LINEAR SYSTEMS
    # Match dimensions of the arrays to be sure to solve the right systems.
    H, fex = xr.broadcast(H, fex, exclude=["radiating_dof", "influenced_dof"])
    H = H.transpose(..., 'radiating_dof', 'influenced_dof')
    fex = fex.transpose(...,  'influenced_dof')

    if wave_direction is not None:  # Legacy behavior for backward compatibility
        H = H.sel(wave_direction=wave_direction)
        fex = fex.sel(wave_direction=wave_direction)

    # Solve and add coordinates
    rao_dims = [d for d in H.dims if d != 'influenced_dof']
    rao_coords = {c: H.coords[c] for c in H.coords if c != 'influenced_dof'}
    rao = xr.DataArray(np.linalg.solve(H.values, fex.values[..., np.newaxis])[..., 0], coords=rao_coords, dims=rao_dims)

    return rao

def compute_rao_from_results(results_data, omega_list, directions, inertia_matrix, hydrostatic_stiffness, 
                             dissipation=None, extra_stiffness=None, vect_dof=None):
    """
    Calcule le RAO à partir d'une liste results_data (sans xarray).

    Parameters
    ----------
    results_data : list
        Liste des dictionnaires résultats (diffraction/radiation) pour chaque omega et direction.
    omega_list : list or array
        Liste des fréquences.
    directions : list or array
        Liste des directions de houle.
    inertia_matrix : ndarray
        Matrice d’inertie (NxN).
    hydrostatic_stiffness : ndarray
        Matrice de raideur hydrostatique (NxN).
    dissipation : ndarray, optional
        Matrice d’amortissement additionnel (e.g. PTO).
    extra_stiffness : ndarray, optional
        Matrice de raideur additionnelle (e.g. amarrage).
    vect_dof : list of str
        Liste des DOFs (ordre cohérent avec les matrices).

    Returns
    -------
    RAO : dict
        RAO[omega][beta] = vector (complex)
    """

    N_dof = len(vect_dof)
    RAO = {}

    for w in omega_list:
        RAO[w] = {}

        # Initialiser matrices radiation
        B = np.zeros((N_dof, N_dof))
        A = np.zeros((N_dof, N_dof))
        # Chercher les résultats radiation pour cette fréquence
        for res in results_data:
            if res['type'] == 'radiation' and np.isclose(res['omega'], w):
                i = vect_dof.index(res['dof'])
                for j, dof_j in enumerate(vect_dof):
                    B[i, j] = res['damping'].get(dof_j, 0.0)
                    A[i, j] = res['added_mass'].get(dof_j, 0.0)

        M = np.array(inertia_matrix)
        K = np.array(hydrostatic_stiffness)
        if extra_stiffness is not None:
            K += extra_stiffness
        if dissipation is not None:
            B += dissipation

        H = -w**2 * (M + A) + 1j * w * B + K

        for beta in directions:
            # Force d'excitation
            F_exc = np.zeros(N_dof, dtype=complex)
            for res in results_data:
                if res['type'] == 'diffraction' and np.isclose(res['omega'], w) and np.isclose(res['beta'], beta):
                    for i, dof in enumerate(vect_dof):
                        F_exc[i] = res['excitation_force'].get(dof, 0.0)

            # Résolution du système
            try:
                RAO_wb = np.linalg.solve(H, F_exc)
            except np.linalg.LinAlgError:
                RAO_wb = np.full(N_dof, np.nan, dtype=complex)

            RAO[w][beta] = RAO_wb

    return RAO
