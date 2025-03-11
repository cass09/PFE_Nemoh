!--------------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!   Contributors list:
!   - G. Delhommeau
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
  MODULE COM_VAR

    IMPLICIT NONE

    ! --- Environment --------------------------------------------------
    ! Volumique mass of the fluid in KG/M**3
    REAL :: RHO
    ! Gravity
    REAL :: G
    ! depth of the domain
    REAL :: Depth
    ! Coordinates where the waves are measured
    REAL :: XEFF,YEFF,ZEFF

    ! --- Geometrical quantities ---------------------------------------
    ! Mesh file
    CHARACTER*80 :: MESHFILE
    INTEGER :: LFILE
    ! No. of points in the surface mesh
    INTEGER :: NP
    ! No of total panels in the surface mesh
    INTEGER :: NFA
    !
    INTEGER :: IMX
    INTEGER :: IXX
    !!!!!!!!!
    !Symmetry: 0 for no symmetry and 1 for symmetry
    INTEGER:: NSYMY
    !
    REAL :: ZER
    ! DIST(I) is the maximum distance (projection on XY plane!) of a panel I from other points
    REAL, DIMENSION(:), ALLOCATABLE :: DIST,TDIS
    !vertices of the panel, size NFA(no of facettes)
    INTEGER, DIMENSION(:), ALLOCATABLE :: M1,M2,M3,M4
    ! Array of cordinates of points of size NP (no of points)
    REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z
    ! Normal vector
    REAL, DIMENSION(:), ALLOCATABLE :: XN,YN,ZN
    ! Array of centre of gravity for each panel
    REAL, DIMENSION(:), ALLOCATABLE :: XG,YG,ZG
     ! Array of Gauss points for each panel
    REAL, DIMENSION(:,:), ALLOCATABLE :: XGA,YGA,ZGA,XJAC
    !Number of Gauss points
    INTEGER:: NG
    ! Array for surface area of the panels
    REAL, DIMENSION(:), ALLOCATABLE :: AIRE

    ! --- Boundary value problem ---------------------------------------
    ! normal velocity array as input
!    REAL, DIMENSION(:), ALLOCATABLE :: NVEL
    ! period
    REAL :: T
    ! Computed and computed potential on original (B) and symmetric boundary(S)
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ZPB,ZPS
    ! Source distribution
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB,ZIGS

    ! --- Solver ------------------------------------------------------
    ! Which solver: (0) direct solver GAUSS ELIMINATION, (1) LU DECOMPOSITION (2): GMRES
    INTEGER:: Indiq_solver,mRestartGMRES,NITERGMRES
    REAL   :: TOLGMRES
    INTEGER,PARAMETER :: ID_DP=0  ! 1 if compiled with -r8 for the double precision, 0 otherwise
    ! Linear complex matrix to be solved
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ZIJ
    ! Storage of inverter matrix to speed up computations (influence coefficients depending only on w)
    COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: AInv,AmatIs
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: Amat
    REAL w_previous
    !
    REAL :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
    ! Variable for storage of Greens function
    REAL :: SP1,SM1,SP2,SM2
    ! Variable for storage of Greens function
    REAL :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2
    REAL :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
    !Values used in interpolation of infinite part of the Greens function
    INTEGER,PARAMETER :: NPINTE=5001
    INTEGER:: NQ
    REAL:: CQ(NPINTE),QQ(NPINTE),AMBDA(31),AR(31)

    ! --- Reading and writing units -----------------------------------
    ! File for visualization of potential
!     INTEGER :: Sav_potential		!AC: unused
    INTEGER :: Switch_Sources   !AC: Always save sources for now (for QTF), updated in INITIALIZATION.f90

  END MODULE COM_VAR
