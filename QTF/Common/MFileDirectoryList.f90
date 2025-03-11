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
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!    List of input/output file
!
!--------------------------------------------------------------------------------------
MODULE MFileDirectoryList

CHARACTER(LEN=*), PARAMETER      ::PreprocDir='QTFPreprocOut'
CHARACTER(LEN=*), PARAMETER      ::LogFILE   ='logfileQTF.txt'
CHARACTER(LEN=*), PARAMETER      ::TotPotFILE='TotalPotentialBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::TotVelFILE='TotalVelocityBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::RadPotFILE='RadiationPotentialBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::RadVelFILE='RadiationVelocityBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::IncPotFILE_FS='IncidentPotentialFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::IncVelFILE_FS='IncidentVelocityFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::TotPotFILE_FS='TotalPotentialFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::TotVelFILE_FS='TotalVelocityFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::RadPotFILE_FS='RadiationPotentialFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::RadVelFILE_FS='RadiationVelocityFreeSurface.bin'
CHARACTER(LEN=*), PARAMETER      ::IncPotFILE_FSdbg='IncidentPotentialFreeSurfaceDbg.dat'
CHARACTER(LEN=*), PARAMETER      ::IncVelFILE_FSdbg='IncidentVelocityFreeSurfaceDbg.dat'
CHARACTER(LEN=*), PARAMETER      ::TotPotFILE_FSdbg='TotalPotentialFreeSurfaceDbg.dat'
CHARACTER(LEN=*), PARAMETER      ::TotVelFILE_FSdbg='TotalVelocityFreeSurfaceDbg.dat'
CHARACTER(LEN=*), PARAMETER      ::RadPotFILE_FSdbg='RadiationPotentialFreeSurfaceDbg.dat'
CHARACTER(LEN=*), PARAMETER      ::RadVelFILE_FSdbg='RadiationVelocityFreeSurfaceDbg.dat'

CHARACTER(LEN=*), PARAMETER      ::  OutQTFDir ='/results/QTF/'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDM ='QTFM_DUOK.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDP ='QTFP_DUOK.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBM='QTFM_HASBO.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBP='QTFP_HASBO.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHFSM='QTFM_HASFS.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHFSP='QTFP_HASFS.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileASYM='QTFM_ASYMP.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileASYP='QTFP_ASYMP.dat'

CHARACTER(LEN=*), PARAMETER      ::  OutFileDM_term ='QTFM_DUOK_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDP_term ='QTFP_DUOK_term_'

CHARACTER(LEN=*), PARAMETER      ::  OutFileHBM_term='QTFM_HASBO_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBP_term='QTFP_HASBO_term_'

CHARACTER(LEN=*), PARAMETER      ::  OutFileHFSM_term='QTFM_HASFS_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHFSP_term='QTFP_HASFS_term_'

CHARACTER(LEN=*), PARAMETER      ::  OutFileASYM_term='QTFM_ASYMP_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileASYP_term='QTFP_ASYMP_term_'



!Output variables
  !DUOK      : Quadratic QTF
  !HASBO     : Potential QTF-Haskind on Body term
  !HASFS     : Potential QTF-Haskind on Free-Surface in Finite domain
  !HASFS_ASYM: Potential QTF-Haskind on Free-Surface in Asymptotic domain
CONTAINS
   SUBROUTINE  make_directory(dirname)
          CHARACTER(LEN=*),       INTENT(IN) ::dirname
          LOGICAL                            ::existdir
          !INQUIRE (DIRECTORY=dirname, EXIST=existdir) !this is Intel-specific
          INQUIRE (FILE=TRIM(dirname)//'/.', EXIST=existdir)
          IF (.NOT.existdir) CALL SYSTEM('mkdir '//dirname)
   END SUBROUTINE


END MODULE
