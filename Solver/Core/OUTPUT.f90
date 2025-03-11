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
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
MODULE OUTPUT

  USE Constants
  USE MMesh

  IMPLICIT NONE

  !PUBLIC :: WRITE_DATA_ON_MESH

  PUBLIC
  INTEGER, PARAMETER :: RAW_OUTPUT=0, TECPLOT_OUTPUT=1
  INTEGER :: output_format = TECPLOT_OUTPUT

CONTAINS
!
!     AC: add write sources for QTF
    SUBROUTINE WRITE_SOURCES(ZIGB,ZIGS,NPanels,filename)
    CHARACTER(LEN=*),           INTENT(IN) :: filename
    INTEGER,                    INTENT(IN) :: NPanels
    COMPLEX,DIMENSION(NPanels), INTENT(IN) :: ZIGB,ZIGS
    INTEGER     :: u,I

    OPEN(NEWUNIT=u, FILE=filename,ACTION='WRITE')
    do I=1,NPanels
      WRITE(u,'(2(X,E14.7))') REAL(ZIGB(I)), IMAG(ZIGB(I))
    end do

    do I=1,NPanels
      WRITE(u,'(2(X,E14.7))') REAL(ZIGS(I)), IMAG(ZIGS(I))
    end do
    CLOSE(u)
    END SUBROUTINE WRITE_SOURCES



  SUBROUTINE WRITE_DATA_ON_MESH(Mesh, cdata, filename)
    ! Save in a file a field of complex values on a mesh.

    TYPE(TMesh),                                   INTENT(IN) :: Mesh
    COMPLEX, DIMENSION(Mesh%NPanels*2**Mesh%ISym), INTENT(IN) :: cdata
    CHARACTER(LEN=*),                              INTENT(IN) :: filename ! Output file

    ! Local variables
    INTEGER            :: u, i, j
    REAL, DIMENSION(3) :: x

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    IF (output_format == TECPLOT_OUTPUT) THEN
      WRITE(u,*) 'VARIABLES="X" "Y" "Z" "abs(p) (Pa)" "angle(p) (rad)"'
      WRITE(u,'(A,I7,A,I7,A)') 'ZONE N=', Mesh%Npoints*2**Mesh%ISym,',E = ', Mesh%Npanels*2**Mesh%ISym,', F=FEPOINT, ET=QUADRILATERAL'
    END IF

    DO i = 1, Mesh%NPanels
      DO j = 1, 4
        WRITE(u, *) Mesh%X(:, Mesh%P(j, i)), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
      END DO
    END DO

    IF (Mesh%ISym == Y_SYMMETRY) THEN
      DO i = Mesh%NPanels+1, 2*Mesh%NPanels
        DO j = 1, 4
          x(:) = Mesh%X(:, Mesh%P(j, i-Mesh%NPanels))
          x(2) = -x(2)
          WRITE(u, *) x(:), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
        END DO
      END DO
    END IF

    IF (output_format == TECPLOT_OUTPUT) THEN
      DO i=1,Mesh%Npanels
        WRITE(u, *) 1+(i-1)*4, 2+(i-1)*4, 3+(i-1)*4, 4+(i-1)*4
      END DO

      IF (Mesh%ISym == Y_SYMMETRY) THEN
        DO i=1,Mesh%Npanels
          WRITE(u, *) 1+(i-1)*4+4*Mesh%NPanels, 2+(i-1)*4+4*Mesh%NPanels, 3+(i-1)*4+4*Mesh%NPanels, 4+(i-1)*4+4*Mesh%NPanels
        END DO
      END IF
    END IF

    CLOSE(u)

  END SUBROUTINE WRITE_DATA_ON_MESH

END MODULE OUTPUT
