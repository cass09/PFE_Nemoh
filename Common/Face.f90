!--------------------------------------------------------------------------------------
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
!   - R. Kurnia
!--------------------------------------------------------------------------------------

MODULE MFace
  ! A single face from a mesh.
  USE CONSTANTS
  USE MMesh
  USE Elementary_functions, ONLY: CROSS_PRODUCT
  IMPLICIT NONE

  PUBLIC        :: Prepare_FaceMesh,VFace_to_Face,New_Face_Extracted_From_Mesh, &
                   Prepare_Waterline
  PRIVATE       :: Face_to_VFACE,PREPARE_GAUSS_QUADRATURE_ON_MESH

  TYPE TDATGQ
    REAL, DIMENSION(16,4):: XGQ
    REAL, DIMENSION(16,4):: YGQ
    REAL, DIMENSION(16,4):: WGQ
  END TYPE
  INTEGER,parameter       :: FLAG_GQ=1
  !INTEGER,parameter       :: NP_GQ=4

    ! A single face extracted from the mesh
  TYPE TFace
    REAL, DIMENSION(3,5)    :: X ! Vertices coordinates, the 5th is the same as the 1st one.
    REAL, DIMENSION(3)      :: XM   ! Center of mass
    REAL, DIMENSION(3)      :: N    ! Normal vector
    REAL                    :: A    ! Area
    REAL                    :: tDis ! Maximal radius of the panel
    INTEGER                 :: NP_GQ! Number of point of the Gauss quadrature
    REAL, ALLOCATABLE       :: dXdXG_WGQ_per_A(:) !a factor dxdXG*WGQ/A for the Gauss Q Integral
    REAL, ALLOCATABLE       :: XM_GQ(:,:)         !Gauss Quad. nodes (x,y,z) in a panel
                                                  !Those allocatable becase we need input NP_GQ from user input
  END TYPE TFace

  ! A vector of Face with size Npanels x 1
  TYPE TVFace
    REAL, DIMENSION(:,:,:),ALLOCATABLE  :: X    ! Vertices coordinates, the 5th is the same as the 1st one.
    REAL, DIMENSION(:,:)  ,ALLOCATABLE  :: XM   ! Center of mass
    REAL, DIMENSION(:,:)  ,ALLOCATABLE  :: N    ! Normal vector
    REAL, DIMENSION(:)    ,ALLOCATABLE  :: A    ! Area
    REAL, DIMENSION(:)    ,ALLOCATABLE  :: tDis ! Maximal radius of the panel
    INTEGER                             :: NP_GQ! Number of point of the Gauss quadrature
    REAL,DIMENSION(:,:)  ,ALLOCATABLE   :: dXdXG_WGQ_per_A  !a factor dxdXG*WGQ/A for the Gauss Q Integral
    REAL,DIMENSION(:,:,:),ALLOCATABLE   :: XM_GQ            !Gauss Quad. nodes (x,y,z) in a panel
  END TYPE

  ! A vector of waterline segments
  TYPE TWLine
    REAL, DIMENSION(:,:),ALLOCATABLE    :: XM             ! Centroid of Waterline segment
    REAL, DIMENSION(:), ALLOCATABLE     :: SegLength      ! Segment length
    INTEGER, DIMENSION(:), ALLOCATABLE  :: IndexPanel     ! Index of panel where the
    INTEGER                             :: NWlineseg      ! Number of waterline segments
  END TYPE

CONTAINS

  SUBROUTINE Prepare_FaceMesh(Mesh,NP_GQ, VFace)

    !INPUT/OUTPUT
    INTEGER     , INTENT(IN)    :: NP_GQ
    TYPE(TMesh) , INTENT(IN)    :: Mesh
    TYPE(TVFace), INTENT(INOUT) :: VFace
    !Local variables
    INTEGER                     :: I,Npanels
    TYPE(TFace)                 :: Face
    TYPE(TDATGQ)                :: DataGQ       !Data Gauss Quadrature
    Npanels=Mesh%Npanels

    Face%NP_GQ=NP_GQ

    Allocate(Face%dXdXG_WGQ_per_A(NP_GQ))
    Allocate(Face%XM_GQ(3,NP_GQ))

    Allocate(VFace%X(Npanels,3,5),VFace%XM(Npanels,3))
    Allocate(VFace%N(Npanels,3),VFace%A(Npanels),VFace%tDis(Npanels))
    Allocate(VFace%dXdXG_WGQ_per_A(Npanels,NP_GQ))
    Allocate(VFace%XM_GQ(Npanels,3,NP_GQ))

    CALL GAUSS_QUADRATURE_DATA(DataGQ)

    DO I=1,Npanels
         CALL NEW_FACE_EXTRACTED_FROM_MESH(Mesh, I, Face)
         CALL PREPARE_GAUSS_QUADRATURE_ON_MESH(Mesh,Face,DATAGQ)
         CALL Face_to_VFace(Face,VFace,I)
    END DO

    DEAllocate(Face%dXdXG_WGQ_per_A)
    DEAllocate(Face%XM_GQ)

    END SUBROUTINE

  SUBROUTINE Face_to_VFace(Face,VFACE,I)
    !INPUT/OUTPUT
    INTEGER,         INTENT(IN)  :: I
    TYPE(TFace),     INTENT(IN)  :: Face
    TYPE(TVFace),    INTENT(INOUT) :: VFace

    VFace%X(I,:,:)              =Face%X(:,:)
    VFace%XM(I,:)               =Face%XM(:)
    VFace%N(I,:)                =Face%N(:)
    VFace%A(I)                  =Face%A
    VFace%tDis(I)               =Face%tDis
    VFace%dXdXG_WGQ_per_A(I,:)  =Face%dXdXG_WGQ_per_A(:)
    VFace%XM_GQ(I,:,:)          =Face%XM_GQ(:,:)
    VFace%NP_GQ                 =Face%NP_GQ
  END SUBROUTINE

  SUBROUTINE VFace_to_Face(VFace,FACE,I)
    !INPUT/OUTPUT
    INTEGER,        INTENT(IN)  :: I
    TYPE(TVFace),   INTENT(IN)  :: VFace
    TYPE(TFace),    INTENT(INOUT) :: Face

    IF(.NOT.ALLOCATED(Face%dXdXG_WGQ_per_A)) THEN
      Allocate(Face%dXdXG_WGQ_per_A(VFace%NP_GQ))
      Allocate(Face%XM_GQ(3,VFace%NP_GQ))
    ENDIF

    Face%X(:,:)              =VFace%X(I,:,:)
    Face%XM(:)               =VFace%XM(I,:)
    Face%N(:)                =VFace%N(I,:)
    Face%A                   =VFace%A(I)
    Face%tDis                =VFace%tDis(I)
    Face%NP_GQ               =VFace%NP_GQ
    Face%dXdXG_WGQ_per_A(:)  =VFace%dXdXG_WGQ_per_A(I,:)
    Face%XM_GQ(:,:)          =VFace%XM_GQ(I,:,:)
  END SUBROUTINE

  SUBROUTINE New_Face_Extracted_From_Mesh(Mesh, I, Face)
    TYPE(TMesh), INTENT(IN)  :: Mesh
    INTEGER,     INTENT(IN)  :: I ! The index of the face in the mesh
    TYPE(TFace), INTENT(INOUT) :: Face

    ! Local variables
    REAL, DIMENSION(3)       :: M0

    Face%X(1:3, 1:4) = Mesh%X(1:3, Mesh%P(1:4, I))
    Face%X(1:3, 5)   = Mesh%X(1:3, Mesh%P(1, I))

    Face%XM(1:3)     = Mesh%XM(1:3, I)
    Face%N(1:3)      = Mesh%N(1:3, I)

    Face%A           = Mesh%A(I)

    ! Compute max radius
    M0(1:3)   = SUM(Face%X(1:3, 1:4), DIM=2)/4 ! Average vertex
    Face%tDis = MAX(                   &
      NORM2(Face%X(1:3, 1) - M0(1:3)), &
      NORM2(Face%X(1:3, 2) - M0(1:3)), &
      NORM2(Face%X(1:3, 3) - M0(1:3)), &
      NORM2(Face%X(1:3, 4) - M0(1:3))  &
      )
  END SUBROUTINE New_Face_Extracted_From_Mesh

   SUBROUTINE PREPARE_GAUSS_QUADRATURE_ON_MESH(Mesh,Face,DATAGQ)
    TYPE(TMesh), INTENT(IN)    :: Mesh
    TYPE(TDATGQ)               :: DataGQ       ! Data Gauss Quadrature
    TYPE(TFace), INTENT(INOUT) :: Face

    ! Local variables
    INTEGER                       :: NP_GQ        ! Number of poits Gauss Quadrature
    INTEGER                       :: IdColumn,L
    REAL, DIMENSION(3)            :: Tangen1,Tangen2,Tangen2c !Unit tangential vector in a panel I
    REAL, DIMENSION(3)            :: Normal          !Unit normal vector
    REAL, DIMENSION(4)            :: XL,YL           !
    REAL                          :: XG1,YG1,AA,BB,CC,DD
    REAL                          :: A,B,C,D

    IF (FLAG_GQ==1) THEN
         NP_GQ=Face%NP_GQ
         IdColumn=NINT(SQRT(REAL(NP_GQ)))

         Tangen1=Face%X(1:3, 3)-Face%X(1:3, 1)
         Tangen2=Face%X(1:3, 4)-Face%X(1:3, 2)
         !Normal=CROSS_PRODUCT(Tangen1,Tangen2)
         !Normal=Normal/sqrt(Normal(1)**2+Normal(2)**2+Normal(3)**2)
         Normal=Face%N(1:3) !Unit normal vector had been calculated before

         Tangen1=Tangen1(1:3)/SQRT(Tangen1(1)**2+Tangen1(2)**2+Tangen1(3)**2) !Unit Tangen
         !Tangen2=Tangen2(1:3)/SQRT(Tangen2(1)**2+Tangen2(2)**2+Tangen2(3)**2)
         Tangen2=CROSS_PRODUCT(Normal,Tangen1)

         XL(:)=Tangen1(1)*(Face%X(1,1:4)-Face%XM(1))&
                   +Tangen1(2)*(Face%X(2,1:4)-Face%XM(2))&
                          +Tangen1(3)*(Face%X(3,1:4)-Face%XM(3))
         YL(:)=-Tangen2(1)*(Face%X(1,1:4)-Face%XM(1))&
                   -Tangen2(2)*(Face%X(2,1:4)-Face%XM(2))&
                           -Tangen2(3)*(Face%X(3,1:4)-Face%XM(3))

         DO L=1,NP_GQ
          AA=.25*(1-DATAGQ%XGQ(L,IdColumn))*(1-DATAGQ%YGQ(L,IdColumn))
          BB=.25*(1-DATAGQ%XGQ(L,IdColumn))*(1+DATAGQ%YGQ(L,IdColumn))
          CC=.25*(1+DATAGQ%XGQ(L,IdColumn))*(1+DATAGQ%YGQ(L,IdColumn))
          DD=.25*(1+DATAGQ%XGQ(L,IdColumn))*(1-DATAGQ%YGQ(L,IdColumn))

          XG1=AA*XL(1)+BB*XL(2)+CC*XL(3)+DD*XL(4)
          YG1=AA*YL(1)+BB*YL(2)+CC*YL(3)+DD*YL(4)

          Face%XM_GQ(1:3,L)=Face%XM(1:3)+Tangen1(1:3)*XG1-Tangen2(1:3)*YG1

          A=(1.-DATAGQ%YGQ(L,IdColumn))*(XL(4)-XL(1))+(1+DATAGQ%YGQ(L,IdColumn))*(XL(3)-XL(2))
          B=(1.-DATAGQ%XGQ(L,IdColumn))*(YL(2)-YL(1))+(1+DATAGQ%XGQ(L,IdColumn))*(YL(3)-YL(4))
          C=(1.-DATAGQ%XGQ(L,IdColumn))*(XL(2)-XL(1))+(1+DATAGQ%XGQ(L,IdColumn))*(XL(3)-XL(4))
          D=(1.-DATAGQ%YGQ(L,IdColumn))*(YL(4)-YL(1))+(1+DATAGQ%YGQ(L,IdColumn))*(YL(3)-YL(2))
          Face%dXdXG_WGQ_per_A(L)=(ABS(A*B-C*D)*DATAGQ%WGQ(L,IdColumn)*.25)/Face%A
         ENDDO
    ELSE
         Face%dXdXG_WGQ_per_A(:)= 1
         Face%XM_GQ(:,1)=Face%XM(:)
    ENDIF
   END SUBROUTINE PREPARE_GAUSS_QUADRATURE_ON_MESH


  SUBROUTINE GAUSS_QUADRATURE_DATA(DATAGQ)
  !This subroutine provide Gauss quadrature nodes and weigthing function data

  !INPUT/OUTPUT
  TYPE(TDATGQ),         INTENT(OUT):: DATAGQ
  !Local
  REAL,DIMENSION(16,4)  :: XGQ,YGQ,WGQ
  INTEGER::I,L

        DATA((XGQ(I,L),I=1,16),L=1,4) /&
         16*0., &
         .57735027, .57735027,-.57735027,-.57735027,12*0.,&
         .77459667, .77459667, .77459667,3*0.,-.77459667,&
        -.77459667,-.77459667,7*0.,&
         .86113631, .86113631, .86113631, .86113631,&
         .33998104, .33998104, .33998104, .33998104,&
        -.33998104,-.33998104,-.33998104,-.33998104,&
        -.86113631,-.86113631,-.86113631,-.86113631/
        DATA((YGQ(I,L),I=1,16),L=1,4)/ &
         16*0., &
        -.57735027,.57735027,-.57735027,.57735027,12*0.,&
        -.77459667,0.,.77459667,-.77459667,0.,.77459667,&
        -.77459667,0.,.77459667,7*0.,&
        -.86113363,-.33998104,.33998104,.86113363,&
        -.86113363,-.33998104,.33998104,.86113363,&
        -.86113363,-.33998104,.33998104,.86113363,&
        -.86113363,-.33998104,.33998104,.86113363/
        DATA((WGQ(I,L),I=1,16),L=1,4)/ &
         1,15*0., &
        .25,.25,.25,.25,12*0.,&
        .07716049,.12345679,.07716049,.12345679,.19753086,&
        .12345679,.07716049,.12345679,.07716049,7*0.,&
        .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1,&
        .56712963E-1,.10632333,.10632333,.56712963E-1,&
        .56712963E-1,.10632333,.10632333,.56712963E-1,&
        .30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1/

  DATAGQ%XGQ=XGQ(:,:)
  DATAGQ%YGQ=YGQ(:,:)
  DATAGQ%WGQ=WGQ(:,:)

  END SUBROUTINE

  SUBROUTINE Prepare_Waterline(VFace,EPS,BodyDiameter,Npanels,WLine)
     !Input/output
     TYPE(TVFace), INTENT(IN)    :: VFace
     REAL,         INTENT(IN)    :: EPS
     REAL,         INTENT(IN)    :: BodyDiameter
     INTEGER,      INTENT(IN)    :: Npanels
     TYPE(TWLine), INTENT(OUT)   :: WLine
     !Local
     REAL,DIMENSION(Npanels,3)     :: XM_WLine
     REAL,DIMENSION(Npanels)       :: SegLength
     INTEGER, DIMENSION(Npanels)   :: IndexPanel
     INTEGER                       :: IlineSeg,I,J
     REAL                          :: Distance,ZER

     IlineSeg=0
     ZER=-EPS*BodyDiameter
     DO I=1,Npanels
        ! restrict to body panels, not the lid panels
        IF(VFace%XM(I,3).LT.ZER) THEN
          DO J=1,4 !Nodes index in the panel-I
            DISTANCE=NORM2(VFace%X(I,:,J+1)-VFace%X(I,:,J))
            IF (DISTANCE.GT.EPS) THEN
               ! restrict to Z nodes on the waterline
               IF (VFace%X(I,3,J+1)+VFace%X(I,3,J).GE.ZER) THEN
                  IlineSeg=IlineSeg+1
                  XM_WLine(IlineSeg,1:2)=(VFace%X(I,1:2,J+1)+VFace%X(I,1:2,J))*0.5 &
                                          +VFace%N(I,1:2)*0.01*DISTANCE
                  XM_WLine(IlineSeg,3)=ZER
                  SegLength(IlineSeg)=DISTANCE
                  IndexPanel(IlineSeg)=I
               ENDIF
            ENDIF
          ENDDO
        ENDIF
     ENDDO
    ALLOCATE(WLine%XM(IlineSeg,3),WLine%SegLength(IlineSeg))
    ALLOCATE(WLine%IndexPanel(IlineSeg))
    WLine%NWlineseg=IlineSeg
    Wline%XM(1:IlineSeg,1:3)=XM_WLine(1:IlineSeg,1:3)
    Wline%SegLength(1:IlineSeg)=SegLength(1:IlineSeg)
    Wline%IndexPanel(1:IlineSeg)=IndexPanel(1:IlineSeg)
  END SUBROUTINE

  ! SUBROUTINE Reflect_Face_Around_XZ_Plane(Face)
  !   ! Change y coordinate into -y
  !   TYPE(TFace), INTENT(INOUT) :: Face

  !   Face%X(2, :) = -Face%X(2, :)
  !   Face%XM(2)   = -Face%XM(2)
  !   Face%N(2)    = -Face%N(2)
  !   Face%X(:, 1:5) = Face%X(:, 5:1:-1) ! Invert order of vertices to be coherent with normal vector.
  ! END SUBROUTINE Reflect_Face_Around_XZ_Plane


  ! SUBROUTINE Reflect_Face_Around_XY_Plane(Face, z0)
  !   ! Change z coordinate into z0 - z
  !   TYPE(TFace), INTENT(INOUT) :: Face
  !   REAL, INTENT(IN)           :: z0

  !   Face%X(3, :) = z0-Face%X(3, :)
  !   Face%XM(3)   = z0-Face%XM(3)
  !   Face%N(3)    =   -Face%N(3)
  !   Face%X(:, 1:5) = Face%X(:, 5:1:-1) ! Invert order of vertices to be coherent with normal vector.
  ! END SUBROUTINE Reflect_Face_Around_XY_Plane

END MODULE MFace
