find_package(LAPACK REQUIRED)

set(POSTPROC_SRC
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/cPackgmres.f
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/zPackgmres.f
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/blas_rot.f
    ${CORE_SRC_DIR}/M_SOLVER.f90
    ${COMMON_SRC_DIR}/Constants.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    ${COMMON_SRC_DIR}/Identification.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/Results.f90
    ${COMMON_SRC_DIR}/Mesh.f90
    ${COMMON_SRC_DIR}/MNemohCal.f90
    MPP_ReadInputFiles.f90
    MPP_Compute_RAOs.f90
    IRF.f90
    Plot_WaveElevation.f90
    Main.f90)

add_executable(postProc ${POSTPROC_SRC})

target_link_libraries(postProc ${LAPACK_LIBRARIES})
