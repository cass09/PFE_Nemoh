find_package(LAPACK REQUIRED)

set(SOLVER_SRC
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/cPackgmres.f
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/zPackgmres.f
    ${EXTERNAL_SRC_DIR}/GMRES_CERFACS/blas_rot.f
    ${CORE_SRC_DIR}/OUTPUT.f90
    ${CORE_SRC_DIR}/M_SOLVER.f90
    ${CORE_SRC_DIR}/GREEN_1.f90
    ${CORE_SRC_DIR}/INITIALIZE_GREEN.f90
    ${CORE_SRC_DIR}/GREEN_2.f90
    ${CORE_SRC_DIR}/SOLVE_BEM_DIRECT.f90
    ${CORE_SRC_DIR}/KOCHIN.f90
    ${CORE_SRC_DIR}/FREESURFACE.f90
    ${CORE_SRC_DIR}/FORCES.f90
    ${COMMON_SRC_DIR}/Constants.f90
    ${COMMON_SRC_DIR}/Logfile.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    ${COMMON_SRC_DIR}/Bodyconditions.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/Mesh.f90
    ${COMMON_SRC_DIR}/Face.f90
    NEMOH.f90)

add_executable(solver ${SOLVER_SRC})

target_link_libraries(solver ${LAPACK_LIBRARIES})
