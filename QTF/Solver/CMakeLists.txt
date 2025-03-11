set(PREPROC_SRC
    ${COMMON_SRC_DIR}/Identification.f90
    ${COMMON_SRC_DIR}/Logfile.f90
    ${COMMON_SRC_DIR}/Constants.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    ${EXTERNAL_SRC_DIR}/linear_interpolation_module.f90
    ${COMMON_SRC_DIR}/MCallInterp.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/MNemohCal.f90
    ${COMMON_SRC_DIR}/Mesh.f90
    ${COMMON_SRC_DIR}/Face.f90
    ${EXTERNAL_SRC_DIR}/MROMBERG.f90
    ${EXTERNAL_SRC_DIR}/MBESSEL.f90
    ${QTF_COMMON_SRC_DIR}/MFileDirectoryList.f90
    ${QTF_COMMON_SRC_DIR}/MReadInputFiles.f90
    MQSolverPreparation.f90
    MQSolverOutputFiles.f90
    MQSolverAsymp.f90
    MQSolver.f90
    Main.f90)

add_executable(QTFsolver ${PREPROC_SRC})