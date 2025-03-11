set(PREPROC_SRC
    ${COMMON_SRC_DIR}/Identification.f90
    ${COMMON_SRC_DIR}/MNemohCal.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/Constants.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    MQpostproc.f90
    Main.f90)

add_executable(QTFpostProc ${PREPROC_SRC})