set(PREPROC_SRC
    ${COMMON_SRC_DIR}/Identification.f90
    ${COMMON_SRC_DIR}/Logfile.f90
    ${COMMON_SRC_DIR}/Constants.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/MNemohCal.f90
    ${COMMON_SRC_DIR}/Mesh.f90
    ${COMMON_SRC_DIR}/Face.f90
    ${CORE_SRC_DIR}/GREEN_1.f90
    ${CORE_SRC_DIR}/INITIALIZE_GREEN.f90
    ${CORE_SRC_DIR}/GREEN_2.f90
    ${QTF_COMMON_SRC_DIR}/MFileDirectoryList.f90
    ${QTF_COMMON_SRC_DIR}/MReadInputFiles.f90
    MInfluenceMatrix.f90
    MQpreprocessor.f90
    Main.f90)

add_executable(QTFpreProc ${PREPROC_SRC})