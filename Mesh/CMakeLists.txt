
set(MESH_SRC
    ${COMMON_SRC_DIR}/Identification.f90
    calCol.f90
    coque.f90
    ExMaillage.f90
    hydre.f90
    Mailleur.f90
    mesh.f90)

add_executable(mesh ${MESH_SRC})

set(HYDROS_SRC
    ${COMMON_SRC_DIR}/Identification.f90
    ${COMMON_SRC_DIR}/Elementary_functions.f90
    ${COMMON_SRC_DIR}/Environment.f90
    ${COMMON_SRC_DIR}/Mesh.f90
    coque.f90
    hydre.f90
    hydrostatic_cal.f90)

add_executable(hydrosCal ${HYDROS_SRC})