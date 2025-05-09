file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/. DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set_number_of_frequencies(${CMAKE_CURRENT_BINARY_DIR} 3 2)

set(TEST_NAME ${TEST_BASE_NAME}_B_FS)

add_test(NAME ${TEST_NAME}_preProc COMMAND preProc .)
add_test(NAME ${TEST_NAME}_solver COMMAND solver .)
add_test(NAME ${TEST_NAME}_postProc COMMAND postProc .)
add_test(NAME ${TEST_NAME}_QTFpreProc COMMAND QTFpreProc .)
add_test(NAME ${TEST_NAME}_QTFsolver COMMAND QTFsolver .)
add_test(NAME ${TEST_NAME}_QTFpostProc COMMAND QTFpostProc .)

set_tests_properties(${TEST_NAME}_preProc PROPERTIES LABELS "NEMOH2;PREPROC")
set_tests_properties(${TEST_NAME}_solver PROPERTIES DEPENDS ${TEST_NAME}_preProc LABELS "NEMOH2;SOLVER")
set_tests_properties(${TEST_NAME}_postProc PROPERTIES DEPENDS ${TEST_NAME}_solver LABELS "NEMOH2;POSTPROC")
set_tests_properties(${TEST_NAME}_QTFpreProc PROPERTIES DEPENDS ${TEST_NAME}_postProc LABELS "NEMOH2;PREPROC;QTF")
set_tests_properties(${TEST_NAME}_QTFsolver PROPERTIES DEPENDS ${TEST_NAME}_QTFpreProc LABELS "NEMOH2;SOLVER;QTF")
set_tests_properties(${TEST_NAME}_QTFpostProc PROPERTIES DEPENDS ${TEST_NAME}_QTFsolver LABELS "NEMOH2;POSTPROC;QTF")
