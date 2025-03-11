function(set_number_of_frequencies NEMOH_TEST_PATH N_FREQ)
    file(READ ${NEMOH_TEST_PATH}/Nemoh.cal FILE_CONTENTS)
    string(REGEX REPLACE "([123][ \t]+)[0-9]+([ \t]+[0-9]+\\.?[0-9]*[ \t]+[0-9]+\\.?[0-9]*[ \t]+! ?Freq type)" "\\1${N_FREQ}\\2" FILE_CONTENTS ${FILE_CONTENTS})
    file(WRITE ${NEMOH_TEST_PATH}/Nemoh.cal ${FILE_CONTENTS})

    if(${ARGC} GREATER 2) # If an additional argument is present
        file(READ ${NEMOH_TEST_PATH}/Nemoh.cal FILE_CONTENTS)
        string(REGEX REPLACE "[0-9]+([ \t]+[0-9]+\\.?[0-9]*[ \t]+[0-9]+\\.?[0-9]*[ \t]+! Number of radial frequencies)" "${ARGV2}\\1" FILE_CONTENTS ${FILE_CONTENTS})
        file(WRITE ${NEMOH_TEST_PATH}/Nemoh.cal ${FILE_CONTENTS})
    endif()
endfunction()

add_subdirectory(1_Cylinder)
add_subdirectory(my_Cylinder)