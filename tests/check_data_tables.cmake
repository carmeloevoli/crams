if(NOT DEFINED CRAMS_DATA_DIR)
  message(FATAL_ERROR "CRAMS_DATA_DIR is not defined")
endif()

set(expected_files
  crams_fragmentation_evoli2019.csv
  crams_fragmentation_evoli2026_st99.csv
  crams_fragmentation_evoli2026_w93.csv
  crams_inelastic_crosec.csv
  crams_inelastic_glauber.csv
  crams_inelastic_tripathi99.csv
)

set(expected_md5_crams_fragmentation_evoli2019_csv 7f52515c55fa4e3a5144355c558c415a)
set(expected_md5_crams_fragmentation_evoli2026_st99_csv 8892f0664678af3ab890f31c7acc8ccf)
set(expected_md5_crams_fragmentation_evoli2026_w93_csv 28813d9b6fc04b878c1d9ddfe27fb73e)
set(expected_md5_crams_inelastic_crosec_csv b5a7b07142ab2376442641728ad60e54)
set(expected_md5_crams_inelastic_glauber_csv 9f6f82d4ffeaee0a453bdb9284c0394f)
set(expected_md5_crams_inelastic_tripathi99_csv 879a4dc440c5a6396f036d481743e84a)

file(GLOB actual_files
  RELATIVE "${CRAMS_DATA_DIR}"
  "${CRAMS_DATA_DIR}/crams_fragmentation_evoli*.csv"
  "${CRAMS_DATA_DIR}/crams_inelastic_*.csv"
)
list(SORT actual_files)

if(NOT actual_files STREQUAL expected_files)
  message(FATAL_ERROR
    "Data table set changed.\n"
    "Expected: ${expected_files}\n"
    "Actual:   ${actual_files}\n"
    "Update tests/check_data_tables.cmake if this change is intentional."
  )
endif()

foreach(file IN LISTS expected_files)
  set(path "${CRAMS_DATA_DIR}/${file}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Missing data table: ${path}")
  endif()

  file(MD5 "${path}" actual_md5)
  string(MAKE_C_IDENTIFIER "expected_md5_${file}" expected_var)

  if(NOT actual_md5 STREQUAL "${${expected_var}}")
    message(FATAL_ERROR
      "MD5 mismatch for ${file}\n"
      "Expected: ${${expected_var}}\n"
      "Actual:   ${actual_md5}\n"
      "Update tests/check_data_tables.cmake if this table change is intentional."
    )
  endif()
endforeach()

message(STATUS "Verified MD5 checksums for ${expected_files}")
