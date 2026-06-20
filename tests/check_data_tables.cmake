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

set(expected_md5_crams_fragmentation_evoli2019_csv b4b4696044903976df05ca383632d929)
set(expected_md5_crams_fragmentation_evoli2026_st99_csv c0d3a47e36ad1baf6533810eab20b3a5)
set(expected_md5_crams_fragmentation_evoli2026_w93_csv 347d6ae779458267b55a3c51b2a69c91)
set(expected_md5_crams_inelastic_crosec_csv d6d72918c4cd5f479164ade08f380f30)
set(expected_md5_crams_inelastic_glauber_csv d575bf2a43ecad46dd405d8a129e2beb)
set(expected_md5_crams_inelastic_tripathi99_csv 8c5d4d51866ae5e8e7fa6700a481f2a5)

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
