###################################################################################
#
# CDMlib - Cartesian Data Management library
#
# Copyright (c) 2013-2017 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################


include_directories(
      ${PROJECT_BINARY_DIR}/include
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/include/inline
      ${TP_INC}
      ${NetCDF_INC}
)

link_directories(
      ${PROJECT_BINARY_DIR}/src
      ${TP_LIB}
      ${NetCDF_LIB}
)



#### Test1 : output uniform Cartesian

add_executable(output_uc ./sample_output_Uni.C)

if(with_MPI)
  target_link_libraries(output_uc -lCDMmpi -lTPmpi)
  set (test_parameters -np 4
                      "output_uc"
                      "${PROJECT_SOURCE_DIR}/Examples/cxx/data/input1.txt"
                      "${PROJECT_SOURCE_DIR}/Examples/cxx/data/input2.txt"
  )
  add_test(NAME TEST_1 COMMAND "mpirun" ${test_parameters}
  )
else()
  target_link_libraries(output_uc -lCDM -lCPM -lTP)
endif()



### Test2 : output non-uniform Cartesian

add_executable(output_nuc ./sample_output_NonUni.C)

if(with_MPI)
  target_link_libraries(output_nuc -lCDMmpi -lTPmpi)
  set (test_parameters -np 4
                      "output_nuc"
                      "${PROJECT_SOURCE_DIR}/Examples/cxx_FDM/data/input1.txt"
                      "${PROJECT_SOURCE_DIR}/Examples/cxx_FDM/data/input2.txt"
  )
  add_test(NAME TEST_2 COMMAND "mpirun" ${test_parameters}
  )
else()
  target_link_libraries(output_nuc -lCDM -lCPM -lTP)
endif()



### Test 3 : input

add_executable(input ./sample_input.C)

if(with_MPI)
  target_link_libraries(input -lCDMmpi -lTPmpi)
else()
  target_link_libraries(input -lCDM -lCPM -lTP)
endif()

if(with_MPI)
  set (test_parameters -np 4
                      "input"
                      "${PROJECT_SOURCE_DIR}/Examples/f90/data/in.jacobi"
  )
  add_test(NAME TEST_3 COMMAND "mpirun" ${test_parameters}
  )
else()
  add_test(TEST_3 input)
endif()
