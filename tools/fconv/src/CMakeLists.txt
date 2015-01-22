message("entered ${PROJECT_SOURCE_DIR}/tools/frm/src")

set(bin_PROGRAMS "fconv")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/tools/fconv/include
  ${PROJECT_SOURCE_DIR}/tools/fconv/include/inline
  ${CPM_CFLAGS}
  ${TP_CFLAGS}
  ${MPI_CFLAGS}
)

set(fconv_SOURCES
  main.C
  InputParam.C
  conv.C
  convMx1.C
  convMxM.C
  convMxN.C
  convOutput.C
  convOutput_SPH.C
  convOutput_BOV.C
  convOutput_AVS.C
  convOutput_VTK.C
  convOutput_PLOT3D.C
)

link_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${CDM_LDFLAGS}
  ${TP_LDFLAGS}
  ${CPM_LDFLAGS}
  ${MPI_LDFLAGS}
)

add_definitions(${BUFFER_TUNING_FLAGS})

add_executable(${bin_PROGRAMS} ${fconv_SOURCES})
target_link_libraries(${bin_PROGRAMS} ${CDM_LIBS} ${TP_LIBS} ${CPM_LIBS}
${MPI_LIBS} ${MPI_LIBS2})

install(TARGETS ${bin_PROGRAMS} RUNTIME DESTINATION
${CMAKE_INSTALL_PREFIX}/bin)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/tools/fconv/include DESTINATION
${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.h")

install(DIRECTORY ${PROJECT_SOURCE_DIR}/tools/fconv/include/inline DESTINATION
${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")
