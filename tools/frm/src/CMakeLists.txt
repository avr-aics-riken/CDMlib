message("entered ${PROJECT_SOURCE_DIR}/tools/frm/src")

set(bin_PROGRAMS "frm")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/tools/frm/include
  ${TP_CFLAGS}
)

set(frm_SOURCES
    Staging.C
    ActiveSubDomain.C
    main.C
)

link_directories(
   ${PROJECT_SOURCE_DIR}/src 
   ${CDM_LDFLAGS}
#   ${TP_DIR}/lib
   ${TP_LDFLAGS}
)

add_definitions(${BUFFER_TUNING_FLAGS})

add_executable(${bin_PROGRAMS} ${frm_SOURCES})
target_link_libraries(${bin_PROGRAMS} ${CDM_LIBS} ${TP_LIBS})

install(TARGETS ${bin_PROGRAMS} RUNTIME DESTINATION
${CMAKE_INSTALL_PREFIX}/bin)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/tools/frm/include DESTINATION
${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.h")
