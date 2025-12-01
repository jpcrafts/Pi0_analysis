# studies/CMakeLists.txt

file(GLOB studies_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

foreach(src ${studies_sources})
    get_filename_component(exe_name "${src}" NAME_WE)
    add_executable("${exe_name}" "${src}")
    target_link_libraries("${exe_name}" PRIVATE project_deps)
endforeach()

# If you want a custom exe name for a particular file:
# add_executable(Pi_0_v5_slim_pairsA
#     Pi_0_v5_slim_pairsA_stable.cxx
# )
# target_link_libraries(Pi_0_v5_slim_pairsA PRIVATE project_deps)
