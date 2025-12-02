# tools/CMakeLists.txt

# Grab all .cxx files in this directory
file(GLOB tools_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

foreach(src ${tools_sources})
    get_filename_component(exe_name "${src}" NAME_WE)
    add_executable("${exe_name}" "${src}")
    target_link_libraries("${exe_name}" PRIVATE project_deps)
endforeach()
