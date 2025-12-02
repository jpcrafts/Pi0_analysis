cmake_minimum_required(VERSION 3.16)

project(Pi0Tools LANGUAGES CXX)

# Make our own cmake/ folder visible
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Centralized dependency handling (ROOT + yaml-cpp, interface target, etc.)
include(Dependencies)

# ==========================
# Top-level standalone .cxx
# ==========================

# All top-level .cxx that should become executables
set(TOP_LEVEL_SOURCES
    Simc_import_take5.cxx
    fit_mgg_minimal.cxx
    reduce_bins.cxx
    make_bins_schema_MxCorr_Mgg.cxx
    fit_mgg_roofit_sb.cxx
    Pi_0_production.cxx
    skimDummy.cxx
)

foreach(src ${TOP_LEVEL_SOURCES})
    get_filename_component(exe_name "${src}" NAME_WE)
    add_executable("${exe_name}" "${src}")
    target_link_libraries("${exe_name}" PRIVATE project_deps)
endforeach()

# If you want specific exe names that differ from filename, you can override,
# e.g.:
# add_executable(Pi_0_v5_slim_pairsA Pi_0_production.cxx)
# target_link_libraries(Pi_0_v5_slim_pairsA PRIVATE project_deps)

# ==========================
# Subdirectories (tools, studies)
# ==========================
add_subdirectory(tools)
add_subdirectory(studies)
