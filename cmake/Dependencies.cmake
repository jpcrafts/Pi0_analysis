# cmake/Dependencies.cmake

include_guard(GLOBAL)

# ==============================
# Global C++ / compile features
# ==============================
# Centralized: everything that links project_deps gets C++17.
add_library(project_deps INTERFACE)
target_compile_features(project_deps INTERFACE cxx_std_17)

# =======================
# ROOT (modern CMake use)
# =======================
# Assumes a proper ROOTConfig.cmake install.
# Add/remove components as needed.
find_package(ROOT CONFIG REQUIRED COMPONENTS
    Core
    Hist
    Tree
    RIO
    Graf
    Gpad
    MathCore
    Matrix
    Physics
    RooFit
    RooFitCore
    TMVA
)

set(PROJECT_ROOT_LIBS
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::RIO
    ROOT::Graf
    ROOT::Gpad
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
    ROOT::RooFit
    ROOT::RooFitCore
    ROOT::TMVA
)

# =================
# yaml-cpp handling
# =================
# Try to find a system yaml-cpp first.
find_package(yaml-cpp CONFIG QUIET)

if(NOT yaml-cpp_FOUND)
    # Fallback: fetch and build yaml-cpp
    include(FetchContent)
    message(STATUS "yaml-cpp not found; fetching with FetchContent")

    FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG yaml-cpp-0.8.0  # Adjust tag if you prefer
    )

    # Disable extras
    set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(yaml-cpp)
endif()

# Normalize the yaml-cpp target name
if(TARGET yaml-cpp::yaml-cpp)
    set(YAML_CPP_LIB yaml-cpp::yaml-cpp)
elseif(TARGET yaml-cpp)
    set(YAML_CPP_LIB yaml-cpp)
else()
    message(FATAL_ERROR "yaml-cpp target not found after find_package/FetchContent")
endif()

# ==========================
# Wire dependencies together
# ==========================

target_link_libraries(project_deps
    INTERFACE
        ${PROJECT_ROOT_LIBS}
        ${YAML_CPP_LIB}
)
