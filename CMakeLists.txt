# From: https://github.com/Mizux/cmake-swig LICENSE: Apache License 2.0

# This file is just an orchestration
cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Apple: Don't modify install_name when touching RPATH.
if(POLICY CMP0068)
  cmake_policy(SET CMP0068 NEW)
endif()

# Define file(GENERATE) behavior for relative paths.
if(POLICY CMP0070)
  cmake_policy(SET CMP0070 NEW)
endif()

# target_sources: use absolute path for INTERFACE_SOURCES.
if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()

# option() honors normal variables. see:
# https://cmake.org/cmake/help/git-stage/policy/CMP0077.html
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

# SWIG: use standard target name.
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()

# SWIG: use SWIG_MODULE_NAME property.
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

# Enable Python find by location
if(POLICY CMP0094)
  cmake_policy(SET CMP0094 NEW)
endif()

project(
  cspy
  VERSION 1.0.4
  LANGUAGES CXX)
message(STATUS "project: ${PROJECT_NAME}")
message(STATUS "version: ${PROJECT_VERSION}")

# Force default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "Choose the type of build, options are: Debug, Release (default), RelWithDebInfo and MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Layout build dir like install dir
include(GNUInstallDirs)
if(UNIX)
  option(BUILD_SHARED_LIBS "Build shared libraries(.so or .dyld)." ON)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  # for multi-config build system (e.g. Xcode, ninja Multi-Config)
  foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
  endforeach()
else()
  option(BUILD_SHARED_LIBS "Build shared libraries (.dll)." OFF)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  # for multi-config builds (e.g. msvc)
  foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG}
        ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
  endforeach()
endif()

# By default only the ortools C++ library is built.
option(BUILD_CXX "Build C++ library" ON)
option(BUILD_PYTHON "Build Python Library" OFF)
option(BUILD_DOTNET "Build Dotnet Library" OFF)
set(LOG_LEVEL
    "WARN"
    CACHE
      STRING
      "Choose the logging level, options are: OFF (default), WARN, INFO, DEBUG, TRACE."
)

message(STATUS "Build C++ library: ${BUILD_CXX}, LOG_LEVEL= ${LOG_LEVEL}")
message(STATUS "Build Python: ${BUILD_PYTHON}")
message(STATUS "Build Dotnet: ${BUILD_DOTNET}")

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
set(LEMON_BUILD_TESTING OFF)
# Remove shared libs for LEMON
set(BUILD_SHARED_LIBS_OLD ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
FetchContent_Declare(
  LEMON
  GIT_REPOSITORY "https://github.com/MultiFlow/LEMON.git"
  GIT_TAG "multiflow/dev")
FetchContent_MakeAvailable(LEMON)
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_OLD})

# Add logger
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
  GIT_TAG "v1.x")
FetchContent_MakeAvailable(spdlog)

enable_testing()
include(cpp)

if(WIN32
   AND (BUILD_PYTHON
        OR BUILD_JAVA
        OR BUILD_DOTNET))
  message(STATUS "Getting SWIG: ...")
  include(swig)
  message(STATUS "Getting SWIG: ...DONE")
endif()

include(python)
include(dotnet)

# Fetch googletest and add C++ tests only if required
if(NOT (BUILD_PYTHON OR BUILD_DOTNET) AND BUILD_TESTING)
  set(FETCHCONTENT_QUIET OFF)
  # Prevent overriding the parent project's compiler/linker settings on Windows
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    GIT_TAG "main")
  FetchContent_MakeAvailable(googletest)
  add_subdirectory(test/cc/)
endif()
