# Common boilerplate
cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prevent in-source builds
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory and run CMake from there")
endif()

# Determine if this is built as a subproject (using add_subdirectory)
# or if it is the master project.
set(MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MASTER_PROJECT ON)
  message(STATUS "Utils is MASTER PROJECT")
  message(STATUS "CMake version: ${CMAKE_VERSION}")
endif()

find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  # # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE UTILS_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE UTILS_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()


project(dominiqs-utils)

# Include zlib as an imported target
find_package(ZLIB REQUIRED)
add_library(ZLIB::Lib STATIC IMPORTED)
set_target_properties(ZLIB::Lib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${ZLIB_INCLUDE_DIRS}
  IMPORTED_LOCATION ${ZLIB_LIBRARIES}
)

# Define libutils
add_library(utils STATIC src/app.cpp src/base64.cpp src/compress.cpp src/it_display.cpp src/maths.cpp src/numbers.cpp 
                         src/path.cpp src/args_parser.cpp src/timer.cpp src/cutpool.cpp src/fileconfig.cpp)

target_include_directories(utils
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(utils PUBLIC fmt::fmt ZLIB::Lib)

add_library(Utils::Lib ALIAS utils)

# Generate version.h
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/utils/version.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/utils/version.h
)

# Add subprojects if this is master project
if (MASTER_PROJECT)
  add_subdirectory(extern/fmt)
endif()
