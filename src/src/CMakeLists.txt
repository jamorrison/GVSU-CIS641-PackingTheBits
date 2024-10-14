# Minimum CMake version
cmake_minimum_required(VERSION 3.21)

# Source and header files
file(GLOB sources *.c)
file(GLOB headers *.h)

# main.c is what executable is compiled from
list(REMOVE_ITEM sources main.c)

# Add directories for included libraries
include_directories(../lib/utils)
include_directories(../lib/htslib/htslib)

# Need to create a core library that gets linked with the executable
add_library(coverage_core ${sources} ${headers})
target_include_directories(coverage_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Build coverage executable
add_executable(coverage main.c)

target_link_libraries(coverage
    coverage_core
    m
    ${LIB_DIR}/utils/libutils.a
    ${LIB_DIR}/htslib/libhts.a
    )

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    target_link_libraries(coverage ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for compilation")
endif(ZLIB_FOUND)

find_package(Threads REQUIRED)
if(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
    include_directories(${Threads_INCLUDE_DIRS})
    target_link_libraries(coverage ${Threads_LIBRARIES})
else()
    message(FATAL_ERROR "pthread not found. Required for compilation")
endif(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
