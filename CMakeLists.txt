cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

SET(PROJECT_NAME fit2x)
PROJECT(${PROJECT_NAME})
ENABLE_LANGUAGE(CXX)
SET(CMAKE_CXX_STANDARD 11)
add_definitions(-D_LIBCPP_ENABLE_CXX17_REMOVED_FEATURES)

## Check for Conda
###########################
IF(DEFINED ENV{CONDA_PREFIX})
    MESSAGE(STATUS "CONDA env seen: --[$ENV{CONDA_PREFIX}]--")
    SET(ENV{PATH} "$ENV{CONDA_PREFIX}/bin;$ENV{PATH}")
    # Define python_command that is used to call python and locate the python libraries
    if (MSVC)
        INCLUDE_DIRECTORIES($ENV{CONDA_PREFIX}/Library/include)
        LINK_DIRECTORIES($ENV{CONDA_PREFIX}/Library/lib)
        SET(python_command $ENV{CONDA_PREFIX}/python.exe)
    ELSE()
        INCLUDE_DIRECTORIES($ENV{CONDA_PREFIX}/include)
        LINK_DIRECTORIES($ENV{CONDA_PREFIX}/lib)
        SET(python_command $ENV{CONDA_PREFIX}/bin/python)
    ENDIF()
ELSE()
    MESSAGE(STATUS "No conda environment defined")
ENDIF()
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

# Windows
##############
if (MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")
    ADD_DEFINITIONS(-DMS_WIN64)
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
    # Requires AVX
    # https://devblogs.microsoft.com/cppblog/simd-extension-to-c-openmp-in-visual-studio/
    # /Oi is for intrinsics
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 -openmp:experimental /Oi /O2")
endif (MSVC)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()
if(LINUX)
    # Requires AVX
    STRING(APPEND CMAKE_CXX_FLAGS " -march=native -Ofast -mavx2 -ffast-math -funroll-loops")
    STRING(APPEND CMAKE_C_FLAGS " -march=native -Ofast -mavx2 -ffast-math -funroll-loops")
endif()

# MACOSX
##############
FIND_PACKAGE(Threads)
if (APPLE)
    # explicit link to libc++ for recent osx versions
    MESSAGE(STATUS "Explicitly link to libc++ on modern osx versions")
    ADD_DEFINITIONS(-stdlib=libc++)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
    # Don't set MACOSX_RPATH by default
    SET(CMAKE_MACOSX_RPATH 0)
    LINK_DIRECTORIES(/usr/local/lib)
    INCLUDE_DIRECTORIES(/usr/local/include)
    # Requires AVX
    # FMA and clang: -ffp-contract=fast
    STRING(APPEND CMAKE_CXX_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops -ffp-contract=fast")
    STRING(APPEND CMAKE_C_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops -ffp-contract=fast")
endif (APPLE)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
OPTION(BUILD_PYTHON_INTERFACE "Build Python interface" ON)

################################################################################
# Source groups
################################################################################
FILE(GLOB SRC_files "${CMAKE_SOURCE_DIR}/src/*.cpp")


# Debug build detection
###########################
IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE("Release build")
    ADD_DEFINITIONS("-DVERBOSE_FIT2X=0")
ELSE (CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("Debug build")
    ADD_DEFINITIONS("-DVERBOSE_FIT2X=0")
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")


## TTTRLIB
LINK_LIBRARIES(tttrlib)

## HDF5 -
# Necessary to fulfil implicit of tttrlib
###########################
FIND_PACKAGE(HDF5 1.10 REQUIRED COMPONENTS C)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
LINK_LIBRARIES(${HDF5_LIBRARIES})
# needed for build on Windows
ADD_DEFINITIONS(-DH5_BUILT_AS_DYNAMIC_LIB)

### Intel MKL
###############
#find_package(MKL REQUIRED)
#if (MKL_FOUND)
#    INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
#    LINK_LIBRARIES(${MKL_CORE_LIBRARY})
#else ()
#    message(WARNING "MKL libs not found")
#endif ()

## Boost
###########################
if (WIN32)
    # use static Boost in Windows
    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_STATIC ON)
    set(Boost_USE_MULTITHREAD)
endif (WIN32)
FIND_PACKAGE(Boost 1.36 REQUIRED COMPONENTS date_time filesystem iostreams)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES})
MESSAGE(${Boost_INCLUDE_DIRS})


## OpenMP
###########################
if (APPLE)
    # dirty hack to make OpenMP work see:
    # https://gitlab.kitware.com/cmake/cmake/issues/18098
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
endif ()
FIND_PACKAGE(OpenMP REQUIRED)
LINK_LIBRARIES(OpenMP::OpenMP_CXX)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
ADD_SUBDIRECTORY(ext)


## Shared / static library target
###########################
# cmake needs unique target name, thus place shared in front of project
ADD_LIBRARY(${PROJECT_NAME}Shared SHARED ${SRC_files})
ADD_LIBRARY(${PROJECT_NAME}Static STATIC ${SRC_files})

# Define which headers are going to be included in the install target
FILE(GLOB HEADERS_FILES "${CMAKE_SOURCE_DIR}/include/*.h")
SET_TARGET_PROPERTIES(${PROJECT_NAME}Shared PROPERTIES PUBLIC_HEADER "${HEADERS_FILES}")
SET_TARGET_PROPERTIES(${PROJECT_NAME}Static PROPERTIES PUBLIC_HEADER "${HEADERS_FILES}")

# By default the output name is the target. To not have different names, i.e.,
# fit2xShard.dll and fit2xStatic.lib the output names are set manually.
SET_TARGET_PROPERTIES(${PROJECT_NAME}Shared PROPERTIES LIBRARY_OUTPUT_NAME "${PROJECT_NAME}")
SET_TARGET_PROPERTIES(${PROJECT_NAME}Shared PROPERTIES OUTPUT_NAME "${PROJECT_NAME}")
SET_TARGET_PROPERTIES(${PROJECT_NAME}Static PROPERTIES LIBRARY_OUTPUT_NAME "${PROJECT_NAME}")
SET_TARGET_PROPERTIES(${PROJECT_NAME}Static PROPERTIES OUTPUT_NAME "${PROJECT_NAME}")

# Add install target
#######################
INSTALL(TARGETS ${PROJECT_NAME}Shared
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME}
)
INSTALL(TARGETS ${PROJECT_NAME}Static
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME}
        )


