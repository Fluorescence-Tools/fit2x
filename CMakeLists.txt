cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)
SET(PROJECT_NAME fit2x)
PROJECT(${PROJECT_NAME})
ENABLE_LANGUAGE(CXX)
SET(CMAKE_CXX_STANDARD 11)
ADD_DEFINITIONS(-D_LIBCPP_ENABLE_CXX17_REMOVED_FEATURES)

################################################################################
###### CMAKE
################################################################################
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
MESSAGE(STATUS "CMAKE_MODULE_PATH='${CMAKE_MODULE_PATH}'")

################################################################################
###### DEFINITIONS
################################################################################
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

################################################################################
# BUILD OPTIONS
################################################################################
OPTION(INSTALL "Install the library and the include files" ON)
OPTION(BUILD_PYTHON_INTERFACE "Build Python interface" ON)
OPTION(BUILD_DOC "Build documentation" OFF)
OPTION(BUILD_R_INTERFACE "Build R interface" OFF)

################################################################################
###### Compiler flags
################################################################################

## Compiler warnings
########################
IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(MSVC)
        ADD_COMPILE_OPTIONS(/Wall)
    else()
        ADD_COMPILE_OPTIONS(-Wall -Wextra -pedantic)
    endif()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")

if(LINUX)
    STRING(APPEND CMAKE_CXX_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops")
    STRING(APPEND CMAKE_C_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops")
endif()
if(APPLE)
    FIND_PACKAGE(Threads)
    # explicit link to libc++ for recent osx versions
    MESSAGE(STATUS "Explicitly link to libc++ on modern osx versions")
    ADD_DEFINITIONS(-stdlib=libc++)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
    # Requires AVX and FMA
    STRING(APPEND CMAKE_CXX_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops")
    STRING(APPEND CMAKE_C_FLAGS " -march=native -Ofast -mavx -ffast-math -funroll-loops")
endif()
if (MSVC)
    ADD_DEFINITIONS(-DMS_WIN64)
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
    # https://devblogs.microsoft.com/cppblog/simd-extension-to-c-openmp-in-visual-studio/
    # /Oi is for intrinsics
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX -openmp:experimental /Oi /O2")
endif (MSVC)

# MACOSX
##############
if (APPLE)
    # Don't set MACOSX_RPATH by default
    # https://github.com/conda/conda/issues/3624
    if(DEFINED ENV{CONDA_PREFIX})
        set(MACOSX_RPATH OFF)
        set(CMAKE_SKIP_RPATH TRUE)
        set(CMAKE_SKIP_INSTALL_RPATH TRUE)
    endif()
endif (APPLE)

# Debug build detection
###########################
IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE("Release build")
    IF(MSVC)
        ADD_COMPILE_OPTIONS(/w)
    ELSE()
        ADD_COMPILE_OPTIONS(-w)
    ENDIF()
ELSE (CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("Debug build")
    ADD_DEFINITIONS("-DVERBOSE_FIT2X=0")
    ## Enable Compiler warnings
    #############################
    IF(MSVC)
        ADD_COMPILE_OPTIONS(/W4)
    ELSE()
        ADD_COMPILE_OPTIONS(-Wall -Wextra -pedantic)
    ENDIF()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")


################################################################################
###### INCLUDES AND LIBRARIES
################################################################################

## CONDA
###########################
FIND_PACKAGE(Conda)

## TTTRLIB
###########################
LINK_LIBRARIES(tttrlib)

## Boost
###########################
FIND_PACKAGE(CondaBoost)
FIND_PACKAGE(Boost 1.36 REQUIRED COMPONENTS filesystem iostreams locale)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_LIBRARIES(Boost::filesystem Boost::iostreams Boost::locale)

## HDF5
# Implicit dependency tttrlib -> HDF5
###########################
FIND_PACKAGE(HDF5 1.10 REQUIRED COMPONENTS C)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
LINK_LIBRARIES(${HDF5_LIBRARIES})
# needed for build on Windows
ADD_DEFINITIONS(-DH5_BUILT_AS_DYNAMIC_LIB)

## OpenMP
###########################
if (APPLE)
    # When using conda use CONDA_PREFIX
    if(DEFINED ENV{COND_PREFIX})
        SET(OMP_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
        SET(OMP_LIB_DIR "$ENV{CONDA_PREFIX}/lib")
    else()
        SET(OMP_INCLUDE_DIR "/usr/local/opt/libomp/include")
        SET(OMP_LIB_DIR "/usr/local/opt/libomp/lib")
    endif()
    # dirty hack to make OpenMP work see:
    # https://gitlab.kitware.com/cmake/cmake/issues/18098
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OMP_INCLUDE_DIR}" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${OMP_INCLUDE_DIR}" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_omp_LIBRARY "${OMP_LIB_DIR}/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
endif ()
FIND_PACKAGE(OpenMP REQUIRED)
LINK_LIBRARIES(OpenMP::OpenMP_CXX)

################################################################################
# Source groups
################################################################################
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
FILE(GLOB SRC_files "${CMAKE_SOURCE_DIR}/src/*.cpp")
ADD_SUBDIRECTORY(ext)

################################################################################
###### Targets
################################################################################

## Shared / static library target
###########################
# cmake needs unique target name, thus place shared in front of project
ADD_LIBRARY(${PROJECT_NAME}Shared SHARED ${SRC_files} include/DecayCurve.h)
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


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    get_property(dirs DIRECTORY . PROPERTY INCLUDE_DIRECTORIES)
    foreach(dir ${dirs})
        message(STATUS "INCLUDE_DIRECTORY='${dir}'")
    endforeach()
    get_property(dirs DIRECTORY . PROPERTY LINK_LIBRARIES)
    foreach(dir ${dirs})
        message(STATUS "LINK_LIBRARIES='${dir}'")
    endforeach()
    get_property(dirs DIRECTORY . PROPERTY LINK_DIRECTORIES)
    foreach(dir ${dirs})
        message(STATUS "LINK_DIRECTORIES='${dir}'")
    endforeach()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")


