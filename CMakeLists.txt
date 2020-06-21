cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

SET(PROJECT_NAME fit2x)
PROJECT(${PROJECT_NAME})
ENABLE_LANGUAGE(CXX)
SET(CMAKE_CXX_STANDARD 11)

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
OPTION(BUILD_PYTHON_INTERFACE "Build Python interface" ON)

################################################################################
# Source groups
################################################################################
FILE(
    GLOB
    SRC_files
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
)

# Optimize compiler flags
###########################
INCLUDE(CMakeHostSystemInformationExtra)
INCLUDE(GetCPUSIMDFeatures)
INCLUDE(CMakeCompilerMachineOption)
# SET(TARGET_ARCHITECTURE skylake)
# CMAKE_COMPILER_MACHINE_OPTION(ARCHITECTURE_FLAG ${TARGET_ARCHITECTURE})
MESSAGE(STATUS "Use compiler option: ${ARCHITECTURE_FLAG}")
ADD_COMPILE_OPTIONS(${ARCHITECTURE_FLAG})


# Debug build detection
###########################
IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE("Release build")
    ADD_DEFINITIONS("-DVERBOSE=0")
ELSE (CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE("Debug build")
    ADD_DEFINITIONS("-DVERBOSE=1")
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")


## Intel MKL
##############
find_package(MKL REQUIRED)
if (MKL_FOUND)
    INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
    LINK_LIBRARIES(${MKL_CORE_LIBRARY})
else ()
    message(WARNING "MKL libs not found")
endif ()


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ap)
ADD_SUBDIRECTORY(ext)

# Add install target
##############
INSTALL(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        PUBLIC_HEADER DESTINATION include
        )

