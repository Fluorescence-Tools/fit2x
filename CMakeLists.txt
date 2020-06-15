cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

project(fit2x C CXX)
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

