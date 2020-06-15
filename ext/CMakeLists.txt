cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 OLD)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
INCLUDE_DIRECTORIES(.)

SET(MODULE_NAME fit2x)

SET_SOURCE_FILES_PROPERTIES(python/fit2x.i PROPERTIES CPLUSPLUS ON)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")


IF(BUILD_PYTHON_INTERFACE)
    FIND_PACKAGE(Python)
    MESSAGE("Building Python interface")
    INCLUDE_DIRECTORIES(${Python_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${Python_NumPy_PATH})
    SWIG_ADD_LIBRARY(
        ${MODULE_NAME}
        LANGUAGE python
        SOURCES python/fit2x.i ${SRC_files}
    )

    EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('LDSHARED'))"
            RESULT_VARIABLE PYTHON_CVPY_PROCESS
            OUTPUT_VARIABLE PYTHON_LDSHARED
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    IF ("${PYTHON_LDSHARED}" MATCHES "dynamic_lookup")
        MESSAGE("Using static linking for Python")
        SET(PYTHON_LINK_LIBRARIES "NO")
    ELSE()
        SET(PYTHON_LINK_LIBRARIES "${PYTHON_LIBRARIES}")
    ENDIF()

    TARGET_LINK_LIBRARIES(
            ${MODULE_NAME}
            ${PYTHON_LIBRARY}
    )
    IF (NOT PYTHON_LINK_LIBRARIES)
        IF (APPLE)
            SET_TARGET_PROPERTIES(${MODULE_NAME} PROPERTIES LINK_FLAGS "-Wl,-flat_namespace,-undefined,dynamic_lookup")
        ENDIF()
    ENDIF()
ENDIF()


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

