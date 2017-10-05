PROJECT( project )
cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
if (POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

#OPTION(useSDL2 "Turn SDL2 on and off" ON)
#OPTION(useFFTW "Turn FTTW on and off" ON)
OPTION(usePython2 "Turn Python2.7 on and off" ON)

MESSAGE("OPTIONS:")
#MESSAGE("useSDL2 = ${useSDL2}")
#MESSAGE("useFFTW = ${useFFTW}")
MESSAGE("usePython2 = ${usePython2}")

file(GLOB_RECURSE LIBS_SRC_PATH "./src/libs/*.cpp")
file(GLOB_RECURSE FUNCTIONS_SRC_PATH "./src/functions/*.cpp")
#file(GLOB_RECURSE CUDA_SRC_PATH "./src/functions/*.cu*")
file(GLOB_RECURSE DEBUG_SRC_PATH "./src/debug/*.cpp")

# sources files
Set (PROJECT_SRCS
    #src/projects/Smoke2D/main.cpp
    #src/projects/Smoke2D/Smoke2D.cpp
    src/projects/EllipsoidFitting/main.cpp
    #src/projects/VesicleChanVese3d/main.cpp
    #src/projects/PMPH/assignment1/cudaInterface.cu
    ${LIBS_SRC_PATH}
    ${FUNCTIONS_SRC_PATH}
 #   ${CUDA_SRC_PATH}
    ${DEBUG_SRC_PATH}
)

MESSAGE("Libs Sources: \n ${LIBS_SRC_PATH}")
MESSAGE("Functions Sources: \n ${FUNCTIONS_SRC_PATH}")
MESSAGE("Cuda Sources: \n ${CUDA_SRC_PATH}")
MESSAGE("Debug Sources: \n ${DEBUG_SRC_PATH}")

#SET(CMAKE_C_COMPILER /usr/bin/g++-4.9.3)
#SET(CMAKE_CXX_COMPILER /usr/bin/g++-4.9.3)

#FIND_PACKAGE(CUDA REQUIRED)
#list(APPEND CUDA_NVCC_FLAGS " -O3 -Xcompiler -fopenmp -std=c++11;") # -arch=sm_20; -rdc=true; -O2;-DVERBOSE;  -lcudadevrt
#SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
ADD_EXECUTABLE( project ${PROJECT_SRCS} )

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)

IF((CMAKE_BUILD_TYPE MATCHES "debug") OR (CMAKE_BUILD_TYPE MATCHES "DEBUG"))
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF((CMAKE_BUILD_TYPE MATCHES "debug") OR (CMAKE_BUILD_TYPE MATCHES "DEBUG"))

IF((CMAKE_BUILD_TYPE MATCHES "release") OR (CMAKE_BUILD_TYPE MATCHES "RELEASE"))
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF((CMAKE_BUILD_TYPE MATCHES "release") OR (CMAKE_BUILD_TYPE MATCHES "RELEASE"))

#enable or disable debugging
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DDEBUG )
ENDIF (CMAKE_BUILD_TYPE MATCHES "Debug")

IF(CMAKE_BUILD_TYPE MATCHES "Release")
    add_definitions(-DNDEBUG )
    add_definitions(-O3)

ENDIF (CMAKE_BUILD_TYPE MATCHES "Release")
set(CMAKE_CXX_FLAGS " -O3 -Wall -Wextra -std=c++11 -fopenmp -Wno-sign-compare")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

find_library(CUDA_CUT_LIBRARY libcudart /opt/cuda/lib64/)

INCLUDE(FindPkgConfig)

#find_package(VTK REQUIRED)
#if (VTK_FOUND)
#        include(${VTK_USE_FILE})
        #include_directories(${VTK_INCLUDE_DIRS})
#        message("VTK found, VTK_USE_FILE: ${VTK_USE_FILE}")
#else()
#        message(ERROR "VTK not found")
#endif()

SET(EIGEN3_INCLUDE_DIR SYSTEM /usr/include/eigen3)
#SET(SDL2_INCLUDE_DIR /usr/include/SDL2/)

#IF (useSDL2)
#    message("test: ${useSDL2}")
#    SET(SDL2_LINKER_FLAGS "-lSDL2")

#    ADD_DEFINITIONS("-DUSE_SDL")

#    PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
#ENDIF (useSDL2)

#IF (useFFTW)
#    SET(FFTW_LINKER_FLAGS "-lfftw3_threads -lfftw3")
#ENDIF (useFFTW)

IF (usePython2)
    SET(PYTHON27_LINKER_FLAGS "-lpython2.7")
    SET(PYTHON27_INCLUDE_DIR SYSTEM /usr/include/python2.7)
ENDIF (usePython2)

INCLUDE_DIRECTORIES (
    /usr/bin
    /usr/include
    /usr/local/include
    ${PYTHON27_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
#    ${SDL2_INCLUDE_DIR}
)
#-lcutil -lcudpp -lcuda -lcudart
SET(CMAKE_EXE_LINKER_FLAGS "-lm -ltiff ${FFTW_LINKER_FLAGS} ${PYTHON27_LINKER_FLAGS}") # ${SDL2_LINKER_FLAGS}

MESSAGE("Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

TARGET_LINK_LIBRARIES( project ${LINK_LIBRARIES} ) #${CUDA_LIBRARY} cuda ${CUDA_LIBRARIES} ${VTK_LIBRARIES}
