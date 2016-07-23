# CMakeLists.txt for KDetSim
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(KDetSim)

# Set ROOT path and try to autodetect it from ENV
set(ROOTSYS $ENV{ROOTSYS} CACHE STRING "ROOT path")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${ROOTSYS}" CACHE PATH "default install path" FORCE)
endif()

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Define includes
include_directories(inc)
set(Project_Includes
    inc/K3D.h
    inc/KDetector.h
    inc/KField.h
    inc/KGeometry.h
    inc/KImplant2D.h
    inc/KImplant3D.h
    inc/KMaterial.h
    inc/KMesh.h
    inc/KPad.h
    inc/KPixel.h
    inc/KStrip.h
    inc/KStruct.h

    inc/fizika.h
    inc/nrutil.h
)

#---Define sources
set(Project_Sources
    src/K3D.cxx
    src/KDetector.cxx
    src/KField.cxx
    src/KGeometry.cxx
    src/KImplant2D.cxx
    src/KImplant3D.cxx
    src/KMaterial.cxx
    src/KMesh.cxx
    src/KMinuit.cxx
    src/KPad.cxx
    src/KPixel.cxx
    src/KStrip.cxx
    src/KStruct.cxx
    src/nrutil.cxx
)

#---ROOT requirements
ROOT_GENERATE_DICTIONARY(
    G__${PROJECT_NAME}

    ${Project_Includes}

    MODULE ${PROJECT_NAME}
    LINKDEF inc/KDetSim_LinkDef.h
)

#---Create a shared library with geneated dictionary
add_library(${PROJECT_NAME} SHARED
    ${Project_Sources}
    G__${PROJECT_NAME}.cxx
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    PUBLIC_HEADER "${Project_Includes}")

target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})

install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION include/${PROJECT_NAME}
)
install(FILES
    ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm
    ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}.rootmap
    DESTINATION lib
)
