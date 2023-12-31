
# CMakeLists.txt for UCorrelator 
cmake_minimum_required(VERSION 2.8.10)
project(UCorrelator)

set(libname "UCorrelator")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


if(DEFINED ENV{PUEO_UTIL_INSTALL_DIR})
  message("PUEO_UTIL_INSTALL_DIR is set to $ENV{PUEO_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  PUEO_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTIL_SHARE_DIR}/cmake/modules/") 

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net)
find_package(FFTW REQUIRED)


link_directories(${UTIL_LIB_DIR})


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
include(${ROOT_USE_FILE})
message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")


include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src/ ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${ROOT_CXX_FLAGS}")

file(GLOB HEADER_FILES
    "src/pueo/*.h"    
)

file(GLOB_RECURSE SOURCE_FILES "src/*.cc")

set(DICTNAME G__${libname})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ZLIB_LIBRARIES} RootFftwWrapper pueoEvent AntarcticaRoot pueoAnalysis ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})


message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")
message("DICTNAME is set to ${DICTNAME}" )

install (FILES ${HEADER_FILES} DESTINATION ${UTIL_INC_DIR}/pueo)
install (TARGETS ${libname}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})
      
if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
  install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
#  install (FILES ${PROJECT_BINARY_DIR}/lib${libname}.rootmap DESTINATION ${UTIL_LIB_DIR})
endif() 

#################################################################################3
### The default build configuration is INSANE. No optimization? What is this 1971? 

message (" Adding new build type") 

set(CMAKE_CXX_FLAGS_DEFAULT 
  "-Os -g -Wall -Wextra" 
  CACHE STRING "c++ Flags used during default eventReaderRoot builds" 
  FORCE ) 

set(CMAKE_C_FLAGS_DEFAULT 
  "-Os -g -Wall"
  CACHE STRING "c Flags used during default eventReaderRoot builds" 
  FORCE ) 

set(CMAKE_EXE_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds" 
  FORCE ) 

set(CMAKE_SHARED_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds" 
  FORCE ) 


mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT) 

if (NOT CMAKE_BUILD_TYPE) 
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE ) 
endif()



if(CMAKE_COMPILER_IS_GNUCXX) 
### someone should do this for clang if they want it to be as fast as possible 
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.1)
  if(NATIVE_ARCH) #define in anitaBuildTool, default on
    add_definitions( -march=native)
  endif()
 endif()
endif() 

option (FORCE_OLD_GPP_ABI "Force old g++ ABI. this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations" OFF ) 
if (FORCE_OLD_GPP_ABI) 
  add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0) 
endif() 
 

option( FFTTOOLS_ENABLE_OPENMP "If FFTtools was compiled with OpenMP support, enable this option to allow multiple FFT's to be computed in parallel") 
if (FFTTOOLS_ENABLE_OPENMP)
  add_definitions(-DFFTTOOLS_COMPILED_WITH_OPENMP)
endif()


option ( ENABLE_OPENMP " Enable OpenMP for multicore usage") 
if(ENABLE_OPENMP) 
  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    add_definitions (-DUCORRELATOR_OPENMP) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else() 
    message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support") 
  endif() 
endif() 

install (DIRECTORY data/ DESTINATION ${UTIL_SHARE_DIR}/UCorrelator)
install (FILES cmake/modules/FindLibConfig.cmake DESTINATION ${UTIL_SHARE_DIR}/cmake/modules)

set (BINSRCDIR drivers) 

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

option ( ENABLE_LIBCONFIG "Enable reading from config file (requires libconfig)" OFF ) 


if (ENABLE_LIBCONFIG)
  find_package(LibConfig) 
  if (LIBCONFIG_FOUND) 
    message( "Found libconfig. Support of reading configuration files enabld") 
    add_definitions( -DENABLE_LIBCONFIG) 
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBCONFIG_EXE_LINKER_FLAGS}")
  else() 
    message("Could not find libonfig. Support of reading config is disabled") 
  endif() 
endif() 


macro(do_driver driver_name) 
  add_executable(${driver_name} ${BINSRCDIR}/${driver_name}.cc) 
  target_link_libraries(${driver_name}  ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} ${libname})
endmacro() 





#do_driver(doWais)
#do_driver(doWaisA4)
#do_driver(doResolutions) 
#do_driver(doSummaryResolutions)
#do_driver(doCorrelationSummaryTree) 
#do_driver(doSimulated) 
#do_driver(doDecimated) 
#do_driver(doDecimatedNoFilter)
#do_driver(doLDB) 
#do_driver(doAll) 
#do_driver(evaluateFiltersA3) 
#do_driver(doGenerateSineSubCache)
#do_driver(eventTxtSaver)
##do_driver(findEventsEv)
#do_driver(makeMinBiasAverage)


