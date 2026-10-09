cmake_minimum_required(VERSION 3.23.0)
project(UCorrelator VERSION 0.0.1)
add_library(${PROJECT_NAME} SHARED)
## ------------- Dependency -------------

find_package(pueoAnalysisFramework CONFIG REQUIRED) 
find_package(OpenMP)
## todo, OpenMP can probably also come from pueoAnalysisFramework but not sure how much that simplifies things.
#  in the meantime, references:
#    1) https://cmake.org/cmake/help/latest/module/FindOpenMP.html#result-variables
#    2) https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE UCORRELATOR_OPENMP)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found, compiling without it.")
endif()

target_sources(${PROJECT_NAME} PRIVATE
  src/AnalysisConfig.cc
  src/Analyzer.cc
  src/AntennaPositions.cc
  src/Baseline.cc
  src/CachedFC.cc
  src/Correlator.cc
  src/EASFitter.cc
  src/Flags.cc
  src/ImageTools.cc
  src/NoiseMachine.cc
  src/PeakFinder.cc
  src/PointingResolutionModel.cc
  src/PolarityMachine.cc
  src/ProbabilityMap.cc
  src/ShapeParameters.cc
  src/SineSubtractCache.cc
  src/SpectrumParameters.cc
  src/TimeDependentAverage.cc
  src/UCFilters.cc
  src/UCKDE.cc
  src/UCMapStatTools.cc
  src/UCorrelatorGUI.cc
  src/Util.cc
  src/WaveformCombiner.cc
)
set(HEADER_FILES
  src/pueo/AnalysisConfig.h
  src/pueo/Analyzer.h
  src/pueo/AntennaPositions.h
  src/pueo/Baseline.h
  src/pueo/CachedFC.h
  src/pueo/Correlator.h
  src/pueo/DeltaT.h
  src/pueo/EASFitter.h
  src/pueo/NoiseMachine.h
  src/pueo/PeakFinder.h
  src/pueo/PointingResolutionModel.h
  src/pueo/PolarityMachine.h
  src/pueo/ProbabilityMap.h
  src/pueo/ShapeParameters.h
  src/pueo/SineSubtractCache.h
  src/pueo/SpectrumParameters.h
  src/pueo/TimeDependentAverage.h
  src/pueo/TrigCache.h
  src/pueo/UCFilters.h
  src/pueo/UCFlags.h
  src/pueo/UCImageTools.h
  src/pueo/UCKDE.h
  src/pueo/UCMapStatTools.h
  src/pueo/UCUtil.h
  src/pueo/UCorrelatorGUI.h
  src/pueo/WaveformCombiner.h
)
target_sources(${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  BASE_DIRS src          # everything in BASE_DIRSc is visible to cmake during the build process,
  FILES ${HEADER_FILES}  # but only FILES will be installed
)

target_link_libraries(${PROJECT_NAME} PUBLIC PUEO::pueoAnalysisFramework)

# CERN ROOT dictionary generation and installation
set(DICTNAME G__${PROJECT_NAME})
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION lib
)

## ------------- Installation -------------

include(GNUInstallDirs)
# the config file has instruction for downstream projects on how to fetch dependencies of this project
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

install(
  TARGETS ${PROJECT_NAME} # just one target, ie this library, for now
  EXPORT ${PROJECT_NAME}Targets # put all targets listed above with this particular "export"
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # install the header files of the target
  LIBRARY DESTINATION lib
)

install(
  EXPORT ${PROJECT_NAME}Targets # auto create and install pueoAnalysisFrameworkTargets.cmake directly to install tree
  NAMESPACE PUEO::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

install(DIRECTORY data/ DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})

## ------------- Options -------------

# Change Log Sep. 25, 2026
# The option FFTTOOLS_ENABLE_OPENMP removed -- this should be derived automatically from RootFftwWrapper

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelwithDebugInfo CACHE STRING "" FORCE)
  # use warning flags in RelWithDebInfo mode, which is set to be the default mode
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:RelWithDebInfo>:-Wall -Wextra -Wno-unused-function>)
endif()

# option (ENABLE_LIBCONFIG "Enable reading from config file (requires libconfig)" OFF)
# if (ENABLE_LIBCONFIG)
#   find_package(LibConfig CONFIG REQUIRED)
#   message( "Found libconfig. Support of reading configuration files enabld") 
#   target_compile_definitions(${PROJECT_NAME} PRIVATE ENABLE_LIBCONFIG) 
# endif() 

