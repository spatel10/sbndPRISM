# Define the subproject name
project(sbnPRISM)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_sbnPRISM_sbnPRISMTools SHARED sbnPRISMTools.cxx)
add_library(sbnanalysis_sbnPRISM_sbnPRISMSelection SHARED sbnPRISMSelection.cxx)

# Define library linking
target_link_libraries(sbnanalysis_sbnPRISM_sbnPRISMSelection sbnanalysis_sbnPRISM_sbnPRISMTools)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_sbnPRISM_sbnPRISMTools DESTINATION lib)
install(TARGETS sbnanalysis_sbnPRISM_sbnPRISMSelection DESTINATION lib)

