# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.10.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.10.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jacob/Desktop/UO/SciVis/differencer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jacob/Desktop/UO/SciVis/differencer

# Include any dependencies generated for this target.
include CMakeFiles/differencer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/differencer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/differencer.dir/flags.make

CMakeFiles/differencer.dir/differencer.cxx.o: CMakeFiles/differencer.dir/flags.make
CMakeFiles/differencer.dir/differencer.cxx.o: differencer.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/jacob/Desktop/UO/SciVis/differencer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/differencer.dir/differencer.cxx.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/differencer.dir/differencer.cxx.o -c /Users/jacob/Desktop/UO/SciVis/differencer/differencer.cxx

CMakeFiles/differencer.dir/differencer.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/differencer.dir/differencer.cxx.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jacob/Desktop/UO/SciVis/differencer/differencer.cxx > CMakeFiles/differencer.dir/differencer.cxx.i

CMakeFiles/differencer.dir/differencer.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/differencer.dir/differencer.cxx.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jacob/Desktop/UO/SciVis/differencer/differencer.cxx -o CMakeFiles/differencer.dir/differencer.cxx.s

CMakeFiles/differencer.dir/differencer.cxx.o.requires:

.PHONY : CMakeFiles/differencer.dir/differencer.cxx.o.requires

CMakeFiles/differencer.dir/differencer.cxx.o.provides: CMakeFiles/differencer.dir/differencer.cxx.o.requires
	$(MAKE) -f CMakeFiles/differencer.dir/build.make CMakeFiles/differencer.dir/differencer.cxx.o.provides.build
.PHONY : CMakeFiles/differencer.dir/differencer.cxx.o.provides

CMakeFiles/differencer.dir/differencer.cxx.o.provides.build: CMakeFiles/differencer.dir/differencer.cxx.o


# Object files for target differencer
differencer_OBJECTS = \
"CMakeFiles/differencer.dir/differencer.cxx.o"

# External object files for target differencer
differencer_EXTERNAL_OBJECTS =

differencer: CMakeFiles/differencer.dir/differencer.cxx.o
differencer: CMakeFiles/differencer.dir/build.make
differencer: /usr/lib/libz.dylib
differencer: /usr/lib/libexpat.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkDomainsChemistryOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersFlowPaths-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersGeneric-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersHyperTree-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersParallelImaging-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersPoints-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersProgrammable-8.1.1.dylib
differencer: /System/Library/Frameworks/Python.framework/Versions/2.7/Python
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkWrappingTools-8.1.a
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersPython-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersSMP-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersSelection-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersTexture-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersTopology-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersVerdict-8.1.1.dylib
differencer: /usr/local/lib/libjpeg.dylib
differencer: /usr/local/lib/libpng.dylib
differencer: /usr/local/lib/libtiff.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkGeovisCore-8.1.1.dylib
differencer: /usr/local/lib/libhdf5.dylib
differencer: /usr/local/lib/libsz.dylib
differencer: /usr/lib/libdl.dylib
differencer: /usr/lib/libm.dylib
differencer: /usr/local/lib/libhdf5_hl.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOAMR-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOEnSight-8.1.1.dylib
differencer: /usr/local/lib/libnetcdf.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOExodus-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOExportOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOImport-8.1.1.dylib
differencer: /usr/lib/libxml2.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOInfovis-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOLSDyna-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOMINC-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOMovie-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOPLY-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOParallel-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOParallelXML-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOSQL-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOTecplotTable-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOVideo-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingMorphological-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingStatistics-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingStencil-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInfovisBoostGraphAlgorithms-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInteractionImage-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingContextOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingFreeTypeFontConfig-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingImage-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingLOD-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingVolumeOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkViewsContext2D-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkViewsInfovis-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkDomainsChemistry-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkWrappingPython27Core-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkverdict-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkproj4-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersAMR-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOExport-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingGL2PSOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkgl2ps-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtklibharu-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkoggtheora-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersParallel-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkexoIIc-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOGeometry-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIONetCDF-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtknetcdfcpp-8.1.1.dylib
differencer: /usr/local/lib/libnetcdf.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkjsoncpp-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkParallelCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOLegacy-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtksqlite-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingOpenGL2-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkglew-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingMath-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkChartsCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingContext2D-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersImaging-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInfovisLayout-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInfovisCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkViewsCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInteractionWidgets-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersHybrid-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingGeneral-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingSources-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersModeling-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingHybrid-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOImage-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkDICOMParser-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkmetaio-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkInteractionStyle-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersExtraction-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersStatistics-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingFourier-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkalglib-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingAnnotation-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingColor-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingVolume-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkImagingCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOXML-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOXMLParser-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkIOCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtklz4-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingLabel-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingFreeType-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkRenderingCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonColor-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersGeometry-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersSources-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersGeneral-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonComputationalGeometry-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkFiltersCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonExecutionModel-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonDataModel-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonMisc-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonSystem-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtksys-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonTransforms-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonMath-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkCommonCore-8.1.1.dylib
differencer: /usr/local/Cellar/vtk/8.1.0/lib/libvtkfreetype-8.1.1.dylib
differencer: /usr/lib/libz.dylib
differencer: CMakeFiles/differencer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/jacob/Desktop/UO/SciVis/differencer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable differencer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/differencer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/differencer.dir/build: differencer

.PHONY : CMakeFiles/differencer.dir/build

CMakeFiles/differencer.dir/requires: CMakeFiles/differencer.dir/differencer.cxx.o.requires

.PHONY : CMakeFiles/differencer.dir/requires

CMakeFiles/differencer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/differencer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/differencer.dir/clean

CMakeFiles/differencer.dir/depend:
	cd /Users/jacob/Desktop/UO/SciVis/differencer && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jacob/Desktop/UO/SciVis/differencer /Users/jacob/Desktop/UO/SciVis/differencer /Users/jacob/Desktop/UO/SciVis/differencer /Users/jacob/Desktop/UO/SciVis/differencer /Users/jacob/Desktop/UO/SciVis/differencer/CMakeFiles/differencer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/differencer.dir/depend

