# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/paulineh/TMA4280/gitLABS4280/project2/poisson

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/paulineh/TMA4280/gitLABS4280/project2/poisson/build

# Include any dependencies generated for this target.
include CMakeFiles/poisson-f.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poisson-f.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poisson-f.dir/flags.make

CMakeFiles/poisson-f.dir/poisson.f90.o: CMakeFiles/poisson-f.dir/flags.make
CMakeFiles/poisson-f.dir/poisson.f90.o: ../poisson.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/paulineh/TMA4280/gitLABS4280/project2/poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/poisson-f.dir/poisson.f90.o"
	/usr/bin/mpif90  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/paulineh/TMA4280/gitLABS4280/project2/poisson/poisson.f90 -o CMakeFiles/poisson-f.dir/poisson.f90.o

CMakeFiles/poisson-f.dir/poisson.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/poisson-f.dir/poisson.f90.i"
	/usr/bin/mpif90  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/paulineh/TMA4280/gitLABS4280/project2/poisson/poisson.f90 > CMakeFiles/poisson-f.dir/poisson.f90.i

CMakeFiles/poisson-f.dir/poisson.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/poisson-f.dir/poisson.f90.s"
	/usr/bin/mpif90  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/paulineh/TMA4280/gitLABS4280/project2/poisson/poisson.f90 -o CMakeFiles/poisson-f.dir/poisson.f90.s

CMakeFiles/poisson-f.dir/poisson.f90.o.requires:

.PHONY : CMakeFiles/poisson-f.dir/poisson.f90.o.requires

CMakeFiles/poisson-f.dir/poisson.f90.o.provides: CMakeFiles/poisson-f.dir/poisson.f90.o.requires
	$(MAKE) -f CMakeFiles/poisson-f.dir/build.make CMakeFiles/poisson-f.dir/poisson.f90.o.provides.build
.PHONY : CMakeFiles/poisson-f.dir/poisson.f90.o.provides

CMakeFiles/poisson-f.dir/poisson.f90.o.provides.build: CMakeFiles/poisson-f.dir/poisson.f90.o


# Object files for target poisson-f
poisson__f_OBJECTS = \
"CMakeFiles/poisson-f.dir/poisson.f90.o"

# External object files for target poisson-f
poisson__f_EXTERNAL_OBJECTS =

poisson-f: CMakeFiles/poisson-f.dir/poisson.f90.o
poisson-f: CMakeFiles/poisson-f.dir/build.make
poisson-f: libcommon.a
poisson-f: /usr/lib/openmpi/lib/libmpi_usempif08.so
poisson-f: /usr/lib/openmpi/lib/libmpi_usempi_ignore_tkr.so
poisson-f: /usr/lib/openmpi/lib/libmpi_mpifh.so
poisson-f: /usr/lib/openmpi/lib/libmpi.so
poisson-f: CMakeFiles/poisson-f.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/paulineh/TMA4280/gitLABS4280/project2/poisson/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable poisson-f"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poisson-f.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poisson-f.dir/build: poisson-f

.PHONY : CMakeFiles/poisson-f.dir/build

CMakeFiles/poisson-f.dir/requires: CMakeFiles/poisson-f.dir/poisson.f90.o.requires

.PHONY : CMakeFiles/poisson-f.dir/requires

CMakeFiles/poisson-f.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poisson-f.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poisson-f.dir/clean

CMakeFiles/poisson-f.dir/depend:
	cd /home/paulineh/TMA4280/gitLABS4280/project2/poisson/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/paulineh/TMA4280/gitLABS4280/project2/poisson /home/paulineh/TMA4280/gitLABS4280/project2/poisson /home/paulineh/TMA4280/gitLABS4280/project2/poisson/build /home/paulineh/TMA4280/gitLABS4280/project2/poisson/build /home/paulineh/TMA4280/gitLABS4280/project2/poisson/build/CMakeFiles/poisson-f.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/poisson-f.dir/depend

