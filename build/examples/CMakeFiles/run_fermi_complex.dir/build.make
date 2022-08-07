# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /software/pexsi_v2.0.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /software/pexsi_v2.0.0/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/run_fermi_complex.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/CMakeFiles/run_fermi_complex.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/run_fermi_complex.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/run_fermi_complex.dir/flags.make

examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o: examples/CMakeFiles/run_fermi_complex.dir/flags.make
examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o: /software/pexsi_v2.0.0/examples/run_fermi_complex.cpp
examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o: examples/CMakeFiles/run_fermi_complex.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/software/pexsi_v2.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o -MF CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o.d -o CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o -c /software/pexsi_v2.0.0/examples/run_fermi_complex.cpp

examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.i"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /software/pexsi_v2.0.0/examples/run_fermi_complex.cpp > CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.i

examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.s"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /software/pexsi_v2.0.0/examples/run_fermi_complex.cpp -o CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.s

# Object files for target run_fermi_complex
run_fermi_complex_OBJECTS = \
"CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o"

# External object files for target run_fermi_complex
run_fermi_complex_EXTERNAL_OBJECTS =

examples/run_fermi_complex: examples/CMakeFiles/run_fermi_complex.dir/run_fermi_complex.cpp.o
examples/run_fermi_complex: examples/CMakeFiles/run_fermi_complex.dir/build.make
examples/run_fermi_complex: src/libpexsi.a
examples/run_fermi_complex: /usr/lib/x86_64-linux-gnu/libmpichcxx.so
examples/run_fermi_complex: /usr/lib/x86_64-linux-gnu/libmpichfort.so
examples/run_fermi_complex: /software/superlu_dist-7.2.0/lib/libsuperlu_dist.a
examples/run_fermi_complex: /usr/lib/x86_64-linux-gnu/libmpich.so
examples/run_fermi_complex: /software/parmetis-4.0.3/lib/libparmetis.a
examples/run_fermi_complex: /software/parmetis-4.0.3/lib/libmetis.a
examples/run_fermi_complex: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/run_fermi_complex: /usr/lib/x86_64-linux-gnu/libblas.so
examples/run_fermi_complex: examples/CMakeFiles/run_fermi_complex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/software/pexsi_v2.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable run_fermi_complex"
	cd /software/pexsi_v2.0.0/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_fermi_complex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/run_fermi_complex.dir/build: examples/run_fermi_complex
.PHONY : examples/CMakeFiles/run_fermi_complex.dir/build

examples/CMakeFiles/run_fermi_complex.dir/clean:
	cd /software/pexsi_v2.0.0/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/run_fermi_complex.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/run_fermi_complex.dir/clean

examples/CMakeFiles/run_fermi_complex.dir/depend:
	cd /software/pexsi_v2.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /software/pexsi_v2.0.0 /software/pexsi_v2.0.0/examples /software/pexsi_v2.0.0/build /software/pexsi_v2.0.0/build/examples /software/pexsi_v2.0.0/build/examples/CMakeFiles/run_fermi_complex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/run_fermi_complex.dir/depend

