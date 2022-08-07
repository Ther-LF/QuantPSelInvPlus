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
include examples/CMakeFiles/driver_ksdft.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/CMakeFiles/driver_ksdft.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/driver_ksdft.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/driver_ksdft.dir/flags.make

examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o: examples/CMakeFiles/driver_ksdft.dir/flags.make
examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o: /software/pexsi_v2.0.0/examples/driver_ksdft.c
examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o: examples/CMakeFiles/driver_ksdft.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/software/pexsi_v2.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o -MF CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o.d -o CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o -c /software/pexsi_v2.0.0/examples/driver_ksdft.c

examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/driver_ksdft.dir/driver_ksdft.c.i"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /software/pexsi_v2.0.0/examples/driver_ksdft.c > CMakeFiles/driver_ksdft.dir/driver_ksdft.c.i

examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/driver_ksdft.dir/driver_ksdft.c.s"
	cd /software/pexsi_v2.0.0/build/examples && /usr/bin/gcc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /software/pexsi_v2.0.0/examples/driver_ksdft.c -o CMakeFiles/driver_ksdft.dir/driver_ksdft.c.s

# Object files for target driver_ksdft
driver_ksdft_OBJECTS = \
"CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o"

# External object files for target driver_ksdft
driver_ksdft_EXTERNAL_OBJECTS =

examples/driver_ksdft: examples/CMakeFiles/driver_ksdft.dir/driver_ksdft.c.o
examples/driver_ksdft: examples/CMakeFiles/driver_ksdft.dir/build.make
examples/driver_ksdft: src/libpexsi.a
examples/driver_ksdft: /usr/lib/x86_64-linux-gnu/libmpichcxx.so
examples/driver_ksdft: /usr/lib/x86_64-linux-gnu/libmpichfort.so
examples/driver_ksdft: /software/superlu_dist-7.2.0/lib/libsuperlu_dist.a
examples/driver_ksdft: /usr/lib/x86_64-linux-gnu/libmpich.so
examples/driver_ksdft: /software/parmetis-4.0.3/lib/libparmetis.a
examples/driver_ksdft: /software/parmetis-4.0.3/lib/libmetis.a
examples/driver_ksdft: /usr/lib/x86_64-linux-gnu/liblapack.so
examples/driver_ksdft: /usr/lib/x86_64-linux-gnu/libblas.so
examples/driver_ksdft: examples/CMakeFiles/driver_ksdft.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/software/pexsi_v2.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable driver_ksdft"
	cd /software/pexsi_v2.0.0/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/driver_ksdft.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/driver_ksdft.dir/build: examples/driver_ksdft
.PHONY : examples/CMakeFiles/driver_ksdft.dir/build

examples/CMakeFiles/driver_ksdft.dir/clean:
	cd /software/pexsi_v2.0.0/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/driver_ksdft.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/driver_ksdft.dir/clean

examples/CMakeFiles/driver_ksdft.dir/depend:
	cd /software/pexsi_v2.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /software/pexsi_v2.0.0 /software/pexsi_v2.0.0/examples /software/pexsi_v2.0.0/build /software/pexsi_v2.0.0/build/examples /software/pexsi_v2.0.0/build/examples/CMakeFiles/driver_ksdft.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/driver_ksdft.dir/depend

