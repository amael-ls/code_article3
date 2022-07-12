
## Phony variables
.PHONY: all clean # To prevent problems if there is a file called 'clean' or 'all'

## General variables
# Compiler
CXX = g++
CXXFLAGS = -std=c++2a -Wall -O2
OUTPUT = ebt

## Directories
# Sources
src_dir_ebt = ./cpp
src_dir_alglib += ./alglib
src_dir = $(src_dir_ebt) $(src_dir_alglib)
include_dirs = $(addprefix -I, $(src_dir))

# Object (*.o)
object_dir_ebt = ./obj_ebt
object_dir_alglib = ./obj_alglib

# Libraries
libs_dir = -L/opt/homebrew/lib

## Input files (headers and scripts)
# EBT files
includes_ebt = $(wildcard *.h++ $(foreach fd, $(src_dir_ebt)/, $(fd)*.h++))
sources_ebt = $(wildcard *.c++ $(foreach fd, $(src_dir_ebt)/, $(fd)*.c++))

# ALGLIB files
includes_alglib = $(wildcard *.h $(foreach fd, $(src_dir_alglib)/, $(fd)*.h))
sources_alglib = $(wildcard *.cpp $(foreach fd, $(src_dir_alglib)/, $(fd)*.cpp))

# Libraries
libs = -ltbb

## Output compilation objects
objects_ebt = $(addprefix $(object_dir_ebt)/, $(notdir $(sources_ebt:.c++=.o)))
objects_alglib = $(addprefix $(object_dir_alglib)/, $(notdir $(sources_alglib:.cpp=.o)))
objects = $(objects_ebt) $(objects_alglib)

all: ebt

ebt: $(objects)
	$(CXX) $(CXXFLAGS) -o $@ $(objects) $(libs_dir) $(libs)

$(object_dir_ebt)/%.o: $(src_dir_ebt)/%.c++ $(includes_ebt) $(includes_alglib)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(include_dirs)

$(object_dir_alglib)/%.o: $(src_dir_alglib)/%.cpp $(includes_alglib)
	mkdir -p $(@D)
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(include_dirs)

clean:
	rm -rf $(object_dir_ebt)/*.o $(object_dir_alglib)/*.o

echoes:
	@echo "Include files for ebt: $(includes_ebt)"
	@echo "Includes files for alglib: $(includes_alglib)"
	@echo "Source files for ebt: $(sources_ebt)"
	@echo "Source files for alglib: $(sources_alglib)"
	@echo "Object files: $(objects)"
	@echo "Include directories: $(include_dirs)"
