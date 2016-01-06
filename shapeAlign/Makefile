# Makefile for shape alignment

SHELL := /bin/bash -e

OBJ_DIR = obj
BIN_DIR = bin

INCLUDES = -Isrc/

CXX = g++-4.9
CXXFLAGS = -Wall -g -D_FILE_OFFSET_BITS=64 -fPIC $(INCLUDES) -fopenmp
LIBS		= -lz -lgsl -lgslcblas

OBJECTS= shapeAlign.o shapeAlignMain.o
SOURCES= shapeAlign.cpp shapeAlign.h ShapeAlignMain.cpp
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

all: $(OBJ_DIR) $(BIN_DIR) $(BUILT_OBJECTS)
	@echo "Building shapeAlign binary."
	@$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/shapeAlign $(OBJ_DIR)/*.o $(LDFLAGS) $(LIBS)
	@echo "Done."
.PHONY: all

$(OBJ_DIR) $(BIN_DIR):
	@mkdir -p $@
	
$(BUILT_OBJECTS): $(FULL_SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) $(CXXFLAGS) -c src/$(*F).cpp -o $@ $(LDFLAGS)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean