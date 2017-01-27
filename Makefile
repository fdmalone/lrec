CC = g++
VPATH = $(PWD)
UTIL = headers/
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = bin/lrec.x
GIT_VERSION := $(shell git describe --dirty --always)
CFLAGS = -c -O3 -DLATTICE=$(LATTICE) -I$(UTIL) -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS = -lgsl -lgslcblas -lm -llapack -lblas -larmadillo

BIN_DIR = bin

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(BIN_DIR)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf $(EXECUTABLE) $(OBJECTS)

$(BIN_DIR):
	mkdir -p $@
