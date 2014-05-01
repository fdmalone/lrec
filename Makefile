CC = g++
VPATH = $(PWD)
UTIL = headers/
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = bin/out
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS = -c -O3 -DLATTICE=$(LATTICE) -I $(UTIL) -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS = -O3 -lgsl -lgslcblas -lm -llapack -lblas
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)
clean:
	rm -rf $(EXECUTABLE) $(OBJECTS)
