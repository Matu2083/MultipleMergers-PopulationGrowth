CC=g++
VERSION = $(shell cat ./version)
CFLAGS=-c -O3 -m64 -std=c++11 -DVERSION=\"${VERSION}\" -I/usr/local/include/ -L/usr/local/lib/
LDFLAGS=-lm -lmpfr -lgmp -lgsl -lgslcblas 
SOURCES=main.cpp myIOParser.cpp myExpectedSFS.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MMC-Phi-MaxLikelihood-LookUp-Lokal.out

all: clean_all $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f *.o
	
.PHONY: clean_all
clean_all:
	rm -f *.o $(EXECUTABLE)

