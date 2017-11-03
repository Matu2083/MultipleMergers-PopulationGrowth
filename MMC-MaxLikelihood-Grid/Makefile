CC=g++
VERSION = $(shell cat ./version)
CFLAGS=-c -O3 -m64 -std=c++11 -DVERSION=\"${VERSION}\"
LDFLAGS=
SOURCES=main.cpp myIOParser.cpp myEstimationMethod.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MMC-MaxLikelihood-Grid.out

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

