CC=g++
VERSION = $(shell cat ./version)
CFLAGS=-c -O3 -std=c++11 -DVERSION=\"${VERSION}\" -DNDEBUG
LDFLAGS=-I/usr/local/include/ -L/usr/local/lib/ -lgsl -lgslcblas 
SOURCES=main.cpp myNode.cpp myTree.cpp myIOParser.cpp myRandomNumbers.cpp mySummaryStatistics.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MMC-Growth.out

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

