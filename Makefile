CC = g++
CFLAGS = -I /usr/include/boost -larmadillo -c -O1 -msse2
LDFLAGS = -I /usr/include/boost -larmadillo -O1 -msse2
SOURCES = shab-pimd.cpp
SOURCES += Paths.cpp
SOURCES += RNG.cpp
SOURCES += rng/sfmt.cpp 
SOURCES += rng/mersenne.cpp 
SOURCES += rng/userintf.cpp
SOURCES += rng/stoc1.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = shab-pimd

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o $(EXECUTABLE)
	
realclean:
	rm -rf ./rng/*.o *.o $(EXECUTABLE)
