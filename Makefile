CC = g++
CFLAGS = -c -I /usr/include/boost -O1 -msse2
LDFLAGS = -I /usr/include/boost -O1 -msse2
LIBFLAGS = -larmadillo
SOURCES = shab-pimd.cpp
SOURCES += RNG.cpp
SOURCES += Paths.cpp
SOURCES += Stats.cpp
SOURCES += rng/sfmt.cpp 
SOURCES += rng/mersenne.cpp 
SOURCES += rng/userintf.cpp
SOURCES += rng/stoc1.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = shab-pimd

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o $(EXECUTABLE)
	
realclean:
	rm -rf ./rng/*.o *.o $(EXECUTABLE)
