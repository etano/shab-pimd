CC = g++
CFLAGS = -larmadillo -I /usr/include/boost -c -O1 -msse2
LDFLAGS = -larmadillo -I /usr/include/boost -O1 -msse2
SOURCES = shab-pimd.cpp
SOURCES += PathsClass.cpp
SOURCES += Observables.cpp
SOURCES += RNG.cpp
SOURCES += Stats.cpp
SOURCES += Staging.cpp
SOURCES += NoseHoover.cpp
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
	rm -rf rng/*.o *.o $(EXECUTABLE)
	

