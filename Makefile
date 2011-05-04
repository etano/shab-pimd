CC = g++
CFLAGS = -c -I /opt/local/include -O1 -msse2
LDFLAGS = -I /opt/local/include -O1 -msse2
LIBFLAGS = -larmadillo -framework Accelerate
SOURCES = shab-pimd.cpp
SOURCES += PathsClass.cpp
SOURCES += Observables.cpp
SOURCES += RNG.cpp
SOURCES += Stats.cpp
SOURCES += Staging.cpp
SOURCES += NormalMode.cpp
SOURCES += NoseHoover.cpp
SOURCES += Langevin.cpp
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
	rm -rf ./rng/*.o *.o $(EXECUTABLE)
	
scrubData:
	rm -rf ./data/figures/*.png ./data/traces/*.dat ./data/output/*.out
