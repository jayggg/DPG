
objects = dpgintegrators.o getcomp.o fluxerr.o enorms.o vertexschwarz.o l2trace.o

uname = $(shell uname -s)
flagshared = 

ifeq ($(uname),Darwin)
	flagshared = -bundle -flat_namespace -undefined suppress 
endif
ifeq ($(uname),Linux)
	flagshared = -shared
endif


%.o : %.cpp 
	ngscxx  -I. -c $? -o $@
#       ngscxx -DDEBUG -I. -c $? -o $@

libDPG.so : $(objects)
	ngscxx $(flagshared) $(objects) -lngsolve -lngfem -o $@

clean:
	rm *.o libDPG.so

all	: libDPG.so
