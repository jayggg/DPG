
VPATH = ./misc:./spaces:./integrators
objects = dpgintegrators.o getcomp.o fluxerr.o enorms.o vertexschwarz.o l2trace.o hcurlintegrators.o periodichcurl.o periodich1.o

%.o : %.cpp
	ngscxx -I. -c $< -o $@
#       ngscxx -DDEBUG -I. -c $? -o $@

libDPG.so : $(objects)
	ngsld -shared $(objects) -lngsolve -lngfem -lngcomp -o $@

clean:
	rm *.o libDPG.so

all	: libDPG.so
