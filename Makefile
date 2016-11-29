
VPATH = ./misc:./spaces:./integrators
objects = dpgintegrators.o getcomp.o fluxerr.o enorms.o  l2trace.o hcurlintegrators.o periodichcurl.o periodich1.o  l2quadpluspace.o l2quadplusfe.o vertexschwarz.o

headers = dpgintegrators.hpp hcurlintegrators.cpp l2quadpluspace.hpp l2quadplusfe.hpp

%.o : %.cpp  $(headers)
	ngscxx -I. -c $< -o $@
#       ngscxx -DDEBUG -I. -c $? -o $@

libDPG.so : $(objects)
	ngsld -shared $(objects) -lngsolve -lngfem -lngcomp -o $@

clean:
	rm *.o libDPG.so

all	: libDPG.so
