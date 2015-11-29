CC_MONO=g++ -Wall -mcmodel=medium

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: fsle3d

readparameters.o: readparameters.cpp readparameters.h date.h
	$(CC) -c readparameters.cpp 

vectorXYZ.o: vectorXYZ.cpp vectorXYZ.h constants.o
	$(CC) -c vectorXYZ.cpp 

constants.o: constants.cpp constants.h
	$(CC) -c constants.cpp 

gridconstruction.o: gridconstruction.cpp gridconstruction.h vectorXYZ.o
	$(CC) -c gridconstruction.cpp

integration.o: integration.cpp integration.h vectorXYZ.o constants.o
	$(CC) -c integration.cpp

vflow.o: vflow.cpp vflow.h vectorXYZ.o constants.o
	$(CC) -c vflow.cpp -lnetcdf_c++

fsle3d.o: fsle3d.cpp
	$(CC) -c fsle3d.cpp 

fsle3d: fsle3d.o readparameters.o gridconstruction.o vectorXYZ.o vflow.o integration.o
	$(CC)  fsle3d.o readparameters.o gridconstruction.o vectorXYZ.o constants.o vflow.o integration.o -o fsle3d -lnetcdf_c++
clean:
	$(RM) *.o 
