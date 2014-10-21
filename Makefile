CC = gcc 
CFLAGS = -O2
OBJ = Lyapunov3D.o Topology.o LocateBoxROMS.o Interpolation.o Integration.o 

LIBS = -L/usr/local/netcdf/lib -lnetcdf -lm

Lyapunov3D: $(OBJ)
	$(CC) -o Lyapunov3D.out $(CFLAGS) $(OBJ) $(LIBS)
clean:
	-rm -f $(OBJ) *~

# *Individual File Dependencies*

Lyapunov3D.o: Lyapunov3D.h

Topology.o: Lyapunov3D.h

LocateBoxROMS.o: Lyapunov3D.h

Interpolation.o: Lyapunov3D.h 

Integration.o: Lyapunov3D.h 
