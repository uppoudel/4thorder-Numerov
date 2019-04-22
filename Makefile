FC=gfortran -c
LD=gfortran
SRC=params.f90 allocatearrays.f90 deallocatearrays.f90 schforthorder.f90
OBJ=params.o allocatearrays.o deallocatearrays.o schforthorder.o
schrodinger:
	$(FC) $(SRC)
	$(LD) $(OBJ) -o schforthorder.x -llapack -lblas
	rm -rf *.o
clean:
	rm -rf *.o *.mod *.x *.dat
        
