out=rpa.out
#FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCC = gfortran
OMP = -fopenmp
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMP) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMP)
clean:
	rm -f *.out *.mod
