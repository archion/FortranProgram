FCC = ifort
#OMPFLAG= -openmp
OMPFLAG= -openmp
#CFLAG=-warn nounused
CFLAG=-warn nounused -traceback
#CFLAG=-warn nounused -check bounds -g -check all -fpe0 -traceback -debug extended
FCCFLAG= -lmy -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95_lp64
#FCCFLAG= -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
#FCCFLAG= -lmy -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95_lp64
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG) $(CFLAG)
clean:
	rm -f *.out *.mod
