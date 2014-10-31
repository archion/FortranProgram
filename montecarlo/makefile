out=vmc_tJ.out
FCC = ifort
OMPFLAG= -openmp
FCCFLAG= -lmy -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.out *.mod
