#out= $(patsubst %.f90,%.out,$(wildcard *.f90))
out = vmc_tJ_ph.out
FCC = ifort
OMPFLAG= -openmp
CFLAG=-warn nounused
FCCFLAG= -lmy -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG) $(CFLAG)
clean:
	rm -f *.out *.mod
