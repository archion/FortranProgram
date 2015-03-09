#out= $(patsubst %.f90,%.out,$(wildcard *.f90))
out = vmc_bdg.out
FCC = ifort
OMPFLAG= -openmp
CFLAG=-warn nounused
FFCFLAG= -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
#FCCFLAG= -lmy -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG) $(CFLAG)
clean:
	rm -f *.out *.mod
