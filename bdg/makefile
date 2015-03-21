out=bdg.out
FCC = ifort
#OMPFLAG= -openmp
#FCCFLAG= -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95 -lmy
FCCFLAG= -lmy -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95_lp64
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.out *.mod
