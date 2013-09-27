out=twogap_ddw.out
#FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
OMPFLAG= -fopenmp
FCC = gfortran
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.out *.mod
