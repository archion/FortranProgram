out=rpa_sc.out
FCC = ifort
OMPFLAG = -openmp
#FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCCFLAG= -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -lpthread -lmkl_lapack95
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.out *.mod
