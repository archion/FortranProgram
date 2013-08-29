out=twogap.out
FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCC = ifort
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ -openmp $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< -openmp
clean:
	rm -f *.out *.mod
