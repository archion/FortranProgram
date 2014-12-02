out=lanczos.out
#FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCCFLAG= -lmy -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lmkl_lapack95
FCC = ifort
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ -openmp $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< -openmp
clean:
	rm -f *.out *.mod
