out=twogap_tJ.out
FCCFLAG= -lmkl_intel -lmkl_sequential -lmkl_core -lpthread
OMPFLAG= -openmp
FCC = ifort
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(OMPFLAG) $(FCCFLAG)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.out *.mod
