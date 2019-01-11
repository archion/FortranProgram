out=libmy.a
source90=$(filter-out $(wildcard lattice*.f90 hamilton*.f90), $(wildcard *.f90)) lattice_final.f90
source77=$(wildcard ./minpack/*.f)
object90=$(source90:.f90=.o)
object77=$(patsubst %.f,%.o,$(notdir $(source77))) $(wildcard ./slatec/*.o)
FCC = ifort
ARFLAG = -rs
OMPFLAG = -qopenmp
#OMPFLAG = -openmp -warn nounused -traceback
all:$(out)
$(out):$(object90) $(object77)
	ar $(ARFLAG) $@ $(object90) $(object77)
$(object77):$(source77)
	$(FCC) -c $(source77) $(OMPFLAG)
$(object90):$(source90)
	$(FCC) -c $(source90) $(OMPFLAG)
#%.o: %.f90
	#$(FCC) -c $< $(OMPFLAG)
#%.o: %.f
	#$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.a *.mod *.o
test:
	@echo $(object77)
