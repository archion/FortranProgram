out=libmy.a
source90=$(wildcard *.f90)
source77=$(wildcard ./minpack/*.f)
object90=$(source90:.f90=.o)
object77=$(patsubst %.f,%.o,$(notdir $(source77)))
FCC = ifort
ARFLAG = -rs
OMPFLAG = -openmp
#OMPFLAG= -warn nounused -traceback
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
