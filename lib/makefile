out=libmy.a
source=$(wildcard *.f90)
object=$(source:.f90=.o)
FCC = ifort
ARFLAG = -rs
OMPFLAG = -openmp
all:$(out)
%.a:$(object)
	ar $(ARFLAG) $@ $(object)
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG)
clean:
	rm -f *.a *.mod *.o
