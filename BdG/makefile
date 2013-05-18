mkllib=/sbin/ifort/mkl/lib/intel64
mklinc=/sbin/ifort/mkl/include
FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
FCC = ifort
TEST: TEST.o
	$(FCC) -o TEST.out TEST.o -I$(mklinc) -I$(mklinc)/intel64/lp64 -L$(mkllib) $(FCCFLAG) -g
TEST.o: TEST.f90
	$(FCC) TEST.f90 -c -g
clean:
	rm TEST.o TEST.out