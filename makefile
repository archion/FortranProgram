ARCH := $(shell getconf LONG_BIT)
FCC = ifort
OMPFLAG= -openmp
CFLAG=-warn nounused
#CFLAG=-warn nounused -traceback
#CFLAG=-warn nounused -check bounds -g -check all -fpe0 -traceback -debug extended
FCCFLAG_64= -lmy -mkl=sequential -lmkl_lapack95_lp64 -lnlopt -lm
FCCFLAG_32= -lmy -mkl=sequential -lmkl_lapack95 -lnlopt -lm
#FCCFLAG_64= -lmy -mkl -lmkl_lapack95_lp64
#FCCFLAG_32= -lmy -mkl -lmkl_lapack95
#FCCFLAG_64= -lmy -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lnlopt -lm
#FCCFLAG_32= -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
#FCCFLAG= -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
FCCFLAG_64= -lmy -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lnlopt -lm
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(FCCFLAG_$(ARCH))
%.o: %.f90
	$(FCC) -c $< $(OMPFLAG) $(CFLAG)
clean:
	rm -f *.out *.mod *.o
