ARCH := $(shell getconf LONG_BIT)
FCC = ifort
#CFLAG= $(AFLAG) -warn nounused -traceback -assume realloc_lhs
CFLAG= $(AFLAG) -warn nounused -traceback
#CFLAG= -O0 -warn nounused -check bounds -g -check all -fpe0 -traceback -debug extended  -fstack-protector  -assume protect_parens  -check -ftrapuv
#LFLAG_64= $(LFLAG) -lmy -mkl=sequential -lmkl_lapack95_lp64 -lnlopt -lm
#LFLAG_32= $(LFLAG) -lmy -mkl=sequential -lmkl_lapack95 -lnlopt -lm
#LFLAG_64= $(LFLAG) -lmy -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lnlopt -lm
#LFLAG_32= -lmy -mkl -lmkl_lapack95 
LFLAG_64= $(LFLAG) -lmy -lnlopt -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lm -ldl
#LFLAG_32= $(LFLAG) -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
#LFLAG_64= $(LFLAG) -lmy -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lnlopt -lm
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(LFLAG_$(ARCH))
%.o: %.f90
	$(FCC) -c $< $(CFLAG)
clean:
	rm -f *.out *.mod *.o
