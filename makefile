ARCH := $(shell getconf LONG_BIT)
FCC = ifort
ifeq ($(DEBUG),A)
	CFLAG= $(AFLAG) -O0 -warn nounused  -check noarg_temp_created -check bounds -g -check all -fpe0 -traceback -debug extended  -fstack-protector  -assume protect_parens  -check -ftrapuv -init=snan,arrays
else ifeq ($(DEBUG),B)
	CFLAG= $(AFLAG) -warn nounused -traceback -init=snan,arrays -fp-speculation=safe
else
	CFLAG= $(AFLAG) -warn nounused
endif
ifeq ($(ARCH),64)
	LFLAG_A= $(LFLAG) -lmy -lnlopt -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lm -ldl
	#LFLAG_A= $(LFLAG) -lmy -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lnlopt -lm
else
	LFLAG_A= $(LFLAG) -lmy -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_lapack95
	#LFLAG_A= $(LFLAG) -lmy -mkl=sequential -lmkl_lapack95 -lnlopt -lm
endif
all:$(out)
%.out: %.o
	$(FCC) $< -o $@ $(LFLAG_A)
%.o: %.f90
	$(FCC) -c $< $(CFLAG)
clean:
	rm -f *.out *.mod *.o
