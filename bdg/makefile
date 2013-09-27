# mkllib=/sbin/ifort/mkl/lib/intel64
# mklinc=/sbin/ifort/mkl/include
# FCC = ifort
# FFT.out: FFT.o
# FFT.o: FFT.f90
	# $(FCC) FFT.f90 -c -g
# clean:
	# rm FFT.o FFT.out
out=bdg_flux.out
#FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCCFLAG= -g -llapack 
FCC = gfortran
#OMPFLAG = -openmp
all : $(out)
%.out : %.o
	$(FCC) $< -o $@ ${OMPFLAG} $(FCCFLAG)
%.o : %.f90
	$(FCC) -c $< ${OMPFLAG}
clean:
	rm -f *.out *.mod
