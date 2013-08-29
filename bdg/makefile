# mkllib=/sbin/ifort/mkl/lib/intel64
# mklinc=/sbin/ifort/mkl/include
# FCC = ifort
# FFT.out: FFT.o
# FFT.o: FFT.f90
	# $(FCC) FFT.f90 -c -g
# clean:
	# rm FFT.o FFT.out
out=BdG_CRS.out
#FCCFLAG= -openmp -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCCFLAG= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
FCC = ifort
all : $(out)
%.out : %.o
	$(FCC) $< -o $@ -openmp $(FCCFLAG)
%.o : %.f90
	$(FCC) -c $< -openmp
clean:
	rm -f *.out *.mod
