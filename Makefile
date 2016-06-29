## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=g++ -std=c++11 -g
#########


## Flags to give the compiler for "release mode"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -I/media/xuejian/WORK/spectra/spectra-0.2.0/spectra-0.2.0/include/ -I/media/xuejian/WORK/spectra/eigen-eigen-07105f7124f9/







obj=main.o Parameter.o OP.o Sub.o QWave.o Super.o DMRGP.o Corr.o
main:$(obj)
	$(CCCOM) -o main $(obj)  $(LIBSPECTRA)
main.o:main.cpp DMRGP.h
	$(CCCOM) -c main.cpp $(LIBSPECTRA)
Parameter.o:Parameter.h Parameter.cpp
	$(CCCOM) -c Parameter.cpp -O2 $(LIBSPECTRA)
OP.o:OP.cpp OP.h Parameter.h
	$(CCCOM) -c OP.cpp -O2 $(LIBSPECTRA)
Sub.o:Sub.cpp Sub.h OP.h  Parameter.h
	$(CCCOM) -c Sub.cpp -O2 $(LIBSPECTRA)
QWave.o:QWave.cpp QWave.h Sub.h
	$(CCCOM) -c QWave.cpp -O2 $(LIBSPECTRA)
Super.o:Super.cpp Super.h 
	$(CCCOM) -c Super.cpp -O2 $(LIBSPECTRA)
DMRGP.o:DMRGP.cpp DMRGP.h SuperEnergy.h Sub.h QWave.h
	$(CCCOM) -c DMRGP.cpp -O2 $(LIBSPECTRA)
Corr.o:Corr.cpp Corr.h Super.h
	$(CCCOM) -c Corr.cpp -O2 $(LIBSPECTRA)
.PHONY:clean
clean:
	rm -f main $(obj)















