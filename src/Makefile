include ../Makefile.config

all: libmatd.a install 

libmatd.a : matd.o
	ar crv $@ matd.o

matd.o : matd.f90
	$(FC) matd.f90 -c $(FFLAGS)

# Preprocess
matd.f90 : matd_.f90 matd.h
	$(FC) matd_.f90 -P -o $@

install : libmatd.a
	cp libmatd.a $(MATD_LIB)
	cp *.mod $(MATD_INCLUDE)

clean :
	rm -f *.mod *.a *.o matd.f90
