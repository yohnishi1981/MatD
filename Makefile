include ./Makefile.config

all: matd ex

# Compile Matd library.
matd:
	$(MAKE) -C src

# Compile examples that use Matd library.
ex:
	$(MAKE) -C examples

.PHONY: cleansrc
cleansrc:
	rm -f ./src/*.mod ./src/*.a ./src/*.o

.PHONY: cleanlib
cleanlib:
	rm -f $(MATD_LIB)/* $(MATD_INCLUDE)/*

.PHONY: cleanex
cleanex:
	rm -f ./examples/bin/*

.PHONY: veryclean
veryclean: cleansrc cleanlib cleanex
