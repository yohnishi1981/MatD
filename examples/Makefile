include ../Makefile.config
SUBDIRS = example*
BINDIR = $(shell pwd)/bin

subdirs:
	for dir in $(SUBDIRS); do \
	  cd $(shell pwd)/$$dir;\
	  $(FC) main.F90 $(FFLAGS) $(DPPFLAGS) -o $(BINDIR)/$$dir -I$(MATD_INCLUDE) $(LIBS) -L$(MATD_LIB) -lmatd;\
	done
