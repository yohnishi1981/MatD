SUBDIRS = example*
BINDIR = $(shell pwd)/bin

subdirs:
	for dir in $(SUBDIRS); do \
	  cd $(shell pwd)/$$dir;\
	  $(FC) main.f90 -o $(BINDIR)/$$dir -I$(MATD_INCLUDE) -L$(MATD_LIB) -lmatd;\
	done
