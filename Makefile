all: matd ex

# The Compiler used in project.
export FC=mpif90

# The directory of Matd library files and include files.
export MATD_DIR := $(shell pwd)
export MATD_LIB = $(MATD_DIR)/lib
export MATD_INCLUDE = $(MATD_DIR)/include
export CFLAGS=

# Compile Matd library.
matd:
	$(MAKE) -C src

# Compile examples that use Matd library.
ex:
	$(MAKE) -C examples
