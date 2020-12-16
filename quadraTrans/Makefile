# Makefile to build class 'quadraTrans~' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
#lib.name = quadraTrans~

# input source file (class name == source file basename)
#class.sources = quadraTrans~.c

# all extra files to be included in binary distribution of the library
#datafiles =

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
#PDLIBBUILDER_DIR=../pd-lib-builder/
#include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder

#--------------------------------------------------
current:
	echo make pd_darwin, pd_linux, pd_nt, pd_irix5, or pd_irix6

clean: ; rm -f -r quadraTrans~.pd_darwin *.o

#----------
pd_darwin: quadraTrans~.pd_darwin

.SUFFIXES: .pd_darwin


DARWININCLUDE =  -I/Applications/Pd-0.51-3.app/Contents/Resources/src -I/usr/local/include

DARWINCFLAGS = -DPD -O3 -Ofast -Wall -W -Wshadow -Wstrict-prototypes \
    -Wno-unused -Wno-parentheses -Wno-switch -Wno-unknown-pragmas -arch i386 -arch x86_64

DARWINLIBS = -L/usr/local/lib -lfftw3f -lm -bundle -undefined suppress -arch i386 -arch x86_64 -flat_namespace

.c.pd_darwin:
	cc $(DARWINCFLAGS) $(DARWININCLUDE) -o $*.o -c $*.c
	cc $(DARWINLIBS)  -o $*.pd_darwin $*.o 
	rm -f $*.o

