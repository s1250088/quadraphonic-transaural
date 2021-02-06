 # Makefile for quadraTrans

PDLIBBUILDER_DIR=.
DOXYGEN=/Applications/Doxygen.app/Contents/Resources/doxygen

cflags =  -I/usr/local/include/ -mmacosx-version-min=10.15

ldflags = -L"/usr/local/lib/" -lfftw3f -lsqlite3 

lib.name = quadraTrans~

class.sources = quadraTrans~.c 

datafiles = quadraTrans~-help.pd

include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder


#-------Documentation-----
doc:
	${DOXYGEN} Doxyfile.cnf
	astyle --style=java --suffix=none --errors-to-stdout *.c 
