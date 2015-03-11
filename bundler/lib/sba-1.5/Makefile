#
# Makefile for Sparse Bundle Adjustment library & demo program
#
include ../../environ.gmake
COMPILER=gcc
CC=$(COMPILER)

INCLUDE_FLAGS=-I../matrix -I../../include

#OPTFLAGS=-O3

OBJS=sba_levmar.o sba_levmar_wrap.o sba_lapack.o sba_crsm.o sba_chkjac.o sba_opencl.o
SRCS=sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c sba_opencl.c
AR=ar
RANLIB=ranlib
MAKE=make

LIB_PATH=-lOpenCL

LIB=libsba.v1.5.a

CFLAGS=-Wall $(OPTFLAGS) $(INCLUDE_FLAGS) $(LIB_PATH)

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) crv $(LIB) $(OBJS)
	$(RANLIB) $(LIB)
	cp $(LIB) ..

sba_levmar.o: sba.h sba_chkjac.h compiler.h sba_opencl.h
sba_levmar_wrap.o: sba.h
sba_lapack.o: sba.h compiler.h
sba_crsm.o: sba.h
sba_chkjac.o: sba.h sba_chkjac.h compiler.h
sba_opencl.o: sba_opencl.h

clean:
	@rm -f *.o

realclean cleanall: clean
	@rm -f libsba.a

depend:
	makedepend -f Makefile $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
