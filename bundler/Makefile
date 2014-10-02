#-----------------------------------------------------------------------------
# Top level makefile for Bundler
#
# Bundler: Structure from Motion for Unordered Photo Collections
# Version: 0.4 04/03/2010
#    http://phototour.cs.washington.edu/bundler/
#-----------------------------------------------------------------------------
# Copyright (c) 2008-2010 University of Washington and Noah Snavely
# All Rights Reserved.
#-----------------------------------------------------------------------------

ANN_TARGET = linux-g++-shared

OS = $(shell uname -o)
ifeq ($(OS), Cygwin)
ANN_TARGET = win32-g++-shared
endif

TOP_DIR = $(shell pwd)
USE_OPENSIFT = 0

default:
# Make libraries
ifeq ($(USE_OPENSIFT), 1)
	cd lib/opensift; $(MAKE); ln -s $(TOP_DIR)/lib/opensift/bin/siftfeat $(TOP_DIR)/bin/siftfeat;  		#Compile and link opensift
endif
	cd lib/5point; $(MAKE)
	cd lib/ann_1.1_char; $(MAKE) $(ANN_TARGET)
	cd lib/imagelib; $(MAKE)
	cd lib/matrix; $(MAKE)
	cd lib/sba-1.5; $(MAKE)
	cd lib/sfm-driver; $(MAKE)
	#cd jpeg/src; $(MAKE)
# Auxiliary libraries
	cd lib/minpack; $(MAKE)
	cd lib/cblas; $(MAKE)
	cd lib/f2c; $(MAKE)
# Main program
	cd src; $(MAKE)


clean:
	cd lib/opensift; $(MAKE) distclean; rm -f $(TOP_DIR)/bin/siftfeat;
	cd lib/5point; $(MAKE) clean
	cd lib/ann_1.1_char; $(MAKE) clean
	rm -f bin/libANN_char.so lib/libANN_char.so lib/ann_1.1_char/lib	#not cleaned up by ann1.1 all 3 are copies
	cd lib/imagelib; $(MAKE) clean
	cd lib/matrix; $(MAKE) clean
	cd lib/sba-1.5; $(MAKE) clean
	rm -f lib/sba-1.5/libsba.v1.5.a						#not cleaned up by sba-1.5
	cd lib/sfm-driver; $(MAKE) clean
	cd lib/minpack; $(MAKE) clean
	cd lib/cblas; $(MAKE) clean
	rm -f lib/cblas/libcblas.a 						#not cleaned up by cblas
	cd lib/f2c; $(MAKE) clean
	cd src; $(MAKE) clean
	rm -f bin/bundler bin/KeyMatchFull bin/Bundle2PMVS bin/Bundle2Vis bin/RadialUndistort bin/Bundle2Ply
	rm -f lib/*.a
