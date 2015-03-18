# Makefile for bundler

#USE_CERES=true
include ../environ.gmake
CC=gcc
CXX=g++
#OPTFLAGS=-O3 -Wall -Wno-unused-result
#OPTFLAGS=-g2
OTHERFLAGS=-std=gnu++0x

OS=$(shell uname -o)

ifeq ($(OS), Cygwin)
BUNDLER=bundler.exe
KEYMATCHFULL=KeyMatchFull.exe
KEYMATCHFAST=key_match_fast.exe
BUNDLE2PMVS=Bundle2PMVS.exe
BUNDLE2VIS=Bundle2Vis.exe
RADIALUNDISTORT=RadialUndistort.exe
OPENCLTEST=opencl_test.exe
else
BUNDLER=bundler
KEYMATCHFULL=KeyMatchFull
KEYMATCHFAST=key_match_fast
BUNDLE2PLY=Bundle2Ply
BUNDLE2PMVS=Bundle2PMVS
BUNDLE2VIS=Bundle2Vis
RADIALUNDISTORT=RadialUndistort
OPENCLTEST=opencl_demo
endif

INCLUDE_PATH=-I../lib/imagelib -I../lib/sfm-driver -I../lib/matrix	\
	-I../lib/5point -I../lib/sba-1.5 -I../lib/ann_1.1_char/include  \
	-I../include -I../lib/tbb/include -I ../lib/flann-1.8.4/src/cpp

LIB_PATH=-L../lib -L../lib/ann_1.1_char/lib -L/usr/lib/x86_64-linux-gnu -L../lib/tbb/lib

BUNDLER_DEFINES=-D__NO_UI__ -D__BUNDLER__ -D__BUNDLER_DISTR__

BUNDLER_OBJS=BaseApp.o BundlerApp.o keys.o Register.o Epipolar.o	\
	Bundle.o BundleFast.o MatchTracks.o Camera.o Geometry.o		\
	ImageData.o SifterUtil.o BaseGeometry.o BundlerGeometry.o	\
	BoundingBox.o BundleAdd.o ComputeTracks.o BruteForceSearch.o	\
	BundleIO.o ProcessBundle.o BundleTwo.o Decompose.o		\
	RelativePose.o Distortion.o TwoFrameModel.o LoadJPEG.o

BUNDLER_BASE_LIBS=-limage -lsfmdrv -lsba.v1.5 -lmatrix -lz -llapack -lblas \
	-lcblas -lminpack -lm -l5point -ljpeg -lANN_char -lgfortran -ltbb  \
	-lOpenCL -lpthread
CERES_LIBS=-lceres -lcholmod -lccolamd -lcolamd -lamd -lcamd -lcxsparse \
	-lsuitesparseconfig -lgomp -lglog -lpthread

# if we are using Ceres, add to options
ifeq ($(USE_CERES),true)
# Edit these two variables as needed
CERES_INCLUDE_PATH=
CERES_LIB_PATH=
BUNDLER_DEFINES+=-D__USE_CERES__
BUNDLER_LIBS=$(CERES_LIBS) $(BUNDLER_BASE_LIBS)
BUNDLER_OBJS+=BundleCeres.o
INCLUDE_PATH+=$(CERES_INCLUDE_PATH)
LIB_PATH+=$(CERES_LIB_PATH)
INCLUDE_PATH+=-I/usr/include/eigen3
else
BUNDLER_LIBS=$(BUNDLER_BASE_LIBS)
endif

CPPFLAGS=$(OPTFLAGS) $(OTHERFLAGS) $(INCLUDE_PATH)

all: $(BUNDLER) $(KEYMATCHFULL) $(KEYMATCHFAST) $(BUNDLE2PLY) $(BUNDLE2PMVS) $(BUNDLE2VIS) $(RADIALUNDISTORT) $(OPENCLTEST)

%.o : %.cpp 
	$(CXX) -c -o $@ $(CPPFLAGS) $(WXFLAGS) $(BUNDLER_DEFINES) $<

$(BUNDLER): $(BUNDLER_OBJS)
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) \
		$(BUNDLER_DEFINES) $(BUNDLER_OBJS) $(BUNDLER_LIBS)
	cp $@ ../bin

$(KEYMATCHFULL): KeyMatchFull.o keys2a.o key_util.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) KeyMatchFull.o keys2a.o key_util.o\
		-lANN_char -lz -ltbb
	cp $@ ../bin
	
$(KEYMATCHFAST): key_match_fast.o keys2a.o key_util.o key_match.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) key_match_fast.o keys2a.o key_util.o key_match.o\
		-lANN_char -lz -ltbb
	cp $@ ../bin

$(BUNDLE2PLY): Bundle2Ply.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) Bundle2Ply.o \
		-lmatrix -llapack -lblas -lcblas -lgfortran -lminpack
	cp $@ ../bin

$(BUNDLE2PMVS): Bundle2PMVS.o LoadJPEG.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) Bundle2PMVS.o LoadJPEG.o \
		-limage -lmatrix -llapack -lblas -lcblas -lgfortran \
		-lminpack -ljpeg
	cp $@ ../bin

$(BUNDLE2VIS): Bundle2Vis.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) Bundle2Vis.o
	cp $@ ../bin

$(RADIALUNDISTORT): RadialUndistort.o LoadJPEG.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) $^ \
		-limage -lmatrix -llapack -lblas -lcblas -lgfortran \
		-lminpack -ljpeg
	cp $@ ../bin

$(OPENCLTEST): opencl_demo.o
	$(CXX) -o $@ $(CPPFLAGS) $(LIB_PATH) opencl_demo.o \
		-lOpenCL
	cp $@ ../bin

clean:
	rm -f *.o *~ $(BUNDLER) $(KEYMATCHFULL) $(BUNDLE2PMVS) \
			$(BUNDLE2VIS) $(RADIALUNDISTORT) $(BUNDLE2PLY) \
			$(OPENCLTEST)
