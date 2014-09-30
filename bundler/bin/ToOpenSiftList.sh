#!/bin/bash
#
# ToOpenSiftList.sh
# Create a script for extracting sift features from a list of images

# Set this variable to your base install path (e.g., /home/foo/bundler)
BIN_PATH=$(dirname $(which $0));

if [ $# -ne 1 ]
then
  echo "Usage: ToOpenSiftList.sh <list_file>"
  exit;
fi

OS=`uname -o`

if [ $OS == "Cygwin" ]
then
    echo "Cygwin not supported by Open Sift"  > /dev/stderr
    exit 1
else
    SIFT=$BIN_PATH/siftfeat
fi

if ! [ -e $SIFT ] ; then
    echo "[ToOpenSiftList] Error: SIFT not found.  Please install SIFT to $BIN_PATH" > /dev/stderr
    exit 1
fi

IMAGE_LIST=$1

awk "{pgm = \$1; key = \$1; sub(\"jpg\$\", \"pgm\", pgm); sub(\"jpg\$\", \"key\", key); print \"/usr/bin/mogrify -format pgm \" \$1 \"; $SIFT -x -o \" key \" \" pgm \"; gzip -f \" key \"; rm \" pgm}" $IMAGE_LIST 
