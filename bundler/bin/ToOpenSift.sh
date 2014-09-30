#!/bin/bash
#
# ToOpenSift.sh
# Create a script for extracting sift features from a set of images

# Set this variable to your base install path (e.g., /home/foo/bundler)
BIN_PATH=$(dirname $(which $0));

IMAGE_DIR="."

if [ $# -eq 1 ]
then
    IMAGE_DIR=$1
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
    echo "[ToOpenSift] Error: SIFT not found.  Please install SIFT to $BIN_PATH" > /dev/stderr
    exit 1
fi

for d in `ls -1 $IMAGE_DIR | egrep "jpg$"`
do 
    pgm_file=$IMAGE_DIR/`echo $d | sed 's/jpg$/pgm/'`
    key_file=$IMAGE_DIR/`echo $d | sed 's/jpg$/key/'`
    echo "mogrify -format pgm $IMAGE_DIR/$d; $SIFT -x -o $key_file $pgm_file; rm $pgm_file; gzip -f $key_file"
done
