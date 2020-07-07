#!/bin/bash
#
# Author: Massimo Torquati
#
# File decompression script leveraging Miniz comp/decomp capabilities.
# This script shows how the decompression of "BIG files" is done.
#
# If the original file was a "BIG file" (i.e. it was split in multiple
# independent files), it opens the file in a temporary directory (using tar),
# decompress each part by using COMPDECOMP and then merge them all in a
# single monolithic file.
# Instead, if the file was not a "BIG file", it just calls COMPDECOMP
# to decompress it.
#

# where the compressor/decompressor executable is
COMPDECOMP=/tmp/compdecomp
# sanity checks
if [ ! -f  $COMPDECOMP ]; then
    echo "Error cannot find the compdecomp executable"
    exit -1
fi
if [ ! -x  $COMPDECOMP ]; then
    echo "Error the compdecomp file is not executable"
    exit -1
fi

# checking input
if [ $# -ne 1 ]; then    
    echo use: $(basename $0) file
    exit -1
fi
infile=$1               
if [ ! -f $infile ]; then    
    echo "Argument is not a file"
    exit -1;   
fi

# check if the file was a "BIG file" compressed in multiple parts
# "BIG compressed files" have a tar magic header 
r=$(file $infile | grep tar)
if [ $? -ne 0 ]; then  # not a "BIG file", simple case
    r=$($COMPDECOMP D $infile >& /dev/null)
    if [ $? -ne 0 ]; then
	echo "Error: cannot decompress the file"
	tmpfile=$infile
	tmpfile+="_decomp"
	rm -f $tmpfile >& /dev/null
	exit -1
    fi
else
    # the file is a TAR composed by multiple compressed parts
    tmpdir=$(mktemp -d)
    if [ $? -ne 0 ]; then  
	echo "Error: cannot create a temporary directory"
	exit -1
    fi
    # opening the tar file
    tar xf $infile -C $tmpdir
    if [ $? -ne 0 ]; then  
	echo "Error: error during untar"
	exit -1
    fi
    # decompressing the single parts, then removing the compressed file
    $COMPDECOMP D $tmpdir >& /dev/null
    if [ $? -ne 0 ]; then  
	echo "Error: error decompressing file parts"
	exit -1
    fi
    # how many parts?
    nblocks=$(find $tmpdir -type f | wc -l)
    
    f1=${infile%.*} # removing the extension  
    f2=${f1##*/}    # this is the filename
    
    # merging the uncompressed parts
    for i in $(seq 1 1 $nblocks); do
	find $tmpdir -type f -name "*.part$i" -exec cat {} >> $tmpdir/$f2 \;
	if [ $? -ne 0 ]; then  
	    echo "Error: error merging file parts"
	    exit -1
	fi    
    done
    mv $tmpdir/$f2 .
    if [ $? -ne 0 ]; then  
	echo "Error: cannot move the uncompressed file in $PWD"
	exit -1
    fi    
    rm -fr $tmpdir >& /dev/null
fi
echo "Done."
echo -n "Remove the $infile file (S/N)?"
read yn
if [ x$yn == x"S" ]; then
    rm -f $infile
fi

exit 0
