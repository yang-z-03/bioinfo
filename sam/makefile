
cc  = gcc
lib = -lm

# debug = -g -O0 -fno-inline
debug = -O3

all: extrumi

shared_sam_c = sam.c
shared_sam_h = sam.h

extrumi: extract-umi.c $(shared_sam_c) $(shared_sam_h)
	$(cc) extract-umi.c $(shared_sam_c) $(shared_sam_h) -o extrumi $(debug)
