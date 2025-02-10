
cc  = gcc
lib = -lm

# debug = -g -O0 -fno-inline
debug = -O3

all: extrumi reptag

shared_sam_c = sam.c
shared_sam_h = sam.h

extrumi: extrumi.c $(shared_sam_c) $(shared_sam_h)
	$(cc) extrumi.c $(shared_sam_c) $(shared_sam_h) -o extrumi $(debug)

reptag: reptag.c $(shared_sam_c) $(shared_sam_h)
	$(cc) reptag.c $(shared_sam_c) $(shared_sam_h) -o reptag $(debug)