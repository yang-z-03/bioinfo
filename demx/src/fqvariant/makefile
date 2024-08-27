
cc = gcc

all: lib test

lib:
	$(cc) -fpic -shared fqvariant.c -o libfqvariant.so

test: lib
	$(cc) libfqvariant.so main.c fqvariant.h -o main