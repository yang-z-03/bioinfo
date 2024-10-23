
#include <stdio.h>

// read all the content of a text sam file into the memory.
// this method will take a great (relatively equal quantity of memory as the
// file on disk) amount of memory. but since the server has that amount of
// memory, i prefer read all at once for better speed.

char* read_all(FILE* f, size_t *fsize);
char* read_all_stdin(size_t *fsize);

// find the next line and returns a pointer pointing to the first character
// of the next line. (this will be the starting point of the sam's name entry)

char* next_line(char* buffer);

// get the length of the current line

int line_length(char* buffer);

// get the length before next tab, or a line-ending.

int field_length(char* buffer);
