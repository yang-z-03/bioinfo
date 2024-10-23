
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>

#include "sam.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

char* read_all(FILE* f, size_t *fsize) {

    // the reading buffer size is initialized to 1 Gb. a sam file usually takes
    // up to ~50 Gbs. so this buffer size is not so big.
  
    #define read_buffer_size (1024 * 1024 * 1024)
    #define allocate_space (1024L * 1024L * 1024L * 100L)

    char* buffer = malloc(sizeof(char) * read_buffer_size);
    size_t content_size = 1;
    
    char *content = malloc(sizeof(char) * allocate_space);
    if (content == NULL) {
        perror("failed to allocate content");
        exit(1);
    }

    content[0] = '\0';

    // note! the fread doesn't append \0 ending to the buffer.
    // we need to append \0 by ourself. otherwise string concatenation
    // will be erronous.
    
    size_t getbyte = fread(buffer, sizeof(char), read_buffer_size - 1, f);
    fprintf(stderr, "get: %ld, last char %c %c. \n", getbyte, buffer[getbyte - 1], buffer[getbyte]);
    
    while(getbyte > 0)
    {
        char *old = content;
        content_size += (getbyte);
        buffer[getbyte] = '\0';
      
        if (content_size > allocate_space)
          content = realloc(content, content_size);
        
        if (content == NULL) {
            perror("failed to reallocate content");
            free(old);
            exit(2);
        }

        strcat(content, buffer);
        getbyte = fread(buffer, sizeof(char), read_buffer_size - 1, f);
        fprintf(stderr, "get: %ld, last char %c %c. \n", getbyte, buffer[getbyte - 1], buffer[getbyte]);
    }

    if (ferror(f)) {
        free(content);
        perror("error reading from file.");
        exit(3);
    }

#ifdef debug

    if (feof(f)) {
        printf("successfully reached end-of-file. \n");
    }

    printf("read sam file: content length %ld \n", content_size);

#endif

    *fsize = content_size;
    return content;
}

char* read_all_stdin(size_t *fsize) {
    return read_all(stdin, fsize);
}

char* next_line(char* buffer) {
    char* next_line = strchr(buffer, '\n');
    if (next_line == NULL) return NULL;
    return next_line + 1;
}

// this line length doesn't contain the \n character.
// it is the exact text length of the line.

int line_length(char* buffer) {
    if (*buffer == '\0') return 0;
    char* next_line = strchr(buffer, '\n');
    if (next_line == NULL) return strlen(buffer);
    return next_line - buffer;
}

int field_length(char* buffer) {
    if (*buffer == '\0') return 0;
    char* next_field = strchr(buffer, '\t');
    if (next_field == NULL) return 0;
    char* next_line = strchr(buffer, '\n');
    if (next_line == NULL) return next_field - buffer;
    return min(next_line - buffer, next_field - buffer);
}