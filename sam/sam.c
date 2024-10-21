
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>

#include "sam.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

char* read_all() {

    // the reading buffer size is initialized to 1 Mb. a sam file usually takes
    // up to ~20 Gbs. so this buffer size is not so big.

    #define read_buffer_size (1024 * 1024)

    char buffer[read_buffer_size];
    size_t content_size = 1;
    
    char *content = malloc(sizeof(char) * read_buffer_size);
    if (content == NULL) {
        perror("failed to allocate content");
        exit(1);
    }

    content[0] = '\0';
    while(fgets(buffer, read_buffer_size, stdin))
    {
        char *old = content;
        content_size += strlen(buffer);
        content = realloc(content, content_size);
        
        if (content == NULL) {
            perror("failed to reallocate content");
            free(old);
            exit(2);
        }

        strcat(content, buffer);
    }

    if (ferror(stdin)) {
        free(content);
        perror("error reading from stdin.");
        exit(3);
    }

#ifdef debug

    if (feof(stdin)) {
        printf("successfully reached end-of-file. \n");
    }

    printf("read sam file: content length %d \n", content_size);

#endif

    return content;
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