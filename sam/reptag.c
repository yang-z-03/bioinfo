
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "sam.h"

int main(int argc, char* argv[]) {

    // read the sam file to the memory
    
    char* content = NULL;
    size_t fsize = 0;
  
    if (strcmp(argv[1], "-") == 0) {
        content = read_all_stdin(&fsize);
    } else {
        FILE* f = fopen(argv[1], "r");
        content = read_all(f, &fsize);
    }

    fprintf(stderr, "read file size: %ld \n", fsize);

    // initialize the output file

    FILE* fout = NULL;
    if (strcmp(argv[1], "-") == 0) fout = stdout;
    else fout = fopen(argv[2], "w");

    if (fout == NULL) {
        perror("failed to initialize output stream.");
        exit(4);
    }

    if (content == NULL) {
        perror("failed to read.");
        exit(1);
    }

    char* begin = content;
  
    // i assume that the name field will not exceed 512 characters.
    // at least in my case it is true.
  
    char name_buffer[512];

    do {

        char firstc = *content;
        if (firstc == '\0') break;
        
        if (firstc == '@') {
            char* nl = next_line(content);
            if (nl == NULL) {
                fprintf(fout, "%s\n", content);
                break;
            } else {
                *(nl - 1) = '\0';
                fprintf(fout, "%s\n", content);
                content = nl;
                continue;
            }
        }
        
        // get the sam records' name
        int should_break = 0;
        while (1) {

            if (*content == 0) { should_break = 1; break; }
            int flen = field_length(content);
            if (flen == 0) { should_break = 1; break; }

            // copy the field to the name buffer
            // now flen must be greater than zero, indicating there is a field.
            strncpy(name_buffer, content, flen);
            name_buffer[flen] = '\0';
                
            if (flen > 6) {
                // replacing condition
                if (strncmp("CB:Z:", name_buffer, 5) == 0) {
                    fprintf(fout, "DB:Z:%s", name_buffer + 5);
                }
                else if (strncmp("DB:Z:", name_buffer, 5) == 0) {
                    fprintf(fout, "CB:Z:%s", name_buffer + 5);
                }
                else fprintf(fout, "%s", name_buffer);
            } else fprintf(fout, "%s", name_buffer);

            if (*(content + flen) == '\0') {
                // end of file
                should_break = 1; break;
            } else if (*(content + flen) == '\t') {
                // goto the next field.
                content += (flen + 1);
                fprintf(fout, "\t");
                continue;
            } else if (*(content + flen) == '\n') {
                // end of line
                content += (flen + 1);
                fprintf(fout, "\n");
                break;
            } else { should_break = 1; break; }
        }
        
        // double check
        if (content == NULL) break;
        if (*content == 0) break;
        if (should_break != 0) break;

    } while (1);

    // fflush(fout);
    fclose(fout);

    return 0;
}