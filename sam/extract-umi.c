
#include <stdio.h>
#include <string.h>

#include "sam.h"

int main() {

    // read the sam file to the memory
    char* content = read_all();

    // i assume that the name field will not exceed 512 characters.
    // at least in my case it is true.
    char name_buffer[512];

    do {

        // get the sam records' name
        int name_length = field_length(content);
        if (name_length == 0) break;

        strncpy(name_buffer, content, name_length);
        name_buffer[name_length] = '\0';

        // the name has a format with: xxxxx_barcode_umi
        char* first_underscore = strchr(name_buffer, '_');
        char* second_underscore = strchr(first_underscore + 1, '_');
        *first_underscore = '\0';
        first_underscore += 1;
        *second_underscore = '\0';
        second_underscore += 1;

        // now the first and second underscore pointer is pointed to the barcode
        // and umi substring, both terminated with \0

        fprintf(stdout, "%s", name_buffer);

        // here, the truncated name is printed, we then print the rest of line

        content += name_length;
        char* next_line_start = next_line(content);
        
        if (next_line_start == NULL) {
            
            // content already contains a \0 at its end.
            fprintf(stdout, "%s", content);

        } else {
            *(next_line_start - 1) = '\0';
            fprintf(stdout, "%s", content);
            content = next_line_start;
        }

        // append the extra tags to the line:
        
        fprintf(
            stdout, "\tCB:Z:%s\tUB:Z:%s\n",
            first_underscore, second_underscore
        );

        if (next_line_start == NULL) break;

    } while (1);

    return 0;
}