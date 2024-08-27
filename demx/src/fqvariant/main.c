
#include "fqvariant.h"
#include <stdio.h>

const char* seqs = "AACGATTG"
                   "AATCTACG"
                   "ATTGATGT"
                   "AGTTATAT"
                   "AAAATTTT";

int main(void) {
    int len = 8; int seg = 1;
    int run = len - seg + 1;
    size_t* results = fqvariant(len, seg, 0, 4, seqs, 5, 8);
    for (int i = 0; i < run; i ++) {
        printf("%d = %d\n", i, results[i]);
    }

    return 0;
}