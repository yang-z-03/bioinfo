
#include <pthread.h>
#include <math.h>
#include <malloc.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

typedef struct { char* data; int rows; int columns; } t_params;
typedef struct { size_t unique; } t_out;

int cmp(char* c1, char* c2, size_t length) {
    for (int x = 0; x < length; x ++) {
        if (c1[x] != c2[x]) return 1;
    }

    return 0;
}

void* unique(void* args) {
    t_params* params = (t_params*) args;
    t_out* out = malloc(sizeof(t_out));
    out->unique = 0;

    char* data = params->data;

    for (size_t x = 0; x < params->columns * (params->rows - 1); x += params->columns) {
        if (data[x] == 0) continue;
        for (size_t y = x + params->columns; y < params->columns * params->rows; y += params->columns) {
            if (data[y] == 0) continue;
            if (cmp(data + x, data + y, params->columns) == 0) data[y] = 0;
        }
    }

    for (size_t x = 0; x < params->columns * (params->rows); x += params->columns) {
        if (data[x] == 0) continue;
        out->unique ++;
    }

    return out;
}

size_t* fqvariant(int lseq, int lsegment, int lscatter, int nthread,
                  const char* data, int rows, int columns) {
    
    int runs = lseq - lsegment + 1;
    int batches = runs / nthread;
    if (runs % nthread != 0) batches += 1;

    pthread_t* threads = (pthread_t*) malloc(runs * sizeof(pthread_t));
    size_t* results = (size_t*) malloc(runs * sizeof(size_t));
    t_params* params = (t_params*) malloc(runs * sizeof(t_params));

    for (int ibatch = 0; ibatch < batches; ibatch ++) {
        int thrstart = nthread * ibatch;
        int thrend = min(nthread * (ibatch + 1), runs);

        for (int ithr = thrstart; ithr < thrend; ithr++) {
            int run = ithr;
            int pstart = max(0, run - lscatter);
            int pend = min(run + lsegment + lscatter, lseq);

            /* now we shall create heap data (only to regions of interest)
               and assign them to threads. */
            
            int srow = rows;
            int scol = pend - pstart;
            char* subset = (char *)malloc(srow * scol * sizeof(char));
            
            for (int r = 0; r < srow; r ++)
                for (int c = 0; c < scol; c ++)
                    subset[r * scol + c] = data[r * columns + pstart + c];
            
            t_params p = { subset, srow, scol };
            params[ithr] = p;
            pthread_create(&threads[ithr], NULL, unique, &params[ithr]);
        }

        /* gather results */

        for (int ithr = thrstart; ithr < thrend; ithr++) {
            pthread_t tid = threads[ithr];
            t_out* tresult = NULL;
            pthread_join(tid, (void**)&tresult);
            results[ithr] = tresult -> unique;
            free(params[ithr].data);
        }
    }

    free(threads);
    free(params);

    return results;
}