#include <R.h>

extern "C" void
es_C_SOM(float *data,
         float *codes,
         float *nhbrdist,
         float *alphasA,
         float *radiiA,
         float *alphasB,
         float *radiiB,
         int *pn,
         int *pdim,
         int *pncodes,
         int *prlen,
         int *dist);

extern "C" void
es_C_BatchSOM(int *pnthreads,
              float *data,
              float *codes,
              float *nhbrdist,
              float *radii,
              int *pn,
              int *pdim,
              int *pncodes,
              int *prlen,
              int *dist);

extern "C" void
es_C_GQTSOM(int *pnthreads,
            float *data,
            int *coords,
            float *codes,
            float *radii,
            int *out_ncodes,
            float *out_codes,
            int *out_coords,
            float *out_emcoords,
            int *pn,
            int *pdim,
            int *pncodes,
            int *prlen,
            int *distf,
            int *nhbr_distf);

extern "C" void
es_C_mapDataToCodes(int *pnthreads,
                    float *data,
                    float *codes,
                    int *pncodes,
                    int *pnd,
                    int *pp,
                    int *nnCodes,
                    float *nnDists,
                    int *dist);
