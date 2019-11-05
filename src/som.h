#include <R.h>

extern "C" void
es_C_SOM(float* data,
         float* codes,
         float* nhbrdist,
         float* alphasA,
         float* radiiA,
         float* alphasB,
         float* radiiB,
         Sint* pn,
         Sint* ppx,
         Sint* pncodes,
         Sint* prlen,
         Sint* dist);

extern "C" void
es_C_BatchSOM(int* pnthreads,
              float* data,
              float* codes,
              float* nhbrdist,
              float* radii,
              Sint* pn,
              Sint* ppx,
              Sint* pncodes,
              Sint* prlen,
              Sint* dist);

extern "C" void
es_C_mapDataToCodes(int* pnthreads,
                    float* data,
                    float* codes,
                    int* pncodes,
                    int* pnd,
                    int* pp,
                    int* nnCodes,
                    float* nnDists,
                    int* dist);
