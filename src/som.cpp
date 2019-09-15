
/* This file is part of EmbedSOM.
 *
 * EmbedSOM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * EmbedSOM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * EmbedSOM. If not, see <https://www.gnu.org/licenses/>.
 *
 * Taken from FlowSOM, Copyright (C) Sofie Van Gassen (2006-)
 * Originally based on code of Ron Wehrens
 *
 * TODO: convert all ints that should be size_t to size_t
 */

#include "som.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R_ext/PrtUtil.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#define RANDIN GetRNGstate ()
#define RANDOUT PutRNGstate ()
#define UNIF unif_rand ()

#include "distfs.h"

extern "C" void es_C_SOM (float* data,
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
                          Sint* dist)
{
	int n = *pn, px = *ppx, ncodes = *pncodes, rlen = *prlen;
	int cd, i, j, k, niter;
	float tmp;
	float (*distf) (float*, float*, int, int, int);

	if (*dist == 1) {
		distf = manh;
	} else if (*dist == 2) {
		distf = sqeucl;
	} else if (*dist == 3) {
		distf = chebyshev;
	} else {
		distf = sqeucl;
	}

	RANDIN;
	niter = rlen * n;

	float thresholdA0 = radiiA[0], alphaA0 = alphasA[0],
	      thresholdADiff = radiiA[1] - radiiA[0],
	      alphaADiff = alphasA[1] - alphasA[0], thresholdB0 = radiiB[0],
	      alphaB0 = alphasB[0], thresholdBDiff = radiiB[1] - radiiB[0],
	      alphaBDiff = alphasB[1] - alphasB[0];

	for (k = 0; k < niter; k++) {
		/* i is a counter over objects in data, cd is a counter over
		   units in the map, and j is a counter over variables */
		i = (int)(n * UNIF); /* Select a random sample */

		/* calculate distances in x and y spaces, and keep track of the
		   nearest code */
		int nearest = 0;
		float nearestd = distf (data + i, codes, px, n, ncodes);
		for (cd = 1; cd < ncodes; cd++) {
			tmp = distf (data + i, codes + cd, px, n, ncodes);
			if (tmp < nearestd) {
				nearest = cd;
				nearestd = tmp;
			}
		}

		float thresholdA = thresholdA0 + k * thresholdADiff / niter,
		      thresholdB = thresholdB0 + k * thresholdBDiff / niter;

		for (cd = 0; cd < ncodes; cd++) {
			float d = nhbrdist[cd + ncodes * nearest];
			if (d > thresholdB) continue;

			float alpha;

			if (d > thresholdA)
				alpha = alphaB0 + k * alphaBDiff / niter;
			else
				alpha = alphaA0 + k * alphaADiff / niter;

			for (j = 0; j < px; j++) {
				tmp = data[i + j * n] - codes[cd + j * ncodes];
				codes[cd + j * ncodes] += tmp * alpha;
			}
		}
	}

	RANDOUT;
}

extern "C" void es_C_mapDataToCodes (float* data,
                                     float* codes,
                                     int* pncodes,
                                     int* pnd,
                                     int* pp,
                                     int* nnCodes,
                                     float* nnDists,
                                     int* dist)
{
	int ncodes = *pncodes, nd = *pnd, p = *pp;
	int i, cd, counter, minid;
	float tmp, mindist;
	float (*distf) (float*, float*, int, int, int);

	if (*dist == 1) {
		distf = manh;
	} else if (*dist == 2) {
		distf = sqeucl;
	} else if (*dist == 3) {
		distf = chebyshev;
	} else {
		distf = sqeucl;
	}

	/* i is a counter over objects in data, cd  is a counter over SOM
	   units, p is the number of columns, nd is the number of datapoints
	   and ncodes is the number of SOM units*/
	counter = -1;
	for (i = 0; i < nd; i++) {
		minid = -1;
		mindist = DBL_MAX;
		for (cd = 0; cd < ncodes; cd++) {
			tmp = distf (&data[i], &codes[cd], p, nd, ncodes);
			// Rprintf("\ndist: %f",tmp2);
			if (tmp < mindist) {
				mindist = tmp;
				minid = cd;
			}
		}
		nnCodes[++counter] = minid + 1;
		nnDists[counter] = (*dist == 2) ? sqrtf (mindist) : mindist;
	}
}
