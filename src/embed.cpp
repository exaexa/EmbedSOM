
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
 */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include <iostream>

// this helps with debugging floating-point overflows and similar nastiness,
// uncomment if needed.
//#define DEBUG_CRASH_ON_FPE

#ifdef DEBUG_CRASH_ON_FPE
#include <fenv.h>
#endif

using namespace std;

// some small numbers first!
static const float min_boost = 0.00001; // lower limit for the parameter

// this is added before normalizing the distances
static const float zero_avoidance = 0.00000000001;

// a tiny epsilon for preventing singularities
static const float koho_gravity = 0.00001;

static inline float sqrf (float n)
{
	return n * n;
}

struct dist_id {
	float dist;
	size_t id;
};

static inline void hswap (dist_id& a, dist_id& b)
{
	dist_id c = a;
	a = b;
	b = c;
}

static void heap_down (dist_id* heap, size_t start, size_t lim)
{
	for (;;) {
		size_t L = 2 * start + 1;
		size_t R = L + 1;
		if (R < lim) {
			float dl = heap[L].dist;
			float dr = heap[R].dist;

			if (dl > dr) {
				if (heap[start].dist >= dl) break;
				hswap (heap[L], heap[start]);
				start = L;
			} else {
				if (heap[start].dist >= dr) break;
				hswap (heap[R], heap[start]);
				start = R;
			}
		} else if (L < lim) {
			if (heap[start].dist < heap[L].dist)
				hswap (heap[L], heap[start]);
			break; // exit safely!
		} else
			break;
	}
}

extern "C" void C_embedSOM (int* psomdim,
                            int* pn,
                            int* pdim,
                            float* pboost,
                            int* pneighbors,
                            float* padjust,
                            int* pncodes,
                            float* points,
                            float* koho,
                            float* emcoords,
                            float* embedding)
{
	size_t topn = *pneighbors;
	const int somdim = *psomdim;
	const size_t n = *pn, indim = *pdim, ncodes = *pncodes;
	float boost = *pboost;

	size_t i, j, k;

#ifdef DEBUG_CRASH_ON_FPE
	feenableexcept (FE_INVALID | FE_OVERFLOW);
#endif

	if (topn > ncodes) topn = ncodes;
	if (boost < min_boost) boost = min_boost;

	vector<dist_id> dists;
	dists.resize (topn);

	float mtx[12]; // only 6 in case of psomdim=2 but who cares

	float* point = points;
	for (size_t ptid = 0; ptid < n; ++ptid, point += indim) {

		// heap-knn
		for (i = 0; i < topn; ++i) {
			float s = 0;
			for (k = 0; k < indim; ++k)
				s += sqrf (point[k] - koho[k + i * indim]);
			dists[i].dist = s;
			dists[i].id = i;
		}

		for (i = 0; i < topn; ++i)
			heap_down (dists.data (), topn - i - 1, topn);

		for (i = topn; i < ncodes; ++i) {
			float s = 0;
			for (k = 0; k < indim; ++k)
				s += sqrf (point[k] - koho[k + i * indim]);
			if (dists[0].dist < s) continue;
			dists[0].dist = s;
			dists[0].id = i;
			heap_down (dists.data (), 0, topn);
		}

		// heapsort the result
		for (i = topn - 1; i > 0; --i) {
			hswap (dists[0], dists[i]);
			heap_down (dists.data (), 0, i);
		}

		// compute scores
		float sum = 0, ssum = 0, min = dists[0].dist;
		for (i = 0; i < topn; ++i) {
			dists[i].dist = sqrtf (dists[i].dist);
			sum += dists[i].dist / (i + 1);
			ssum += 1 / float(i + 1);
			if (dists[i].dist < min) min = dists[i].dist;
		}

		sum = -ssum / (zero_avoidance + sum * boost);

		for (i = 0; i < topn; ++i)
			dists[i].dist = expf ((dists[i].dist - min) * sum);

		// prepare the matrix for 2x2 linear eqn
		if (somdim == 2)
			for (i = 0; i < 6; ++i)
				mtx[i] = 0; // it's stored by columns!
		if (somdim == 3)
			for (i = 0; i < 12; ++i) mtx[i] = 0;

		for (i = 0; i < topn; ++i) {
			// add a really tiny influence of the point to prevent
			// singularities
			size_t idx = dists[i].id;
			float ix, iy, iz;
			if (somdim == 2) {
				ix = emcoords[2 * idx + 0];
				iy = emcoords[2 * idx + 1];
			}
			if (somdim == 3) {
				ix = emcoords[3 * idx + 0];
				iy = emcoords[3 * idx + 1];
				iz = emcoords[3 * idx + 2];
			}
			float pi = dists[i].dist;
			float gs = koho_gravity * dists[i].dist;
			if (somdim == 2) {
				mtx[0] += gs;
				mtx[3] += gs;
				mtx[4] += gs * ix;
				mtx[5] += gs * iy;
			}
			if (somdim == 3) {
				mtx[0] += gs;
				mtx[4] += gs;
				mtx[8] += gs;
				mtx[9] += gs * ix;
				mtx[10] += gs * iy;
				mtx[11] += gs * iz;
			}

			for (j = i + 1; j < topn; ++j) {

				size_t jdx = dists[j].id;
				float jx, jy, jz;
				if (somdim == 2) {
					jx = emcoords[2 * jdx + 0];
					jy = emcoords[2 * jdx + 1];
				}
				if (somdim == 3) {
					jx = emcoords[3 * jdx + 0];
					jy = emcoords[3 * jdx + 1];
					jz = emcoords[3 * jdx + 2];
				}
				float pj = dists[j].dist;

				float scalar = 0, sqdist = 0;
				for (k = 0; k < indim; ++k) {
					float tmp = koho[k + indim * jdx] -
					            koho[k + indim * idx];
					sqdist += tmp * tmp;
					scalar += tmp * (point[k] -
					                 koho[k + indim * idx]);
				}

				if (scalar != 0) {
					if (sqdist == 0)
						continue;
					else
						scalar /= sqdist;
				}

				if (somdim == 2) {
					const float hx = jx - ix;
					const float hy = jy - iy;
					const float hpxy = hx * hx + hy * hy;
					const float ihpxy = 1 / hpxy;

					const float s =
					  pi * pj / powf (hpxy, *padjust);

					const float diag = s * hx * hy * ihpxy;
					const float rhsc =
					  s * (scalar +
					       (hx * ix + hy * iy) * ihpxy);

					mtx[0] += s * hx * hx * ihpxy;
					mtx[1] += diag;
					mtx[2] += diag;
					mtx[3] += s * hy * hy * ihpxy;
					mtx[4] += hx * rhsc;
					mtx[5] += hy * rhsc;
				}

				if (somdim == 3) {
					const float hx = jx - ix;
					const float hy = jy - iy;
					const float hz = jz - iz;
					const float hpxyz =
					  hx * hx + hy * hy + hz * hz;
					const float s =
					  pi * pj / powf (hpxyz, *padjust);
					const float ihpxyz = 1 / hpxyz;
					const float sihpxyz = s * ihpxyz;

					const float rhsc =
					  s * (scalar +
					       (hx * ix + hy * iy + hz * iz) *
					         ihpxyz);

					mtx[0] += sihpxyz * hx * hx;
					mtx[1] += sihpxyz * hx * hy;
					mtx[2] += sihpxyz * hx * hz;
					mtx[3] += sihpxyz * hy * hx;
					mtx[4] += sihpxyz * hy * hy;
					mtx[5] += sihpxyz * hy * hz;
					mtx[6] += sihpxyz * hz * hx;
					mtx[7] += sihpxyz * hz * hy;
					mtx[8] += sihpxyz * hz * hz;
					mtx[9] += hx * rhsc;
					mtx[10] += hy * rhsc;
					mtx[11] += hz * rhsc;
				}
			}
		}

		// cramer (output is stored R-style by columns)
		if (somdim == 2) {
			float det = mtx[0] * mtx[3] - mtx[1] * mtx[2];
			embedding[ptid] =
			  (mtx[4] * mtx[3] - mtx[5] * mtx[2]) / det;
			embedding[ptid + n] =
			  (mtx[0] * mtx[5] - mtx[1] * mtx[4]) / det;
		}
		if (somdim == 3) {
			float det =
			  mtx[0] * mtx[4] * mtx[8] + mtx[1] * mtx[5] * mtx[6] +
			  mtx[2] * mtx[3] * mtx[7] - mtx[0] * mtx[5] * mtx[7] -
			  mtx[1] * mtx[3] * mtx[8] - mtx[2] * mtx[4] * mtx[6];
			embedding[ptid] = (mtx[9] * mtx[4] * mtx[8] +
			                   mtx[10] * mtx[5] * mtx[6] +
			                   mtx[11] * mtx[3] * mtx[7] -
			                   mtx[9] * mtx[5] * mtx[7] -
			                   mtx[10] * mtx[3] * mtx[8] -
			                   mtx[11] * mtx[4] * mtx[6]) /
			                  det;
			embedding[ptid + n] = (mtx[0] * mtx[10] * mtx[8] +
			                       mtx[1] * mtx[11] * mtx[6] +
			                       mtx[2] * mtx[9] * mtx[7] -
			                       mtx[0] * mtx[11] * mtx[7] -
			                       mtx[1] * mtx[9] * mtx[8] -
			                       mtx[2] * mtx[10] * mtx[6]) /
			                      det;
			embedding[ptid + 2 * n] = (mtx[0] * mtx[4] * mtx[11] +
			                           mtx[1] * mtx[5] * mtx[9] +
			                           mtx[2] * mtx[3] * mtx[10] -
			                           mtx[0] * mtx[5] * mtx[10] -
			                           mtx[1] * mtx[3] * mtx[11] -
			                           mtx[2] * mtx[4] * mtx[9]) /
			                          det;
		}
	}

#ifdef DEBUG_CRASH_ON_FPE
	fedisableexcept (FE_INVALID | FE_OVERFLOW);
#endif
}

#include "som.h"
#include <R.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef cMethods[] = {
	{ "C_embedSOM", (DL_FUNC)&C_embedSOM, 11 },
	{ "es_C_SOM", (DL_FUNC)&es_C_SOM, 12 },
	{ "es_C_mapDataToCodes", (DL_FUNC)&es_C_mapDataToCodes, 8 },
	{ NULL, NULL, 0 }
};

void R_init_EmbedSOM (DllInfo* info)
{
	R_registerRoutines (info, cMethods, NULL, NULL, NULL);
}
