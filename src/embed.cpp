
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
static const float min_boost = 0.001; // lower limit for the parameter

// this is added before calculating log-distances
static const float zero_avoidance = 0.00000000001;

// a tiny epsilon for preventing singularities
static const float koho_gravity = 0.00001;

// maximal exponent allowed in scores
static const float score_exp_norm = 15;
static const float inv_neg_score_norm = -1 / score_exp_norm;

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

static inline float logscore_normalize (float a)
{
	/* normalize out extreme values. We've only got 8 bits of exponent
	 * precision in floats, so after adding some safety margin we should
	 * fit between 10^-10 and 10^10. That translates to around exp(23) in
	 * base-e as a limit. Therefore we run this through logsig to get it at
	 * least between -20 and 20. */
	return score_exp_norm / (1 + expf (inv_neg_score_norm * a));
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

extern "C" void C_embedSOM (int* pn,
                            int* pdim,
                            float* pboost,
                            int* pneighbors,
                            float* padjust,
                            int* pxdim,
                            int* pydim,
                            float* points,
                            float* koho,
                            float* embedding)
{
	size_t n = *pn, indim = *pdim, topn = *pneighbors, xdim = *pxdim,
	       ydim = *pydim;
	float boost = *pboost;

	size_t i, j, k;

#ifdef DEBUG_CRASH_ON_FPE
	feenableexcept (FE_INVALID | FE_OVERFLOW);
#endif

	if (topn > xdim * ydim) topn = xdim * ydim;
	if (boost < min_boost) boost = min_boost;

	vector<dist_id> dists;
	dists.resize (topn);

	const float score_mod = -float(indim > 1 ? indim - 1 : 1) / (2 * boost);

	float mtx[6];

	float* point = points;
	for (size_t ptid = 0; ptid < n; ++ptid, point += indim) {

		for (i = 0; i < topn; ++i) {
			float s = 0;
			for (k = 0; k < indim; ++k)
				s += sqrf (point[k] - koho[k + i * indim]);
			dists[i].dist = s;
			dists[i].id = i;
		}

		for (i = 0; i < topn; ++i)
			heap_down (dists.data (), topn - i - 1, topn);

		for (i = topn; i < xdim * ydim; ++i) {
			float s = 0;
			for (k = 0; k < indim; ++k)
				s += sqrf (point[k] - koho[k + i * indim]);
			if (dists[0].dist < s) continue;
			dists[0].dist = s;
			dists[0].id = i;
			heap_down (dists.data (), 0, topn);
		}

		for (i = 0; i < topn; ++i)
			dists[i].dist = logf (zero_avoidance + dists[i].dist);

		float sum = 0;
		for (i = 0; i < topn; ++i) sum += dists[i].dist;
		sum /= topn;

		for (i = 0; i < topn; ++i)
			dists[i].dist = expf (logscore_normalize (
			  score_mod * (dists[i].dist - sum)));

		// prepare the matrix for 2x2 linear eqn
		for (i = 0; i < 6; ++i) mtx[i] = 0; // it's stored by columns!

		for (i = 0; i < topn; ++i) {
			// add a really tiny influence of the point to prevent
			// singularities
			size_t idx = dists[i].id;
			float ix = idx % xdim, iy = idx / xdim;
			float pi = dists[i].dist;
			float gs = koho_gravity * dists[i].dist;
			mtx[0] += gs;
			mtx[3] += gs;
			mtx[4] += gs * ix;
			mtx[5] += gs * iy;

			for (j = i + 1; j < topn; ++j) {

				size_t jdx = dists[j].id;
				float jx = jdx % xdim, jy = jdx / xdim;
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

				const float hx = jx - ix;
				const float hy = jy - iy;
				const float hpxy = hx * hx + hy * hy;
				const float ihpxy = 1 / hpxy;

				const float s = pi * pj / powf (hpxy, *padjust);

				const float diag = s * hx * hy * ihpxy;
				const float rhsc =
				  s * (scalar + (hx * ix + hy * iy) * ihpxy);

				mtx[0] += s * hx * hx * ihpxy;
				mtx[1] += diag;
				mtx[2] += diag;
				mtx[3] += s * hy * hy * ihpxy;
				mtx[4] += hx * rhsc;
				mtx[5] += hy * rhsc;
			}
		}

		// cramer
		float det = mtx[0] * mtx[3] - mtx[1] * mtx[2];
		// output is stored R-style by columns
		embedding[ptid] = (mtx[4] * mtx[3] - mtx[5] * mtx[2]) / det;
		embedding[ptid + n] = (mtx[0] * mtx[5] - mtx[1] * mtx[4]) / det;
	}

#ifdef DEBUG_CRASH_ON_FPE
	fedisableexcept (FE_INVALID | FE_OVERFLOW);
#endif
}

#include <R.h>
#include <R_ext/Rdynload.h>

static const R_CMethodDef cMethods[] = {
	{ "C_embedSOM", (DL_FUNC)&C_embedSOM, 10 },
	{ NULL, NULL, 0 }
};

void R_init_EmbedSOM (DllInfo* info)
{
	R_registerRoutines (info, cMethods, NULL, NULL, NULL);
}
