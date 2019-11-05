
/* This file is part of EmbedSOM.
 *
 * Copyright (C) 2018-2019 Mirek Kratochvil <exa.exa@gmail.com>
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
#include <thread>
#include <vector>

// this helps with debugging floating-point overflows and similar nastiness,
// uncomment if needed.
//#define DEBUG_CRASH_ON_FPE

#include "distfs.h"

#ifdef DEBUG_CRASH_ON_FPE
#include <fenv.h>
#endif

using namespace std;

// some small numbers first!
static const float min_boost = 1e-5; // lower limit for the parameter

// this affects how steeply the score decreases for near-farthest codes
static const float max_avoidance = 10;

// this is added before normalizing the distances
static const float zero_avoidance = 1e-10;

// a tiny epsilon for preventing singularities
static const float koho_gravity = 1e-5;

static inline float
sqrf(float n)
{
	return n * n;
}

struct dist_id
{
	float dist;
	size_t id;
};

static inline void
hswap(dist_id& a, dist_id& b)
{
	dist_id c = a;
	a = b;
	b = c;
}

static void
heap_down(dist_id* heap, size_t start, size_t lim)
{
	for (;;) {
		size_t L = 2 * start + 1;
		size_t R = L + 1;
		if (R < lim) {
			float dl = heap[L].dist;
			float dr = heap[R].dist;

			if (dl > dr) {
				if (heap[start].dist >= dl)
					break;
				hswap(heap[L], heap[start]);
				start = L;
			} else {
				if (heap[start].dist >= dr)
					break;
				hswap(heap[R], heap[start]);
				start = R;
			}
		} else if (L < lim) {
			if (heap[start].dist < heap[L].dist)
				hswap(heap[L], heap[start]);
			break; // exit safely!
		} else
			break;
	}
}

template<class distf, int embed_dim, bool threaded>
static void
embedsom(const size_t threads,
         const size_t n,
         const size_t kohos,
         const size_t dim,
         const float boost,
         const size_t topn,
         const float adjust,
         const float* points,
         const float* koho,
         const float* emcoords,
         float* embedding)
{
	if (threaded) {
		vector<thread> ts(threads);
		for (size_t i = 0; i < threads; ++i)
			ts[i] = thread(
			  [&](size_t thread_id) {
				  size_t dbegin = thread_id * n / threads,
				         dend = (thread_id + 1) * n / threads;
				  const float* d = points + dbegin * dim;
				  float* e = embedding + dbegin * embed_dim;
				  size_t nd = dend - dbegin;

				  embedsom<distf, embed_dim, false>(1,
				                                    nd,
				                                    kohos,
				                                    dim,
				                                    boost,
				                                    topn,
				                                    adjust,
				                                    d,
				                                    koho,
				                                    emcoords,
				                                    e);
			  },
			  i);
		for (auto& t : ts)
			t.join();
		return;
	}

	size_t i, j, k;
	const size_t topnn = topn < kohos ? topn + 1 : topn;

#ifdef DEBUG_CRASH_ON_FPE
	feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

	vector<dist_id> dists;
	dists.resize(topnn);

	float mtx[embed_dim * (1 + embed_dim)];

	for (size_t ptid = 0; ptid < n; ++ptid) {
		const float* point = points + dim * ptid;

		// heap-knn
		for (i = 0; i < topnn; ++i) {
			dists[i].dist = distf::comp(point, koho + i * dim, dim);
			dists[i].id = i;
		}

		for (i = 0; i < topnn; ++i)
			heap_down(dists.data(), topnn - i - 1, topnn);

		for (i = topnn; i < kohos; ++i) {
			float s = distf::comp(point, koho + i * dim, dim);
			if (dists[0].dist < s)
				continue;
			dists[0].dist = s;
			dists[0].id = i;
			heap_down(dists.data(), 0, topnn);
		}

		// heapsort the result
		for (i = topnn - 1; i > 0; --i) {
			hswap(dists[0], dists[i]);
			heap_down(dists.data(), 0, i);
		}

		// compute distance distribution (assume normal b/c why not)
		float mean = 0, sd = 0, wsum = 0;
		for (i = 0; i < topnn; ++i) {
			const float tmp = distf::back(dists[i].dist);
			const float w = 1 / float(i + 1);
			mean += tmp * w;
			sd += tmp * tmp * w;
			wsum += w;
			dists[i].dist = tmp;
		}

		mean /= wsum;
		sd = boost / sqrtf(sd / wsum - sqrf(mean));
		const float nmax = max_avoidance / dists[topnn - 1].dist;

		// get the score!
		for (i = 0; i < topn; ++i)
			if (topn < topnn)
				dists[i].dist =
				  expf((mean - dists[i].dist) * sd) *
				  (1 -
				   expf(dists[i].dist * nmax - max_avoidance));
			else
				dists[i].dist =
				  expf((mean - dists[i].dist) * sd);

		// prepare the eqn matrix
		for (i = 0; i < embed_dim * (1 + embed_dim); ++i)
			mtx[i] = 0;

		for (i = 0; i < topn; ++i) {
			// add a really tiny influence of the point to prevent
			// singularities
			size_t idx = dists[i].id;
			float ix, iy, iz;
			if (embed_dim == 2) {
				ix = emcoords[2 * idx + 0];
				iy = emcoords[2 * idx + 1];
			}
			if (embed_dim == 3) {
				ix = emcoords[3 * idx + 0];
				iy = emcoords[3 * idx + 1];
				iz = emcoords[3 * idx + 2];
			}
			float pi = dists[i].dist;
			float gs = koho_gravity * dists[i].dist;
			if (embed_dim == 2) {
				mtx[0] += gs;
				mtx[3] += gs;
				mtx[4] += gs * ix;
				mtx[5] += gs * iy;
			}
			if (embed_dim == 3) {
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
				if (embed_dim == 2) {
					jx = emcoords[2 * jdx + 0];
					jy = emcoords[2 * jdx + 1];
				}
				if (embed_dim == 3) {
					jx = emcoords[3 * jdx + 0];
					jy = emcoords[3 * jdx + 1];
					jz = emcoords[3 * jdx + 2];
				}
				float pj = dists[j].dist;

#ifndef USE_INTRINS
				float scalar = 0, sqdist = 0;
				for (k = 0; k < dim; ++k) {
					float tmp = koho[k + dim * jdx] -
					            koho[k + dim * idx];
					sqdist += tmp * tmp;
					scalar += tmp * (point[k] -
					                 koho[k + dim * idx]);
				}
#else
				float scalar, sqdist;
				{
					const float *ki = koho + dim * idx,
					            *kj = koho + dim * jdx,
					            *pp = point, *ke = ki + dim,
					            *kie = ke - 3;

					__m128 sca = _mm_setzero_ps(),
					       sqd = _mm_setzero_ps();
					for (; ki < kie;
					     ki += 4, kj += 4, pp += 4) {
						__m128 ti = _mm_loadu_ps(ki);
						__m128 tmp = _mm_sub_ps(
						  _mm_loadu_ps(kj), ti);
						sqd = _mm_add_ps(
						  sqd, _mm_mul_ps(tmp, tmp));
						sca = _mm_add_ps(
						  sca,
						  _mm_mul_ps(
						    tmp,
						    _mm_sub_ps(_mm_loadu_ps(pp),
						               ti)));
					}
					scalar =
					  sca[0] + sca[1] + sca[2] + sca[3];
					sqdist =
					  sqd[0] + sqd[1] + sqd[2] + sqd[3];
					for (; ki < ke; ++ki, ++kj, ++pp) {
						float tmp = *kj - *ki;
						sqdist += tmp * tmp;
						scalar += tmp * (*pp - *ki);
					}
				}
#endif

				if (sqdist == 0)
					continue;
				else
					scalar /= sqdist;

				if (embed_dim == 2) {
					const float hx = jx - ix;
					const float hy = jy - iy;
					const float hpxy = hx * hx + hy * hy;
					if (hpxy < zero_avoidance)
						continue;
					const float ihpxy = 1 / hpxy;

					const float s =
					  pi * pj * powf(1 + hpxy, -adjust) *
					  expf(-sqrf(scalar - .5));
					const float sihpxy = s * ihpxy;

					const float diag = s * hx * hy * ihpxy;
					const float rhsc =
					  s * (scalar +
					       (hx * ix + hy * iy) * ihpxy);

					mtx[0] += hx * hx * sihpxy;
					mtx[1] += hx * hy * sihpxy;
					mtx[2] += hy * hx * sihpxy;
					mtx[3] += hy * hy * sihpxy;
					mtx[4] += hx * rhsc;
					mtx[5] += hy * rhsc;
				}

				if (embed_dim == 3) {
					const float hx = jx - ix;
					const float hy = jy - iy;
					const float hz = jz - iz;
					const float hpxyz =
					  hx * hx + hy * hy + hz * hz;
					if (hpxyz < zero_avoidance)
						continue;
					const float ihpxyz = 1 / hpxyz;

					const float s =
					  pi * pj * powf(1 + hpxyz, -adjust) *
					  expf(-sqrf(scalar - .5));
					const float sihpxyz = s * ihpxyz;

					const float rhsc =
					  s * (scalar +
					       (hx * ix + hy * iy + hz * iz) *
					         ihpxyz);

					mtx[0] += hx * hx * sihpxyz;
					mtx[1] += hx * hy * sihpxyz;
					mtx[2] += hx * hz * sihpxyz;
					mtx[3] += hy * hx * sihpxyz;
					mtx[4] += hy * hy * sihpxyz;
					mtx[5] += hy * hz * sihpxyz;
					mtx[6] += hz * hx * sihpxyz;
					mtx[7] += hz * hy * sihpxyz;
					mtx[8] += hz * hz * sihpxyz;
					mtx[9] += hx * rhsc;
					mtx[10] += hy * rhsc;
					mtx[11] += hz * rhsc;
				}
			}
		}

		// cramer (output is stored R-style by columns)
		if (embed_dim == 2) {
			float det = mtx[0] * mtx[3] - mtx[1] * mtx[2];
			embedding[ptid * 2] =
			  (mtx[4] * mtx[3] - mtx[5] * mtx[2]) / det;
			embedding[ptid * 2 + 1] =
			  (mtx[0] * mtx[5] - mtx[1] * mtx[4]) / det;
		}
		if (embed_dim == 3) {
			float det =
			  mtx[0] * mtx[4] * mtx[8] + mtx[1] * mtx[5] * mtx[6] +
			  mtx[2] * mtx[3] * mtx[7] - mtx[0] * mtx[5] * mtx[7] -
			  mtx[1] * mtx[3] * mtx[8] - mtx[2] * mtx[4] * mtx[6];
			embedding[ptid * 3] = (mtx[9] * mtx[4] * mtx[8] +
			                       mtx[10] * mtx[5] * mtx[6] +
			                       mtx[11] * mtx[3] * mtx[7] -
			                       mtx[9] * mtx[5] * mtx[7] -
			                       mtx[10] * mtx[3] * mtx[8] -
			                       mtx[11] * mtx[4] * mtx[6]) /
			                      det;
			embedding[ptid * 3 + 1] = (mtx[0] * mtx[10] * mtx[8] +
			                           mtx[1] * mtx[11] * mtx[6] +
			                           mtx[2] * mtx[9] * mtx[7] -
			                           mtx[0] * mtx[11] * mtx[7] -
			                           mtx[1] * mtx[9] * mtx[8] -
			                           mtx[2] * mtx[10] * mtx[6]) /
			                          det;
			embedding[ptid * 3 + 2] = (mtx[0] * mtx[4] * mtx[11] +
			                           mtx[1] * mtx[5] * mtx[9] +
			                           mtx[2] * mtx[3] * mtx[10] -
			                           mtx[0] * mtx[5] * mtx[10] -
			                           mtx[1] * mtx[3] * mtx[11] -
			                           mtx[2] * mtx[4] * mtx[9]) /
			                          det;
		}
	}

#ifdef DEBUG_CRASH_ON_FPE
	fedisableexcept(FE_INVALID | FE_OVERFLOW);
#endif
}

extern "C" void
C_embedSOM(int* pnthreads,
           int* pedim,
           int* pn,
           int* pkohos,
           int* pdim,
           int* pdist,
           float* pboost,
           int* pneighbors,
           float* padjust,
           float* points,
           float* koho,
           float* emcoords,
           float* embedding)
{
	int embeddim = *pedim;
	size_t n = *pn, dim = *pdim, kohos = *pkohos;

	int threads = *pnthreads;

	if (threads < 0)
		threads = 1;
	if (threads == 0)
		threads = thread::hardware_concurrency();

	auto emf = threads == 1 ? embedsom<distfs::sqeucl, 2, false>
	                        : embedsom<distfs::sqeucl, 2, true>;
	if (threads == 1) {
		if (embeddim == 2) {
			if (*pdist == 1)
				emf = embedsom<distfs::manh, 2, false>;
			if (*pdist == 3)
				emf = embedsom<distfs::chebyshev, 2, false>;
		} else if (embeddim == 3) {
			emf = embedsom<distfs::sqeucl, 3, false>;
			if (*pdist == 1)
				emf = embedsom<distfs::manh, 3, false>;
			if (*pdist == 3)
				emf = embedsom<distfs::chebyshev, 3, false>;
		} else
			return; // waat.
	} else {
		if (embeddim == 2) {
			if (*pdist == 1)
				emf = embedsom<distfs::manh, 2, true>;
			if (*pdist == 3)
				emf = embedsom<distfs::chebyshev, 2, true>;
		} else if (embeddim == 3) {
			emf = embedsom<distfs::sqeucl, 3, true>;
			if (*pdist == 1)
				emf = embedsom<distfs::manh, 3, true>;
			if (*pdist == 3)
				emf = embedsom<distfs::chebyshev, 3, true>;
		} else
			return; // waat.
	}

	size_t topn = *pneighbors;
	if (topn > kohos)
		topn = kohos;
	float boost = *pboost;
	if (boost < min_boost)
		boost = min_boost;
	emf(threads,
	    n,
	    kohos,
	    dim,
	    *pboost,
	    topn,
	    *padjust,
	    points,
	    koho,
	    emcoords,
	    embedding);
}

#include "som.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

static const R_CMethodDef cMethods[] = {
	{ "C_embedSOM", (DL_FUNC)&C_embedSOM, 13 },
	{ "es_C_SOM", (DL_FUNC)&es_C_SOM, 12 },
	{ "es_C_mapDataToCodes", (DL_FUNC)&es_C_mapDataToCodes, 8 },
	{ NULL, NULL, 0 }
};

void
R_init_EmbedSOM(DllInfo* info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}
