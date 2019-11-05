
/* This file is part of EmbedSOM.
 *
 * Copyright (C) 2018-2019 Mirek Kratochvil <exa.exa@gmail.com>
 *
 * Based on code from FlowSOM,
 * Copyright (C) 2016-2019 Sofie Van Gassen et al.
 *
 * Originally based on code of Ron Wehrens
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

#include "som.h"

#include <thread>
#include <vector>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R_ext/PrtUtil.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#include "distfs.h"
#include "use_intrins.h"

template<class distf>
static void
som(size_t n,
    size_t k,
    size_t dim,
    size_t rlen,
    const float* points,
    float* koho,
    const float* nhbrdist,
    float alphasA[2],
    float radiiA[2],
    float alphasB[2],
    float radiiB[2])
{
	size_t niter = rlen * n;

	float thresholdA0 = radiiA[0], alphaA0 = alphasA[0],
	      thresholdADiff = radiiA[1] - radiiA[0],
	      alphaADiff = alphasA[1] - alphasA[0], thresholdB0 = radiiB[0],
	      alphaB0 = alphasB[0], thresholdBDiff = radiiB[1] - radiiB[0],
	      alphaBDiff = alphasB[1] - alphasB[0];

	RANDIN;

	for (size_t iter = 0; iter < niter; ++iter) {
		size_t point = size_t(n * UNIF);
		if (point >= n)
			point = 0; // careful there.

		float riter = iter / (float)niter;

		size_t nearest = 0;
		{
			float nearestd =
			  distf::comp(points + dim * point, koho, dim);
			for (size_t i = 1; i < k; ++i) {
				float tmp = distf::comp(
				  points + dim * point, koho + dim * i, dim);
				if (tmp < nearestd) {
					nearest = i;
					nearestd = tmp;
				}
			}
		}

		float thresholdA = thresholdA0 + riter * thresholdADiff,
		      thresholdB = thresholdB0 + riter * thresholdBDiff,
		      alphaA = alphaA0 + riter * alphaADiff,
		      alphaB = alphaB0 + riter * alphaBDiff;

		for (size_t i = 0; i < k; ++i) {
			float d = nhbrdist[i + k * nearest];

			float alpha;

			if (d > thresholdA) {
				if (d > thresholdB)
					continue;
				alpha = alphaB;
			} else
				alpha = alphaA;

#ifndef USE_INTRINS
			for (size_t j = 0; j < dim; ++j)
				koho[j + i * dim] +=
				  alpha *
				  (points[j + point * dim] - koho[j + i * dim]);
#else
			float *ki = koho + (i * dim), *ke = ki + dim,
			      *kie = ke - 3;
			const float* pi = points + (point * dim);
			__m128 al = _mm_set1_ps(alpha);
			__m128 k = _mm_loadu_ps(ki);
			for (; ki < kie; ki += 4, pi += 4) {
				__m128 kn = _mm_loadu_ps(ki + 4);
				_mm_storeu_ps(
				  ki,
				  _mm_add_ps(
				    k,
				    _mm_mul_ps(
				      al, _mm_sub_ps(_mm_loadu_ps(pi), k))));
				k = kn;
			}
			for (; ki < ke; ++ki, ++pi)
				*ki += alpha * (*pi - *ki);
#endif
		}
	}

	RANDOUT;
}

template<class distf, bool threaded>
static void
mapNNs(size_t threads,
       size_t n,
       size_t k,
       size_t dim,
       const float* points,
       const float* koho,
       int* mapping,
       float* dists)
{
	if (threaded) {
		std::vector<std::thread> ts(threads);
		for (size_t i = 0; i < threads; ++i)
			ts[i] = std::thread(
			  [&](size_t thread_id) {
				  size_t dbegin = thread_id * n / threads,
				         dend = (thread_id + 1) * n / threads;
				  const float* d = points + dbegin * dim;
				  int* m = mapping + dbegin;
				  float* dd = dists + dbegin;
				  size_t nd = dend - dbegin;

				  mapNNs<distf, false>(
				    1, nd, k, dim, d, koho, m, dd);
			  },
			  i);
		for (auto& t : ts)
			t.join();
		return;
	}

	for (size_t point = 0; point < n; ++point) {
		size_t nearest = 0;
		float nearestd = distf::comp(points + dim * point, koho, dim);
		for (size_t i = 1; i < k; ++i) {
			float tmp = distf::comp(
			  points + dim * point, koho + dim * i, dim);
			if (tmp < nearestd) {
				nearest = i;
				nearestd = tmp;
			}
		}

		mapping[point] = nearest;
		dists[point] = distf::back(nearestd);
	}
}

extern "C" void
es_C_SOM(float* points,
         float* koho,
         float* nhbrdist,
         float* alphasA,
         float* radiiA,
         float* alphasB,
         float* radiiB,
         Sint* pn,
         Sint* pdim,
         Sint* pkohos,
         Sint* prlen,
         Sint* dist)
{
	int n = *pn, dim = *pdim, kohos = *pkohos, rlen = *prlen;

	auto somf = som<distfs::sqeucl>;
	if (*dist == 1)
		somf = som<distfs::manh>;
	else if (*dist == 3)
		somf = som<distfs::chebyshev>;

	somf(n,
	     kohos,
	     dim,
	     rlen,
	     points,
	     koho,
	     nhbrdist,
	     alphasA,
	     radiiA,
	     alphasB,
	     radiiB);
}

extern "C" void
es_C_mapDataToCodes(int* pnthreads,
                    float* points,
                    float* koho,
                    int* pn,
                    int* pdim,
                    int* pkohos,
                    int* mapping,
                    float* dists,
                    int* dist)
{
	int n = *pn, dim = *pdim, kohos = *pkohos;

	int threads = *pnthreads;

	if (threads < 0)
		threads = 1;
	if (threads == 0)
		threads = std::thread::hardware_concurrency();

	auto mapf = threads == 1 ? mapNNs<distfs::sqeucl, false>
	                         : mapNNs<distfs::sqeucl, true>;
	if (*dist == 1)
		mapf = threads == 1 ? mapNNs<distfs::manh, false>
		                    : mapNNs<distfs::manh, true>;
	else if (*dist == 3)
		mapf = threads == 1 ? mapNNs<distfs::chebyshev, false>
		                    : mapNNs<distfs::chebyshev, true>;

	mapf(threads, n, kohos, dim, points, koho, mapping, dists);
}
