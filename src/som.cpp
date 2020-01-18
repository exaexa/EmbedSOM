
/* This file is part of EmbedSOM.
 *
 * Copyright (C) 2018-2020 Mirek Kratochvil <exa.exa@gmail.com>
 *
 * Parts of the code are based on FlowSOM,
 * Copyright (C) 2016-2019 Sofie Van Gassen et al.
 * (These were originally based on code of Ron Wehrens.)
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

#include <complex>
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
    const float *points,
    float *koho,
    const float *nhbrdist,
    const float alphasA[2],
    const float radiiA[2],
    const float alphasB[2],
    const float radiiB[2])
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
			const float *pi = points + (point * dim);
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

constexpr float min_radius = 1e-10;

template<class distf, bool threaded>
static void
bsom(size_t threads,
     size_t n,
     size_t kohos,
     size_t dim,
     size_t epochs,
     const float *points,
     float *koho,
     const float *nhbrdist,
     const float *radii)
{
	std::vector<std::thread> ts(threads);
	std::vector<std::vector<float>> taccs, tws;
	taccs.resize(threads);
	for (auto &acc : taccs)
		acc.resize(kohos * dim);
	tws.resize(threads);
	for (auto &ws : tws)
		ws.resize(kohos);
	std::vector<float> weights(kohos);
	std::vector<float> oldkoho(kohos * dim);

	for (size_t epoch = 0; epoch < epochs; ++epoch) {
		auto reduce_func = [&](size_t thread_id) {
			std::vector<float> &ws = tws[thread_id];
			std::vector<float> &acc = taccs[thread_id];
			size_t dbegin = thread_id * n / threads,
			       dend = (thread_id + 1) * n / threads;
			const float *d = points + dbegin * dim;
			size_t nd = dend - dbegin;
			for (auto &x : acc)
				x = 0;
			for (auto &x : ws)
				x = 0;

			for (size_t i = 0; i < nd; ++i) {
				size_t closest = 0;
				float closestd =
				  distf::comp(d + i * dim, koho, dim);
				for (size_t j = 1; j < kohos; ++j) {
					float tmp = distf::comp(
					  d + i * dim, koho + j * dim, dim);
					if (tmp < closestd) {
						closest = j;
						closestd = tmp;
					}
				}

				ws[closest] += 1;
				for (size_t k = 0; k < dim; ++k)
					acc[closest * dim + k] +=
					  d[i * dim + k];
			}
		};

		if (threaded) {
			for (size_t i = 0; i < threads; ++i)
				ts[i] = std::thread(reduce_func, i);
			for (auto &t : ts)
				t.join();
			for (size_t i = 1; i < threads; ++i)
				for (size_t j = 0; j < kohos * dim; ++j)
					taccs[0][j] += taccs[i][j];
			for (size_t i = 1; i < threads; ++i)
				for (size_t j = 0; j < kohos; ++j)
					tws[0][j] += tws[i][j];

		} else
			reduce_func(0);

		std::copy_n(koho, kohos * dim, oldkoho.data());
		std::fill_n(koho, kohos * dim, 0);
		std::fill_n(weights.data(), kohos, 0);

		const float invSqSigma =
		  -powf(std::max(min_radius, radii[epoch]), -2);

		for (size_t si = 0; si < kohos; ++si)
			for (size_t di = 0; di < kohos; ++di) {
				float w = exp(sqrf(nhbrdist[si * kohos + di]) *
				              invSqSigma);
				for (size_t k = 0; k < dim; ++k)
					koho[di * dim + k] +=
					  taccs[0][si * dim + k] * w;
				weights[di] += tws[0][si] * w;
			}

		for (size_t i = 0; i < kohos; ++i)
			if (weights[i] > 0)
				for (size_t k = 0; k < dim; ++k)
					koho[i * dim + k] /= weights[i];
			else
				std::copy_n(oldkoho.data() + i * dim,
				            dim,
				            koho + i * dim);
	}
}

struct kohoid_t
{
	unsigned level;
	int x;
	int y;

	using comp = std::complex<float>;

	inline comp grid_pos() const
	{
		return comp(x, y) / float(1 << level) - comp(1, 1) +
		       comp(.5, .5) / float(1 << level);
	}
	inline void sink()
	{
		level++;
		x *= 2;
		y *= 2;
	}
	inline void offset(unsigned xo, unsigned yo)
	{
		x += xo;
		y += yo;
	}
};

inline static float
inv_level_radius(float dim, float level)
{
	return powf(4, level / dim);
}

template<class distf, class nhbr_distf, bool threaded>
static void
gqtsom(size_t threads,
       size_t n,
       size_t in_kohos,
       size_t dim,
       size_t epochs,
       const float *points,
       const int *coords,
       const float *codes,
       const float *radii,
       int *out_kohos,
       float *out_koho,
       int *out_coords,
       float *out_emcoords)
{
	size_t target_kohos = *out_kohos;

	std::vector<float> koho(in_kohos * dim, 0);
	std::copy(codes, codes + in_kohos * dim, koho.begin());
	std::vector<kohoid_t> kohoid(in_kohos);
	for (size_t i = 0; i < in_kohos; ++i) {
		kohoid[i].level = coords[i * 3 + 0];
		kohoid[i].x = coords[i * 3 + 1];
		kohoid[i].y = coords[i * 3 + 2];
	}

	std::vector<std::thread> ts;
	if (threaded)
		ts.resize(threads);

	std::vector<std::vector<float>> taccs, tws;
	taccs.resize(threads);
	tws.resize(threads);

	for (size_t epoch = 0; epoch < epochs; ++epoch) {

		float epochRadius = std::max(min_radius, radii[epoch]);

		auto get_gauss_factor_dstlevel =
		  [&](kohoid_t s, kohoid_t d, float radius) -> float {
			auto spos = s.grid_pos();
			auto dpos = d.grid_pos();
			/* the type conversion here seems a bit harsh
			 * but the standard says that it should work */
			return exp(-sqrf(nhbr_distf::back(nhbr_distf::comp(
			                   (float *)&spos, (float *)&dpos, 2)) *
			                 (1 << d.level) / radius));
		};

		const size_t kohos = kohoid.size();

		// find nearest neighbors for SOM update
		auto collect_kohos = [&](size_t thread_id) {
			std::vector<float> &ws = tws[thread_id];
			std::vector<float> &acc = taccs[thread_id];
			size_t dbegin = thread_id * n / threads,
			       dend = (thread_id + 1) * n / threads;
			const float *d = points + dbegin * dim;
			size_t nd = dend - dbegin;
			acc.resize(kohos * dim);
			ws.resize(kohos);

			for (auto &x : acc)
				x = 0;
			for (auto &x : ws)
				x = 0;

			for (size_t i = 0; i < nd; ++i) {
				size_t closest = 0;
				float closestd =
				  distf::comp(d + i * dim, koho.data(), dim) *
				  inv_level_radius(dim, kohoid[0].level);
				for (size_t j = 1; j < kohos; ++j) {
					float tmp =
					  distf::comp(d + i * dim,
					              koho.data() + j * dim,
					              dim) *
					  inv_level_radius(dim,
					                   kohoid[j].level);
					if (tmp < closestd) {
						closest = j;
						closestd = tmp;
					}
				}

				ws[closest] += 1;
				for (size_t k = 0; k < dim; ++k)
					acc[closest * dim + k] +=
					  d[i * dim + k];
			}
		};

		if (threaded) {
			for (size_t i = 0; i < threads; ++i)
				ts[i] = std::thread(collect_kohos, i);
			for (size_t i = 0; i < threads; ++i)
				ts[i].join();

			// collect
			for (size_t i = 1; i < threads; ++i)
				for (size_t j = 0; j < kohos * dim; ++j)
					taccs[0][j] += taccs[i][j];
			for (size_t i = 1; i < threads; ++i)
				for (size_t j = 0; j < kohos; ++j)
					tws[0][j] += tws[i][j];
		} else
			collect_kohos(0);

		std::vector<float> diffs(kohos, 0), weights(kohos, 0), oldkoho;
		oldkoho.swap(koho);
		koho.resize(kohos * dim, 0);

		// gaussify
		for (size_t si = 0; si < kohos; ++si)
			for (size_t di = 0; di < kohos; ++di) {
				auto factor = get_gauss_factor_dstlevel(
				  kohoid[si], kohoid[di], epochRadius);
				for (size_t k = 0; k < dim; ++k)
					koho[di * dim + k] +=
					  taccs[0][si * dim + k] * factor;
				weights[di] += tws[0][si] * factor;
			}

		// normalize
		for (size_t i = 0; i < kohos; ++i)
			if (weights[i] > 0) {
				for (size_t k = 0; k < dim; ++k)
					koho[i * dim + k] /= weights[i];
				diffs[i] = distf::comp(koho.data() + i * dim,
				                       oldkoho.data() + i * dim,
				                       dim) *
				           weights[i];
			} else
				std::copy_n(oldkoho.data() + i * dim,
				            dim,
				            koho.data() + i * dim);

		// do not spawn new kohos in the last epoch
		if (epoch + 1 >= epochs)
			break;

		// find kohos with largest movement
		std::vector<std::pair<float, size_t>> sqes(kohoid.size());
		for (size_t i = 0; i < kohos; ++i)
			sqes[i] =
			  std::make_pair(diffs[i] / (1 + kohoid[i].level), i);

		size_t target_nodes =
		  ((epoch)*target_kohos + (epochs - epoch - 2) * in_kohos) /
		  (epochs - 2);
		if (target_nodes <= kohos)
			continue;
		if (target_nodes > kohos * 4)
			target_nodes = kohos * 4;
		size_t to_expand = (target_nodes - kohos) / 3;
		/* NB: because we only have target_nodes-sized output space,
		 * this needs to hold here:
		 *
		 *   3*to_expand+kohoid.size() <= target_nodes
		 */
		std::partial_sort(sqes.begin(),
		                  sqes.begin() + to_expand,
		                  sqes.end(),
		                  std::greater<std::pair<float, size_t>>());

		// expand
		koho.reserve(kohos + to_expand * 3 * dim);
		kohoid.reserve(kohos + to_expand * 3);
		for (size_t expanding = 0; expanding < to_expand; ++expanding) {
			const size_t i = sqes[expanding].second;

			kohoid_t a[4];
			std::vector<float> tmpkoho(4 * dim, 0);

			a[0] = kohoid[i];
			a[0].sink();
			a[1] = a[2] = a[3] = a[0];
			a[1].offset(1, 0);
			a[2].offset(0, 1);
			a[3].offset(1, 1);

			auto gauss_guess = [&](kohoid_t k, float *out) {
				float w = 0;
				std::fill_n(out, dim, 0);
				for (size_t ki = 0; ki < kohoid.size(); ++ki) {
					auto factor = get_gauss_factor_dstlevel(
					  kohoid[ki], k, 1);
					for (size_t d = 0; d < dim; ++d)
						out[d] +=
						  koho[d + dim * ki] * factor;
					w += factor;
				}
				for (size_t d = 0; d < dim; ++d)
					if (w > 0)
						out[d] /= w;
					else
						out[d] = koho[d + dim * i];
			};

			for (size_t j = 0; j < 4; ++j)
				gauss_guess(a[j], tmpkoho.data() + j * dim);

			kohoid[i] = a[0];
			kohoid.push_back(a[1]);
			kohoid.push_back(a[2]);
			kohoid.push_back(a[3]);

			std::copy(tmpkoho.begin(),
			          tmpkoho.begin() + dim,
			          koho.begin() + i * dim);
			for (size_t j = 1; j < 4; ++j)
				koho.insert(koho.end(),
				            tmpkoho.begin() + j * dim,
				            tmpkoho.begin() + (j + 1) * dim);
		}
	}

	size_t kohos = kohoid.size();
	if (kohos > target_kohos)
		kohos = target_kohos; // this should not happen but let's be
		                      // safe right?
	*out_kohos = kohos;
	for (size_t i = 0; i < kohos; ++i) {
		for (size_t k = 0; k < dim; ++k)
			out_koho[i * dim + k] = koho[i * dim + k];
		out_coords[i * 3 + 0] = kohoid[i].level;
		out_coords[i * 3 + 1] = kohoid[i].x;
		out_coords[i * 3 + 2] = kohoid[i].y;
		auto p = kohoid[i].grid_pos();
		out_emcoords[i * 2 + 0] = p.real();
		out_emcoords[i * 2 + 1] = p.imag();
	}
}

template<class distf, bool threaded>
static void
mapNNs(size_t threads,
       size_t n,
       size_t k,
       size_t dim,
       const float *points,
       const float *koho,
       int *mapping,
       float *dists)
{
	if (threaded) {
		std::vector<std::thread> ts(threads);
		for (size_t i = 0; i < threads; ++i)
			ts[i] = std::thread(
			  [&](size_t thread_id) {
				  size_t dbegin = thread_id * n / threads,
				         dend = (thread_id + 1) * n / threads;
				  const float *d = points + dbegin * dim;
				  int *m = mapping + dbegin;
				  float *dd = dists + dbegin;
				  size_t nd = dend - dbegin;

				  mapNNs<distf, false>(
				    1, nd, k, dim, d, koho, m, dd);
			  },
			  i);
		for (auto &t : ts)
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

/*
 * R interfaces
 */

extern "C" void
es_C_SOM(float *points,
         float *koho,
         float *nhbrdist,
         float *alphasA,
         float *radiiA,
         float *alphasB,
         float *radiiB,
         Sint *pn,
         Sint *pdim,
         Sint *pkohos,
         Sint *prlen,
         Sint *dist)
{
	size_t n = *pn, dim = *pdim, kohos = *pkohos, rlen = *prlen;

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
es_C_BatchSOM(int *pnthreads,
              float *points,
              float *koho,
              float *nhbrdist,
              float *radii,
              Sint *pn,
              Sint *pdim,
              Sint *pkohos,
              Sint *prlen,
              Sint *dist)
{

	int n = *pn, dim = *pdim, kohos = *pkohos, rlen = *prlen;

	int threads = *pnthreads;

	if (threads < 0)
		threads = 1;
	if (threads == 0)
		threads = std::thread::hardware_concurrency();

	auto somf = threads == 1 ? bsom<distfs::sqeucl, false>
	                         : bsom<distfs::sqeucl, true>;
	if (*dist == 1)
		somf = threads == 1 ? bsom<distfs::manh, false>
		                    : bsom<distfs::manh, true>;
	else if (*dist == 3)
		somf = threads == 1 ? bsom<distfs::chebyshev, false>
		                    : bsom<distfs::chebyshev, true>;

	somf(threads, n, kohos, dim, rlen, points, koho, nhbrdist, radii);
}

extern "C" void
es_C_GQTSOM(int *pnthreads,
            float *points,
            int *coords,
            float *codes,
            float *radii,
            int *out_ncodes,
            float *out_codes,
            int *out_coords,
            float *out_emcoords,
            int *pn,
            int *pdim,
            int *pkohos,
            int *prlen,
            int *distf,
            int *nhbr_distf)
{
	size_t n = *pn, dim = *pdim, kohos = *pkohos, rlen = *prlen,
	       threads = *pnthreads;

	if (threads < 0)
		threads = 1;
	if (threads == 0)
		threads = std::thread::hardware_concurrency();

	if (kohos < 2)
		return; // needed for computability of QE

	auto somf = threads == 1 ? gqtsom<distfs::sqeucl, distfs::sqeucl, false>
	                         : gqtsom<distfs::sqeucl, distfs::sqeucl, true>;
/* it would be great if C++ supported some nice way of data promotion, like with
 * autogenerated switch or so. */
#define choose_somf(DI, DF, NDI, NDF)                                          \
	if (*distf == DI && *nhbr_distf == NDI)                                \
	somf = threads == 1 ? gqtsom<DF, NDF, false> : gqtsom<DF, NDF, true>

	choose_somf(1, distfs::manh, 1, distfs::manh);
	choose_somf(1, distfs::manh, 2, distfs::sqeucl);
	choose_somf(1, distfs::manh, 3, distfs::chebyshev);
	choose_somf(2, distfs::sqeucl, 1, distfs::manh);
	choose_somf(2, distfs::sqeucl, 2, distfs::sqeucl);
	choose_somf(2, distfs::sqeucl, 3, distfs::chebyshev);
	choose_somf(3, distfs::chebyshev, 1, distfs::manh);
	choose_somf(3, distfs::chebyshev, 2, distfs::sqeucl);
	choose_somf(3, distfs::chebyshev, 3, distfs::chebyshev);

#undef choose_somf

	somf(threads,
	     n,
	     kohos,
	     dim,
	     rlen,
	     points,
	     coords,
	     codes,
	     radii,
	     out_ncodes,
	     out_codes,
	     out_coords,
	     out_emcoords);
}

extern "C" void
es_C_mapDataToCodes(int *pnthreads,
                    float *points,
                    float *koho,
                    int *pn,
                    int *pdim,
                    int *pkohos,
                    int *mapping,
                    float *dists,
                    int *dist)
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
