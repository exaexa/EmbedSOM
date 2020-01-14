
/* This file is part of EmbedSOM.
 *
 * Copyright (C) 2018-2020 Mirek Kratochvil <exa.exa@gmail.com>
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

#include "distfs.h"

using namespace std;

// some small numbers first!
static const float min_boost = 1e-5; // lower limit for the parameter

// this affects how steeply the score decreases for near-farthest codes
static const float max_avoidance = 10;

// this is added before normalizing the distances
static const float zero_avoidance = 1e-10;

// a tiny epsilon for preventing singularities
static const float koho_gravity = 1e-5;

/*
 * KNN computation (mainly heap)
 */
struct dist_id
{
	float dist;
	size_t id;
};

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
				swap(heap[L], heap[start]);
				start = L;
			} else {
				if (heap[start].dist >= dr)
					break;
				swap(heap[R], heap[start]);
				start = R;
			}
		} else if (L < lim) {
			if (heap[start].dist < heap[L].dist)
				swap(heap[L], heap[start]);
			break; // exit safely!
		} else
			break;
	}
}

template<class distf>
void
knn(const float* point,
    const float* koho,
    size_t kohos,
    size_t dim,
    size_t topnn,
    vector<dist_id>& dists)
{
	size_t i;

	// push first topnn kohos
	for (i = 0; i < topnn; ++i) {
		dists[i].dist = distf::comp(point, koho + i * dim, dim);
		dists[i].id = i;
	}

	// make a heap
	for (i = 0; i < topnn; ++i)
		heap_down(dists.data(), topnn - i - 1, topnn);

	// insert the rest
	for (i = topnn; i < kohos; ++i) {
		float s = distf::comp(point, koho + i * dim, dim);
		if (dists[0].dist < s)
			continue;
		dists[0].dist = s;
		dists[0].id = i;
		heap_down(dists.data(), 0, topnn);
	}

	// heapsort the NNs
	for (i = topnn - 1; i > 0; --i) {
		swap(dists[0], dists[i]);
		heap_down(dists.data(), 0, i);
	}
}

/*
 * Projection- and fitting-related helpers
 */

template<int embed_dim>
void
add_gravity(const float* emcoords, float score, float* mtx)
{
	float gs = score * koho_gravity;
	if (embed_dim == 2) {
		mtx[0] += gs;
		mtx[3] += gs;
		mtx[4] += gs * emcoords[0];
		mtx[5] += gs * emcoords[1];
	}
	if (embed_dim == 3) {
		mtx[0] += gs;
		mtx[4] += gs;
		mtx[8] += gs;
		mtx[9] += gs * emcoords[0];
		mtx[10] += gs * emcoords[1];
		mtx[11] += gs * emcoords[2];
	}
}

template<int embed_dim>
inline static float
dotp_ec(const float* a, const float* b)
{
	float r = 0;
	for (size_t i = 0; i < embed_dim; ++i)
		r += a[i] * b[i];
	return r;
}

template<int embed_dim>
void
add_approximation(float score_i,
                  float score_j,
                  const float* iec,
                  const float* jec,
                  float scalar_proj,
                  float adjust,
                  float* mtx)
{
	float h[embed_dim], hp = 0;
	for (size_t i = 0; i < embed_dim; ++i)
		hp += sqrf(h[i] = jec[i] - iec[i]);
	if (hp < zero_avoidance)
		return;

	const float s = score_i * score_j * powf(1 + hp, -adjust) *
	                expf(-sqrf(scalar_proj - .5));
	const float sihp = s / hp;
	const float rhsc = s * (scalar_proj + dotp_ec<embed_dim>(h, iec) / hp);

	if (embed_dim == 2) {

		mtx[0] += h[0] * h[0] * sihp;
		mtx[1] += h[0] * h[1] * sihp;
		mtx[2] += h[1] * h[0] * sihp;
		mtx[3] += h[1] * h[1] * sihp;
		mtx[4] += h[0] * rhsc;
		mtx[5] += h[1] * rhsc;
	}

	if (embed_dim == 3) {
		mtx[0] += h[0] * h[0] * sihp;
		mtx[1] += h[0] * h[1] * sihp;
		mtx[2] += h[0] * h[2] * sihp;
		mtx[3] += h[1] * h[0] * sihp;
		mtx[4] += h[1] * h[1] * sihp;
		mtx[5] += h[1] * h[2] * sihp;
		mtx[6] += h[2] * h[0] * sihp;
		mtx[7] += h[2] * h[1] * sihp;
		mtx[8] += h[2] * h[2] * sihp;
		mtx[9] += h[0] * rhsc;
		mtx[10] += h[1] * rhsc;
		mtx[11] += h[2] * rhsc;
	}
}

template<int embed_dim>
void
solve_lin_eq(const float* mtx, float* embedding)
{
	if (embed_dim == 2) {
		float det = mtx[0] * mtx[3] - mtx[1] * mtx[2];
		embedding[0] = (mtx[4] * mtx[3] - mtx[5] * mtx[2]) / det;
		embedding[1] = (mtx[0] * mtx[5] - mtx[1] * mtx[4]) / det;
	}
	if (embed_dim == 3) {
		// this looks ugly
		float det =
		  mtx[0] * mtx[4] * mtx[8] + mtx[1] * mtx[5] * mtx[6] +
		  mtx[2] * mtx[3] * mtx[7] - mtx[0] * mtx[5] * mtx[7] -
		  mtx[1] * mtx[3] * mtx[8] - mtx[2] * mtx[4] * mtx[6];
		embedding[0] =
		  (mtx[9] * mtx[4] * mtx[8] + mtx[10] * mtx[5] * mtx[6] +
		   mtx[11] * mtx[3] * mtx[7] - mtx[9] * mtx[5] * mtx[7] -
		   mtx[10] * mtx[3] * mtx[8] - mtx[11] * mtx[4] * mtx[6]) /
		  det;
		embedding[1] =
		  (mtx[0] * mtx[10] * mtx[8] + mtx[1] * mtx[11] * mtx[6] +
		   mtx[2] * mtx[9] * mtx[7] - mtx[0] * mtx[11] * mtx[7] -
		   mtx[1] * mtx[9] * mtx[8] - mtx[2] * mtx[10] * mtx[6]) /
		  det;
		embedding[2] =
		  (mtx[0] * mtx[4] * mtx[11] + mtx[1] * mtx[5] * mtx[9] +
		   mtx[2] * mtx[3] * mtx[10] - mtx[0] * mtx[5] * mtx[10] -
		   mtx[1] * mtx[3] * mtx[11] - mtx[2] * mtx[4] * mtx[9]) /
		  det;
	}
}

/*
 * EmbedSOM function for a single point
 */
template<class distf, int embed_dim>
static void
embedsom_point(const size_t kohos,
               const size_t dim,
               const float boost,
               const size_t topn,
               const float adjust,
               const float* point,
               const float* koho,
               const float* emcoords,
               float* embedding,
               vector<dist_id>& dists)
{
	const size_t topnn = topn < kohos ? topn + 1 : topn;

	knn<distf>(point, koho, kohos, dim, topnn, dists);

	// compute the distance distribution for the scores
	float mean = 0, sd = 0, wsum = 0;
	for (size_t i = 0; i < topnn; ++i) {
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

	// convert the stuff to scores
	for (size_t i = 0; i < topn; ++i)
		if (topn < topnn)
			dists[i].dist =
			  expf((mean - dists[i].dist) * sd) *
			  (1 - expf(dists[i].dist * nmax - max_avoidance));
		else
			dists[i].dist = expf((mean - dists[i].dist) * sd);

	// create the empty equation matrix
	float mtx[embed_dim * (1 + embed_dim)];
	fill(mtx, mtx + embed_dim * (1 + embed_dim), 0);

	// for all points in the neighborhood
	for (size_t i = 0; i < topn; ++i) {
		size_t idx = dists[i].id;
		float score_i = dists[i].dist; // score of 'i'

		float iec[embed_dim]; // emcoords for point 'i'
		copy_n(emcoords + embed_dim * idx, embed_dim, iec);

		/* this adds a really tiny influence of the point to
		 * prevent singularities */
		add_gravity<embed_dim>(iec, score_i, mtx);

		// for all combinations of point 'i' with points in the
		// neighborhood
		for (size_t j = i + 1; j < topn; ++j) {

			size_t jdx = dists[j].id;
			float score_j = dists[j].dist; // score of 'j'

			float jec[embed_dim]; // emcoords for point 'j'
			copy_n(emcoords + embed_dim * jdx, embed_dim, jec);

			float scalar, sqdist;
			distf::proj(koho + dim * idx,
			            koho + dim * jdx,
			            point,
			            dim,
			            scalar,
			            sqdist);

			if (sqdist == 0)
				continue;
			else
				scalar /= sqdist;

			add_approximation<embed_dim>(
			  score_i, score_j, iec, jec, scalar, adjust, mtx);
		}
	}

	solve_lin_eq<embed_dim>(mtx, embedding);
}

/*
 * EmbedSOM batch for multiple points
 */
template<class distf, int embed_dim>
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
	// spawn more batches in threads if parallelization is required
	if (threads > 1) {
		vector<thread> ts(threads);
		for (size_t i = 0; i < threads; ++i)
			ts[i] = thread(
			  [&](size_t thread_id) {
				  size_t dbegin = thread_id * n / threads,
				         dend = (thread_id + 1) * n / threads;
				  const float* d = points + dbegin * dim;
				  float* e = embedding + dbegin * embed_dim;
				  size_t nd = dend - dbegin;

				  embedsom<distf, embed_dim>(1,
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

	// single-thread version
	const size_t topnn = topn < kohos ? topn + 1 : topn;

	vector<dist_id> dists;
	dists.resize(topnn);

	for (size_t i = 0; i < n; ++i)
		embedsom_point<distf, embed_dim>(kohos,
		                                 dim,
		                                 boost,
		                                 topn,
		                                 adjust,
		                                 points + dim * i,
		                                 koho,
		                                 emcoords,
		                                 embedding + embed_dim * i,
		                                 dists);
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

	auto emf = embedsom<distfs::sqeucl, 2>;
	if (embeddim == 2) {
		if (*pdist == 1)
			emf = embedsom<distfs::manh, 2>;
		if (*pdist == 3)
			emf = embedsom<distfs::chebyshev, 2>;
	} else if (embeddim == 3) {
		emf = embedsom<distfs::sqeucl, 3>;
		if (*pdist == 1)
			emf = embedsom<distfs::manh, 3>;
		if (*pdist == 3)
			emf = embedsom<distfs::chebyshev, 3>;
	} else
		return; // waat.

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
	{ "es_C_BatchSOM", (DL_FUNC)&es_C_BatchSOM, 10 },
	{ "es_C_mapDataToCodes", (DL_FUNC)&es_C_mapDataToCodes, 9 },
	{ NULL, NULL, 0 }
};

void
R_init_EmbedSOM(DllInfo* info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}
