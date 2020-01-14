
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

#include "use_intrins.h"
#include <algorithm>
#include <cmath>

static inline float
sqrf(float n)
{
	return n * n;
}

namespace distfs {
using std::abs;
using std::max;

struct sqeucl
{
	inline static float back(float x) { return sqrt(x); }
	inline static float comp(const float *p1,
	                         const float *p2,
	                         const size_t dim)
	{
#ifndef USE_INTRINS
		float sqdist = 0;
		for (size_t i = 0; i < dim; ++i) {
			float tmp = p1[i] - p2[i];
			sqdist += tmp * tmp;
		}
		return sqdist;
#else
		const float *p1e = p1 + dim, *p1ie = p1e - 3;

		__m128 s = _mm_setzero_ps();
		for (; p1 < p1ie; p1 += 4, p2 += 4) {
			__m128 tmp =
			  _mm_sub_ps(_mm_loadu_ps(p1), _mm_loadu_ps(p2));
			s = _mm_add_ps(_mm_mul_ps(tmp, tmp), s);
		}
		float sqdist = s[0] + s[1] + s[2] + s[3];
		for (; p1 < p1e; ++p1, ++p2) {
			float tmp = *p1 - *p2;
			sqdist += tmp * tmp;
		}
		return sqdist;
#endif
	}

	inline static void proj(const float *la,
	                        const float *lb,
	                        const float *p,
	                        const size_t dim,
	                        float &o_scalar,
	                        float &o_sqdist)
	{

#ifndef USE_INTRINS
		float scalar = 0, sqdist = 0;
		for (size_t k = 0; k < dim; ++k) {
			float tmp = lb[k] - la[k];
			sqdist += tmp * tmp;
			scalar += tmp * (p[k] - la[k]);
		}
#else
		float scalar, sqdist;
		{
			const float *ke = la + dim, *kie = ke - 3;

			__m128 sca = _mm_setzero_ps(), sqd = _mm_setzero_ps();
			for (; la < kie; la += 4, lb += 4, p += 4) {
				__m128 ti = _mm_loadu_ps(la);
				__m128 tmp = _mm_sub_ps(_mm_loadu_ps(lb), ti);
				sqd = _mm_add_ps(sqd, _mm_mul_ps(tmp, tmp));
				sca = _mm_add_ps(
				  sca,
				  _mm_mul_ps(tmp,
				             _mm_sub_ps(_mm_loadu_ps(p), ti)));
			}
			scalar = sca[0] + sca[1] + sca[2] + sca[3];
			sqdist = sqd[0] + sqd[1] + sqd[2] + sqd[3];
			for (; la < ke; ++la, ++lb, ++p) {
				float tmp = *lb - *la;
				sqdist += tmp * tmp;
				scalar += tmp * (*p - *la);
			}
		}
#endif
		o_scalar = scalar;
		o_sqdist = sqdist;
	}
};

#ifdef USE_INTRINS
inline static __m128
abs_mask(void)
{
	__m128i minus1 = _mm_set1_epi32(-1);
	return _mm_castsi128_ps(_mm_srli_epi32(minus1, 1));
}
inline static __m128
vec_abs(__m128 v)
{
	return _mm_and_ps(abs_mask(), v);
}
#endif

struct manh
{
	inline static float back(float x) { return x; }
	inline static float comp(const float *p1,
	                         const float *p2,
	                         const size_t dim)
	{
#ifndef USE_INTRINS
		float mdist = 0;
		for (size_t i = 0; i < dim; ++i) {
			mdist += abs(p1[i] - p2[i]);
		}
		return mdist;
#else
		const float *p1e = p1 + dim, *p1ie = p1e - 3;

		__m128 s = _mm_setzero_ps();
		for (; p1 < p1ie; p1 += 4, p2 += 4) {
			s = _mm_add_ps(s,
			               vec_abs(_mm_sub_ps(_mm_loadu_ps(p1),
			                                  _mm_loadu_ps(p2))));
		}
		float mdist = s[0] + s[1] + s[2] + s[3];
		for (; p1 < p1e; ++p1, ++p2) {
			mdist += abs(*p1 - *p2);
		}
		return mdist;
#endif
	}

	inline static void proj(const float *la,
	                        const float *lb,
	                        const float *p,
	                        const size_t dim,
	                        float &o_scalar,
	                        float &o_sqdist)
	{
		// TODO verify if this is right, TODO SIMD
		float max = 0;
		size_t maxdim = 0;
		for (size_t i = 0; i < dim; ++i) {
			float tmp = abs(lb[i] - la[i]);
			if (tmp > max) {
				max = tmp;
				maxdim = i;
			}
		}
		o_scalar = p[maxdim] - la[maxdim];
		o_sqdist = lb[maxdim] - la[maxdim];
	}
};

struct chebyshev
{
	inline static float back(float x) { return x; }
	inline static float comp(const float *p1,
	                         const float *p2,
	                         const size_t dim)
	{
#ifndef USE_INTRINS
		float cdist = 0;
		for (size_t i = 0; i < dim; ++i) {
			cdist = max(cdist, abs(p1[i] - p2[i]));
		}
		return cdist;
#else
		const float *p1e = p1 + dim, *p1ie = p1e - 3;

		__m128 s = _mm_setzero_ps();
		for (; p1 < p1ie; p1 += 4, p2 += 4) {
			s = _mm_max_ps(s,
			               vec_abs(_mm_sub_ps(_mm_loadu_ps(p1),
			                                  _mm_loadu_ps(p2))));
		}
		float cdist = max(s[0], max(s[1], max(s[2], s[3])));
		for (; p1 < p1e; ++p1, ++p2) {
			cdist = max(cdist, abs(*p1 - *p2));
		}
		return cdist;
#endif
	}

	inline static void proj(const float *la,
	                        const float *lb,
	                        const float *p,
	                        const size_t dim,
	                        float &o_scalar,
	                        float &o_sqdist)
	{
		// TODO FIXME
		o_scalar = 0;
		o_sqdist = 0;
	}
};
};
