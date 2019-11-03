
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

namespace distfs {
using std::abs;
using std::max;

struct sqeucl
{
	inline static float back(float x) { return sqrt(x); }
	inline static float comp(const float* p1,
	                         const float* p2,
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
	inline static float comp(const float* p1,
	                         const float* p2,
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
};

struct chebyshev
{
	inline static float back(float x) { return x; }
	inline static float comp(const float* p1,
	                         const float* p2,
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
};
};
