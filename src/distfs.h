
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
 */

static float sqeucl (float* p1, float* p2, int px, int n, int ncodes)
{
	int j;
	float tmp;

	float xdist = 0.0;
	for (j = 0; j < px; j++) {
		tmp = p1[j * n] - p2[j * ncodes];
		xdist += tmp * tmp;
	}
	return xdist;
}

static float manh (float* p1, float* p2, int px, int n, int ncodes)
{
	int j;
	float xdist = 0.0, tmp;
	for (j = 0; j < px; j++) {
		tmp = p1[j * n] - p2[j * ncodes];
		xdist += abs (tmp);
	}
	return xdist;
}

static float chebyshev (float* p1, float* p2, int px, int n, int ncodes)
{
	int j;
	float xdist = 0.0, tmp;
	for (j = 0; j < px; j++) {
		tmp = abs(p1[j * n] - p2[j * ncodes]);
		if (tmp > xdist) xdist = tmp;
	}
	return xdist;
}
