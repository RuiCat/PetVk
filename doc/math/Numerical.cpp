/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Numerical.h
 * Desc:    Numerical library
 * Version: 1.11
 * Author:  Alexander Zapryagaev <frustum@unigine.com>
 *
 * This file is part of the Unigine engine (http://unigine.com/).
 *
 * Your use and or redistribution of this software in source and / or
 * binary form, with or without modification, is subject to: (i) your
 * ongoing acceptance of and compliance with the terms and conditions of
 * the Unigine License Agreement; and (ii) your inclusion of this notice
 * in any version of this software that you use or redistribute.
 * A copy of the Unigine License Agreement is available by contacting
 * Unigine Corp. at http://unigine.com/
 */

#include "Numerical.h"

/******************************************************************************\
*
* Bisection solver
*
\******************************************************************************/

/*
 */
static float bi_polynom(const float *c,int size,float x) {
	float ret = 0.0f;
	if(size) {
		ret = c[0];
		float value = x;
		for(int i = 1; i < size; i++) {
			ret += c[i] * value;
			value *= x;
		}
	}
	return ret;
}

/*
 */
static int bi_bisection(float &x,const float *c,int size,float x0,float x1,float iepsilon) {
	float y0 = bi_polynom(c,size,x0);
	if(Math::abs(y0) < EPSILON) {
		x = x0;
		return 1;
	}
	float y1 = bi_polynom(c,size,x1);
	if(Math::abs(y1) < EPSILON) {
		x = x1;
		return 1;
	}
	if(y0 * y1 > 0.0f) {
		return 0;
	}
	int num_iterations = Math::ftoi(Math::log((x1 - x0) * iepsilon) / LOG2 + 0.5f);
	for(int i = 0; i < num_iterations; i++) {
		x = (x0 + x1) * 0.5f;
		float y = bi_polynom(c,size,x);
		float v = y0 * y;
		if(v < 0.0f) {
			x1 = x;
			y1 = y;
		} else if(v > 0.0f) {
			x0 = x;
			y0 = y;
		} else {
			break;
		}
	}
	return 1;
}

/*
 */
static int bi_solve(float *ret,const float *c,int size,float x0,float x1,float iepsilon) {
	float *r = ret;
	float roots[32];
	float derivative[32];
	if(size == 1) {
		if(bi_bisection(*r,c,size,x0,x1,iepsilon)) r++;
	} else {
		for(int i = 1; i < size; i++) {
			derivative[i - 1] = c[i] * Math::itof(i);
		}
		int num = bi_solve(roots,derivative,size - 1,x0,x1,iepsilon);
		if(num > 0) {
			if(bi_bisection(*r,c,size,x0,roots[0],iepsilon)) r++;
			for(int i = 0; i < num - 1; i++) {
				if(bi_bisection(*r,c,size,roots[i],roots[i + 1],iepsilon)) r++;
			}
			if(bi_bisection(*r,c,size,roots[num - 1],x1,iepsilon)) r++;
		} else {
			if(bi_bisection(*r,c,size,x0,x1,iepsilon)) r++;
		}
	}
	return (int)(r - ret);
}

/*
 */
int biSolve(float *ret,const float *c,int size,float iepsilon) {
	float range = 0.0f;
	int degree = size - 1;
	for(int i = degree; i >= 0; i--) {
		if(Math::abs(c[i]) > EPSILON) break;
	}
	float idegree = 1.0f / c[degree];
	for(int i = 0; i < degree; i++) {
		float r = Math::abs(c[i]) * idegree;
		if(range < r) range = r;
	}
	range += 1.0f;
	return bi_solve(ret,c,size,-range,range,iepsilon);
}

/******************************************************************************\
*
* LU solver
*
\******************************************************************************/

/*
 */
int luDecompose(float *m,int size,int *index) {
	int ret = 1;
	for(int i = 0; i < size; i++) {
		index[i] = i;
	}
	float *diagonal = m;
	for(int i = 1; i < size; i++) {
		float *d0 = diagonal + 1;
		float *d1 = diagonal + size;
		for(int j = i; j < size; j++) {
			float f = *d0;
			*d0++ = *d1;
			*d1 = f;
			d1 += size;
		}
		diagonal += size + 1;
	}
	diagonal = m;
	for(int i = 0; i < size; i++) {
		int num = i;
		float max = Math::abs(*diagonal);
		const float *s = diagonal + size;
		for(int j = i + 1; j < size; j++) {
			float f = Math::abs(*s);
			s += size;
			if(max < f) {
				max = f;
				num = j;
			}
		}
		if(Math::abs(max) < EPSILON) {
			return 0;
		}
		if(i != num) {
			ret = -ret;
			int j = index[i];
			index[i] = index[num];
			index[num] = j;
			float *d0 = m + size * i;
			float *d1 = m + size * num;
			for(int j = 0; j < size; j++) {
				float d = *d0;
				*d0 = *d1;
				*d1 = d;
				d0++;
				d1++;
			}
		}
		float f = 1.0f / *diagonal;
		float *d = diagonal + size;
		for(int j = i + 1; j < size; j++) {
			*d *= f;
			d += size;
		}
		d = diagonal + size;
		int to = size - ((size - i - 1) & 3);
		for(int j = i + 1; j < to; j += 4) {
			Simd::eliminate(d + 1,diagonal + 1,d,size,size - i - 1);
			d += size << 2;
		}
		for(int j = to; j < size; j++) {
			float f = *d++;
			Simd::mad(d,diagonal + 1,-f,d,size - i - 1);
			d += size - 1;
		}
		diagonal += size + 1;
	}
	return ret;
}

/*
 */
float luDeterminant(const float *m,int size,int ret) {
	const float *s = m;
	float det = Math::itof(ret);
	for(int i = 0; i < size; i++) {
		det *= s[size * i + i];
	}
	return det;
}

/*
 */
void luSolve3(float *ret,const float *m,const float *b,const int *index) {
	float *d = ret;
	const float *s = m;
	d[0] = b[index[0]];
	d[1] = b[index[1]] - s[3] * d[0];
	d[2] = b[index[2]] - s[6] * d[0] - s[7] * d[1];
	d[2] = (d[2]) / s[8];
	d[1] = (d[1] - s[5] * d[2]) / s[4];
	d[0] = (d[0] - s[1] * d[1] - s[2] * d[2]) / s[0];
}

void luSolve4(float *ret,const float *m,const float *b,const int *index) {
	float *d = ret;
	const float *s = m;
	d[0] = b[index[0]];
	d[1] = b[index[1]] - s[4] * d[0];
	d[2] = b[index[2]] - s[8] * d[0] - s[9] * d[1];
	d[3] = b[index[3]] - s[12] * d[0] - s[13] * d[1] - s[14] * d[2];
	d[3] = (d[3]) / s[15];
	d[2] = (d[2] - s[11] * d[3]) / s[10];
	d[1] = (d[1] - s[16] * d[2] - s[17] * d[3]) / s[5];
	d[0] = (d[0] - s[1]  * d[1] - s[2]  * d[2] - s[3] * d[3]) / s[0];
}

void luSolve5(float *ret,const float *m,const float *b,const int *index) {
	float *d = ret;
	const float *s = m;
	d[0] = b[index[0]];
	d[1] = b[index[1]] - s[5] * d[0];
	d[2] = b[index[2]] - s[10] * d[0] - s[11] * d[1];
	d[3] = b[index[3]] - s[15] * d[0] - s[16] * d[1] - s[17] * d[2];
	d[4] = b[index[4]] - s[20] * d[0] - s[21] * d[1] - s[22] * d[2] - s[23] * d[3];
	d[4] = (d[4]) / s[24];
	d[3] = (d[3] - s[19] * d[4]) / s[18];
	d[2] = (d[2] - s[13] * d[3] - s[14] * d[4]) / s[12];
	d[1] = (d[1] - s[7]  * d[2] - s[8]  * d[3] - s[9] * d[4]) / s[6];
	d[0] = (d[0] - s[1]  * d[1] - s[2]  * d[2] - s[3] * d[3] - s[4] * d[4]) / s[0];
}

void luSolve6(float *ret,const float *m,const float *b,const int *index) {
	float *d = ret;
	const float *s = m;
	d[0] = b[index[0]];
	d[1] = b[index[1]] - s[6] * d[0];
	d[2] = b[index[2]] - s[12] * d[0] - s[13] * d[1];
	d[3] = b[index[3]] - s[18] * d[0] - s[19] * d[1] - s[20] * d[2];
	d[4] = b[index[4]] - s[24] * d[0] - s[25] * d[1] - s[26] * d[2] - s[27] * d[3];
	d[5] = b[index[5]] - s[30] * d[0] - s[31] * d[1] - s[32] * d[2] - s[33] * d[3] - s[34] * d[4];
	d[5] = (d[5]) / s[35];
	d[4] = (d[4] - s[29] * d[5]) / s[28];
	d[3] = (d[3] - s[22] * d[4] - s[23] * d[5]) / s[21];
	d[2] = (d[2] - s[15] * d[3] - s[16] * d[4] - s[17] * d[5]) / s[14];
	d[1] = (d[1] - s[8]  * d[2] - s[9]  * d[3] - s[10] * d[4] - s[11] * d[5]) / s[7];
	d[0] = (d[0] - s[1]  * d[1] - s[2]  * d[2] - s[3]  * d[3] - s[4]  * d[4] - s[5] * d[5]) / s[0];
}

void luSolve(float *ret,const float *m,const float *b,int size,const int *index) {
	float sum;
	float *d = ret;
	const float *s = m;
	for(int i = 0; i < size; i++) {
		Simd::dot(sum,s,d,i);
		d[i] = b[index[i]] - sum;
		s += size;
	}
	s = m + size * size - 1;
	for(int i = size - 1; i >= 0; i--) {
		Simd::dot(sum,s + 1,d + i + 1,size - i - 1);
		d[i] = (d[i] - sum) / *s;
		s -= size + 1;
	}
}
