/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    SHBasis.cpp
 * Desc:    Spherical harmonics orthogonal basis
 * Version: 1.08
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

#include "MathLib.h"
#include "SHBasis.h"

/*
 */
SHBasis::SHBasis() {
	
}

SHBasis::~SHBasis() {
	
}

/*
 */
double SHBasis::factorial(int x) const {
	static const double table[16] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0,
		40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0 };
	if(x < 16) return table[x];
	double ret = table[15];
	for(int i = 16; i <= x; i++) ret *= (double)i;
	return ret;
}

/*
 */
double SHBasis::K(int l,int m) const {
	return Math::sqrt((2.0 * l + 1.0) * factorial(l - m) / (4.0 * 3.14159265358979323846 * factorial(l + m)));
}

double SHBasis::P(int l,int m,double x) const {
	double pmm = 1.0;
	if(m > 0) {
		double fact = 1.0;
		double somx2 = Math::sqrt((1.0 - x) * (1.0 + x));
		for(int i = 0; i < m; i++) {
			pmm *= -fact * somx2;
			fact += 2.0;
		}
	}
	if(l == m) return pmm;
	double pmmp1 = x * (2.0 * m + 1.0) * pmm;
	if(l == m + 1) return pmmp1;
	double pll = 0.0;
	for(int i = m + 2; i <= l; i++) {
		pll = ((2.0 * i - 1.0) * x * pmmp1 - (i + m - 1.0) * pmm) / (i - m);
		pmm = pmmp1;
		pmmp1 = pll;
	}
	return pll;
}

/*
 */
double SHBasis::get(int l,int m,const float *dir) const {
	double length = Math::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
	return get(l,m,Math::atan2(-dir[1],-dir[0]),Math::acos(dir[2] / length));
}

double SHBasis::get(int l,int m,double phi,double theta) const {
	const double sqrt2 = Math::sqrt(2.0);
	if(m > 0) return sqrt2 * K(l,m) * Math::cos(phi * m) * P(l,m,Math::cos(theta));
	if(m < 0) return sqrt2 * K(l,-m) * Math::sin(-phi * m) * P(l,-m,Math::cos(theta));
	return K(l,0) * P(l,m,Math::cos(theta));
}
