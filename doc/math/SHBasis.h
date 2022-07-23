/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    SHBasis.h
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

#ifndef __SH_BASIS_H__
#define __SH_BASIS_H__

#include "Base.h"

/*
 */
class SHBasis {
		
	public:
		
		SHBasis();
		~SHBasis();
		
		double factorial(int x) const;
		
		double get(int l,int m,const float *dir) const;
		double get(int l,int m,double phi,double theta) const;
		
	private:
		
		double K(int l,int m) const;
		double P(int l,int m,double x) const;
};

#endif /* __SH_BASIS_H__ */
