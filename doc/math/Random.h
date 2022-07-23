/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Random.h
 * Desc:    Random number generator
 * Version: 1.05
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

#ifndef __RANDOM_H__
#define __RANDOM_H__

#include "Base.h"

/*
 */
class Random {
		
	public:
		
		Random();
		Random(unsigned int seed);
		~Random();
		
		void setSeed(unsigned int seed) const;
		unsigned int getSeed() const;
		
		enum {
			A = 1664525,
			C = 1013904223,
			MAX_RANDOM = 0x7fffffff,
		};
		
		INLINE unsigned int get() const {
			seed = (unsigned int)((unsigned long long)seed * A + C) & MAX_RANDOM;
			return seed;
		}
		
		int getInt(int from,int to) const;
		float getFloat(float from,float to) const;
		double getDouble(double from,double to) const;
		
	private:
		
		mutable unsigned int seed;
};

#endif /* __RANDOM_H__ */
