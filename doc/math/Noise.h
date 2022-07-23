/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Noise.h
 * Desc:    Perlin noise
 * Version: 1.02
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

#ifndef __NOISE_H__
#define __NOISE_H__

#include "Base.h"

/*
 */
class Noise {
		
	public:
		
		Noise();
		Noise(unsigned int seed);
		~Noise();
		
		void setSeed(unsigned int seed);
		unsigned int getSeed() const;
		
		float get1(float x) const;
		float get2(float x,float y) const;
		float get3(float x,float y,float z) const;
		
		float getTurbulence1(float x,int frequency) const;
		float getTurbulence2(float x,float y,int frequency) const;
		float getTurbulence3(float x,float y,float z,int frequency) const;
		
		float getTileable1(float x,float width) const;
		float getTileable2(float x,float y,float width,float height) const;
		float getTileable3(float x,float y,float z,float width,float height,float depth) const;
		
		float getTileableTurbulence1(float x,float width,int frequency) const;
		float getTileableTurbulence2(float x,float y,float width,float height,int frequency) const;
		float getTileableTurbulence3(float x,float y,float z,float width,float height,float depth,int frequency) const;
		
	private:
		
		enum {
			A = 1664525,
			C = 1013904223,
			MAX_RANDOM = 0x7fffffff,
			SAMPLES = 256,
		};
		
		unsigned int get_random_int();
		float get_random_float();
		
		unsigned int seed;
		
		int permutation[SAMPLES * 2 + 2];
		float gradient1[SAMPLES * 2 + 2][1];
		float gradient2[SAMPLES * 2 + 2][2];
		float gradient3[SAMPLES * 2 + 2][3];
};

#endif /* __NOISE_H__ */
