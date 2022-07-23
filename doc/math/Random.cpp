/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Random.cpp
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

#include "Random.h"

/*
 */
Random::Random() {
	setSeed((unsigned int)time(NULL));
}

Random::Random(unsigned int seed) {
	setSeed(seed);
}

Random::~Random() {

}

/*
 */
void Random::setSeed(unsigned int s) const {
	seed = s;
}

unsigned int Random::getSeed() const {
	return seed;
}

/*
 */
int Random::getInt(int from,int to) const {
	
	int range = to - from;
	if(range <= 1) return from;
	
	unsigned int rand = get();
	
	if(range < 0xffff) {
		rand >>= 16;
		unsigned int rand_max = ((MAX_RANDOM >> 16) / range) * range;
		while(rand > rand_max) rand = get() >> 16;
	} else {
		unsigned int rand_max = (MAX_RANDOM / range) * range;
		while(rand > rand_max) rand = get();
	}
	
	return from + rand % range;
}

/*
 */
float Random::getFloat(float from,float to) const {
	
	union { unsigned int i; float f; } rand = { (0x3f800000 | (get() & 0x007fffff)) };
	
	return from + (rand.f - 1.0f) * (to - from);
}

/*
 */
double Random::getDouble(double from,double to) const {
	
	return from + (double)get() / MAX_RANDOM * (to - from);
}
