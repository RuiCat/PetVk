/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Noise.cpp
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

#include "Noise.h"

/*
 */
Noise::Noise() {
	setSeed((unsigned int)time(NULL));
}

Noise::Noise(unsigned int seed) {
	setSeed(seed);
}

Noise::~Noise() {
	
}

/*
 */
unsigned int Noise::get_random_int() {
	seed = (unsigned int)((unsigned long long)seed * A + C) & MAX_RANDOM;
	return seed;
}

float Noise::get_random_float() {
	return (float)get_random_int() / MAX_RANDOM * 2.0f - 1.0f;
}

/*
 */
static void normalize2(float *v) {
	float ilength = 1.0f / sqrtf(v[0] * v[0] + v[1] * v[1]);
	v[0] = v[0] * ilength;
	v[1] = v[1] * ilength;
}

static void normalize3(float *v) {
	float ilength = 1.0f / sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] = v[0] * ilength;
	v[1] = v[1] * ilength;
	v[2] = v[2] * ilength;
}

/*
 */
void Noise::setSeed(unsigned int s) {
	
	seed = s;
	if(seed < 1) seed = 1;
	
	// make random arrays
	for(int i = 0; i < SAMPLES; i++) {
		permutation[i] = i;
		gradient1[i][0] = get_random_float();
		gradient2[i][0] = get_random_float();
		gradient2[i][1] = get_random_float();
		gradient3[i][0] = get_random_float();
		gradient3[i][1] = get_random_float();
		gradient3[i][2] = get_random_float();
		normalize2(gradient2[i]);
		normalize3(gradient3[i]);
	}
	
	// permutation array
	for(int i = 0; i < SAMPLES; i++) {
		int j = (get_random_int() >> 16) % SAMPLES;
		int k = permutation[i];
		permutation[i] = permutation[j];
		permutation[j] = k;
	}
	
	// duplicate randoms
	for(int i = 0; i < SAMPLES + 2; i++) {
		permutation[SAMPLES + i] = permutation[i];
		gradient1[SAMPLES + i][0] = gradient1[i][0];
		gradient2[SAMPLES + i][0] = gradient2[i][0];
		gradient2[SAMPLES + i][1] = gradient2[i][1];
		gradient3[SAMPLES + i][0] = gradient3[i][0];
		gradient3[SAMPLES + i][1] = gradient3[i][1];
		gradient3[SAMPLES + i][2] = gradient3[i][2];
	}
}

unsigned int Noise::getSeed() const {
	return seed;
}

/*
 */
static INLINE void setup(float x,int &b0,int &b1,float &r0,float &r1) {
	b0 = (int)(x + 4096.0f) & 255;
	b1 = (b0 + 1) & 255;
	r0 = (x + 4096.0f) - (int)(x + 4096.0f);
	r1 = r0 - 1.0f;
}

static INLINE float fade(float x) {
	return x * x * x * (x * (x * 6.0f - 15.0f) + 10.0f);
}

static INLINE float dot2(const float *v,float x,float y) {
	return v[0] * x + v[1] * y;
}

static INLINE float dot3(const float *v,float x,float y,float z) {
	return v[0] * x + v[1] * y + v[2] * z;
}

static INLINE float lerp(float v0,float v1,float k) {
	return v0 + (v1 - v0) * k;
}

/*
 */
float Noise::get1(float x) const {
	
	int bx0,bx1;
	float rx0,rx1;
	
	setup(x,bx0,bx1,rx0,rx1);
	
	float sx = fade(rx0);
	
	float a = gradient1[permutation[bx0]][0] * rx0;
	float b = gradient1[permutation[bx1]][0] * rx1;
	
	return lerp(a,b,sx);
}

float Noise::get2(float x,float y) const {
	
	int bx0,bx1,by0,by1;
	float rx0,rx1,ry0,ry1;
	
	setup(x,bx0,bx1,rx0,rx1);
	setup(y,by0,by1,ry0,ry1);
	
	float sx = fade(rx0);
	float sy = fade(ry0);
	
	int b00 = permutation[permutation[bx0] + by0];
	int b10 = permutation[permutation[bx1] + by0];
	int b01 = permutation[permutation[bx0] + by1];
	int b11 = permutation[permutation[bx1] + by1];
	
	float a = lerp(dot2(gradient2[b00],rx0,ry0),dot2(gradient2[b10],rx1,ry0),sx);
	float b = lerp(dot2(gradient2[b01],rx0,ry1),dot2(gradient2[b11],rx1,ry1),sx);
	
	return lerp(a,b,sy);
}

float Noise::get3(float x,float y,float z) const {
	
	int bx0,bx1,by0,by1,bz0,bz1;
	float rx0,rx1,ry0,ry1,rz0,rz1;
	
	setup(x,bx0,bx1,rx0,rx1);
	setup(y,by0,by1,ry0,ry1);
	setup(z,bz0,bz1,rz0,rz1);
	
	float sx = fade(rx0);
	float sy = fade(ry0);
	float sz = fade(rz0);
	
	int b00 = permutation[permutation[bx0] + by0];
	int b10 = permutation[permutation[bx1] + by0];
	int b01 = permutation[permutation[bx0] + by1];
	int b11 = permutation[permutation[bx1] + by1];
	
	float a0 = lerp(dot3(gradient3[b00 + bz0],rx0,ry0,rz0),dot3(gradient3[b10 + bz0],rx1,ry0,rz0),sx);
	float b0 = lerp(dot3(gradient3[b01 + bz0],rx0,ry1,rz0),dot3(gradient3[b11 + bz0],rx1,ry1,rz0),sx);
	
	float a1 = lerp(dot3(gradient3[b00 + bz1],rx0,ry0,rz1),dot3(gradient3[b10 + bz1],rx1,ry0,rz1),sx);
	float b1 = lerp(dot3(gradient3[b01 + bz1],rx0,ry1,rz1),dot3(gradient3[b11 + bz1],rx1,ry1,rz1),sx);
	
	float c0 = lerp(a0,b0,sy);
	float c1 = lerp(a1,b1,sy);
	
	return lerp(c0,c1,sz);
}

/*
 */
float Noise::getTurbulence1(float x,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += get1(x * f) / f;
	}
	return res;
}

float Noise::getTurbulence2(float x,float y,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += get2(x * f,y * f) / f;
	}
	return res;
}

float Noise::getTurbulence3(float x,float y,float z,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += get3(x * f,y * f,z * f) / f;
	}
	return res;
}

/*
 */
float Noise::getTileable1(float x,float width) const {
	x -= (int)(x / width) * width;
	if(x < 0.0f) x += width;
	return (get1(x) * (width - x) + get1(x - width) * x) / width;
}

float Noise::getTileable2(float x,float y,float width,float height) const {
	x -= (int)(x / width) * width;
	y -= (int)(y / height) * height;
	if(x < 0.0f) x += width;
	if(y < 0.0f) y += height;
	return (get2(x,y) * (width - x) * (height - y) +
			get2(x - width,y) * x * (height - y) +
			get2(x,y - height) * (width - x) * y +
			get2(x - width,y - height) * x * y) / (width * height);
}

float Noise::getTileable3(float x,float y,float z,float width,float height,float depth) const {
	x -= (int)(x / width) * width;
	y -= (int)(y / height) * height;
	z -= (int)(z / depth) * depth;
	if(x < 0.0f) x += width;
	if(y < 0.0f) y += height;
	if(z < 0.0f) z += depth;
	return (get3(x,y,z) * (width - x) * (height - y) * (depth - z) +
			get3(x - width,y,z) * x * (height - y) * (depth - z) +
			get3(x,y - height,z) * (width - x) * y * (depth - z) +
			get3(x - width,y - height,z) * x * y * (depth - z) +
			get3(x,y,z - depth) * (width - x) * (height - y) * z +
			get3(x - width,y,z - depth) * x * (height - y) * z +
			get3(x,y - height,z - depth) * (width - x) * y * z +
			get3(x - width,y - height,z - depth) * x * y * z) / (width * height * depth);
}

/*
 */
float Noise::getTileableTurbulence1(float x,float width,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += getTileable1(x * f,width * f) / f;
	}
	return res;
}

float Noise::getTileableTurbulence2(float x,float y,float width,float height,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += getTileable2(x * f,y * f,width * f,height * f) / f;
	}
	return res;
}

float Noise::getTileableTurbulence3(float x,float y,float z,float width,float height,float depth,int frequency) const {
	float res = 0.0f;
	for(int f = frequency; f >= 1; f >>= 1) {
		res += getTileable3(x * f,y * f,z * f,width * f,height * f,depth * f) / f;
	}
	return res;
}
