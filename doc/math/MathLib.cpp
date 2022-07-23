/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    MathLib.cpp
 * Desc:    Math library
 * Version: 1.75
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

#include "Log.h"
#include "MathLib.h"

/******************************************************************************\
*
* MathLib
*
\******************************************************************************/

/*
 */
Math::Math() {
	
}

/*
 */
float Math::exp2Fast(float v) {
	int i = ftoi(v - 0.5f);
	v = v - itof(i);
	return IntFloat((i + 127) << 23).f * (((((0.0018775767f * v + 0.0089893397f) * v + 0.055826318f) * v + 0.24015361f) * v + 0.69315308f) * v + 0.99999994f);
}

/*
 */
float Math::log2Fast(float v) {
	int i = IntFloat(v).i;
	int e = ((i >> 23) & 0xff) - 127;
	v = IntFloat((i & 0x007fffff) | 0x3f800000).f;
	return itof(e) + (v - 1.0f) * (((((-0.034436006f * v + 0.31821337f) * v - 1.2315303f) * v + 2.5988452f) * v - 3.3241990f) * v + 3.1157899f);
}

/*
 */
float Math::sinFast(float a) {
	float a2 = a * a;
	return a * (((((-2.39e-08f * a2 + 2.7526e-06f) * a2 - 1.98409e-04f) * a2 + 0.0083333315f) * a2 - 0.1666666664f) * a2 + 1.0f);
}

/*
 */
float Math::cosFast(float a) {
	float a2 = a * a;
	return (((((-2.605e-07f * a2 + 2.47609e-05f) * a2 - 0.0013888397f) * a2 + 0.0416666418f) * a2 - 0.4999999963f) * a2 + 1.0f);
}

/*
 */
float Math::asinFast(float v) {
    return PI05 - (((((((-0.0012624911f * v + 0.0066700901f) * v - 0.0170881256f) * v + 0.0308918810f) * v - 0.0501743046f) * v + 0.0889789874f) * v - 0.2145988016f) * v + 1.5707963050f) * Math::sqrt(1.0f - v);
}

/*
 */
float Math::acosFast(float v) {
	return (((((((-0.0012624911f * v + 0.0066700901f) * v - 0.0170881256f) * v + 0.0308918810f) * v - 0.0501743046f) * v + 0.0889789874f) * v - 0.2145988016f) * v + 1.5707963050f) * Math::sqrt(1.0f - v);
}

/*
 */
void Math::sincos(float angle,float &s,float &c) {
	#if defined(_WIN32) && defined(ARCH_X86)
		__asm {
			fld angle
			fsincos
			mov edx, s
			mov ecx, c
			fstp dword ptr [ecx]
			fstp dword ptr [edx]
		}
	#elif defined(_LINUX) && !defined(ARCH_ARM)
		double res_0,res_1;
		asm volatile("fsincos" : "=t"(res_0), "=u"(res_1) : "0"(angle));
		s = (float)res_1;
		c = (float)res_0;
	#else
		s = Math::sin(angle);
		c = Math::cos(angle);
	#endif
}

void Math::sincos(double angle,double &s,double &c) {
	#if defined(_WIN32) && defined(ARCH_X86)
		__asm {
			fld angle
			fsincos
			mov edx, s
			mov ecx, c
			fstp qword ptr [ecx]
			fstp qword ptr [edx]
		}
	#elif defined(_LINUX) && !defined(ARCH_ARM)
		asm volatile("fsincos" : "=t"(c), "=u"(s) : "0"(angle));
	#else
		s = Math::sin(angle);
		c = Math::cos(angle);
	#endif
}

/*
 */
template <class Type>
static INLINE Type bezier_func(const Type *curve,float k0) {
	float k1 = 1.0f - k0;
	float k00 = k0 * k0;
	float k11 = k1 * k1;
	return curve[0] * (k11 * k1) + 3.0f * (curve[1] * (k11 * k0) + curve[2] * (k00 * k1)) + curve[3] * (k00 * k0);
}

float Math::bezier(const float *times,const float *values,float time) {
	float k0 = 0.0f;
	float k1 = 1.0f;
	float t1 = bezier_func(times,k1) - time;
	if(abs(t1) < EPSILON) return bezier_func(values,k1);
	float t0 = bezier_func(times,k0) - time;
	if(abs(t0) < EPSILON || t0 * t1 > 0.0f) return bezier_func(values,k0);
	for(int i = 0; i < 20; i++) {
		float k = (k0 + k1) * 0.5f;
		float t = bezier_func(times,k) - time;
		float v = t0 * t;
		if(v < 0.0f) {
			k1 = k;
			t1 = t;
		} else if(v > 0.0f) {
			k0 = k;
			t0 = t;
		} else {
			return bezier_func(values,k);
		}
	}
	return bezier_func(values,k0);
}

double Math::bezier(const float *times,const double *values,float time) {
	float k0 = 0.0f;
	float k1 = 1.0f;
	float t1 = bezier_func(times,k1) - time;
	if(abs(t1) < EPSILON) return bezier_func(values,k1);
	float t0 = bezier_func(times,k0) - time;
	if(abs(t0) < EPSILON || t0 * t1 > 0.0f) return bezier_func(values,k0);
	for(int i = 0; i < 20; i++) {
		float k = (k0 + k1) * 0.5f;
		float t = bezier_func(times,k) - time;
		float v = t0 * t;
		if(v < 0.0f) {
			k1 = k;
			t1 = t;
		} else if(v > 0.0f) {
			k0 = k;
			t0 = t;
		} else {
			return bezier_func(values,k);
		}
	}
	return bezier_func(values,k0);
}

/******************************************************************************\
*
* Memory
*
\******************************************************************************/

/*
 */
#ifndef _WEBGL

/*
 */
static INLINE void math_memset_64(void *dest,unsigned int src,size_t size) {
	
	ASSERT_ALIGNED16(dest);
	
	#if defined(USE_SSE) || defined(USE_ALTIVEC)
		
		ATTRIBUTE_ALIGNED16(unsigned int data[4]);
		data[0] = src;
		data[1] = src;
		data[2] = src;
		data[3] = src;
		
		#ifdef USE_SSE
			__m128 res = _mm_load_ps((const float*)data);
		#elif USE_ALTIVEC
			vec_uchar16 res = vec_ld(0,(const unsigned char*)data);
		#endif
		
	#endif
	
	unsigned int *d = (unsigned int*)dest;
	for(size_t i = size >> 6; i > 0; i--) {
		
		#ifdef USE_SSE
			
			_mm_stream_ps((float*)d + 0,res);
			_mm_stream_ps((float*)d + 4,res);
			_mm_stream_ps((float*)d + 8,res);
			_mm_stream_ps((float*)d + 12,res);
			
		#elif USE_ALTIVEC
			
			vec_st(res,0,(unsigned char*)d);
			vec_st(res,16,(unsigned char*)d);
			vec_st(res,32,(unsigned char*)d);
			vec_st(res,48,(unsigned char*)d);
			
		#else
			
			d[0] = src;
			d[1] = src;
			d[2] = src;
			d[3] = src;
			d[4] = src;
			d[5] = src;
			d[6] = src;
			d[7] = src;
			d[8] = src;
			d[9] = src;
			d[10] = src;
			d[11] = src;
			d[12] = src;
			d[13] = src;
			d[14] = src;
			d[15] = src;
			
		#endif
		
		d += 16;
	}
}

/*
 */
void Math::memset(void *dest,int c,size_t size) {
	
	unsigned char *d = (unsigned char*)dest;
	
	unsigned int src = (c) | (c << 8) | (c << 16) | (c << 24);
	
	if(size & ~63) {
		
		while(IS_ALIGNED16(d) == 0) {
			*d++ = c;
			size--;
		}
		
		if(size & ~63) {
			math_memset_64(d,src,size);
			d += (size & ~63);
			size &= 63;
		}
	}
	
	if(size & ~15) {
		
		for(size_t i = size >> 4; i > 0; i--) {
			*(unsigned int*)(d + 0) = src;
			*(unsigned int*)(d + 4) = src;
			*(unsigned int*)(d + 8) = src;
			*(unsigned int*)(d + 12) = src;
			d += 16;
		}
		
		size &= 15;
	}
	
	if(size) {
		
		for(size_t i = size; i > 0; i--) {
			*d++ = c;
		}
	}
}

/*
 */
static INLINE void math_memcpy_64a(void *dest,const void *src,size_t size) {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(dest);
	
	unsigned int *d = (unsigned int*)dest;
	const unsigned int *s = (const unsigned int*)src;
	for(size_t i = size >> 6; i > 0; i--) {
		
		#ifdef USE_SSE
			
			_mm_prefetch((const char*)s + 512,_MM_HINT_NTA);
			
			__m128 res_0 = _mm_load_ps((const float*)s + 0);
			__m128 res_1 = _mm_load_ps((const float*)s + 4);
			__m128 res_2 = _mm_load_ps((const float*)s + 8);
			__m128 res_3 = _mm_load_ps((const float*)s + 12);
			
			_mm_stream_ps((float*)d + 0,res_0);
			_mm_stream_ps((float*)d + 4,res_1);
			_mm_stream_ps((float*)d + 8,res_2);
			_mm_stream_ps((float*)d + 12,res_3);
			
		#elif USE_ALTIVEC
			
			__dcbt((const unsigned char*)s + 512);
			
			vec_uchar16 res_0 = vec_ld(0,(const unsigned char*)s);
			vec_uchar16 res_1 = vec_ld(16,(const unsigned char*)s);
			vec_uchar16 res_2 = vec_ld(32,(const unsigned char*)s);
			vec_uchar16 res_3 = vec_ld(48,(const unsigned char*)s);
			
			vec_st(res_0,0,(unsigned char*)d);
			vec_st(res_1,16,(unsigned char*)d);
			vec_st(res_2,32,(unsigned char*)d);
			vec_st(res_3,48,(unsigned char*)d);
			
		#elif defined(USE_NEON) && (defined(_LINUX) || defined(_ANDROID))
			
			asm volatile(
				"pld  [ %r1, #512 ]			\n"
				"vldm   %r1, { d0 - d7 }	\n"
				"vstm   %r0, { d0 - d7 }	\n"
				: : "r"(d), "r"(s)
				: "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7"
			);
			
		#elif defined(ARCH_ARM) && (defined(_LINUX) || defined(_ANDROID))
			
			asm volatile(
				"mov   r0, %r0				\n"
				"mov   r1, %r1				\n"
				"pld [ r1, #512 ]			\n"
				"ldm   r1!, { r2 - r9 }		\n"
				"stm   r0!, { r2 - r9 }		\n"
				"ldm   r1,  { r2 - r9 }		\n"
				"stm   r0,  { r2 - r9 }		\n"
				: : "r"(d), "r"(s)
				: "r0", "r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9"
			);
			
		#else
			
			d[0] = s[0];
			d[1] = s[1];
			d[2] = s[2];
			d[3] = s[3];
			d[4] = s[4];
			d[5] = s[5];
			d[6] = s[6];
			d[7] = s[7];
			d[8] = s[8];
			d[9] = s[9];
			d[10] = s[10];
			d[11] = s[11];
			d[12] = s[12];
			d[13] = s[13];
			d[14] = s[14];
			d[15] = s[15];
			
		#endif
		
		s += 16;
		d += 16;
	}
}

static INLINE void math_memcpy_64u(void *dest,const void *src,size_t size) {
	
	ASSERT_ALIGNED16(dest);
	
	#ifdef USE_ALTIVEC
		vec_uchar16 mask = vec_lvsl(0,(const unsigned char*)src);
	#endif
	
	unsigned int *d = (unsigned int*)dest;
	const unsigned int *s = (const unsigned int*)src;
	for(size_t i = size >> 6; i > 0; i--) {
		
		#ifdef USE_SSE
			
			_mm_prefetch((const char*)s + 512,_MM_HINT_NTA);
			
			__m128 res_0 = _mm_loadu_ps((const float*)s + 0);
			__m128 res_1 = _mm_loadu_ps((const float*)s + 4);
			__m128 res_2 = _mm_loadu_ps((const float*)s + 8);
			__m128 res_3 = _mm_loadu_ps((const float*)s + 12);
			
			_mm_stream_ps((float*)d + 0,res_0);
			_mm_stream_ps((float*)d + 4,res_1);
			_mm_stream_ps((float*)d + 8,res_2);
			_mm_stream_ps((float*)d + 12,res_3);
			
		#elif USE_ALTIVEC
			
			__dcbt((const unsigned char*)s + 512);
			
			vec_uchar16 res_0 = vec_ld(0,(const unsigned char*)s);
			vec_uchar16 res_1 = vec_ld(15,(const unsigned char*)s);
			vec_uchar16 res_2 = vec_ld(31,(const unsigned char*)s);
			vec_uchar16 res_3 = vec_ld(47,(const unsigned char*)s);
			vec_uchar16 res_4 = vec_ld(63,(const unsigned char*)s);
			
			vec_st(vec_perm(res_0,res_1,mask),0,(unsigned char*)d);
			vec_st(vec_perm(res_1,res_2,mask),16,(unsigned char*)d);
			vec_st(vec_perm(res_2,res_3,mask),32,(unsigned char*)d);
			vec_st(vec_perm(res_3,res_4,mask),48,(unsigned char*)d);
			
		#else
			
			d[0] = s[0];
			d[1] = s[1];
			d[2] = s[2];
			d[3] = s[3];
			d[4] = s[4];
			d[5] = s[5];
			d[6] = s[6];
			d[7] = s[7];
			d[8] = s[8];
			d[9] = s[9];
			d[10] = s[10];
			d[11] = s[11];
			d[12] = s[12];
			d[13] = s[13];
			d[14] = s[14];
			d[15] = s[15];
			
		#endif
		
		s += 16;
		d += 16;
	}
}

/*
 */
void Math::memcpy(void *dest,const void *src,size_t size) {
	
	unsigned char *d = (unsigned char*)dest;
	const unsigned char *s = (const unsigned char*)src;
	
	if(size & ~63) {
		
		while(IS_ALIGNED16(d) == 0) {
			*d++ = *s++;
			size--;
		}
		
		if(size & ~63) {
			if(IS_ALIGNED16(s)) math_memcpy_64a(d,s,size);
			else math_memcpy_64u(d,s,size);
			s += (size & ~63);
			d += (size & ~63);
			size &= 63;
		}
	}
	
	if(size & ~15) {
		
		while(IS_ALIGNED4(d) == 0) {
			*d++ = *s++;
			size--;
		}
		
		for(size_t i = size >> 4; i > 0; i--) {
			*(unsigned int*)(d + 0) = *(const unsigned int*)(s + 0);
			*(unsigned int*)(d + 4) = *(const unsigned int*)(s + 4);
			*(unsigned int*)(d + 8) = *(const unsigned int*)(s + 8);
			*(unsigned int*)(d + 12) = *(const unsigned int*)(s + 12);
			s += 16;
			d += 16;
		}
		
		size &= 15;
	}
	
	if(size & ~3) {
		
		for(size_t i = size >> 2; i > 0; i--) {
			*(unsigned int*)d = *(const unsigned int*)s;
			s += 4;
			d += 4;
		}
		
		size &= 3;
	}
	
	if(size) {
		
		for(size_t i = size; i > 0; i--) {
			*d++ = *s++;
		}
	}
}

/*
 */
int Math::memcmp(const void *src_0,const void *src_1,size_t size) {
	
	const unsigned char *s0 = (const unsigned char*)src_0;
	const unsigned char *s1 = (const unsigned char*)src_1;
	
	if(size & ~3) {
		
		for(size_t i = size >> 2; i > 0; i--) {
			unsigned int i0 = *(unsigned int*)s0;
			unsigned int i1 = *(unsigned int*)s1;
			if(i0 != i1) {
				for(int j = 0; j < 4; j++) {
					unsigned char c0 = *s0++;
					unsigned char c1 = *s1++;
					if(c0 < c1) return -1;
					if(c0 > c1) return 1;
				}
			}
			s0 += 4;
			s1 += 4;
		}
		
		size &= 3;
	}
	
	if(size) {
		
		for(size_t i = size; i > 0; i--) {
			unsigned char c0 = *s0++;
			unsigned char c1 = *s1++;
			if(c0 < c1) return -1;
			if(c0 > c1) return 1;
		}
	}
	
	return 0;
}

#endif /* _WEBGL */

/******************************************************************************\
*
* vec2
*
\******************************************************************************/

/*
 */
const vec2 vec2_zero(0.0f);
const vec2 vec2_one(1.0f);
const vec2 vec2_epsilon(EPSILON);
const vec2 vec2_infinity(INFINITY);

/*
 */
vec2 min(const vec2 &v0,const vec2 &v1) {
	vec2 ret;
	#ifdef USE_NEON
		ret.vec = vmin_f32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
	#endif
	return ret;
}

vec2 max(const vec2 &v0,const vec2 &v1) {
	vec2 ret;
	#ifdef USE_NEON
		ret.vec = vmax_f32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
	#endif
	return ret;
}

vec2 clamp(const vec2 &v,const vec2 &v0,const vec2 &v1) {
	vec2 ret;
	#ifdef USE_NEON
		ret.vec = vmin_f32(vmax_f32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
	#endif
	return ret;
}

vec2 saturate(const vec2 &v) {
	vec2 ret;
	#ifdef USE_NEON
		ret.vec = vmin_f32(vmax_f32(v.vec,vec2_zero.vec),vec2_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
	#endif
	return ret;
}

vec2 lerp(const vec2 &v0,const vec2 &v1,float k) {
	vec2 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* vec3
*
\******************************************************************************/

/*
 */
const vec3 vec3_zero(0.0f);
const vec3 vec3_one(1.0f);
const vec3 vec3_epsilon(EPSILON);
const vec3 vec3_infinity(INFINITY);

/*
 */
vec3 min(const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(v0.vec,v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(v0.vec,v1.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
	#endif
	return ret;
}

vec3 max(const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	#ifdef USE_SSE
		ret.vec = _mm_max_ps(v0.vec,v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_max(v0.vec,v1.vec);
	#elif USE_NEON
		ret.vec = vmaxq_f32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
	#endif
	return ret;
}

vec3 clamp(const vec3 &v,const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(_mm_max_ps(v.vec,v0.vec),v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(vec_max(v.vec,v0.vec),v1.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(vmaxq_f32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
	#endif
	return ret;
}

vec3 saturate(const vec3 &v) {
	vec3 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(_mm_max_ps(v.vec,vec3_zero.vec),vec3_one.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(vec_max(v.vec,vec3_zero.vec),vec3_one.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(vmaxq_f32(v.vec,vec3_zero.vec),vec3_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
		ret.z = saturate(v.z);
	#endif
	return ret;
}

vec3 lerp(const vec3 &v0,const vec3 &v1,float k) {
	vec3 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* vec4
*
\******************************************************************************/

/*
 */
const vec4 vec4_zero(0.0f);
const vec4 vec4_one(1.0f);
const vec4 vec4_epsilon(EPSILON);
const vec4 vec4_infinity(INFINITY);

/*
 */
vec4 min(const vec4 &v0,const vec4 &v1) {
	vec4 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(v0.vec,v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(v0.vec,v1.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
		ret.w = min(v0.w,v1.w);
	#endif
	return ret;
}

vec4 max(const vec4 &v0,const vec4 &v1) {
	vec4 ret;
	#ifdef USE_SSE
		ret.vec = _mm_max_ps(v0.vec,v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_max(v0.vec,v1.vec);
	#elif USE_NEON
		ret.vec = vmaxq_f32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
		ret.w = max(v0.w,v1.w);
	#endif
	return ret;
}

vec4 clamp(const vec4 &v,const vec4 &v0,const vec4 &v1) {
	vec4 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(_mm_max_ps(v.vec,v0.vec),v1.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(vec_max(v.vec,v0.vec),v1.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(vmaxq_f32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
		ret.w = clamp(v.w,v0.w,v1.w);
	#endif
	return ret;
}

vec4 saturate(const vec4 &v) {
	vec4 ret;
	#ifdef USE_SSE
		ret.vec = _mm_min_ps(_mm_max_ps(v.vec,vec4_zero.vec),vec4_one.vec);
	#elif USE_ALTIVEC
		ret.vec = vec_min(vec_max(v.vec,vec4_zero.vec),vec4_one.vec);
	#elif USE_NEON
		ret.vec = vminq_f32(vmaxq_f32(v.vec,vec4_zero.vec),vec4_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
		ret.z = saturate(v.z);
		ret.w = saturate(v.w);
	#endif
	return ret;
}

vec4 lerp(const vec4 &v0,const vec4 &v1,float k) {
	vec4 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* dvec2
*
\******************************************************************************/

/*
 */
const dvec2 dvec2_zero(0.0);
const dvec2 dvec2_one(1.0);
const dvec2 dvec2_epsilon(EPSILON);
const dvec2 dvec2_infinity(INFINITY);

/*
 */
dvec2 min(const dvec2 &v0,const dvec2 &v1) {
	dvec2 ret;
	#ifdef USE_SSE2
		ret.vec = _mm_min_pd(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
	#endif
	return ret;
}

dvec2 max(const dvec2 &v0,const dvec2 &v1) {
	dvec2 ret;
	#ifdef USE_SSE2
		ret.vec = _mm_max_pd(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
	#endif
	return ret;
}

dvec2 clamp(const dvec2 &v,const dvec2 &v0,const dvec2 &v1) {
	dvec2 ret;
	#ifdef USE_SSE2
		ret.vec = _mm_min_pd(_mm_max_pd(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
	#endif
	return ret;
}

dvec2 saturate(const dvec2 &v) {
	dvec2 ret;
	#ifdef USE_SSE2
		ret.vec = _mm_min_pd(_mm_max_pd(v.vec,dvec2_zero.vec),dvec2_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
	#endif
	return ret;
}

dvec2 lerp(const dvec2 &v0,const dvec2 &v1,double k) {
	dvec2 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* dvec3
*
\******************************************************************************/

/*
 */
const dvec3 dvec3_zero(0.0);
const dvec3 dvec3_one(1.0);
const dvec3 dvec3_epsilon(EPSILON);
const dvec3 dvec3_infinity(INFINITY);

/*
 */
dvec3 min(const dvec3 &v0,const dvec3 &v1) {
	dvec3 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(v0.vec0,v1.vec0);
		ret.vec1 = _mm_min_sd(v0.vec1,v1.vec1);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
	#endif
	return ret;
}

dvec3 max(const dvec3 &v0,const dvec3 &v1) {
	dvec3 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_max_pd(v0.vec0,v1.vec0);
		ret.vec1 = _mm_max_sd(v0.vec1,v1.vec1);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
	#endif
	return ret;
}

dvec3 clamp(const dvec3 &v,const dvec3 &v0,const dvec3 &v1) {
	dvec3 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(_mm_max_pd(v.vec0,v0.vec0),v1.vec0);
		ret.vec1 = _mm_min_sd(_mm_max_sd(v.vec1,v0.vec1),v1.vec1);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
	#endif
	return ret;
}

dvec3 saturate(const dvec3 &v) {
	dvec3 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(_mm_max_pd(v.vec0,dvec2_zero.vec),dvec2_one.vec);
		ret.vec1 = _mm_min_sd(_mm_max_sd(v.vec1,dvec2_zero.vec),dvec2_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
		ret.z = saturate(v.z);
	#endif
	return ret;
}

dvec3 lerp(const dvec3 &v0,const dvec3 &v1,double k) {
	dvec3 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* dvec4
*
\******************************************************************************/

/*
 */
const dvec4 dvec4_zero(0.0);
const dvec4 dvec4_one(1.0);
const dvec4 dvec4_epsilon(EPSILON);
const dvec4 dvec4_infinity(INFINITY);

/*
 */
dvec4 min(const dvec4 &v0,const dvec4 &v1) {
	dvec4 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(v0.vec0,v1.vec0);
		ret.vec1 = _mm_min_pd(v0.vec1,v1.vec1);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
		ret.w = min(v0.w,v1.w);
	#endif
	return ret;
}

dvec4 max(const dvec4 &v0,const dvec4 &v1) {
	dvec4 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_max_pd(v0.vec0,v1.vec0);
		ret.vec1 = _mm_max_pd(v0.vec1,v1.vec1);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
		ret.w = max(v0.w,v1.w);
	#endif
	return ret;
}

dvec4 clamp(const dvec4 &v,const dvec4 &v0,const dvec4 &v1) {
	dvec4 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(_mm_max_pd(v.vec0,v0.vec0),v1.vec0);
		ret.vec1 = _mm_min_pd(_mm_max_pd(v.vec1,v0.vec1),v1.vec1);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
		ret.w = clamp(v.w,v0.w,v1.w);
	#endif
	return ret;
}

dvec4 saturate(const dvec4 &v) {
	dvec4 ret;
	#ifdef USE_SSE2
		ret.vec0 = _mm_min_pd(_mm_max_pd(v.vec0,dvec2_zero.vec),dvec2_one.vec);
		ret.vec1 = _mm_min_pd(_mm_max_pd(v.vec1,dvec2_zero.vec),dvec2_one.vec);
	#else
		ret.x = saturate(v.x);
		ret.y = saturate(v.y);
		ret.z = saturate(v.z);
		ret.w = saturate(v.w);
	#endif
	return ret;
}

dvec4 lerp(const dvec4 &v0,const dvec4 &v1,double k) {
	dvec4 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* hvec2
*
\******************************************************************************/

/*
 */
const hvec2 hvec2_zero(0.0f);
const hvec2 hvec2_one(1.0f);

/******************************************************************************\
*
* hvec3
*
\******************************************************************************/

/*
 */
const hvec3 hvec3_zero(0.0f);
const hvec3 hvec3_one(1.0f);

/******************************************************************************\
*
* hvec4
*
\******************************************************************************/

/*
 */
const hvec4 hvec4_zero(0.0f);
const hvec4 hvec4_one(1.0f);

/******************************************************************************\
*
* ivec2
*
\******************************************************************************/

/*
 */
const ivec2 ivec2_zero(0);
const ivec2 ivec2_one(1);

/*
 */
ivec2 min(const ivec2 &v0,const ivec2 &v1) {
	ivec2 ret;
	#ifdef USE_NEON
		ret.vec = vmin_s32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
	#endif
	return ret;
}

ivec2 max(const ivec2 &v0,const ivec2 &v1) {
	ivec2 ret;
	#ifdef USE_NEON
		ret.vec = vmax_s32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
	#endif
	return ret;
}

ivec2 clamp(const ivec2 &v,const ivec2 &v0,const ivec2 &v1) {
	ivec2 ret;
	#ifdef USE_NEON
		ret.vec = vmin_s32(vmax_s32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
	#endif
	return ret;
}

ivec2 lerp(const ivec2 &v0,const ivec2 &v1,int k) {
	ivec2 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* ivec3
*
\******************************************************************************/

/*
 */
const ivec3 ivec3_zero(0);
const ivec3 ivec3_one(1);

/*
 */
ivec3 min(const ivec3 &v0,const ivec3 &v1) {
	ivec3 ret;
	#ifdef USE_NEON
		ret.vec = vminq_s32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
	#endif
	return ret;
}

ivec3 max(const ivec3 &v0,const ivec3 &v1) {
	ivec3 ret;
	#ifdef USE_NEON
		ret.vec = vmaxq_s32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
	#endif
	return ret;
}

ivec3 clamp(const ivec3 &v,const ivec3 &v0,const ivec3 &v1) {
	ivec3 ret;
	#ifdef USE_NEON
		ret.vec = vminq_s32(vmaxq_s32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
	#endif
	return ret;
}

ivec3 lerp(const ivec3 &v0,const ivec3 &v1,int k) {
	ivec3 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* ivec4
*
\******************************************************************************/

/*
 */
const ivec4 ivec4_zero(0);
const ivec4 ivec4_one(1);

/*
 */
ivec4 min(const ivec4 &v0,const ivec4 &v1) {
	ivec4 ret;
	#ifdef USE_NEON
		ret.vec = vminq_s32(v0.vec,v1.vec);
	#else
		ret.x = min(v0.x,v1.x);
		ret.y = min(v0.y,v1.y);
		ret.z = min(v0.z,v1.z);
		ret.w = min(v0.w,v1.w);
	#endif
	return ret;
}

ivec4 max(const ivec4 &v0,const ivec4 &v1) {
	ivec4 ret;
	#ifdef USE_NEON
		ret.vec = vmaxq_s32(v0.vec,v1.vec);
	#else
		ret.x = max(v0.x,v1.x);
		ret.y = max(v0.y,v1.y);
		ret.z = max(v0.z,v1.z);
		ret.w = max(v0.w,v1.w);
	#endif
	return ret;
}

ivec4 clamp(const ivec4 &v,const ivec4 &v0,const ivec4 &v1) {
	ivec4 ret;
	#ifdef USE_NEON
		ret.vec = vminq_s32(vmaxq_s32(v.vec,v0.vec),v1.vec);
	#else
		ret.x = clamp(v.x,v0.x,v1.x);
		ret.y = clamp(v.y,v0.y,v1.y);
		ret.z = clamp(v.z,v0.z,v1.z);
		ret.w = clamp(v.w,v0.w,v1.w);
	#endif
	return ret;
}

ivec4 lerp(const ivec4 &v0,const ivec4 &v1,int k) {
	ivec4 ret;
	return lerp(ret,v0,v1,k);
}

/******************************************************************************\
*
* bvec4
*
\******************************************************************************/

/*
 */
const bvec4 bvec4_zero((unsigned char)0);
const bvec4 bvec4_one((unsigned char)255);

/*
 */
bvec4 min(const bvec4 &v0,const bvec4 &v1) {
	bvec4 ret;
	ret.x = min(v0.x,v1.x);
	ret.y = min(v0.y,v1.y);
	ret.z = min(v0.z,v1.z);
	ret.w = min(v0.w,v1.w);
	return ret;
}

bvec4 max(const bvec4 &v0,const bvec4 &v1) {
	bvec4 ret;
	ret.x = max(v0.x,v1.x);
	ret.y = max(v0.y,v1.y);
	ret.z = max(v0.z,v1.z);
	ret.w = max(v0.w,v1.w);
	return ret;
}

bvec4 clamp(const bvec4 &v,const bvec4 &v0,const bvec4 &v1) {
	bvec4 ret;
	ret.x = clamp(v.x,v0.x,v1.x);
	ret.y = clamp(v.y,v0.y,v1.y);
	ret.z = clamp(v.z,v0.z,v1.z);
	ret.w = clamp(v.w,v0.w,v1.w);
	return ret;
}

/******************************************************************************\
*
* mat2
*
\******************************************************************************/

/*
 */
static const float mat2_identity_values[4] = {
	1.0f, 0.0f,
	0.0f, 1.0f,
};

/*
 */
const mat2 mat2_zero(0.0f);
const mat2 mat2_one(1.0f);
const mat2 mat2_identity(mat2_identity_values);

/*
 */
mat2::mat2(float v) {
	m00 = v; m01 = v;
	m10 = v; m11 = v;
}

mat2::mat2(const mat2 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

mat2::mat2(const mat3 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

mat2::mat2(const mat4 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

mat2::mat2(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01;
	m10 = (float)m.m10; m11 = (float)m.m11;
}

mat2::mat2(const float *m) {
	m00 = m[0]; m01 = m[2];
	m10 = m[1]; m11 = m[3];
}

/*
 */
mat2 mat2::operator-() const {
	mat2 ret;
	ret.m00 = -m00; ret.m01 = -m01;
	ret.m10 = -m10; ret.m11 = -m11;
	return ret;
}

mat2 &mat2::operator*=(float v) {
	return mul(*this,*this,v);
}

mat2 &mat2::operator*=(const mat2 &m) {
	return mul(*this,mat2(*this),m);
}

mat2 &mat2::operator+=(const mat2 &m) {
	return add(*this,*this,m);
}

mat2 &mat2::operator-=(const mat2 &m) {
	return sub(*this,*this,m);
}

/*
 */
void mat2::setRow(int row,const vec2 &v) {
	assert((unsigned int)row < 2 && "mat2::setRow(): bad row");
	mat[row + 0] = v.x;
	mat[row + 2] = v.y;
}

vec2 mat2::getRow(int row) const {
	assert((unsigned int)row < 2 && "mat2::getRow(): bad row");
	return vec2(mat[row + 0],mat[row + 2]);
}

/*
 */
void mat2::setColumn(int column,const vec2 &v) {
	assert((unsigned int)column < 2 && "mat2::setColumn(): bad column");
	column *= 2;
	mat[column + 0] = v.x;
	mat[column + 1] = v.y;
}

vec2 mat2::getColumn(int column) const {
	assert((unsigned int)column < 2 && "mat2::getColumn(): bad column");
	return vec2(mat + column * 2);
}

/*
 */
void mat2::set(const mat2 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

void mat2::set(const mat3 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

void mat2::set(const mat4 &m) {
	m00 = m.m00; m01 = m.m01;
	m10 = m.m10; m11 = m.m11;
}

void mat2::set(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01;
	m10 = (float)m.m10; m11 = (float)m.m11;
}

void mat2::set(const float *m) {
	m00 = m[0]; m01 = m[2];
	m10 = m[1]; m11 = m[3];
}

void mat2::get(float *m) const {
	m[0] = m00; m[2] = m01;
	m[1] = m10; m[3] = m11;
}

/*
 */
void mat2::setZero() {
	m00 = 0.0f; m01 = 0.0f;
	m10 = 0.0f; m11 = 0.0f;
}

void mat2::setIdentity() {
	m00 = 1.0f; m01 = 0.0f;
	m10 = 0.0f; m11 = 1.0f;
}

void mat2::setRotate(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c; m01 = -s;
	m10 = s; m11 = c;
}

void mat2::setScale(const vec2 &v) {
	m00 = v.x;  m01 = 0.0f;
	m10 = 0.0f; m11 = v.y;
}

/*
 */
int operator==(const mat2 &m0,const mat2 &m1) {
	return compare(m0,m1);
}

int operator!=(const mat2 &m0,const mat2 &m1) {
	return !compare(m0,m1);
}

mat2 operator*(const mat2 &m,float v) {
	mat2 ret;
	return mul(ret,m,v);
}

vec2 operator*(const mat2 &m,const vec2 &v) {
	vec2 ret;
	return mul(ret,m,v);
}

vec2 operator*(const vec2 &v,const mat2 &m) {
	vec2 ret;
	return mul(ret,v,m);
}

dvec2 operator*(const mat2 &m,const dvec2 &v) {
	dvec2 ret;
	return mul(ret,m,v);
}

dvec2 operator*(const dvec2 &v,const mat2 &m) {
	dvec2 ret;
	return mul(ret,v,m);
}

mat2 operator*(const mat2 &m0,const mat2 &m1) {
	mat2 ret;
	return mul(ret,m0,m1);
}

mat2 operator+(const mat2 &m0,const mat2 &m1) {
	mat2 ret;
	return add(ret,m0,m1);
}

mat2 operator-(const mat2 &m0,const mat2 &m1) {
	mat2 ret;
	return sub(ret,m0,m1);
}

/*
 */
int compare(const mat2 &m0,const mat2 &m1) {
	return (compare(m0.m00,m1.m00) && compare(m0.m10,m1.m10)) &&
		compare(m0.m01,m1.m01) && compare(m0.m11,m1.m11);
}

int compare(const mat2 &m0,const mat2 &m1,float epsilon) {
	return (compare(m0.m00,m1.m00,epsilon) && compare(m0.m10,m1.m10,epsilon)) &&
		compare(m0.m01,m1.m01,epsilon) && compare(m0.m11,m1.m11,epsilon);
}

float trace(const mat2 &m) {
	return m.m00 + m.m11;
}

float determinant(const mat2 &m) {
	return m.m00 * m.m11 - m.m10 * m.m01;
}

mat2 &mul(mat2 &ret,const mat2 &m,float v) {
	ret.m00 = m.m00 * v; ret.m01 = m.m01 * v;
	ret.m10 = m.m10 * v; ret.m11 = m.m11 * v;
	return ret;
}

vec2 &mul(vec2 &ret,const mat2 &m,const vec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y;
	ret.y = m.m10 * v.x + m.m11 * v.y;
	return ret;
}

vec2 &mul(vec2 &ret,const vec2 &v,const mat2 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

dvec2 &mul(dvec2 &ret,const mat2 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y;
	ret.y = m.m10 * v.x + m.m11 * v.y;
	return ret;
}

dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat2 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

mat2 &mul(mat2 &ret,const mat2 &m0,const mat2 &m1) {
	ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10;
	ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10;
	ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11;
	ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11;
	return ret;
}

mat2 &add(mat2 &ret,const mat2 &m0,const mat2 &m1) {
	ret.m00 = m0.m00 + m1.m00; ret.m01 = m0.m01 + m1.m01;
	ret.m10 = m0.m10 + m1.m10; ret.m11 = m0.m11 + m1.m11;
	return ret;
}

mat2 &sub(mat2 &ret,const mat2 &m0,const mat2 &m1) {
	ret.m00 = m0.m00 - m1.m00; ret.m01 = m0.m01 - m1.m01;
	ret.m10 = m0.m10 - m1.m10; ret.m11 = m0.m11 - m1.m11;
	return ret;
}

mat2 &transpose(mat2 &ret,const mat2 &m) {
	ret.m00 = m.m00; ret.m01 = m.m10;
	ret.m10 = m.m01; ret.m11 = m.m11;
	return ret;
}

mat2 &inverse(mat2 &ret,const mat2 &m) {
	float idet = Math::rcp(determinant(m));
	ret.m00 =  m.m11 * idet;
	ret.m10 = -m.m10 * idet;
	ret.m01 = -m.m01 * idet;
	ret.m11 =  m.m00 * idet;
	return ret;
}

mat2 &inverse(mat2 &ret,const mat2 &m,float det) {
	float idet = Math::rcp(det);
	ret.m00 =  m.m11 * idet;
	ret.m10 = -m.m10 * idet;
	ret.m01 = -m.m01 * idet;
	ret.m11 =  m.m00 * idet;
	return ret;
}

/*
 */
mat2 transpose(const mat2 &m) {
	mat2 ret;
	return transpose(ret,m);
}

mat2 inverse(const mat2 &m) {
	mat2 ret;
	return inverse(ret,m);
}

mat2 inverse(const mat2 &m,float det) {
	mat2 ret;
	return inverse(ret,m,det);
}

/******************************************************************************\
*
* mat3
*
\******************************************************************************/

/*
 */
static const float mat3_identity_values[12] = {
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f,
};

/*
 */
const mat3 mat3_zero(0.0f);
const mat3 mat3_one(1.0f);
const mat3 mat3_identity(mat3_identity_values);

/*
 */
mat3::mat3(float v) {
	m00 = v;    m01 = v;    m02 = v;
	m10 = v;    m11 = v;    m12 = v;
	m20 = v;    m21 = v;    m22 = v;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f;
}

mat3::mat3(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = 0.0f;
	m20 = 0.0f;  m21 = 0.0f;  m22 = 1.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;
}

mat3::mat3(const mat3 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
		m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;
	#endif
}

mat3::mat3(const mat4 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
		m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;
	#endif
}

mat3::mat3(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01; m02 = (float)m.m02;
	m10 = (float)m.m10; m11 = (float)m.m11; m12 = (float)m.m12;
	m20 = (float)m.m20; m21 = (float)m.m21; m22 = (float)m.m22;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;
}

mat3::mat3(const quat &q) {
	mat3 m = q.getMat3();
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
		m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;
	#endif
}

mat3::mat3(const float *m) {
	m00 = m[0]; m01 = m[4]; m02 = m[8];
	m10 = m[1]; m11 = m[5]; m12 = m[9];
	m20 = m[2]; m21 = m[6]; m22 = m[10];
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f;
}

/*
 */
mat3 mat3::operator-() const {
	mat3 ret;
	ret.m00 = -m00; ret.m01 = -m01; ret.m02 = -m02;
	ret.m10 = -m10; ret.m11 = -m11; ret.m12 = -m12;
	ret.m20 = -m20; ret.m21 = -m21; ret.m22 = -m22;
	return ret;
}

mat3 &mat3::operator*=(float v) {
	return mul(*this,*this,v);
}

mat3 &mat3::operator*=(const mat3 &m) {
	return mul(*this,mat3(*this),m);
}

mat3 &mat3::operator+=(const mat3 &m) {
	return add(*this,*this,m);
}

mat3 &mat3::operator-=(const mat3 &m) {
	return sub(*this,*this,m);
}

/*
 */
void mat3::setRow(int row,const vec3 &v) {
	assert((unsigned int)row < 3 && "mat3::setRow(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
}

vec3 mat3::getRow(int row) const {
	assert((unsigned int)row < 3 && "mat3::getRow(): bad row");
	return vec3(mat[row + 0],mat[row + 4],mat[row + 8]);
}

/*
 */
void mat3::setColumn(int column,const vec3 &v) {
	assert((unsigned int)column < 3 && "mat3::setColumn(): bad column");
	column *= 4;
	mat[column + 0] = v.x;
	mat[column + 1] = v.y;
	mat[column + 2] = v.z;
}

vec3 mat3::getColumn(int column) const {
	assert((unsigned int)column < 3 && "mat3::getColumn(): bad column");
	return vec3(mat + column * 4);
}

/*
 */
void mat3::setDiagonal(const vec3 &v) {
	m00 = v.x;  m01 = 0.0f; m02 = 0.0f;
	m10 = 0.0f; m11 = v.y;  m12 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = v.z;
}

vec3 mat3::getDiagonal() const {
	return vec3(m00,m11,m22);
}

/*
 */
void mat3::set(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = 0.0f;
	m20 = 0.0f;  m21 = 0.0f;  m22 = 1.0f;
}

void mat3::set(const mat3 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
	#endif
}

void mat3::set(const mat4 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
	#endif
}

void mat3::set(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01; m02 = (float)m.m02;
	m10 = (float)m.m10; m11 = (float)m.m11; m12 = (float)m.m12;
	m20 = (float)m.m20; m21 = (float)m.m21; m22 = (float)m.m22;
}

void mat3::set(const quat &q) {
	mat3 m = q.getMat3();
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02;
		m10 = m.m10; m11 = m.m11; m12 = m.m12;
		m20 = m.m20; m21 = m.m21; m22 = m.m22;
	#endif
}

void mat3::set(const float *m) {
	m00 = m[0]; m01 = m[4]; m02 = m[8];
	m10 = m[1]; m11 = m[5]; m12 = m[9];
	m20 = m[2]; m21 = m[6]; m22 = m[10];
}

void mat3::get(float *m) const {
	m[0] = m00; m[4] = m01; m[8] = m02;
	m[1] = m10; m[5] = m11; m[9] = m12;
	m[2] = m20; m[6] = m21; m[10] = m22;
}

/*
 */
void mat3::setZero() {
	m00 = 0.0f; m01 = 0.0f; m02 = 0.0f;
	m10 = 0.0f; m11 = 0.0f; m12 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 0.0f;
}

void mat3::setIdentity() {
	m00 = 1.0f; m01 = 0.0f; m02 = 0.0f;
	m10 = 0.0f; m11 = 1.0f; m12 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 1.0f;
}

void mat3::setSkewSymmetric(const vec3 &v) {
	m00 = 0.0f; m01 = -v.z; m02 =  v.y;
	m10 =  v.z; m11 =  0.0; m12 = -v.x;
	m20 = -v.y; m21 =  v.x; m22 =  0.0;
}

void mat3::setRotate(const vec3 &axis,float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	vec3 v = normalize(axis);
	float xx = v.x * v.x;
	float yy = v.y * v.y;
	float zz = v.z * v.z;
	float xy = v.x * v.y;
	float yz = v.y * v.z;
	float zx = v.z * v.x;
	float xs = v.x * s;
	float ys = v.y * s;
	float zs = v.z * s;
	m00 = (1.0f - c) * xx + c;  m01 = (1.0f - c) * xy - zs; m02 = (1.0f - c) * zx + ys;
	m10 = (1.0f - c) * xy + zs; m11 = (1.0f - c) * yy + c;  m12 = (1.0f - c) * yz - xs;
	m20 = (1.0f - c) * zx - ys; m21 = (1.0f - c) * yz + xs; m22 = (1.0f - c) * zz + c;
}

void mat3::setRotateX(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = 1.0f; m01 = 0.0f; m02 = 0.0f;
	m10 = 0.0f; m11 = c;    m12 = -s;
	m20 = 0.0f; m21 = s;    m22 = c;
}

void mat3::setRotateY(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;    m01 = 0.0f; m02 = s;
	m10 = 0.0f; m11 = 1.0f; m12 = 0.0f;
	m20 = -s;   m21 = 0.0f; m22 = c;
}

void mat3::setRotateZ(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;    m01 = -s;   m02 = 0.0f;
	m10 = s;    m11 = c;    m12 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 1.0f;
}

void mat3::setScale(const vec3 &v) {
	m00 = v.x;  m01 = 0.0f; m02 = 0.0f;
	m10 = 0.0f; m11 = v.y;  m12 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = v.z;
}

/*
 */
quat mat3::getQuat() const {
	quat ret;
	float s,trace = m00 + m11 + m22;
	if(trace > 0.0f) {
		s = Math::sqrt(trace + 1.0f);
		ret.w = s * 0.5f;
		s = Math::rcp(s) * 0.5f;
		ret.x = (m21 - m12) * s;
		ret.y = (m02 - m20) * s;
		ret.z = (m10 - m01) * s;
	} else {
		int i = 0;
		if(m11 > m00) i = 1;
		if(m22 > mat[4 * i + i]) i = 2;
		switch(i) {
			case 0:
				s = Math::sqrt(m00 - m11 - m22 + 1.0f);
				ret.x = s * 0.5f;
				s = Math::rcp(s) * 0.5f;
				ret.y = (m10 + m01) * s;
				ret.z = (m02 + m20) * s;
				ret.w = (m21 - m12) * s;
				break;
			case 1:
				s = Math::sqrt(m11 - m00 - m22 + 1.0f);
				ret.y = s * 0.5f;
				s = Math::rcp(s) * 0.5f;
				ret.z = (m21 + m12) * s;
				ret.x = (m10 + m01) * s;
				ret.w = (m02 - m20) * s;
				break;
			case 2:
				s = Math::sqrt(m22 - m00 - m11 + 1.0f);
				ret.z = s * 0.5f;
				s = Math::rcp(s) * 0.5f;
				ret.x = (m02 + m20) * s;
				ret.y = (m21 + m12) * s;
				ret.w = (m10 - m01) * s;
				break;
		}
	}
	return ret;
}

/*
 */
int operator==(const mat3 &m0,const mat3 &m1) {
	return compare(m0,m1);
}

int operator!=(const mat3 &m0,const mat3 &m1) {
	return !compare(m0,m1);
}

mat3 operator*(const mat3 &m,float v) {
	mat3 ret;
	return mul(ret,m,v);
}

vec2 operator*(const mat3 &m,const vec2 &v) {
	vec2 ret;
	return mul(ret,m,v);
}

vec2 operator*(const vec2 &v,const mat3 &m) {
	vec2 ret;
	return mul(ret,v,m);
}

vec3 operator*(const mat3 &m,const vec3 &v) {
	vec3 ret;
	return mul(ret,m,v);
}

vec3 operator*(const vec3 &v,const mat3 &m) {
	vec3 ret;
	return mul(ret,v,m);
}

dvec2 operator*(const mat3 &m,const dvec2 &v) {
	dvec2 ret;
	return mul(ret,m,v);
}

dvec2 operator*(const dvec2 &v,const mat3 &m) {
	dvec2 ret;
	return mul(ret,v,m);
}

dvec3 operator*(const mat3 &m,const dvec3 &v) {
	dvec3 ret;
	return mul(ret,m,v);
}

dvec3 operator*(const dvec3 &v,const mat3 &m) {
	dvec3 ret;
	return mul(ret,v,m);
}

mat3 operator*(const mat3 &m0,const mat3 &m1) {
	mat3 ret;
	return mul(ret,m0,m1);
}

mat3 operator+(const mat3 &m0,const mat3 &m1) {
	mat3 ret;
	return add(ret,m0,m1);
}

mat3 operator-(const mat3 &m0,const mat3 &m1) {
	mat3 ret;
	return sub(ret,m0,m1);
}

/*
 */
int compare(const mat3 &m0,const mat3 &m1) {
	return (compare(m0.m00,m1.m00) && compare(m0.m10,m1.m10) && compare(m0.m20,m1.m20) &&
		compare(m0.m01,m1.m01) && compare(m0.m11,m1.m11) && compare(m0.m21,m1.m21) &&
		compare(m0.m02,m1.m02) && compare(m0.m12,m1.m12) && compare(m0.m22,m1.m22));
}

int compare(const mat3 &m0,const mat3 &m1,float epsilon) {
	return (compare(m0.m00,m1.m00,epsilon) && compare(m0.m10,m1.m10,epsilon) && compare(m0.m20,m1.m20,epsilon) &&
		compare(m0.m01,m1.m01,epsilon) && compare(m0.m11,m1.m11,epsilon) && compare(m0.m21,m1.m21,epsilon) &&
		compare(m0.m02,m1.m02,epsilon) && compare(m0.m12,m1.m12,epsilon) && compare(m0.m22,m1.m22,epsilon));
}

float trace(const mat3 &m) {
	return m.m00 + m.m11 + m.m22;
}

float determinant(const mat3 &m) {
	float det = 0.0f;
	det =  m.m00 * (m.m11 * m.m22 - m.m12 * m.m21);
	det -= m.m01 * (m.m10 * m.m22 - m.m12 * m.m20);
	det += m.m02 * (m.m10 * m.m21 - m.m11 * m.m20);
	return det;
}

mat3 &mul(mat3 &ret,const mat3 &m,float v) {
	#ifdef USE_SSE
		__m128 temp = _mm_set1_ps(v);
		ret.col0 = _mm_mul_ps(m.col0,temp);
		ret.col1 = _mm_mul_ps(m.col1,temp);
		ret.col2 = _mm_mul_ps(m.col2,temp);
	#elif USE_ALTIVEC
		vec_float4 temp = vec_splats(v);
		vec_float4 zero = vec_splats(0.0f);
		ret.col0 = vec_madd(m.col0,temp,zero);
		ret.col1 = vec_madd(m.col1,temp,zero);
		ret.col2 = vec_madd(m.col2,temp,zero);
	#elif USE_NEON
		ret.col0 = vmulq_n_f32(m.col0,v);
		ret.col1 = vmulq_n_f32(m.col1,v);
		ret.col2 = vmulq_n_f32(m.col2,v);
	#else
		ret.m00 = m.m00 * v; ret.m01 = m.m01 * v; ret.m02 = m.m02 * v;
		ret.m10 = m.m10 * v; ret.m11 = m.m11 * v; ret.m12 = m.m12 * v;
		ret.m20 = m.m20 * v; ret.m21 = m.m21 * v; ret.m22 = m.m22 * v;
	#endif
	return ret;
}

vec2 &mul(vec2 &ret,const mat3 &m,const vec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12;
	return ret;
}

vec2 &mul(vec2 &ret,const vec2 &v,const mat3 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21;
	return ret;
}

vec3 &mul(vec3 &ret,const mat3 &m,const vec3 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,W));
		ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,W),zero);
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
		ret.vec = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmulq_lane_f32(m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		ret.vec = vmlaq_lane_f32(res_1,m.col2,high,0);
	#else
		ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
		ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
		ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	#endif
	return ret;
}

vec3 &mul(vec3 &ret,const vec3 &v,const mat3 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dvec2 &mul(dvec2 &ret,const mat3 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12;
	return ret;
}

dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat3 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21;
	return ret;
}

dvec3 &mul(dvec3 &ret,const mat3 &m,const dvec3 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	return ret;
}

dvec3 &mul(dvec3 &ret,const dvec3 &v,const mat3 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

mat3 &mul(mat3 &ret,const mat3 &m,const vec3 &v) {
	ret.m00 =  m.m01 * v.z - m.m02 * v.y;
	ret.m10 =  m.m11 * v.z - m.m12 * v.y;
	ret.m20 =  m.m21 * v.z - m.m22 * v.y;
	ret.m01 = -m.m00 * v.z + m.m02 * v.x;
	ret.m11 = -m.m10 * v.z + m.m12 * v.x;
	ret.m21 = -m.m20 * v.z + m.m22 * v.x;
	ret.m02 =  m.m00 * v.y - m.m01 * v.x;
	ret.m12 =  m.m10 * v.y - m.m11 * v.x;
	ret.m22 =  m.m20 * v.y - m.m21 * v.x;
	return ret;
}

mat3 &mul(mat3 &ret,const vec3 &v,const mat3 &m) {
	ret.m00 = -v.z * m.m10 + v.y * m.m20;
	ret.m10 =  v.z * m.m00 - v.x * m.m20;
	ret.m20 = -v.y * m.m00 + v.x * m.m10;
	ret.m01 = -v.z * m.m11 + v.y * m.m21;
	ret.m11 =  v.z * m.m01 - v.x * m.m21;
	ret.m21 = -v.y * m.m01 + v.x * m.m11;
	ret.m02 = -v.z * m.m12 + v.y * m.m22;
	ret.m12 =  v.z * m.m02 - v.x * m.m22;
	ret.m22 = -v.y * m.m02 + v.x * m.m12;
	return ret;
}

mat3 &mul(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col0,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col0,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col0,Z,Z,Z,W));
		ret.col0 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col1,X,X,X,W));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col1,Y,Y,Y,W));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col1,Z,Z,Z,W));
		ret.col1 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col2,X,X,X,W));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col2,Y,Y,Y,W));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col2,Z,Z,Z,W));
		ret.col2 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 xxxw = VEC_PERM2(X,X,X,W);
		vec_uchar16 yyyw = VEC_PERM2(Y,Y,Y,W);
		vec_uchar16 zzzw = VEC_PERM2(Z,Z,Z,W);
		ret.col0 = vec_madd(m0.col0,vec_perm(m1.col0,m1.col0,xxxw),zero);
		ret.col1 = vec_madd(m0.col0,vec_perm(m1.col1,m1.col1,xxxw),zero);
		ret.col2 = vec_madd(m0.col0,vec_perm(m1.col2,m1.col2,xxxw),zero);
		ret.col0 = vec_madd(m0.col1,vec_perm(m1.col0,m1.col0,yyyw),ret.col0);
		ret.col1 = vec_madd(m0.col1,vec_perm(m1.col1,m1.col1,yyyw),ret.col1);
		ret.col2 = vec_madd(m0.col1,vec_perm(m1.col2,m1.col2,yyyw),ret.col2);
		ret.col0 = vec_madd(m0.col2,vec_perm(m1.col0,m1.col0,zzzw),ret.col0);
		ret.col1 = vec_madd(m0.col2,vec_perm(m1.col1,m1.col1,zzzw),ret.col1);
		ret.col2 = vec_madd(m0.col2,vec_perm(m1.col2,m1.col2,zzzw),ret.col2);
	#elif USE_NEON
		float32x2_t low_0 = vget_low_f32(m1.col0);
		float32x2_t low_1 = vget_low_f32(m1.col1);
		float32x2_t low_2 = vget_low_f32(m1.col2);
		ret.col0 = vmulq_lane_f32(m0.col0,low_0,0);
		ret.col1 = vmulq_lane_f32(m0.col0,low_1,0);
		ret.col2 = vmulq_lane_f32(m0.col0,low_2,0);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col1,low_0,1);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col1,low_1,1);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col1,low_2,1);
		float32x2_t high_0 = vget_high_f32(m1.col0);
		float32x2_t high_1 = vget_high_f32(m1.col1);
		float32x2_t high_2 = vget_high_f32(m1.col2);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col2,high_0,0);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col2,high_1,0);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col2,high_2,0);
	#else
		ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20;
		ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20;
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21;
		ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21;
		ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22;
		ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
	#endif
	return ret;
}

mat3 &add(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	#ifdef USE_SSE
		ret.col0 = _mm_add_ps(m0.col0,m1.col0);
		ret.col1 = _mm_add_ps(m0.col1,m1.col1);
		ret.col2 = _mm_add_ps(m0.col2,m1.col2);
	#elif USE_ALTIVEC
		ret.col0 = vec_add(m0.col0,m1.col0);
		ret.col1 = vec_add(m0.col1,m1.col1);
		ret.col2 = vec_add(m0.col2,m1.col2);
	#elif USE_NEON
		ret.col0 = vaddq_f32(m0.col0,m1.col0);
		ret.col1 = vaddq_f32(m0.col1,m1.col1);
		ret.col2 = vaddq_f32(m0.col2,m1.col2);
	#else
		ret.m00 = m0.m00 + m1.m00; ret.m01 = m0.m01 + m1.m01; ret.m02 = m0.m02 + m1.m02;
		ret.m10 = m0.m10 + m1.m10; ret.m11 = m0.m11 + m1.m11; ret.m12 = m0.m12 + m1.m12;
		ret.m20 = m0.m20 + m1.m20; ret.m21 = m0.m21 + m1.m21; ret.m22 = m0.m22 + m1.m22;
	#endif
	return ret;
}

mat3 &sub(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	#ifdef USE_SSE
		ret.col0 = _mm_sub_ps(m0.col0,m1.col0);
		ret.col1 = _mm_sub_ps(m0.col1,m1.col1);
		ret.col2 = _mm_sub_ps(m0.col2,m1.col2);
	#elif USE_ALTIVEC
		ret.col0 = vec_sub(m0.col0,m1.col0);
		ret.col1 = vec_sub(m0.col1,m1.col1);
		ret.col2 = vec_sub(m0.col2,m1.col2);
	#elif USE_NEON
		ret.col0 = vsubq_f32(m0.col0,m1.col0);
		ret.col1 = vsubq_f32(m0.col1,m1.col1);
		ret.col2 = vsubq_f32(m0.col2,m1.col2);
	#else
		ret.m00 = m0.m00 - m1.m00; ret.m01 = m0.m01 - m1.m01; ret.m02 = m0.m02 - m1.m02;
		ret.m10 = m0.m10 - m1.m10; ret.m11 = m0.m11 - m1.m11; ret.m12 = m0.m12 - m1.m12;
		ret.m20 = m0.m20 - m1.m20; ret.m21 = m0.m21 - m1.m21; ret.m22 = m0.m22 - m1.m22;
	#endif
	return ret;
}

mat3 &orthonormalize(mat3 &ret,const mat3 &m) {
	#ifdef USE_SSE
		__m128 x_yzx = _MM_SWIZZLE(m.col0,Y,Z,X,W);
		__m128 x_zxy = _MM_SWIZZLE(m.col0,Z,X,Y,W);
		__m128 y_yzx = _MM_SWIZZLE(m.col1,Y,Z,X,W);
		__m128 y_zxy = _MM_SWIZZLE(m.col1,Z,X,Y,W);
		__m128 z = _mm_sub_ps(_mm_mul_ps(x_yzx,y_zxy),_mm_mul_ps(x_zxy,y_yzx));
		__m128 z_yzx = _MM_SWIZZLE(z,Y,Z,X,W);
		__m128 z_zxy = _MM_SWIZZLE(z,Z,X,Y,W);
		__m128 y = _mm_sub_ps(_mm_mul_ps(z_yzx,x_zxy),_mm_mul_ps(z_zxy,x_yzx));
		__m128 col_0 = _mm_mul_ps(m.col0,m.col0);
		__m128 col_1 = _mm_mul_ps(y,y);
		__m128 col_2 = _mm_mul_ps(z,z);
		__m128 res_0 = _mm_shuffle_ps(col_0,col_1,_MM_PERM2(X,Y,X,Y));
		__m128 res_1 = _mm_shuffle_ps(col_0,col_1,_MM_PERM2(Z,W,Z,W));
		__m128 res_2 = _mm_shuffle_ps(col_2,col_2,_MM_PERM2(X,Y,X,Y));
		__m128 res_3 = _mm_shuffle_ps(col_2,col_2,_MM_PERM2(Z,W,Z,W));
		__m128 row_0 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(X,Z,X,Z));
		__m128 row_1 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(Y,W,Y,W));
		__m128 row_2 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(X,Z,X,Z));
		__m128 ilength = _mm_rsqrt_ps_nr(_mm_add_ps(_mm_add_ps(row_0,row_1),row_2));
		ret.col0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(ilength,X,X,X,X));
		ret.col1 = _mm_mul_ps(y,_MM_SWIZZLE(ilength,Y,Y,Y,Y));
		ret.col2 = _mm_mul_ps(z,_MM_SWIZZLE(ilength,Z,Z,Z,Z));
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 yzxw = VEC_PERM2(Y,Z,X,W);
		vec_uchar16 zxyw = VEC_PERM2(Z,X,Y,W);
		vec_float4 x_yzx = vec_perm(m.col0,m.col0,yzxw);
		vec_float4 x_zxy = vec_perm(m.col0,m.col0,zxyw);
		vec_float4 y_yzx = vec_perm(m.col1,m.col1,yzxw);
		vec_float4 y_zxy = vec_perm(m.col1,m.col1,zxyw);
		vec_float4 z = vec_sub(vec_madd(x_yzx,y_zxy,zero),vec_madd(x_zxy,y_yzx,zero));
		vec_float4 z_yzx = vec_perm(z,z,yzxw);
		vec_float4 z_zxy = vec_perm(z,z,zxyw);
		vec_float4 y = vec_sub(vec_madd(z_yzx,x_zxy,zero),vec_madd(z_zxy,x_yzx,zero));
		vec_float4 col_0 = vec_madd(m.col0,m.col0,zero);
		vec_float4 col_1 = vec_madd(y,y,zero);
		vec_float4 col_2 = vec_madd(z,z,zero);
		vec_float4 res_0 = vec_perm(col_0,col_1,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_1 = vec_perm(col_0,col_1,VEC_PERM2(Z,W,Z,W));
		vec_float4 res_2 = vec_perm(col_2,col_2,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_3 = vec_perm(col_2,col_2,VEC_PERM2(Z,W,Z,W));
		vec_float4 row_0 = vec_perm(res_0,res_2,VEC_PERM2(X,Z,X,Z));
		vec_float4 row_1 = vec_perm(res_0,res_2,VEC_PERM2(Y,W,Y,W));
		vec_float4 row_2 = vec_perm(res_1,res_3,VEC_PERM2(X,Z,X,Z));
		vec_float4 ilength = vec_rsqrt_nr(vec_add(vec_add(row_0,row_1),row_2));
		ret.col0 = vec_madd(m.col0,VEC_SWIZZLE(ilength,X,X,X,X),zero);
		ret.col1 = vec_madd(y,VEC_SWIZZLE(ilength,Y,Y,Y,Y),zero);
		ret.col2 = vec_madd(z,VEC_SWIZZLE(ilength,Z,Z,Z,Z),zero);
	#elif USE_NEON
		float32x4_t z = vcrossq_f32(m.col0,m.col1);
		float32x4_t y = vcrossq_f32(z,m.col0);
		ret.col0 = vnormalize3q_f32(m.col0);
		ret.col1 = vnormalize3q_f32(y);
		ret.col2 = vnormalize3q_f32(z);
	#else
		vec3 x = vec3(m.m00,m.m10,m.m20);
		vec3 y = vec3(m.m01,m.m11,m.m21);
		vec3 z = cross(x,y);
		cross(y,z,x);
		x.normalize();
		y.normalize();
		z.normalize();
		ret.m00 = x.x; ret.m01 = y.x; ret.m02 = z.x;
		ret.m10 = x.y; ret.m11 = y.y; ret.m12 = z.y;
		ret.m20 = x.z; ret.m21 = y.z; ret.m22 = z.z;
	#endif
	return ret;
}

mat3 &transpose(mat3 &ret,const mat3 &m) {
	ret.m00 = m.m00; ret.m01 = m.m10; ret.m02 = m.m20;
	ret.m10 = m.m01; ret.m11 = m.m11; ret.m12 = m.m21;
	ret.m20 = m.m02; ret.m21 = m.m12; ret.m22 = m.m22;
	return ret;
}

mat3 &inverse(mat3 &ret,const mat3 &m) {
	float idet = Math::rcp(determinant(m));
	ret.m00 =  (m.m11 * m.m22 - m.m12 * m.m21) * idet;
	ret.m10 = -(m.m10 * m.m22 - m.m12 * m.m20) * idet;
	ret.m20 =  (m.m10 * m.m21 - m.m11 * m.m20) * idet;
	ret.m01 = -(m.m01 * m.m22 - m.m02 * m.m21) * idet;
	ret.m11 =  (m.m00 * m.m22 - m.m02 * m.m20) * idet;
	ret.m21 = -(m.m00 * m.m21 - m.m01 * m.m20) * idet;
	ret.m02 =  (m.m01 * m.m12 - m.m02 * m.m11) * idet;
	ret.m12 = -(m.m00 * m.m12 - m.m02 * m.m10) * idet;
	ret.m22 =  (m.m00 * m.m11 - m.m01 * m.m10) * idet;
	return ret;
}

mat3 &inverse(mat3 &ret,const mat3 &m,float det) {
	float idet = Math::rcp(det);
	ret.m00 =  (m.m11 * m.m22 - m.m12 * m.m21) * idet;
	ret.m10 = -(m.m10 * m.m22 - m.m12 * m.m20) * idet;
	ret.m20 =  (m.m10 * m.m21 - m.m11 * m.m20) * idet;
	ret.m01 = -(m.m01 * m.m22 - m.m02 * m.m21) * idet;
	ret.m11 =  (m.m00 * m.m22 - m.m02 * m.m20) * idet;
	ret.m21 = -(m.m00 * m.m21 - m.m01 * m.m20) * idet;
	ret.m02 =  (m.m01 * m.m12 - m.m02 * m.m11) * idet;
	ret.m12 = -(m.m00 * m.m12 - m.m02 * m.m10) * idet;
	ret.m22 =  (m.m00 * m.m11 - m.m01 * m.m10) * idet;
	return ret;
}

/*
 */
mat3 orthonormalize(const mat3 &m) {
	mat3 ret;
	return orthonormalize(ret,m);
}

mat3 transpose(const mat3 &m) {
	mat3 ret;
	return transpose(ret,m);
}

mat3 inverse(const mat3 &m) {
	mat3 ret;
	return inverse(ret,m);
}

mat3 inverse(const mat3 &m,float det) {
	mat3 ret;
	return inverse(ret,m,det);
}

/*
 */
mat3 rotate3(const vec3 &axis,float angle) {
	mat3 ret;
	ret.setRotate(axis,angle);
	return ret;
}

mat3 rotate3(float x,float y,float z,float angle) {
	return rotate3(vec3(x,y,z),angle);
}

mat3 rotate3(const quat &q) {
	return q.getMat3();
}

mat3 rotateX3(float angle) {
	mat3 ret;
	ret.setRotateX(angle);
	return ret;
}

mat3 rotateY3(float angle) {
	mat3 ret;
	ret.setRotateY(angle);
	return ret;
}

mat3 rotateZ3(float angle) {
	mat3 ret;
	ret.setRotateZ(angle);
	return ret;
}

mat3 scale3(const vec3 &v) {
	mat3 ret;
	ret.setScale(v);
	return ret;
}

mat3 scale3(float x,float y,float z) {
	return scale3(vec3(x,y,z));
}

/*
 */
mat3 jacobi(const mat3 &m,mat3 &v) {
	mat3 j,ret = m;
	v.setIdentity();
	float prev = INFINITY;
	const int num_iterations = 50;
	for(int n = 0; n < num_iterations; n++) {
		int p = 0;
		int q = 1;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				if(i == j) continue;
				if(Math::abs(ret.get(i,j)) > Math::abs(ret.get(p,q))) {
					p = i;
					q = j;
				}
			}
		}
		float c = 1.0f;
		float s = 0.0f;
		if(Math::abs(ret.get(p,q)) > EPSILON) {
			float r = (ret.get(q,q) - ret.get(p,p)) / (ret.get(p,q) * 2.0f);
			if(r >= 0.0f) r = Math::rcp(r + Math::sqrt(1.0f + r * r));
			else r = -Math::rcp(-r + Math::sqrt(1.0f + r * r));
			c = Math::rsqrt(1.0f + r * r);
			s = r * c;
		}
		j.setIdentity();
		j.set(p,p,c);
		j.set(p,q,s);
		j.set(q,p,-s);
		j.set(q,q,c);
		v = v * j;
		ret = transpose(j) * ret * j;
		float sum = 0.0f;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				if(i == j) continue;
				sum += ret.get(i,j) * ret.get(i,j);
			}
		}
		if(prev <= sum) break;
		prev = sum;
	}
	return ret;
}

/******************************************************************************\
*
* mat4
*
\******************************************************************************/

/*
 */
static const float mat4_identity_values[16] = {
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 1.0f,
};

/*
 */
const mat4 mat4_zero(0.0f);
const mat4 mat4_one(1.0f);
const mat4 mat4_identity(mat4_identity_values);

/*
 */
mat4::mat4(float v) {
	m00 = v; m01 = v; m02 = v; m03 = v;
	m10 = v; m11 = v; m12 = v; m13 = v;
	m20 = v; m21 = v; m22 = v; m23 = v;
	m30 = v; m31 = v; m32 = v; m33 = v;
}

mat4::mat4(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0f; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f;  m21 = 0.0f;  m22 = 1.0f; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f; m33 = 1.0f;
}

mat4::mat4(const mat3 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0f;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

mat4::mat4(const mat4 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
		col3 = m.col3;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
		m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
		m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
		m30 = m.m30; m31 = m.m31; m32 = m.m32; m33 = m.m33;
	#endif
}

mat4::mat4(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01; m02 = (float)m.m02; m03 = (float)m.m03;
	m10 = (float)m.m10; m11 = (float)m.m11; m12 = (float)m.m12; m13 = (float)m.m13;
	m20 = (float)m.m20; m21 = (float)m.m21; m22 = (float)m.m22; m23 = (float)m.m23;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

mat4::mat4(const quat &q) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0f;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

mat4::mat4(const float *m) {
	m00 = m[0]; m01 = m[4]; m02 = m[8];  m03 = m[12];
	m10 = m[1]; m11 = m[5]; m12 = m[9];  m13 = m[13];
	m20 = m[2]; m21 = m[6]; m22 = m[10]; m23 = m[14];
	m30 = m[3]; m31 = m[7]; m32 = m[11]; m33 = m[15];
}

mat4::mat4(const mat3 &m,const vec3 &v) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

mat4::mat4(const quat &q,const vec3 &v) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

/*
 */
mat4 mat4::operator-() const {
	mat4 ret;
	ret.m00 = -m00; ret.m01 = -m01; ret.m02 = -m02; ret.m03 = -m03;
	ret.m10 = -m10; ret.m11 = -m11; ret.m12 = -m12; ret.m13 = -m13;
	ret.m20 = -m20; ret.m21 = -m21; ret.m22 = -m22; ret.m23 = -m23;
	ret.m30 = -m30; ret.m31 = -m31; ret.m32 = -m32; ret.m33 = -m33;
	return ret;
}

mat4 &mat4::operator*=(float v) {
	return mul(*this,*this,v);
}

mat4 &mat4::operator*=(const mat4 &m) {
	return mul(*this,mat4(*this),m);
}

mat4 &mat4::operator+=(const mat4 &m) {
	return add(*this,*this,m);
}

mat4 &mat4::operator-=(const mat4 &m) {
	return sub(*this,*this,m);
}

/*
 */
void mat4::set(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0f; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f;  m21 = 0.0f;  m22 = 1.0f; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f; m33 = 1.0f;
}

void mat4::set(const mat3 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0f;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

void mat4::set(const mat4 &m) {
	#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
		col0 = m.col0;
		col1 = m.col1;
		col2 = m.col2;
		col3 = m.col3;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
		m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
		m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
		m30 = m.m30; m31 = m.m31; m32 = m.m32; m33 = m.m33;
	#endif
}

void mat4::set(const dmat4 &m) {
	m00 = (float)m.m00; m01 = (float)m.m01; m02 = (float)m.m02; m03 = (float)m.m03;
	m10 = (float)m.m10; m11 = (float)m.m11; m12 = (float)m.m12; m13 = (float)m.m13;
	m20 = (float)m.m20; m21 = (float)m.m21; m22 = (float)m.m22; m23 = (float)m.m23;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

void mat4::set(const quat &q) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0f;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

void mat4::set(const float *m) {
	m00 = m[0]; m01 = m[4]; m02 = m[8];  m03 = m[12];
	m10 = m[1]; m11 = m[5]; m12 = m[9];  m13 = m[13];
	m20 = m[2]; m21 = m[6]; m22 = m[10]; m23 = m[14];
	m30 = m[3]; m31 = m[7]; m32 = m[11]; m33 = m[15];
}

void mat4::set(const mat3 &m,const vec3 &v) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

void mat4::set(const quat &q,const vec3 &v) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

void mat4::get(float *m) const {
	m[0] = m00; m[4] = m01; m[8] = m02;  m[12] = m03;
	m[1] = m10; m[5] = m11; m[9] = m12;  m[13] = m13;
	m[2] = m20; m[6] = m21; m[10] = m22; m[14] = m23;
	m[3] = m30; m[7] = m31; m[11] = m32; m[15] = m33;
}

/*
 */
void mat4::setRow(int row,const vec4 &v) {
	assert((unsigned int)row < 4 && "mat4::setRow(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
	mat[row + 12] = v.w;
}

void mat4::setRow3(int row,const vec3 &v) {
	assert((unsigned int)row < 4 && "mat4::setRow3(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
}

vec4 mat4::getRow(int row) const {
	assert((unsigned int)row < 4 && "mat4::getRow(): bad row");
	return vec4(mat[row + 0],mat[row + 4],mat[row + 8],mat[row + 12]);
}

vec3 mat4::getRow3(int row) const {
	assert((unsigned int)row < 4 && "mat4::getRow3(): bad row");
	return vec3(mat[row + 0],mat[row + 4],mat[row + 8]);
}

/*
 */
void mat4::setColumn(int column,const vec4 &v) {
	assert((unsigned int)column < 4 && "mat4::setColumn(): bad column");
	column *= 4;
	mat[column + 0] = v.x;
	mat[column + 1] = v.y;
	mat[column + 2] = v.z;
	mat[column + 3] = v.w;
}

void mat4::setColumn3(int column,const vec3 &v) {
	assert((unsigned int)column < 4 && "mat4::setColumn3(): bad column");
	column *= 4;
	mat[column + 0] = v.x;
	mat[column + 1] = v.y;
	mat[column + 2] = v.z;
}

vec4 mat4::getColumn(int column) const {
	assert((unsigned int)column < 4 && "mat4::getColumn(): bad column");
	return vec4(mat + column * 4);
}

vec3 mat4::getColumn3(int column) const {
	assert((unsigned int)column < 4 && "mat4::getColumn3(): bad column");
	return vec3(mat + column * 4);
}

/*
 */
void mat4::setDiagonal(const vec4 &v) {
	m00 = v.x;  m01 = 0.0f; m02 = 0.0f; m03 = 0.0f;
	m10 = 0.0f; m11 = v.y;  m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = v.z;  m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = v.w;
}

vec4 mat4::getDiagonal() const {
	return vec4(m00,m11,m22,m33);
}

/*
 */
void mat4::setZero() {
	m00 = 0.0f; m01 = 0.0f; m02 = 0.0f; m03 = 0.0f;
	m10 = 0.0f; m11 = 0.0f; m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 0.0f; m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 0.0f;
}

void mat4::setIdentity() {
	m00 = 1.0f; m01 = 0.0f; m02 = 0.0f; m03 = 0.0f;
	m10 = 0.0f; m11 = 1.0f; m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 1.0f; m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

void mat4::setTranslate(const vec3 &v) {
	m00 = 1.0f; m01 = 0.0f; m02 = 0.0f; m03 = v.x;
	m10 = 0.0f; m11 = 1.0f; m12 = 0.0f; m13 = v.y;
	m20 = 0.0f; m21 = 0.0f; m22 = 1.0f; m23 = v.z;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

void mat4::setRotate(const vec3 &axis,float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	vec3 v = normalize(axis);
	float xx = v.x * v.x;
	float yy = v.y * v.y;
	float zz = v.z * v.z;
	float xy = v.x * v.y;
	float yz = v.y * v.z;
	float zx = v.z * v.x;
	float xs = v.x * s;
	float ys = v.y * s;
	float zs = v.z * s;
	m00 = (1.0f - c) * xx + c;  m01 = (1.0f - c) * xy - zs; m02 = (1.0f - c) * zx + ys; m03 = 0.0f;
	m10 = (1.0f - c) * xy + zs; m11 = (1.0f - c) * yy + c;  m12 = (1.0f - c) * yz - xs; m13 = 0.0f;
	m20 = (1.0f - c) * zx - ys; m21 = (1.0f - c) * yz + xs; m22 = (1.0f - c) * zz + c;  m23 = 0.0f;
	m30 = 0.0f;                 m31 = 0.0f;                 m32 = 0.0f;                 m33 = 1.0f;
}

void mat4::setRotateX(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = 1.0f; m01 = 0.0f; m02 = 0.0f; m03 = 0.0f;
	m10 = 0.0f; m11 = c;    m12 = -s;   m13 = 0.0f;
	m20 = 0.0f; m21 = s;    m22 = c;    m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

void mat4::setRotateY(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;    m01 = 0.0f; m02 = s;    m03 = 0.0f;
	m10 = 0.0f; m11 = 1.0f; m12 = 0.0f; m13 = 0.0f;
	m20 = -s;   m21 = 0.0f; m22 = c;    m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

void mat4::setRotateZ(float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;    m01 = -s;   m02 = 0.0f; m03 = 0.0f;
	m10 = s;    m11 = c;    m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = 1.0f; m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

void mat4::setScale(const vec3 &v) {
	m00 = v.x;  m01 = 0.0f; m02 = 0.0f; m03 = 0.0f;
	m10 = 0.0f; m11 = v.y;  m12 = 0.0f; m13 = 0.0f;
	m20 = 0.0f; m21 = 0.0f; m22 = v.z;  m23 = 0.0f;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f; m33 = 1.0f;
}

/*
 */
int operator==(const mat4 &m0,const mat4 &m1) {
	return compare(m0,m1);
}

int operator!=(const mat4 &m0,const mat4 &m1) {
	return !compare(m0,m1);
}

mat4 operator*(const mat4 &m,float v) {
	mat4 ret;
	return mul(ret,m,v);
}

vec2 operator*(const mat4 &m,const vec2 &v) {
	vec2 ret;
	return mul(ret,m,v);
}

vec2 operator*(const vec2 &v,const mat4 &m) {
	vec2 ret;
	return mul(ret,v,m);
}

vec3 operator*(const mat4 &m,const vec3 &v) {
	vec3 ret;
	return mul(ret,m,v);
}

vec3 operator*(const vec3 &v,const mat4 &m) {
	vec3 ret;
	return mul(ret,v,m);
}

vec4 operator*(const mat4 &m,const vec4 &v) {
	vec4 ret;
	return mul(ret,m,v);
}

vec4 operator*(const vec4 &v,const mat4 &m) {
	vec4 ret;
	return mul(ret,v,m);
}

dvec2 operator*(const mat4 &m,const dvec2 &v) {
	dvec2 ret;
	return mul(ret,m,v);
}

dvec2 operator*(const dvec2 &v,const mat4 &m) {
	dvec2 ret;
	return mul(ret,v,m);
}

dvec3 operator*(const mat4 &m,const dvec3 &v) {
	dvec3 ret;
	return mul(ret,m,v);
}

dvec3 operator*(const dvec3 &v,const mat4 &m) {
	dvec3 ret;
	return mul(ret,v,m);
}

dvec4 operator*(const mat4 &m,const dvec4 &v) {
	dvec4 ret;
	return mul(ret,m,v);
}

dvec4 operator*(const dvec4 &v,const mat4 &m) {
	dvec4 ret;
	return mul(ret,v,m);
}

mat4 operator*(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return mul(ret,m0,m1);
}

mat4 operator+(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return add(ret,m0,m1);
}

mat4 operator-(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return sub(ret,m0,m1);
}

/*
 */
int compare(const mat4 &m0,const mat4 &m1) {
	return (compare(m0.m00,m1.m00) && compare(m0.m10,m1.m10) && compare(m0.m20,m1.m20) && compare(m0.m30,m1.m30) &&
		compare(m0.m01,m1.m01) && compare(m0.m11,m1.m11) && compare(m0.m21,m1.m21) && compare(m0.m31,m1.m31) &&
		compare(m0.m02,m1.m02) && compare(m0.m12,m1.m12) && compare(m0.m22,m1.m22) && compare(m0.m32,m1.m32) &&
		compare(m0.m03,m1.m03) && compare(m0.m13,m1.m13) && compare(m0.m23,m1.m23) && compare(m0.m33,m1.m33));
}

int compare(const mat4 &m0,const mat4 &m1,float epsilon) {
	return (compare(m0.m00,m1.m00,epsilon) && compare(m0.m10,m1.m10,epsilon) && compare(m0.m20,m1.m20,epsilon) && compare(m0.m30,m1.m30,epsilon) &&
		compare(m0.m01,m1.m01,epsilon) && compare(m0.m11,m1.m11,epsilon) && compare(m0.m21,m1.m21,epsilon) && compare(m0.m31,m1.m31,epsilon) &&
		compare(m0.m02,m1.m02,epsilon) && compare(m0.m12,m1.m12,epsilon) && compare(m0.m22,m1.m22,epsilon) && compare(m0.m32,m1.m32,epsilon) &&
		compare(m0.m03,m1.m03,epsilon) && compare(m0.m13,m1.m13,epsilon) && compare(m0.m23,m1.m23,epsilon) && compare(m0.m33,m1.m33,epsilon));
}

float trace(const mat4 &m) {
	return m.m00 + m.m11 + m.m22 + m.m33;
}

float determinant(const mat4 &m) {
	float det = 0.0f;
	det =  m.m00 * (m.m11 * (m.m22 * m.m33 - m.m23 * m.m32) - m.m12 * (m.m21 * m.m33 - m.m23 * m.m31) + m.m13 * (m.m21 * m.m32 - m.m22 * m.m31));
	det -= m.m01 * (m.m10 * (m.m22 * m.m33 - m.m23 * m.m32) - m.m12 * (m.m20 * m.m33 - m.m23 * m.m30) + m.m13 * (m.m20 * m.m32 - m.m22 * m.m30));
	det += m.m02 * (m.m10 * (m.m21 * m.m33 - m.m23 * m.m31) - m.m11 * (m.m20 * m.m33 - m.m23 * m.m30) + m.m13 * (m.m20 * m.m31 - m.m21 * m.m30));
	det -= m.m03 * (m.m10 * (m.m21 * m.m32 - m.m22 * m.m31) - m.m11 * (m.m20 * m.m32 - m.m22 * m.m30) + m.m12 * (m.m20 * m.m31 - m.m21 * m.m30));
	return det;
}

float determinant3(const mat4 &m) {
	float det = 0.0f;
	det =  m.m00 * (m.m11 * m.m22 - m.m12 * m.m21);
	det -= m.m01 * (m.m10 * m.m22 - m.m12 * m.m20);
	det += m.m02 * (m.m10 * m.m21 - m.m11 * m.m20);
	return det;
}

mat4 &mul(mat4 &ret,const mat4 &m,float v) {
	#ifdef USE_SSE
		__m128 temp = _mm_set1_ps(v);
		ret.col0 = _mm_mul_ps(m.col0,temp);
		ret.col1 = _mm_mul_ps(m.col1,temp);
		ret.col2 = _mm_mul_ps(m.col2,temp);
		ret.col3 = _mm_mul_ps(m.col3,temp);
	#elif USE_ALTIVEC
		vec_float4 temp = vec_splats(v);
		vec_float4 zero = vec_splats(0.0f);
		ret.col0 = vec_madd(m.col0,temp,zero);
		ret.col1 = vec_madd(m.col1,temp,zero);
		ret.col2 = vec_madd(m.col2,temp,zero);
		ret.col3 = vec_madd(m.col3,temp,zero);
	#elif USE_NEON
		ret.col0 = vmulq_n_f32(m.col0,v);
		ret.col1 = vmulq_n_f32(m.col1,v);
		ret.col2 = vmulq_n_f32(m.col2,v);
		ret.col3 = vmulq_n_f32(m.col3,v);
	#else
		ret.m00 = m.m00 * v; ret.m01 = m.m01 * v; ret.m02 = m.m02 * v; ret.m03 = m.m03 * v;
		ret.m10 = m.m10 * v; ret.m11 = m.m11 * v; ret.m12 = m.m12 * v; ret.m13 = m.m13 * v;
		ret.m20 = m.m20 * v; ret.m21 = m.m21 * v; ret.m22 = m.m22 * v; ret.m23 = m.m23 * v;
		ret.m30 = m.m30 * v; ret.m31 = m.m31 * v; ret.m32 = m.m32 * v; ret.m33 = m.m33 * v;
	#endif
	return ret;
}

vec2 &mul(vec2 &ret,const mat4 &m,const vec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m03;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m13;
	return ret;
}

vec2 &mul(vec2 &ret,const vec2 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m30;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m31;
	return ret;
}

vec3 &mul(vec3 &ret,const mat4 &m,const vec3 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,W));
		ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,m.col3));
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,W),m.col3);
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
		ret.vec = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmlaq_lane_f32(m.col3,m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		ret.vec = vmlaq_lane_f32(res_1,m.col2,high,0);
	#else
		ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
		ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
		ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
	#endif
	return ret;
}

vec3 &mul(vec3 &ret,const vec3 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z + m.m30;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z + m.m31;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z + m.m32;
	return ret;
}

vec4 &mul(vec4 &ret,const mat4 &m,const vec4 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(m.col3,_MM_SWIZZLE(v.vec,W,W,W,W));
		ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,X),vec_splats(0.0f));
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,Y),res_0);
		vec_float4 res_2 = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,Z),res_1);
		ret.vec = vec_madd(m.col3,VEC_SWIZZLE(v.vec,W,W,W,W),res_2);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmulq_lane_f32(m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,m.col2,high,0);
		ret.vec = vmlaq_lane_f32(res_2,m.col3,high,1);
	#else
		ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w;
		ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w;
		ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w;
		ret.w = m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33 * v.w;
	#endif
	return ret;
}

vec4 &mul(vec4 &ret,const vec4 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z + m.m30 * v.w;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z + m.m31 * v.w;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z + m.m32 * v.w;
	ret.w = m.m03 * v.x + m.m13 * v.y + m.m23 * v.z + m.m33 * v.w;
	return ret;
}

dvec2 &mul(dvec2 &ret,const mat4 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m03;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m13;
	return ret;
}

dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m30;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m31;
	return ret;
}

dvec3 &mul(dvec3 &ret,const mat4 &m,const dvec3 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
	return ret;
}

dvec3 &mul(dvec3 &ret,const dvec3 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z + m.m30;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z + m.m31;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z + m.m32;
	return ret;
}

dvec4 &mul(dvec4 &ret,const mat4 &m,const dvec4 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w;
	ret.w = m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33 * v.w;
	return ret;
}

dvec4 &mul(dvec4 &ret,const dvec4 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z + m.m30 * v.w;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z + m.m31 * v.w;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z + m.m32 * v.w;
	ret.w = m.m03 * v.x + m.m13 * v.y + m.m23 * v.z + m.m33 * v.w;
	return ret;
}

vec2 &mul3(vec2 &ret,const mat4 &m,const vec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y;
	ret.y = m.m10 * v.x + m.m11 * v.y;
	return ret;
}

vec2 &mul3(vec2 &ret,const vec2 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

vec3 &mul3(vec3 &ret,const mat4 &m,const vec3 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,W));
		ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,W),vec_splats(0.0f));
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
		ret.vec = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmulq_lane_f32(m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		ret.vec = vmlaq_lane_f32(res_1,m.col2,high,0);
	#else
		ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
		ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
		ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	#endif
	return ret;
}

vec3 &mul3(vec3 &ret,const vec3 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

vec4 &mul3(vec4 &ret,const mat4 &m,const vec4 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,W));
		ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,W),vec_splats(0.0f));
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
		ret.vec = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmulq_lane_f32(m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		ret.vec = vmlaq_lane_f32(res_1,m.col2,high,0);
	#else
		ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
		ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
		ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	#endif
	return ret;
}

vec4 &mul3(vec4 &ret,const vec4 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dvec2 &mul3(dvec2 &ret,const mat4 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y;
	ret.y = m.m10 * v.x + m.m11 * v.y;
	return ret;
}

dvec2 &mul3(dvec2 &ret,const dvec2 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

dvec3 &mul3(dvec3 &ret,const mat4 &m,const dvec3 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	return ret;
}

dvec3 &mul3(dvec3 &ret,const dvec3 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dvec4 &mul3(dvec4 &ret,const mat4 &m,const dvec4 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	return ret;
}

dvec4 &mul3(dvec4 &ret,const dvec4 &v,const mat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

mat4 &mul(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col0,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col0,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col0,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(m0.col3,_MM_SWIZZLE(m1.col0,W,W,W,W));
		ret.col0 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col1,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col1,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col1,Z,Z,Z,Z));
		res_3 = _mm_mul_ps(m0.col3,_MM_SWIZZLE(m1.col1,W,W,W,W));
		ret.col1 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col2,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col2,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col2,Z,Z,Z,Z));
		res_3 = _mm_mul_ps(m0.col3,_MM_SWIZZLE(m1.col2,W,W,W,W));
		ret.col2 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col3,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col3,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col3,Z,Z,Z,Z));
		res_3 = _mm_mul_ps(m0.col3,_MM_SWIZZLE(m1.col3,W,W,W,W));
		ret.col3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
		vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
		vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
		vec_uchar16 wwww = VEC_PERM2(W,W,W,W);
		ret.col0 = vec_madd(m0.col0,vec_perm(m1.col0,m1.col0,xxxx),zero);
		ret.col1 = vec_madd(m0.col0,vec_perm(m1.col1,m1.col1,xxxx),zero);
		ret.col2 = vec_madd(m0.col0,vec_perm(m1.col2,m1.col2,xxxx),zero);
		ret.col3 = vec_madd(m0.col0,vec_perm(m1.col3,m1.col3,xxxx),zero);
		ret.col0 = vec_madd(m0.col1,vec_perm(m1.col0,m1.col0,yyyy),ret.col0);
		ret.col1 = vec_madd(m0.col1,vec_perm(m1.col1,m1.col1,yyyy),ret.col1);
		ret.col2 = vec_madd(m0.col1,vec_perm(m1.col2,m1.col2,yyyy),ret.col2);
		ret.col3 = vec_madd(m0.col1,vec_perm(m1.col3,m1.col3,yyyy),ret.col3);
		ret.col0 = vec_madd(m0.col2,vec_perm(m1.col0,m1.col0,zzzz),ret.col0);
		ret.col1 = vec_madd(m0.col2,vec_perm(m1.col1,m1.col1,zzzz),ret.col1);
		ret.col2 = vec_madd(m0.col2,vec_perm(m1.col2,m1.col2,zzzz),ret.col2);
		ret.col3 = vec_madd(m0.col2,vec_perm(m1.col3,m1.col3,zzzz),ret.col3);
		ret.col0 = vec_madd(m0.col3,vec_perm(m1.col0,m1.col0,wwww),ret.col0);
		ret.col1 = vec_madd(m0.col3,vec_perm(m1.col1,m1.col1,wwww),ret.col1);
		ret.col2 = vec_madd(m0.col3,vec_perm(m1.col2,m1.col2,wwww),ret.col2);
		ret.col3 = vec_madd(m0.col3,vec_perm(m1.col3,m1.col3,wwww),ret.col3);
	#elif USE_NEON
		float32x2_t low_0 = vget_low_f32(m1.col0);
		float32x2_t low_1 = vget_low_f32(m1.col1);
		float32x2_t low_2 = vget_low_f32(m1.col2);
		float32x2_t low_3 = vget_low_f32(m1.col3);
		ret.col0 = vmulq_lane_f32(m0.col0,low_0,0);
		ret.col1 = vmulq_lane_f32(m0.col0,low_1,0);
		ret.col2 = vmulq_lane_f32(m0.col0,low_2,0);
		ret.col3 = vmulq_lane_f32(m0.col0,low_3,0);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col1,low_0,1);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col1,low_1,1);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col1,low_2,1);
		ret.col3 = vmlaq_lane_f32(ret.col3,m0.col1,low_3,1);
		float32x2_t high_0 = vget_high_f32(m1.col0);
		float32x2_t high_1 = vget_high_f32(m1.col1);
		float32x2_t high_2 = vget_high_f32(m1.col2);
		float32x2_t high_3 = vget_high_f32(m1.col3);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col2,high_0,0);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col2,high_1,0);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col2,high_2,0);
		ret.col3 = vmlaq_lane_f32(ret.col3,m0.col2,high_3,0);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col3,high_0,1);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col3,high_1,1);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col3,high_2,1);
		ret.col3 = vmlaq_lane_f32(ret.col3,m0.col3,high_3,1);
	#else
		ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20 + m0.m03 * m1.m30;
		ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20 + m0.m13 * m1.m30;
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20 + m0.m23 * m1.m30;
		ret.m30 = m0.m30 * m1.m00 + m0.m31 * m1.m10 + m0.m32 * m1.m20 + m0.m33 * m1.m30;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21 + m0.m03 * m1.m31;
		ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21 + m0.m13 * m1.m31;
		ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21 + m0.m23 * m1.m31;
		ret.m31 = m0.m30 * m1.m01 + m0.m31 * m1.m11 + m0.m32 * m1.m21 + m0.m33 * m1.m31;
		ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22 + m0.m03 * m1.m32;
		ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22 + m0.m13 * m1.m32;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22 + m0.m23 * m1.m32;
		ret.m32 = m0.m30 * m1.m02 + m0.m31 * m1.m12 + m0.m32 * m1.m22 + m0.m33 * m1.m32;
		ret.m03 = m0.m00 * m1.m03 + m0.m01 * m1.m13 + m0.m02 * m1.m23 + m0.m03 * m1.m33;
		ret.m13 = m0.m10 * m1.m03 + m0.m11 * m1.m13 + m0.m12 * m1.m23 + m0.m13 * m1.m33;
		ret.m23 = m0.m20 * m1.m03 + m0.m21 * m1.m13 + m0.m22 * m1.m23 + m0.m23 * m1.m33;
		ret.m33 = m0.m30 * m1.m03 + m0.m31 * m1.m13 + m0.m32 * m1.m23 + m0.m33 * m1.m33;
	#endif
	return ret;
}

mat4 &mul4(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col0,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col0,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col0,Z,Z,Z,Z));
		ret.col0 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col1,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col1,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col1,Z,Z,Z,Z));
		ret.col1 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col2,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col2,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col2,Z,Z,Z,Z));
		ret.col2 = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		res_0 = _mm_mul_ps(m0.col0,_MM_SWIZZLE(m1.col3,X,X,X,X));
		res_1 = _mm_mul_ps(m0.col1,_MM_SWIZZLE(m1.col3,Y,Y,Y,Y));
		res_2 = _mm_mul_ps(m0.col2,_MM_SWIZZLE(m1.col3,Z,Z,Z,Z));
		ret.col3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,m0.col3));
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
		vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
		vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
		ret.col0 = vec_madd(m0.col0,vec_perm(m1.col0,m1.col0,xxxx),zero);
		ret.col1 = vec_madd(m0.col0,vec_perm(m1.col1,m1.col1,xxxx),zero);
		ret.col2 = vec_madd(m0.col0,vec_perm(m1.col2,m1.col2,xxxx),zero);
		ret.col3 = vec_madd(m0.col0,vec_perm(m1.col3,m1.col3,xxxx),m0.col3);
		ret.col0 = vec_madd(m0.col1,vec_perm(m1.col0,m1.col0,yyyy),ret.col0);
		ret.col1 = vec_madd(m0.col1,vec_perm(m1.col1,m1.col1,yyyy),ret.col1);
		ret.col2 = vec_madd(m0.col1,vec_perm(m1.col2,m1.col2,yyyy),ret.col2);
		ret.col3 = vec_madd(m0.col1,vec_perm(m1.col3,m1.col3,yyyy),ret.col3);
		ret.col0 = vec_madd(m0.col2,vec_perm(m1.col0,m1.col0,zzzz),ret.col0);
		ret.col1 = vec_madd(m0.col2,vec_perm(m1.col1,m1.col1,zzzz),ret.col1);
		ret.col2 = vec_madd(m0.col2,vec_perm(m1.col2,m1.col2,zzzz),ret.col2);
		ret.col3 = vec_madd(m0.col2,vec_perm(m1.col3,m1.col3,zzzz),ret.col3);
	#elif USE_NEON
		float32x2_t low_0 = vget_low_f32(m1.col0);
		float32x2_t low_1 = vget_low_f32(m1.col1);
		float32x2_t low_2 = vget_low_f32(m1.col2);
		float32x2_t low_3 = vget_low_f32(m1.col3);
		ret.col0 = vmulq_lane_f32(m0.col0,low_0,0);
		ret.col1 = vmulq_lane_f32(m0.col0,low_1,0);
		ret.col2 = vmulq_lane_f32(m0.col0,low_2,0);
		ret.col3 = vmlaq_lane_f32(m0.col3,m0.col0,low_3,0);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col1,low_0,1);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col1,low_1,1);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col1,low_2,1);
		ret.col3 = vmlaq_lane_f32(ret.col3,m0.col1,low_3,1);
		float32x2_t high_0 = vget_high_f32(m1.col0);
		float32x2_t high_1 = vget_high_f32(m1.col1);
		float32x2_t high_2 = vget_high_f32(m1.col2);
		float32x2_t high_3 = vget_high_f32(m1.col3);
		ret.col0 = vmlaq_lane_f32(ret.col0,m0.col2,high_0,0);
		ret.col1 = vmlaq_lane_f32(ret.col1,m0.col2,high_1,0);
		ret.col2 = vmlaq_lane_f32(ret.col2,m0.col2,high_2,0);
		ret.col3 = vmlaq_lane_f32(ret.col3,m0.col2,high_3,0);
	#else
		ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20;
		ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20;
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m30 = 0.0f;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21;
		ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21;
		ret.m31 = 0.0f;
		ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22;
		ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
		ret.m32 = 0.0f;
		ret.m03 = m0.m00 * m1.m03 + m0.m01 * m1.m13 + m0.m02 * m1.m23 + m0.m03;
		ret.m13 = m0.m10 * m1.m03 + m0.m11 * m1.m13 + m0.m12 * m1.m23 + m0.m13;
		ret.m23 = m0.m20 * m1.m03 + m0.m21 * m1.m13 + m0.m22 * m1.m23 + m0.m23;
		ret.m33 = 1.0f;
	#endif
	return ret;
}

mat4 &mul3(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20;
	ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20;
	ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
	ret.m30 = 0.0f;
	ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
	ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21;
	ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21;
	ret.m31 = 0.0f;
	ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22;
	ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22;
	ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
	ret.m32 = 0.0f;
	ret.m03 = 0.0f;
	ret.m13 = 0.0f;
	ret.m23 = 0.0f;
	ret.m33 = 1.0f;
	return ret;
}

mat4 &mult(mat4 &ret,const mat4 &m,const vec3 &v) {
	#ifdef USE_SSE
		ret.col0 = m.col0;
		ret.col1 = m.col1;
		ret.col2 = m.col2;
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,Z));
		ret.col3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,m.col3));
	#elif USE_ALTIVEC
		ret.col0 = m.col0;
		ret.col1 = m.col1;
		ret.col2 = m.col2;
		ret.col3 = vec_madd(m.col0,vec_perm(v.vec,v.vec,VEC_PERM2(X,X,X,X)),m.col3);
		ret.col3 = vec_madd(m.col1,vec_perm(v.vec,v.vec,VEC_PERM2(Y,Y,Y,Y)),ret.col3);
		ret.col3 = vec_madd(m.col2,vec_perm(v.vec,v.vec,VEC_PERM2(Z,Z,Z,Z)),ret.col3);
	#elif USE_NEON
		ret.col0 = m.col0;
		ret.col1 = m.col1;
		ret.col2 = m.col2;
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		ret.col3 = vmlaq_lane_f32(m.col3,m.col0,low,0);
		ret.col3 = vmlaq_lane_f32(ret.col3,m.col1,low,1);
		ret.col3 = vmlaq_lane_f32(ret.col3,m.col2,high,0);
	#else
		ret.m00 = m.m00; ret.m01 = m.m01; ret.m02 = m.m02;
		ret.m10 = m.m10; ret.m11 = m.m11; ret.m12 = m.m12;
		ret.m20 = m.m20; ret.m21 = m.m21; ret.m22 = m.m22;
		ret.m30 = m.m30; ret.m31 = m.m31; ret.m32 = m.m32;
		ret.m03 = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
		ret.m13 = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
		ret.m23 = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
		ret.m33 = m.m33;
	#endif
	return ret;
}

vec3 &proj(vec3 &ret,const mat4 &m,const vec3 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,Z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,m.col3));
		ret.vec = _mm_div_ps(res_3,_MM_SWIZZLE(res_3,W,W,W,W));
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,X),m.col3);
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,Y),res_0);
		vec_float4 res_2 = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,Z),res_1);
		ret.vec = vec_madd(res_2,vec_rcp_nr(VEC_SWIZZLE(res_2,W,W,W,W)),vec_splats(0.0f));
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmlaq_lane_f32(m.col3,m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,m.col2,high,0);
		ret.vec = vmulq_lane_f32(res_2,vrcp_nr_f32(vget_high_f32(res_2)),1);
	#else
		float iw = Math::rcp(m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33);
		ret.x = (m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03) * iw;
		ret.y = (m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13) * iw;
		ret.z = (m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23) * iw;
	#endif
	return ret;
}

vec4 &proj(vec4 &ret,const mat4 &m,const vec4 &v) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(v.vec,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(m.col1,_MM_SWIZZLE(v.vec,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(m.col2,_MM_SWIZZLE(v.vec,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(m.col3,_MM_SWIZZLE(v.vec,W,W,W,W));
		__m128 res_4 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		ret.vec = _mm_div_ps(res_4,_MM_SWIZZLE(res_4,W,W,W,W));
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 res_0 = vec_madd(m.col0,VEC_SWIZZLE(v.vec,X,X,X,X),zero);
		vec_float4 res_1 = vec_madd(m.col1,VEC_SWIZZLE(v.vec,Y,Y,Y,Y),res_0);
		vec_float4 res_2 = vec_madd(m.col2,VEC_SWIZZLE(v.vec,Z,Z,Z,Z),res_1);
		vec_float4 res_3 = vec_madd(m.col3,VEC_SWIZZLE(v.vec,W,W,W,W),res_2);
		ret.vec = vec_madd(res_3,vec_rcp_nr(VEC_SWIZZLE(res_3,W,W,W,W)),zero);
	#elif USE_NEON
		float32x2_t low = vget_low_f32(v.vec);
		float32x2_t high = vget_high_f32(v.vec);
		float32x4_t res_0 = vmulq_lane_f32(m.col0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,m.col1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,m.col2,high,0);
		float32x4_t res_3 = vmlaq_lane_f32(res_2,m.col3,high,1);
		ret.vec = vmulq_lane_f32(res_3,vrcp_nr_f32(vget_high_f32(res_3)),1);
	#else
		float iw = Math::rcp(m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33 * v.w);
		ret.x = (m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w) * iw;
		ret.y = (m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w) * iw;
		ret.z = (m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w) * iw;
		ret.w = 1.0f;
	#endif
	return ret;
}

dvec3 &proj(dvec3 &ret,const mat4 &m,const dvec3 &v) {
	double iw = Math::rcp(m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33);
	ret.x = (m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03) * iw;
	ret.y = (m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13) * iw;
	ret.z = (m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23) * iw;
	return ret;
}

dvec4 &proj(dvec4 &ret,const mat4 &m,const dvec4 &v) {
	double iw = Math::rcp(m.m30 * v.x + m.m31 * v.y + m.m32 * v.z + m.m33 * v.w);
	ret.x = (m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w) * iw;
	ret.y = (m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w) * iw;
	ret.z = (m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w) * iw;
	ret.w = 1.0f;
	return ret;
}

mat4 &add(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	#ifdef USE_SSE
		ret.col0 = _mm_add_ps(m0.col0,m1.col0);
		ret.col1 = _mm_add_ps(m0.col1,m1.col1);
		ret.col2 = _mm_add_ps(m0.col2,m1.col2);
		ret.col3 = _mm_add_ps(m0.col3,m1.col3);
	#elif USE_ALTIVEC
		ret.col0 = vec_add(m0.col0,m1.col0);
		ret.col1 = vec_add(m0.col1,m1.col1);
		ret.col2 = vec_add(m0.col2,m1.col2);
		ret.col3 = vec_add(m0.col3,m1.col3);
	#elif USE_NEON
		ret.col0 = vaddq_f32(m0.col0,m1.col0);
		ret.col1 = vaddq_f32(m0.col1,m1.col1);
		ret.col2 = vaddq_f32(m0.col2,m1.col2);
		ret.col3 = vaddq_f32(m0.col3,m1.col3);
	#else
		ret.m00 = m0.m00 + m1.m00; ret.m01 = m0.m01 + m1.m01; ret.m02 = m0.m02 + m1.m02; ret.m03 = m0.m03 + m1.m03;
		ret.m10 = m0.m10 + m1.m10; ret.m11 = m0.m11 + m1.m11; ret.m12 = m0.m12 + m1.m12; ret.m13 = m0.m13 + m1.m13;
		ret.m20 = m0.m20 + m1.m20; ret.m21 = m0.m21 + m1.m21; ret.m22 = m0.m22 + m1.m22; ret.m23 = m0.m23 + m1.m23;
		ret.m30 = m0.m30 + m1.m30; ret.m31 = m0.m31 + m1.m31; ret.m32 = m0.m32 + m1.m32; ret.m33 = m0.m33 + m1.m33;
	#endif
	return ret;
}

mat4 &sub(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	#ifdef USE_SSE
		ret.col0 = _mm_sub_ps(m0.col0,m1.col0);
		ret.col1 = _mm_sub_ps(m0.col1,m1.col1);
		ret.col2 = _mm_sub_ps(m0.col2,m1.col2);
		ret.col3 = _mm_sub_ps(m0.col3,m1.col3);
	#elif USE_ALTIVEC
		ret.col0 = vec_sub(m0.col0,m1.col0);
		ret.col1 = vec_sub(m0.col1,m1.col1);
		ret.col2 = vec_sub(m0.col2,m1.col2);
		ret.col3 = vec_sub(m0.col3,m1.col3);
	#elif USE_NEON
		ret.col0 = vsubq_f32(m0.col0,m1.col0);
		ret.col1 = vsubq_f32(m0.col1,m1.col1);
		ret.col2 = vsubq_f32(m0.col2,m1.col2);
		ret.col3 = vsubq_f32(m0.col3,m1.col3);
	#else
		ret.m00 = m0.m00 - m1.m00; ret.m01 = m0.m01 - m1.m01; ret.m02 = m0.m02 - m1.m02; ret.m03 = m0.m03 - m1.m03;
		ret.m10 = m0.m10 - m1.m10; ret.m11 = m0.m11 - m1.m11; ret.m12 = m0.m12 - m1.m12; ret.m13 = m0.m13 - m1.m13;
		ret.m20 = m0.m20 - m1.m20; ret.m21 = m0.m21 - m1.m21; ret.m22 = m0.m22 - m1.m22; ret.m23 = m0.m23 - m1.m23;
		ret.m30 = m0.m30 - m1.m30; ret.m31 = m0.m31 - m1.m31; ret.m32 = m0.m32 - m1.m32; ret.m33 = m0.m33 - m1.m33;
	#endif
	return ret;
}

mat4 &orthonormalize(mat4 &ret,const mat4 &m) {
	#ifdef USE_SSE
		__m128 x_yzx = _MM_SWIZZLE(m.col0,Y,Z,X,W);
		__m128 x_zxy = _MM_SWIZZLE(m.col0,Z,X,Y,W);
		__m128 y_yzx = _MM_SWIZZLE(m.col1,Y,Z,X,W);
		__m128 y_zxy = _MM_SWIZZLE(m.col1,Z,X,Y,W);
		__m128 z = _mm_sub_ps(_mm_mul_ps(x_yzx,y_zxy),_mm_mul_ps(x_zxy,y_yzx));
		__m128 z_yzx = _MM_SWIZZLE(z,Y,Z,X,W);
		__m128 z_zxy = _MM_SWIZZLE(z,Z,X,Y,W);
		__m128 y = _mm_sub_ps(_mm_mul_ps(z_yzx,x_zxy),_mm_mul_ps(z_zxy,x_yzx));
		__m128 col_0 = _mm_mul_ps(m.col0,m.col0);
		__m128 col_1 = _mm_mul_ps(y,y);
		__m128 col_2 = _mm_mul_ps(z,z);
		__m128 res_0 = _mm_shuffle_ps(col_0,col_1,_MM_PERM2(X,Y,X,Y));
		__m128 res_1 = _mm_shuffle_ps(col_0,col_1,_MM_PERM2(Z,W,Z,W));
		__m128 res_2 = _mm_shuffle_ps(col_2,col_2,_MM_PERM2(X,Y,X,Y));
		__m128 res_3 = _mm_shuffle_ps(col_2,col_2,_MM_PERM2(Z,W,Z,W));
		__m128 row_0 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(X,Z,X,Z));
		__m128 row_1 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(Y,W,Y,W));
		__m128 row_2 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(X,Z,X,Z));
		__m128 ilength = _mm_rsqrt_ps_nr(_mm_add_ps(_mm_add_ps(row_0,row_1),row_2));
		ret.col0 = _mm_mul_ps(m.col0,_MM_SWIZZLE(ilength,X,X,X,X));
		ret.col1 = _mm_mul_ps(y,_MM_SWIZZLE(ilength,Y,Y,Y,Y));
		ret.col2 = _mm_mul_ps(z,_MM_SWIZZLE(ilength,Z,Z,Z,Z));
		ret.col3 = m.col3;
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 yzxw = VEC_PERM2(Y,Z,X,W);
		vec_uchar16 zxyw = VEC_PERM2(Z,X,Y,W);
		vec_float4 x_yzx = vec_perm(m.col0,m.col0,yzxw);
		vec_float4 x_zxy = vec_perm(m.col0,m.col0,zxyw);
		vec_float4 y_yzx = vec_perm(m.col1,m.col1,yzxw);
		vec_float4 y_zxy = vec_perm(m.col1,m.col1,zxyw);
		vec_float4 z = vec_sub(vec_madd(x_yzx,y_zxy,zero),vec_madd(x_zxy,y_yzx,zero));
		vec_float4 z_yzx = vec_perm(z,z,yzxw);
		vec_float4 z_zxy = vec_perm(z,z,zxyw);
		vec_float4 y = vec_sub(vec_madd(z_yzx,x_zxy,zero),vec_madd(z_zxy,x_yzx,zero));
		vec_float4 col_0 = vec_madd(m.col0,m.col0,zero);
		vec_float4 col_1 = vec_madd(y,y,zero);
		vec_float4 col_2 = vec_madd(z,z,zero);
		vec_float4 res_0 = vec_perm(col_0,col_1,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_1 = vec_perm(col_0,col_1,VEC_PERM2(Z,W,Z,W));
		vec_float4 res_2 = vec_perm(col_2,col_2,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_3 = vec_perm(col_2,col_2,VEC_PERM2(Z,W,Z,W));
		vec_float4 row_0 = vec_perm(res_0,res_2,VEC_PERM2(X,Z,X,Z));
		vec_float4 row_1 = vec_perm(res_0,res_2,VEC_PERM2(Y,W,Y,W));
		vec_float4 row_2 = vec_perm(res_1,res_3,VEC_PERM2(X,Z,X,Z));
		vec_float4 ilength = vec_rsqrt_nr(vec_add(vec_add(row_0,row_1),row_2));
		ret.col0 = vec_madd(m.col0,VEC_SWIZZLE(ilength,X,X,X,X),zero);
		ret.col1 = vec_madd(y,VEC_SWIZZLE(ilength,Y,Y,Y,Y),zero);
		ret.col2 = vec_madd(z,VEC_SWIZZLE(ilength,Z,Z,Z,Z),zero);
		ret.col3 = m.col3;
	#elif USE_NEON
		float32x4_t z = vcrossq_f32(m.col0,m.col1);
		float32x4_t y = vcrossq_f32(z,m.col0);
		ret.col0 = vnormalize3q_f32(m.col0);
		ret.col1 = vnormalize3q_f32(y);
		ret.col2 = vnormalize3q_f32(z);
		ret.col3 = m.col3;
	#else
		vec3 x = vec3(m.m00,m.m10,m.m20);
		vec3 y = vec3(m.m01,m.m11,m.m21);
		vec3 z = cross(x,y);
		cross(y,z,x);
		x.normalize();
		y.normalize();
		z.normalize();
		ret.m00 = x.x;   ret.m01 = y.x;   ret.m02 = z.x;   ret.m03 = m.m03;
		ret.m10 = x.y;   ret.m11 = y.y;   ret.m12 = z.y;   ret.m13 = m.m13;
		ret.m20 = x.z;   ret.m21 = y.z;   ret.m22 = z.z;   ret.m23 = m.m23;
		ret.m30 = m.m30; ret.m31 = m.m31; ret.m32 = m.m32; ret.m33 = m.m33;
	#endif
	return ret;
}

mat4 &rotation(mat4 &ret,const mat4 &m) {
	ret.m00 = m.m00; ret.m01 = m.m01; ret.m02 = m.m02; ret.m03 = 0.0f;
	ret.m10 = m.m10; ret.m11 = m.m11; ret.m12 = m.m12; ret.m13 = 0.0f;
	ret.m20 = m.m20; ret.m21 = m.m21; ret.m22 = m.m22; ret.m23 = 0.0f;
	ret.m30 = 0.0f;  ret.m31 = 0.0f;  ret.m32 = 0.0f;  ret.m33 = 1.0f;
	return ret;
}

mat4 &transpose(mat4 &ret,const mat4 &m) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_shuffle_ps(m.col0,m.col1,_MM_PERM2(X,Y,X,Y));
		__m128 res_1 = _mm_shuffle_ps(m.col0,m.col1,_MM_PERM2(Z,W,Z,W));
		__m128 res_2 = _mm_shuffle_ps(m.col2,m.col3,_MM_PERM2(X,Y,X,Y));
		__m128 res_3 = _mm_shuffle_ps(m.col2,m.col3,_MM_PERM2(Z,W,Z,W));
		ret.col0 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(X,Z,X,Z));
		ret.col1 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(Y,W,Y,W));
		ret.col2 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(X,Z,X,Z));
		ret.col3 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(Y,W,Y,W));
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_perm(m.col0,m.col1,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_1 = vec_perm(m.col0,m.col1,VEC_PERM2(Z,W,Z,W));
		vec_float4 res_2 = vec_perm(m.col2,m.col3,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_3 = vec_perm(m.col2,m.col3,VEC_PERM2(Z,W,Z,W));
		ret.col0 = vec_perm(res_0,res_2,VEC_PERM2(X,Z,X,Z));
		ret.col1 = vec_perm(res_0,res_2,VEC_PERM2(Y,W,Y,W));
		ret.col2 = vec_perm(res_1,res_3,VEC_PERM2(X,Z,X,Z));
		ret.col3 = vec_perm(res_1,res_3,VEC_PERM2(Y,W,Y,W));
	#else
		ret.m00 = m.m00; ret.m01 = m.m10; ret.m02 = m.m20; ret.m03 = m.m30;
		ret.m10 = m.m01; ret.m11 = m.m11; ret.m12 = m.m21; ret.m13 = m.m31;
		ret.m20 = m.m02; ret.m21 = m.m12; ret.m22 = m.m22; ret.m23 = m.m32;
		ret.m30 = m.m03; ret.m31 = m.m13; ret.m32 = m.m23; ret.m33 = m.m33;
	#endif
	return ret;
}

mat4 &transpose3(mat4 &ret,const mat4 &m) {
	ret.m00 = m.m00; ret.m01 = m.m10; ret.m02 = m.m20; ret.m03 = m.m03;
	ret.m10 = m.m01; ret.m11 = m.m11; ret.m12 = m.m21; ret.m13 = m.m13;
	ret.m20 = m.m02; ret.m21 = m.m12; ret.m22 = m.m22; ret.m23 = m.m23;
	ret.m30 = m.m30; ret.m31 = m.m31; ret.m32 = m.m32; ret.m33 = m.m33;
	return ret;
}

mat4 &inverse(mat4 &ret,const mat4 &m) {
	#ifdef USE_SSE
		__m128 res_0 = _mm_shuffle_ps(m.col0,m.col1,_MM_PERM2(X,Y,X,Y));
		__m128 res_1 = _mm_shuffle_ps(m.col0,m.col1,_MM_PERM2(Z,W,Z,W));
		__m128 res_2 = _mm_shuffle_ps(m.col2,m.col3,_MM_PERM2(X,Y,X,Y));
		__m128 res_3 = _mm_shuffle_ps(m.col2,m.col3,_MM_PERM2(Z,W,Z,W));
		__m128 row_0 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(X,Z,X,Z));
		__m128 row_1 = _mm_shuffle_ps(res_2,res_0,_MM_PERM2(Y,W,Y,W));
		__m128 row_2 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(X,Z,X,Z));
		__m128 row_3 = _mm_shuffle_ps(res_3,res_1,_MM_PERM2(Y,W,Y,W));
		__m128 temp = _mm_mul_ps(row_2,row_3);
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		res_0 = _mm_mul_ps(row_1,temp);
		res_1 = _mm_mul_ps(row_0,temp);
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_0 = _mm_sub_ps(_mm_mul_ps(row_1,temp),res_0);
		res_1 = _mm_sub_ps(_mm_mul_ps(row_0,temp),res_1);
		res_1 = _MM_SWIZZLE(res_1,Z,W,X,Y);
		temp = _mm_mul_ps(row_1,row_2);
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		res_0 = _mm_add_ps(_mm_mul_ps(row_3,temp),res_0);
		res_3 = _mm_mul_ps(row_0,temp);
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_0 = _mm_sub_ps(res_0,_mm_mul_ps(row_3,temp));
		res_3 = _mm_sub_ps(_mm_mul_ps(row_0,temp),res_3);
		res_3 = _MM_SWIZZLE(res_3,Z,W,X,Y);
		temp = _mm_mul_ps(row_3,_MM_SWIZZLE(row_1,Z,W,X,Y));
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		row_2 = _MM_SWIZZLE(row_2,Z,W,X,Y);
		res_0 = _mm_add_ps(_mm_mul_ps(row_2,temp),res_0);
		res_2 = _mm_mul_ps(row_0,temp);
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_0 = _mm_sub_ps(res_0,_mm_mul_ps(row_2,temp));
		res_2 = _mm_sub_ps(_mm_mul_ps(row_0,temp),res_2);
		res_2 = _MM_SWIZZLE(res_2,Z,W,X,Y);
		temp = _mm_mul_ps(row_0,row_1);
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		res_2 = _mm_add_ps(_mm_mul_ps(row_3,temp),res_2);
		res_3 = _mm_sub_ps(_mm_mul_ps(row_2,temp),res_3);
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_2 = _mm_sub_ps(_mm_mul_ps(row_3,temp),res_2);
		res_3 = _mm_sub_ps(res_3,_mm_mul_ps(row_2,temp));
		temp = _mm_mul_ps(row_0,row_3);
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		res_1 = _mm_sub_ps(res_1,_mm_mul_ps(row_2,temp));
		res_2 = _mm_add_ps(_mm_mul_ps(row_1,temp),res_2);
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_1 = _mm_add_ps(_mm_mul_ps(row_2,temp),res_1);
		res_2 = _mm_sub_ps(res_2,_mm_mul_ps(row_1,temp));
		temp = _mm_mul_ps(row_0,row_2);
		temp = _MM_SWIZZLE(temp,Y,X,W,Z);
		res_1 = _mm_add_ps(_mm_mul_ps(row_3,temp),res_1);
		res_3 = _mm_sub_ps(res_3,_mm_mul_ps(row_1,temp));
		temp = _MM_SWIZZLE(temp,Z,W,X,Y);
		res_1 = _mm_sub_ps(res_1,_mm_mul_ps(row_3,temp));
		res_3 = _mm_add_ps(_mm_mul_ps(row_1,temp),res_3);
		__m128 det = _mm_mul_ps(row_0,res_0);
		det = _mm_add_ps(det,_MM_SWIZZLE(det,Y,X,W,Z));
		det = _mm_add_ss(det,_MM_SWIZZLE(det,Z,W,X,Y));
		temp = _MM_SWIZZLE(_mm_rcp_ss_nr(det),X,X,X,X);
		ret.col0 = _mm_mul_ps(res_0,temp);
		ret.col1 = _mm_mul_ps(res_1,temp);
		ret.col2 = _mm_mul_ps(res_2,temp);
		ret.col3 = _mm_mul_ps(res_3,temp);
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_uchar16 yxwz = VEC_PERM2(Y,X,W,Z);
		vec_uchar16 zwxy = VEC_PERM2(Z,W,X,Y);
		vec_float4 res_0 = vec_perm(m.col0,m.col1,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_1 = vec_perm(m.col0,m.col1,VEC_PERM2(Z,W,Z,W));
		vec_float4 res_2 = vec_perm(m.col2,m.col3,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_3 = vec_perm(m.col2,m.col3,VEC_PERM2(Z,W,Z,W));
		vec_float4 row_0 = vec_perm(res_0,res_2,VEC_PERM2(X,Z,X,Z));
		vec_float4 row_1 = vec_perm(res_2,res_0,VEC_PERM2(Y,W,Y,W));
		vec_float4 row_2 = vec_perm(res_1,res_3,VEC_PERM2(X,Z,X,Z));
		vec_float4 row_3 = vec_perm(res_3,res_1,VEC_PERM2(Y,W,Y,W));
		vec_float4 temp = vec_madd(row_2,row_3,zero);
		temp = vec_perm(temp,temp,yxwz);
		res_0 = vec_madd(row_1,temp,zero);
		res_1 = vec_madd(row_0,temp,zero);
		temp = vec_perm(temp,temp,zwxy);
		res_0 = vec_sub(vec_madd(row_1,temp,zero),res_0);
		res_1 = vec_sub(vec_madd(row_0,temp,zero),res_1);
		res_1 = vec_perm(res_1,res_1,zwxy);
		temp = vec_madd(row_1,row_2,zero);
		temp = vec_perm(temp,temp,yxwz);
		res_0 = vec_madd(row_3,temp,res_0);
		res_3 = vec_madd(row_0,temp,zero);
		temp = vec_perm(temp,temp,zwxy);
		res_0 = vec_nmsub(row_3,temp,res_0);
		res_3 = vec_sub(vec_madd(row_0,temp,zero),res_3);
		res_3 = vec_perm(res_3,res_3,zwxy);
		temp = vec_madd(vec_perm(row_1,row_1,zwxy),row_3,zero);
		temp = vec_perm(temp,temp,yxwz);
		row_2 = vec_perm(row_2,row_2,zwxy);
		res_0 = vec_madd(row_2,temp,res_0);
		res_2 = vec_madd(row_0,temp,zero);
		temp = vec_perm(temp,temp,zwxy);
		res_0 = vec_nmsub(row_2,temp,res_0);
		res_2 = vec_sub(vec_madd(row_0,temp,zero),res_2);
		res_2 = vec_perm(res_2,res_2,zwxy);
		temp = vec_madd(row_0,row_1,zero);
		temp = vec_perm(temp,temp,yxwz);
		res_2 = vec_madd(row_3,temp,res_2);
		res_3 = vec_sub(vec_madd(row_2,temp,zero),res_3);
		temp = vec_perm(temp,temp,zwxy);
		res_2 = vec_sub(vec_madd(row_3,temp,zero),res_2);
		res_3 = vec_nmsub(row_2,temp,res_3);
		temp = vec_madd(row_0,row_3,zero);
		temp = vec_perm(temp,temp,yxwz);
		res_1 = vec_nmsub(row_2,temp,res_1);
		res_2 = vec_madd(row_1,temp,res_2);
		temp = vec_perm(temp,temp,zwxy);
		res_1 = vec_madd(row_2,temp,res_1);
		res_2 = vec_nmsub(row_1,temp,res_2);
		temp = vec_madd(row_0,row_2,zero);
		temp = vec_perm(temp,temp,yxwz);
		res_1 = vec_madd(row_3,temp,res_1);
		res_3 = vec_nmsub(row_1,temp,res_3);
		temp = vec_perm(temp,temp,zwxy);
		res_1 = vec_nmsub(row_3,temp,res_1);
		res_3 = vec_madd(row_1,temp,res_3);
		vec_float4 det = vec_madd(row_0,res_0,zero);
		det = vec_add(det,vec_sld(det,det,8));
		det = vec_add(det,vec_sld(det,det,4));
		temp = vec_rcp_nr(det);
		ret.col0 = vec_madd(res_0,temp,zero);
		ret.col1 = vec_madd(res_1,temp,zero);
		ret.col2 = vec_madd(res_2,temp,zero);
		ret.col3 = vec_madd(res_3,temp,zero);
	#elif defined(USE_NEON) && (defined(_WINRT) || defined(_IOS))
		float32x4x4_t transpose = vld4q_f32((const float32_t*)m.mat);
		float32x4_t row_0 = transpose.val[0];
		float32x4_t row_1 = vextq_f32(transpose.val[1],transpose.val[1],2);
		float32x4_t row_2 = transpose.val[2];
		float32x4_t row_3 = vextq_f32(transpose.val[3],transpose.val[3],2);
		float32x4_t temp = vmulq_f32(row_2,row_3);
		temp = vrev64q_f32(temp);
		float32x4_t res_0 = vmulq_f32(row_1,temp);
		float32x4_t res_1 = vmulq_f32(row_0,temp);
		temp = vextq_f32(temp,temp,2);
		res_0 = vnegq_f32(vmlsq_f32(res_0,row_1,temp));
		res_1 = vnegq_f32(vmlsq_f32(res_1,row_0,temp));
		res_1 = vextq_f32(res_1,res_1,2);
		temp = vmulq_f32(row_1,row_2);
		temp = vrev64q_f32(temp);
		res_0 = vmlaq_f32(res_0,row_3,temp);
		float32x4_t res_3 = vmulq_f32(row_0,temp);
		temp = vextq_f32(temp,temp,2);
		res_0 = vmlsq_f32(res_0,row_3,temp);
		res_3 = vnegq_f32(vmlsq_f32(res_3,row_0,temp));
		res_3 = vextq_f32(res_3,res_3,2);
		temp = vmulq_f32(row_3,vextq_f32(row_1,row_1,2));
		temp = vrev64q_f32(temp);
		row_2 = vextq_f32(row_2,row_2,2);
		res_0 = vmlaq_f32(res_0,row_2,temp);
		float32x4_t res_2 = vmulq_f32(row_0,temp);
		temp = vextq_f32(temp,temp,2);
		res_0 = vmlsq_f32(res_0,row_2,temp);
		res_2 = vnegq_f32(vmlsq_f32(res_2,row_0,temp));
		res_2 = vextq_f32(res_2,res_2,2);
		temp = vmulq_f32(row_0,row_1);
		temp = vrev64q_f32(temp);
		res_2 = vmlaq_f32(res_2,row_3,temp);
		res_3 = vnegq_f32(vmlsq_f32(res_3,row_2,temp));
		temp = vextq_f32(temp,temp,2);
		res_2 = vnegq_f32(vmlsq_f32(res_2,row_3,temp));
		res_3 = vmlsq_f32(res_3,row_2,temp);
		temp = vmulq_f32(row_0,row_3);
		temp = vrev64q_f32(temp);
		res_1 = vmlsq_f32(res_1,row_2,temp);
		res_2 = vmlaq_f32(res_2,row_1,temp);
		temp = vextq_f32(temp,temp,2);
		res_1 = vmlaq_f32(res_1,row_2,temp);
		res_2 = vmlsq_f32(res_2,row_1,temp);
		temp = vmulq_f32(row_0,row_2);
		temp = vrev64q_f32(temp);
		res_1 = vmlaq_f32(res_1,row_3,temp);
		res_3 = vmlsq_f32(res_3,row_1,temp);
		temp = vextq_f32(temp,temp,2);
		res_1 = vmlsq_f32(res_1,row_3,temp);
		res_3 = vmlaq_f32(res_3,row_1,temp);
		float32x4_t det = vmulq_f32(row_0,res_0);
		det = vaddq_f32(det,vextq_f32(det,det,2));
		det = vaddq_f32(det,vextq_f32(det,det,1));
		float32x2_t idet = vrcp_nr_f32(vget_low_f32(det));
		ret.col0 = vmulq_lane_f32(res_0,idet,0);
		ret.col1 = vmulq_lane_f32(res_1,idet,0);
		ret.col2 = vmulq_lane_f32(res_2,idet,0);
		ret.col3 = vmulq_lane_f32(res_3,idet,0);
	#elif USE_NEON
		asm volatile(
			"vld4.32  { d0, d2, d4, d6 }, [%r1]!	\n"
			"vld4.32  { d1, d3, d5, d7 }, [%r1]		\n"
			"vext.32    q1, q1, q1, #2				\n"
			"vext.32    q3, q3, q3, #2				\n"
			"vmul.f32   q8, q2, q3					\n"
			"vrev64.32  q8, q8						\n"
			"vmul.f32   q4, q1, q8					\n"
			"vmul.f32   q5, q0, q8					\n"
			"vext.32    q8, q8, q8, #2				\n"
			"vmls.f32   q4, q1, q8					\n"
			"vneg.f32   q4, q4						\n"
			"vmls.f32   q5, q0, q8					\n"
			"vneg.f32   q5, q5						\n"
			"vext.32    q5, q5, q5, #2				\n"
			"vmul.f32   q8, q1, q2					\n"
			"vrev64.32  q8, q8						\n"
			"vmla.f32   q4, q3, q8					\n"
			"vmul.f32   q7, q0, q8					\n"
			"vext.f32   q8, q8, q8, #2				\n"
			"vmls.f32   q4, q3, q8					\n"
			"vmls.f32   q7, q0, q8					\n"
			"vneg.f32   q7, q7						\n"
			"vext.32    q7, q7, q7, #2				\n"
			"vext.32    q8, q1, q1, #2				\n"
			"vmul.f32   q8, q3, q8					\n"
			"vrev64.32  q8, q8						\n"
			"vext.32    q2, q2, q2, #2				\n"
			"vmla.f32   q4, q2, q8					\n"
			"vmul.f32   q6, q0, q8					\n"
			"vext.f32   q8, q8, q8, #2				\n"
			"vmls.f32   q4, q2, q8					\n"
			"vmls.f32   q6, q0, q8					\n"
			"vneg.f32   q6, q6						\n"
			"vext.32    q6, q6, q6, #2				\n"
			"vmul.f32   q8, q0, q1					\n"
			"vrev64.32  q8, q8						\n"
			"vmla.f32   q6, q3, q8					\n"
			"vmls.f32   q7, q2, q8					\n"
			"vneg.f32   q7, q7						\n"
			"vext.f32   q8, q8, q8, #2				\n"
			"vmls.f32   q6, q3, q8					\n"
			"vneg.f32   q6, q6						\n"
			"vmls.f32   q7, q2, q8					\n"
			"vmul.f32   q8, q0, q3					\n"
			"vrev64.32  q8, q8						\n"
			"vmls.f32   q5, q2, q8					\n"
			"vmla.f32   q6, q1, q8					\n"
			"vext.f32   q8, q8, q8, #2				\n"
			"vmla.f32   q5, q2, q8					\n"
			"vmls.f32   q6, q1, q8					\n"
			"vmul.f32   q8, q0, q2					\n"
			"vrev64.32  q8, q8						\n"
			"vmla.f32   q5, q3, q8					\n"
			"vmls.f32   q7, q1, q8					\n"
			"vext.f32   q8, q8, q8, #2				\n"
			"vmls.f32   q5, q3, q8					\n"
			"vmla.f32   q7, q1, q8					\n"
			"vmul.f32   q0, q0, q4					\n"
			"vext.32    q1, q0, q0, #2				\n"
			"vadd.f32   q0, q0, q1					\n"
			"vext.32    q1, q0, q0, #1				\n"
			"vadd.f32   q0, q0, q1					\n"
			"vrecpe.f32 q1, q0						\n"
			"vrecps.f32 q2, q1, q0					\n"
			"vmul.f32   q1, q1, q2					\n"
			"vrecps.f32 q2, q1, q0					\n"
			"vmul.f32   q8, q1, q2					\n"
			"vmul.f32   q0, q4, q8					\n"
			"vmul.f32   q1, q5, q8					\n"
			"vmul.f32   q2, q6, q8					\n"
			"vmul.f32   q3, q7, q8					\n"
			"vstm      %r0, { q0 - q3 }				\n"
			: : "r"(ret.mat), "r"(m.mat)
			: "q0", "q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8"
		);
	#else
		float temp[12];
		temp[0]  = m.m22 * m.m33;
		temp[1]  = m.m23 * m.m32;
		temp[2]  = m.m21 * m.m33;
		temp[3]  = m.m23 * m.m31;
		temp[4]  = m.m21 * m.m32;
		temp[5]  = m.m22 * m.m31;
		temp[6]  = m.m20 * m.m33;
		temp[7]  = m.m23 * m.m30;
		temp[8]  = m.m20 * m.m32;
		temp[9]  = m.m22 * m.m30;
		temp[10] = m.m20 * m.m31;
		temp[11] = m.m21 * m.m30;
		ret.m00  = temp[0] * m.m11 + temp[3] * m.m12 + temp[4]  * m.m13;
		ret.m00 -= temp[1] * m.m11 + temp[2] * m.m12 + temp[5]  * m.m13;
		ret.m10  = temp[1] * m.m10 + temp[6] * m.m12 + temp[9]  * m.m13;
		ret.m10 -= temp[0] * m.m10 + temp[7] * m.m12 + temp[8]  * m.m13;
		ret.m20  = temp[2] * m.m10 + temp[7] * m.m11 + temp[10] * m.m13;
		ret.m20 -= temp[3] * m.m10 + temp[6] * m.m11 + temp[11] * m.m13;
		ret.m30  = temp[5] * m.m10 + temp[8] * m.m11 + temp[11] * m.m12;
		ret.m30 -= temp[4] * m.m10 + temp[9] * m.m11 + temp[10] * m.m12;
		ret.m01  = temp[1] * m.m01 + temp[2] * m.m02 + temp[5]  * m.m03;
		ret.m01 -= temp[0] * m.m01 + temp[3] * m.m02 + temp[4]  * m.m03;
		ret.m11  = temp[0] * m.m00 + temp[7] * m.m02 + temp[8]  * m.m03;
		ret.m11 -= temp[1] * m.m00 + temp[6] * m.m02 + temp[9]  * m.m03;
		ret.m21  = temp[3] * m.m00 + temp[6] * m.m01 + temp[11] * m.m03;
		ret.m21 -= temp[2] * m.m00 + temp[7] * m.m01 + temp[10] * m.m03;
		ret.m31  = temp[4] * m.m00 + temp[9] * m.m01 + temp[10] * m.m02;
		ret.m31 -= temp[5] * m.m00 + temp[8] * m.m01 + temp[11] * m.m02;
		temp[0]  = m.m02 * m.m13;
		temp[1]  = m.m03 * m.m12;
		temp[2]  = m.m01 * m.m13;
		temp[3]  = m.m03 * m.m11;
		temp[4]  = m.m01 * m.m12;
		temp[5]  = m.m02 * m.m11;
		temp[6]  = m.m00 * m.m13;
		temp[7]  = m.m03 * m.m10;
		temp[8]  = m.m00 * m.m12;
		temp[9]  = m.m02 * m.m10;
		temp[10] = m.m00 * m.m11;
		temp[11] = m.m01 * m.m10;
		ret.m02  = temp[0]  * m.m31 + temp[3]  * m.m32 + temp[4]  * m.m33;
		ret.m02 -= temp[1]  * m.m31 + temp[2]  * m.m32 + temp[5]  * m.m33;
		ret.m12  = temp[1]  * m.m30 + temp[6]  * m.m32 + temp[9]  * m.m33;
		ret.m12 -= temp[0]  * m.m30 + temp[7]  * m.m32 + temp[8]  * m.m33;
		ret.m22  = temp[2]  * m.m30 + temp[7]  * m.m31 + temp[10] * m.m33;
		ret.m22 -= temp[3]  * m.m30 + temp[6]  * m.m31 + temp[11] * m.m33;
		ret.m32  = temp[5]  * m.m30 + temp[8]  * m.m31 + temp[11] * m.m32;
		ret.m32 -= temp[4]  * m.m30 + temp[9]  * m.m31 + temp[10] * m.m32;
		ret.m03  = temp[2]  * m.m22 + temp[5]  * m.m23 + temp[1]  * m.m21;
		ret.m03 -= temp[4]  * m.m23 + temp[0]  * m.m21 + temp[3]  * m.m22;
		ret.m13  = temp[8]  * m.m23 + temp[0]  * m.m20 + temp[7]  * m.m22;
		ret.m13 -= temp[6]  * m.m22 + temp[9]  * m.m23 + temp[1]  * m.m20;
		ret.m23  = temp[6]  * m.m21 + temp[11] * m.m23 + temp[3]  * m.m20;
		ret.m23 -= temp[10] * m.m23 + temp[2]  * m.m20 + temp[7]  * m.m21;
		ret.m33  = temp[10] * m.m22 + temp[4]  * m.m20 + temp[9]  * m.m21;
		ret.m33 -= temp[8]  * m.m21 + temp[11] * m.m22 + temp[5]  * m.m20;
		float idet = Math::rcp(m.m00 * ret.m00 + m.m01 * ret.m10 + m.m02 * ret.m20 + m.m03 * ret.m30);
		ret.m00 *= idet;
		ret.m10 *= idet;
		ret.m20 *= idet;
		ret.m30 *= idet;
		ret.m01 *= idet;
		ret.m11 *= idet;
		ret.m21 *= idet;
		ret.m31 *= idet;
		ret.m02 *= idet;
		ret.m12 *= idet;
		ret.m22 *= idet;
		ret.m32 *= idet;
		ret.m03 *= idet;
		ret.m13 *= idet;
		ret.m23 *= idet;
		ret.m33 *= idet;
	#endif
	return ret;
}

mat4 &inverse4(mat4 &ret,const mat4 &m) {
	float idet = Math::rcp(determinant3(m));
	ret.m00 =  (m.m11 * m.m22 - m.m12 * m.m21) * idet;
	ret.m10 = -(m.m10 * m.m22 - m.m12 * m.m20) * idet;
	ret.m20 =  (m.m10 * m.m21 - m.m11 * m.m20) * idet;
	ret.m30 = 0.0f;
	ret.m01 = -(m.m01 * m.m22 - m.m02 * m.m21) * idet;
	ret.m11 =  (m.m00 * m.m22 - m.m02 * m.m20) * idet;
	ret.m21 = -(m.m00 * m.m21 - m.m01 * m.m20) * idet;
	ret.m31 = 0.0f;
	ret.m02 =  (m.m01 * m.m12 - m.m02 * m.m11) * idet;
	ret.m12 = -(m.m00 * m.m12 - m.m02 * m.m10) * idet;
	ret.m22 =  (m.m00 * m.m11 - m.m01 * m.m10) * idet;
	ret.m32 = 0.0f;
	ret.m03 = -(ret.m00 * m.m03 + ret.m01 * m.m13 + ret.m02 * m.m23);
	ret.m13 = -(ret.m10 * m.m03 + ret.m11 * m.m13 + ret.m12 * m.m23);
	ret.m23 = -(ret.m20 * m.m03 + ret.m21 * m.m13 + ret.m22 * m.m23);
	ret.m33 = 1.0f;
	return ret;
}

mat4 &lerp(mat4 &ret,const mat4 &m0,const mat4 &m1,float k) {
	vec3 positions[3];
	quat rotations[3];
	vec3 scales[3];
	decomposeTransform(m0,positions[0],rotations[0],scales[0]);
	decomposeTransform(m1,positions[1],rotations[1],scales[1]);
	lerp(positions[2],positions[0],positions[1],k);
	slerp(rotations[2],rotations[0],rotations[1],k);
	lerp(scales[2],scales[0],scales[1],k);
	return composeTransform(ret,positions[2],rotations[2],scales[2]);
}

/*
 */
mat4 orthonormalize(const mat4 &m) {
	mat4 ret;
	return orthonormalize(ret,m);
}

mat4 rotation(const mat4 &m) {
	mat4 ret;
	return rotation(ret,m);
}

mat4 transpose(const mat4 &m) {
	mat4 ret;
	return transpose(ret,m);
}

mat4 transpose3(const mat4 &m) {
	mat4 ret;
	return transpose3(ret,m);
}

mat4 inverse(const mat4 &m) {
	mat4 ret;
	return inverse(ret,m);
}

mat4 inverse4(const mat4 &m) {
	mat4 ret;
	return inverse4(ret,m);
}

mat4 lerp(const mat4 &m0,const mat4 &m1,float k) {
	mat4 ret;
	return lerp(ret,m0,m1,k);
}

/*
 */
mat4 translate(const vec3 &v) {
	mat4 ret;
	ret.setTranslate(v);
	return ret;
}

mat4 translate(float x,float y,float z) {
	mat4 ret;
	ret.setTranslate(vec3(x,y,z));
	return ret;
}

mat4 rotate(const vec3 &axis,float angle) {
	mat4 ret;
	ret.setRotate(axis,angle);
	return ret;
}

mat4 rotate(float x,float y,float z,float angle) {
	return rotate(vec3(x,y,z),angle);
}

mat4 rotate(const quat &q) {
	return mat4(q.getMat3());
}

mat4 rotateX(float angle) {
	mat4 ret;
	ret.setRotateX(angle);
	return ret;
}

mat4 rotateY(float angle) {
	mat4 ret;
	ret.setRotateY(angle);
	return ret;
}

mat4 rotateZ(float angle) {
	mat4 ret;
	ret.setRotateZ(angle);
	return ret;
}

mat4 scale(const vec3 &v) {
	mat4 ret;
	ret.setScale(v);
	return ret;
}

mat4 scale(float x,float y,float z) {
	return scale(vec3(x,y,z));
}

/*
 */
mat4 reflect(const vec4 &plane) {
	mat4 ret;
	float x = plane.x;
	float y = plane.y;
	float z = plane.z;
	float x2 = x * 2.0f;
	float y2 = y * 2.0f;
	float z2 = z * 2.0f;
	ret.m00 = 1.0f - x * x2; ret.m01 = -y * x2;       ret.m02 = -z * x2;       ret.m03 = -plane.w * x2;
	ret.m10 = -x * y2;       ret.m11 = 1.0f - y * y2; ret.m12 = -z * y2;       ret.m13 = -plane.w * y2;
	ret.m20 = -x * z2;       ret.m21 = -y * z2;       ret.m22 = 1.0f - z * z2; ret.m23 = -plane.w * z2;
	ret.m30 = 0.0f;          ret.m31 = 0.0f;          ret.m32 = 0.0f;          ret.m33 = 1.0f;
	return ret;
}

mat4 ortho(float l,float r,float b,float t,float n,float f) {
	mat4 ret;
	float rl = r - l;
	float tb = t - b;
	float fn = f - n;
	ret.m00 = 2.0f / rl; ret.m01 = 0.0f;      ret.m02 = 0.0f;       ret.m03 = -(r + l) / rl;
	ret.m10 = 0.0f;      ret.m11 = 2.0f / tb; ret.m12 = 0.0f;       ret.m13 = -(t + b) / tb;
	ret.m20 = 0.0f;      ret.m21 = 0.0f;      ret.m22 = -2.0f / fn; ret.m23 = -(f + n) / fn;
	ret.m30 = 0.0f;      ret.m31 = 0.0f;      ret.m32 = 0.0f;       ret.m33 = 1.0f;
	return ret;
}

mat4 frustum(float l,float r,float b,float t,float n,float f) {
	mat4 ret;
	float rl = r - l;
	float tb = t - b;
	float fn = f - n;
	ret.m00 = 2.0f * n / rl; ret.m01 = 0.0f;          ret.m02 = (r + l) / rl;  ret.m03 = 0.0f;
	ret.m10 = 0.0f;          ret.m11 = 2.0f * n / tb; ret.m12 = (t + b) / tb;  ret.m13 = 0.0f;
	ret.m20 = 0.0f;          ret.m21 = 0.0f;          ret.m22 = -(f + n) / fn; ret.m23 = -2.0f * f * n / fn;
	ret.m30 = 0.0f;          ret.m31 = 0.0f;          ret.m32 = -1.0f;         ret.m33 = 0.0f;
	return ret;
}

mat4 perspective(float fov,float aspect,float n,float f) {
	mat4 ret;
	float h = 1.0f;
	float w = aspect;
	if(!compare(fov,90.0f)) {
		h = Math::tan(fov * DEG2RAD * 0.5f);
		w = h * aspect;
	}
	float fn = f - n;
	ret.m00 = 1.0f / w; ret.m01 = 0.0f;     ret.m02 = 0.0f;          ret.m03 = 0.0f;
	ret.m10 = 0.0f;     ret.m11 = 1.0f / h; ret.m12 = 0.0f;          ret.m13 = 0.0f;
	ret.m20 = 0.0f;     ret.m21 = 0.0f;     ret.m22 = -(f + n) / fn; ret.m23 = -2.0f * f * n / fn;
	ret.m30 = 0.0f;     ret.m31 = 0.0f;     ret.m32 = -1.0f;         ret.m33 = 0.0f;
	return ret;
}

mat4 setTo(const vec3 &position,const vec3 &direction,const vec3 &up) {
	mat4 ret;
	vec3 z = normalize(position - direction);
	vec3 x = normalize(cross(up,z));
	vec3 y = normalize(cross(z,x));
	ret.m00 = x.x;  ret.m01 = y.x;  ret.m02 = z.x;  ret.m03 = position.x;
	ret.m10 = x.y;  ret.m11 = y.y;  ret.m12 = z.y;  ret.m13 = position.y;
	ret.m20 = x.z;  ret.m21 = y.z;  ret.m22 = z.z;  ret.m23 = position.z;
	ret.m30 = 0.0f; ret.m31 = 0.0f; ret.m32 = 0.0f; ret.m33 = 1.0f;
	return ret;
}

mat4 lookAt(const vec3 &position,const vec3 &direction,const vec3 &up) {
	mat4 ret,m0,m1;
	vec3 z = normalize(position - direction);
	vec3 x = normalize(cross(up,z));
	vec3 y = normalize(cross(z,x));
	m0.m00 = x.x;  m0.m01 = x.y;  m0.m02 = x.z;  m0.m03 = 0.0f;
	m0.m10 = y.x;  m0.m11 = y.y;  m0.m12 = y.z;  m0.m13 = 0.0f;
	m0.m20 = z.x;  m0.m21 = z.y;  m0.m22 = z.z;  m0.m23 = 0.0f;
	m0.m30 = 0.0f; m0.m31 = 0.0f; m0.m32 = 0.0f; m0.m33 = 1.0f;
	m1.setTranslate(-position);
	return mul(ret,m0,m1);
}

mat4 obliqueProjection(const mat4 &projection,const vec4 &plane) {
	mat4 ret,transform;
	vec4 v = plane * inverse(projection);
	int z = Math::abs(Math::ftoi(v.z));
	if(z) v /= Math::itof(z);
	v.w -= 1.0f;
	if(v.z < 0.0f) v = -v;
	transform.setIdentity();
	transform.setRow(2,v);
	return mul(ret,transform,projection);
}

mat4 symmetryProjection(const mat4 &projection) {
	mat4 ret;
	float r = (projection.m02 + 1.0f) / projection.m00;
	float t = (projection.m12 + 1.0f) / projection.m11;
	float l = r - 2.0f / projection.m00;
	float b = t - 2.0f / projection.m11;
	float w = max(Math::abs(r),Math::abs(l));
	float h = max(Math::abs(t),Math::abs(b));
	ret.m00 = 1.0f / w; ret.m01 = 0.0f;     ret.m02 = 0.0f;           ret.m03 = 0.0f;
	ret.m10 = 0.0f;     ret.m11 = 1.0f / h; ret.m12 = 0.0f;           ret.m13 = 0.0f;
	ret.m20 = 0.0f;     ret.m21 = 0.0f;     ret.m22 = projection.m22; ret.m23 = projection.m23;
	ret.m30 = 0.0f;     ret.m31 = 0.0f;     ret.m32 = -1.0f;          ret.m33 = 0.0f;
	return ret;
}

mat4 cubeTransform(int face) {
	static const mat4 transforms[6] = {
		rotateY(90.0f) * rotateX(180.0f),
		rotateY(-90.0f) * rotateX(180.0f),
		rotateX(-90.0f),
		rotateX(90.0f),
		rotateX(180.0f),
		rotateZ(180.0f),
	};
	assert((unsigned int)face <= 5 && "cubeTransform(): bad face number");
	return transforms[face];
}

/*
 */
void decomposeTransform(const mat4 &m,vec4 &position,quat &rot) {
	mat3 rotate,scale;
	mat3 rotation = mat3(m);
	orthonormalize(rotate,rotation);
	mul(scale,transpose(rotate),rotation);
	position.x = m.m03;
	position.y = m.m13;
	position.z = m.m23;
	position.w = (scale.m00 + scale.m11 + scale.m22) * (1.0f / 3.0f);
	rot = rotate.getQuat();
}

mat4 &composeTransform(mat4 &ret,const vec4 &position,const quat &rot) {
	float x2 = (rot.x + rot.x) * position.w;
	float y2 = (rot.y + rot.y) * position.w;
	float z2 = (rot.z + rot.z) * position.w;
	float xx2 = rot.x * x2;
	float yy2 = rot.y * y2;
	float zz2 = rot.z * z2;
	float zx2 = rot.z * x2;
	float xy2 = rot.x * y2;
	float yz2 = rot.y * z2;
	float wx2 = rot.w * x2;
	float wy2 = rot.w * y2;
	float wz2 = rot.w * z2;
	ret.m00 = position.w - yy2 - zz2;
	ret.m10 = xy2 + wz2;
	ret.m20 = zx2 - wy2;
	ret.m30 = 0.0f;
	ret.m01 = xy2 - wz2;
	ret.m11 = position.w - xx2 - zz2;
	ret.m21 = yz2 + wx2;
	ret.m31 = 0.0f;
	ret.m02 = zx2 + wy2;
	ret.m12 = yz2 - wx2;
	ret.m22 = position.w - xx2 - yy2;
	ret.m32 = 0.0f;
	ret.m03 = position.x;
	ret.m13 = position.y;
	ret.m23 = position.z;
	ret.m33 = 1.0f;
	return ret;
}

void decomposeTransform(const mat4 &m,vec3 &position,quat &rot,vec3 &s) {
	mat3 rotate,scale;
	mat3 rotation = mat3(m);
	orthonormalize(rotate,rotation);
	mul(scale,transpose(rotate),rotation);
	position.x = m.m03;
	position.y = m.m13;
	position.z = m.m23;
	rot = rotate.getQuat();
	s.x = scale.m00;
	s.y = scale.m11;
	s.z = scale.m22;
}

mat4 &composeTransform(mat4 &ret,const vec3 &position,const quat &rot,const vec3 &s) {
	mat3 rotation,scale;
	scale.setDiagonal(s);
	ret.set(mul(rotation,rot.getMat3(),scale),position);
	return ret;
}

void decomposeProjection(const mat4 &projection,float &znear,float &zfar) {
	if(compare(projection.m32,-1.0f)) {
		znear = projection.m23 / (projection.m22 - 1.0f);
		zfar = projection.m23 / (projection.m22 + 1.0f);
	} else {
		znear = (projection.m23 + 1.0f) / projection.m22;
		zfar = (projection.m23 - 1.0f) / projection.m22;
	}
}

/*
 */
const mat4 &hardwareProjectionGL(const mat4 &projection,int width,int height) {
	return projection;
}

const mat4 &hardwareProjectionD3D9(const mat4 &projection,int width,int height) {
	static mat4 ret,temp;
	static const mat4 offset = translate(0.0f,0.0f,0.5f) * scale(1.0f,1.0f,0.5f);
	return mul(ret,translate(-1.0f / width,1.0f / height,0.0f),mul(temp,offset,projection));
}

const mat4 &hardwareProjectionD3D10(const mat4 &projection,int width,int height) {
	static mat4 ret;
	static const mat4 offset = translate(0.0f,0.0f,0.5f) * scale(1.0f,1.0f,0.5f);
	return mul(ret,offset,projection);
}

const mat4 &(*hardwareProjection)(const mat4 &projection,int width,int height) = hardwareProjectionGL;

/*
 */
void Math::setGL() {
	hardwareProjection = hardwareProjectionGL;
}

void Math::setD3D9() {
	hardwareProjection = hardwareProjectionD3D9;
}

void Math::setD3D10() {
	hardwareProjection = hardwareProjectionD3D10;
}

/******************************************************************************\
*
* dmat4
*
\******************************************************************************/

/*
 */
static const double dmat4_identity_values[12] = {
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 0.0,
};

/*
 */
const dmat4 dmat4_zero(0.0);
const dmat4 dmat4_one(1.0);
const dmat4 dmat4_identity(dmat4_identity_values);

/*
 */
dmat4::dmat4(double v) {
	m00 = v; m01 = v; m02 = v; m03 = v;
	m10 = v; m11 = v; m12 = v; m13 = v;
	m20 = v; m21 = v; m22 = v; m23 = v;
}

dmat4::dmat4(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = 0.0; m13 = 0.0;
	m20 = 0.0;   m21 = 0.0;   m22 = 1.0; m23 = 0.0;
}

dmat4::dmat4(const mat3 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0;
}

dmat4::dmat4(const mat4 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
}

dmat4::dmat4(const dmat4 &m) {
	#ifdef USE_SSE2
		col0 = m.col0; col1 = m.col1;
		col2 = m.col2; col3 = m.col3;
		col4 = m.col4; col5 = m.col5;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
		m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
		m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
	#endif
}

dmat4::dmat4(const quat &q) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0;
}

dmat4::dmat4(const double *m) {
	m00 = m[0]; m01 = m[3]; m02 = m[6]; m03 = m[9];
	m10 = m[1]; m11 = m[4]; m12 = m[7]; m13 = m[10];
	m20 = m[2]; m21 = m[5]; m22 = m[8]; m23 = m[11];
}

dmat4::dmat4(const mat3 &m,const dvec3 &v) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
}

dmat4::dmat4(const quat &q,const dvec3 &v) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
}

/*
 */
dmat4 dmat4::operator-() const {
	dmat4 ret;
	ret.m00 = -m00; ret.m01 = -m01; ret.m02 = -m02; ret.m03 = -m03;
	ret.m10 = -m10; ret.m11 = -m11; ret.m12 = -m12; ret.m13 = -m13;
	ret.m20 = -m20; ret.m21 = -m21; ret.m22 = -m22; ret.m23 = -m23;
	return ret;
}

dmat4 &dmat4::operator*=(double v) {
	return mul(*this,*this,v);
}

dmat4 &dmat4::operator*=(const dmat4 &m) {
	return mul(*this,dmat4(*this),m);
}

dmat4 &dmat4::operator+=(const dmat4 &m) {
	return add(*this,*this,m);
}

dmat4 &dmat4::operator-=(const dmat4 &m) {
	return sub(*this,*this,m);
}

/*
 */
void dmat4::set(const mat2 &m) {
	m00 = m.m00; m01 = m.m01; m02 = 0.0; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = 0.0; m13 = 0.0;
	m20 = 0.0;   m21 = 0.0;   m22 = 1.0; m23 = 0.0;
}

void dmat4::set(const mat3 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0;
}

void dmat4::set(const mat4 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
}

void dmat4::set(const dmat4 &m) {
	#ifdef USE_SSE2
		col0 = m.col0; col1 = m.col1;
		col2 = m.col2; col3 = m.col3;
		col4 = m.col4; col5 = m.col5;
	#else
		m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = m.m03;
		m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = m.m13;
		m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = m.m23;
	#endif
}

void dmat4::set(const quat &q) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0;
}

void dmat4::set(const double *m) {
	m00 = m[0]; m01 = m[3]; m02 = m[6]; m03 = m[9];
	m10 = m[1]; m11 = m[4]; m12 = m[7]; m13 = m[10];
	m20 = m[2]; m21 = m[5]; m22 = m[8]; m23 = m[11];
}

void dmat4::set(const mat3 &m,const dvec3 &v) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
}

void dmat4::set(const quat &q,const dvec3 &v) {
	mat3 m = q.getMat3();
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = v.x;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = v.y;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = v.z;
}

void dmat4::get(double *m) const {
	m[0] = m00; m[3] = m01; m[6] = m02; m[9] = m03;
	m[1] = m10; m[4] = m11; m[7] = m12; m[10] = m13;
	m[2] = m20; m[5] = m21; m[8] = m22; m[11] = m23;
}

/*
 */
void dmat4::setRow(int row,const dvec4 &v) {
	assert((unsigned int)row < 4 && "dmat4::setRow(): bad row");
	switch(row) {
		case 0: m00 = v.x; m01 = v.y; m02 = v.z; m03 = v.w; break;
		case 1: m10 = v.x; m11 = v.y; m12 = v.z; m13 = v.w; break;
		case 2: m20 = v.x; m21 = v.y; m22 = v.z; m23 = v.w; break;
	}
}

void dmat4::setRow3(int row,const dvec3 &v) {
	assert((unsigned int)row < 4 && "dmat4::setRow3(): bad row");
	switch(row) {
		case 0: m00 = v.x; m01 = v.y; m02 = v.z; break;
		case 1: m10 = v.x; m11 = v.y; m12 = v.z; break;
		case 2: m20 = v.x; m21 = v.y; m22 = v.z; break;
	}
}

dvec4 dmat4::getRow(int row) const {
	assert((unsigned int)row < 4 && "dmat4::getRow(): bad row");
	switch(row) {
		case 0: return dvec4(m00,m01,m02,m03);
		case 1: return dvec4(m10,m11,m12,m13);
		case 2: return dvec4(m20,m21,m22,m23);
		case 3: return dvec4(0.0,0.0,0.0,1.0);
	}
	return dvec4_zero;
}

dvec3 dmat4::getRow3(int row) const {
	assert((unsigned int)row < 4 && "dmat4::getRow3(): bad row");
	switch(row) {
		case 0: return dvec3(m00,m01,m02);
		case 1: return dvec3(m10,m11,m12);
		case 2: return dvec3(m20,m21,m22);
		case 3: return dvec3(0.0,0.0,0.0);
	}
	return dvec3_zero;
}

/*
 */
void dmat4::setColumn(int column,const dvec4 &v) {
	assert((unsigned int)column < 4 && "dmat4::setColumn(): bad column");
	switch(column) {
		case 0: m00 = v.x; m10 = v.y; m20 = v.z; break;
		case 1: m01 = v.x; m11 = v.y; m21 = v.z; break;
		case 2: m02 = v.x; m12 = v.y; m22 = v.z; break;
		case 3: m03 = v.x; m13 = v.y; m23 = v.z; break;
	}
}

void dmat4::setColumn3(int column,const dvec3 &v) {
	assert((unsigned int)column < 4 && "dmat4::setColumn3(): bad column");
	switch(column) {
		case 0: m00 = v.x; m10 = v.y; m20 = v.z; break;
		case 1: m01 = v.x; m11 = v.y; m21 = v.z; break;
		case 2: m02 = v.x; m12 = v.y; m22 = v.z; break;
		case 3: m03 = v.x; m13 = v.y; m23 = v.z; break;
	}
}

dvec4 dmat4::getColumn(int column) const {
	assert((unsigned int)column < 4 && "dmat4::getColumn(): bad column");
	switch(column) {
		case 0: return dvec4(m00,m10,m20,0.0);
		case 1: return dvec4(m01,m11,m21,0.0);
		case 2: return dvec4(m02,m12,m22,0.0);
		case 3: return dvec4(m03,m13,m23,1.0);
	}
	return dvec4_zero;
}

dvec3 dmat4::getColumn3(int column) const {
	assert((unsigned int)column < 4 && "dmat4::getColumn3(): bad column");
	switch(column) {
		case 0: return dvec3(m00,m10,m20);
		case 1: return dvec3(m01,m11,m21);
		case 2: return dvec3(m02,m12,m22);
		case 3: return dvec3(m03,m13,m23);
	}
	return dvec3_zero;
}

/*
 */
void dmat4::setZero() {
	m00 = 0.0; m01 = 0.0; m02 = 0.0; m03 = 0.0;
	m10 = 0.0; m11 = 0.0; m12 = 0.0; m13 = 0.0;
	m20 = 0.0; m21 = 0.0; m22 = 0.0; m23 = 0.0;
}

void dmat4::setIdentity() {
	m00 = 1.0; m01 = 0.0; m02 = 0.0; m03 = 0.0;
	m10 = 0.0; m11 = 1.0; m12 = 0.0; m13 = 0.0;
	m20 = 0.0; m21 = 0.0; m22 = 1.0; m23 = 0.0;
}

void dmat4::setTranslate(const dvec3 &v) {
	m00 = 1.0; m01 = 0.0; m02 = 0.0; m03 = v.x;
	m10 = 0.0; m11 = 1.0; m12 = 0.0; m13 = v.y;
	m20 = 0.0; m21 = 0.0; m22 = 1.0; m23 = v.z;
}

void dmat4::setRotate(const dvec3 &axis,double angle) {
	double s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	dvec3 v = normalize(axis);
	double xx = v.x * v.x;
	double yy = v.y * v.y;
	double zz = v.z * v.z;
	double xy = v.x * v.y;
	double yz = v.y * v.z;
	double zx = v.z * v.x;
	double xs = v.x * s;
	double ys = v.y * s;
	double zs = v.z * s;
	m00 = (1.0 - c) * xx + c;  m01 = (1.0 - c) * xy - zs; m02 = (1.0 - c) * zx + ys; m03 = 0.0;
	m10 = (1.0 - c) * xy + zs; m11 = (1.0 - c) * yy + c;  m12 = (1.0 - c) * yz - xs; m13 = 0.0;
	m20 = (1.0 - c) * zx - ys; m21 = (1.0 - c) * yz + xs; m22 = (1.0 - c) * zz + c;  m23 = 0.0;
}

void dmat4::setRotateX(double angle) {
	double s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = 1.0; m01 = 0.0; m02 = 0.0; m03 = 0.0;
	m10 = 0.0; m11 = c;   m12 = -s;  m13 = 0.0;
	m20 = 0.0; m21 = s;   m22 = c;   m23 = 0.0;
}

void dmat4::setRotateY(double angle) {
	double s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;   m01 = 0.0; m02 = s;   m03 = 0.0;
	m10 = 0.0; m11 = 1.0; m12 = 0.0; m13 = 0.0;
	m20 = -s;  m21 = 0.0; m22 = c;   m23 = 0.0;
}

void dmat4::setRotateZ(double angle) {
	double s,c;
	Math::sincos(angle * DEG2RAD,s,c);
	m00 = c;   m01 = -s;  m02 = 0.0; m03 = 0.0;
	m10 = s;   m11 = c;   m12 = 0.0; m13 = 0.0;
	m20 = 0.0; m21 = 0.0; m22 = 1.0; m23 = 0.0;
}

void dmat4::setScale(const dvec3 &v) {
	m00 = v.x; m01 = 0.0; m02 = 0.0; m03 = 0.0;
	m10 = 0.0; m11 = v.y; m12 = 0.0; m13 = 0.0;
	m20 = 0.0; m21 = 0.0; m22 = v.z; m23 = 0.0;
}

/*
 */
int operator==(const dmat4 &m0,const dmat4 &m1) {
	return compare(m0,m1);
}

int operator!=(const dmat4 &m0,const dmat4 &m1) {
	return !compare(m0,m1);
}

dmat4 operator*(const dmat4 &m,double v) {
	dmat4 ret;
	return mul(ret,m,v);
}

vec2 operator*(const dmat4 &m,const vec2 &v) {
	vec2 ret;
	return mul(ret,m,v);
}

vec2 operator*(const vec2 &v,const dmat4 &m) {
	vec2 ret;
	return mul(ret,v,m);
}

vec3 operator*(const dmat4 &m,const vec3 &v) {
	vec3 ret;
	return mul(ret,m,v);
}

vec3 operator*(const vec3 &v,const dmat4 &m) {
	vec3 ret;
	return mul(ret,v,m);
}

vec4 operator*(const dmat4 &m,const vec4 &v) {
	vec4 ret;
	return mul(ret,m,v);
}

vec4 operator*(const vec4 &v,const dmat4 &m) {
	vec4 ret;
	return mul(ret,v,m);
}

dvec2 operator*(const dmat4 &m,const dvec2 &v) {
	dvec2 ret;
	return mul(ret,m,v);
}

dvec2 operator*(const dvec2 &v,const dmat4 &m) {
	dvec2 ret;
	return mul(ret,v,m);
}

dvec3 operator*(const dmat4 &m,const dvec3 &v) {
	dvec3 ret;
	return mul(ret,m,v);
}

dvec3 operator*(const dvec3 &v,const dmat4 &m) {
	dvec3 ret;
	return mul(ret,v,m);
}

dvec4 operator*(const dmat4 &m,const dvec4 &v) {
	dvec4 ret;
	return mul(ret,m,v);
}

dvec4 operator*(const dvec4 &v,const dmat4 &m) {
	dvec4 ret;
	return mul(ret,v,m);
}

dmat4 operator*(const dmat4 &m0,const dmat4 &m1) {
	dmat4 ret;
	return mul(ret,m0,m1);
}

dmat4 operator+(const dmat4 &m0,const dmat4 &m1) {
	dmat4 ret;
	return add(ret,m0,m1);
}

dmat4 operator-(const dmat4 &m0,const dmat4 &m1) {
	dmat4 ret;
	return sub(ret,m0,m1);
}

/*
 */
int compare(const dmat4 &m0,const dmat4 &m1) {
	return (compare(m0.m00,m1.m00) && compare(m0.m10,m1.m10) && compare(m0.m20,m1.m20) &&
		compare(m0.m01,m1.m01) && compare(m0.m11,m1.m11) && compare(m0.m21,m1.m21) &&
		compare(m0.m02,m1.m02) && compare(m0.m12,m1.m12) && compare(m0.m22,m1.m22) &&
		compare(m0.m03,m1.m03) && compare(m0.m13,m1.m13) && compare(m0.m23,m1.m23));
}

int compare(const dmat4 &m0,const dmat4 &m1,double epsilon) {
	return (compare(m0.m00,m1.m00,epsilon) && compare(m0.m10,m1.m10,epsilon) && compare(m0.m20,m1.m20,epsilon) &&
		compare(m0.m01,m1.m01,epsilon) && compare(m0.m11,m1.m11,epsilon) && compare(m0.m21,m1.m21,epsilon) &&
		compare(m0.m02,m1.m02,epsilon) && compare(m0.m12,m1.m12,epsilon) && compare(m0.m22,m1.m22,epsilon) &&
		compare(m0.m03,m1.m03,epsilon) && compare(m0.m13,m1.m13,epsilon) && compare(m0.m23,m1.m23,epsilon));
}

double determinant(const dmat4 &m) {
	double det = 0.0;
	det =  m.m00 * (m.m11 * m.m22 - m.m12 * m.m21);
	det -= m.m01 * (m.m10 * m.m22 - m.m12 * m.m20);
	det += m.m02 * (m.m10 * m.m21 - m.m11 * m.m20);
	return det;
}

dmat4 &mul(dmat4 &ret,const dmat4 &m,double v) {
	#ifdef USE_SSE2
		__m128d temp = _mm_set1_pd(v);
		ret.col0 = _mm_mul_pd(m.col0,temp); ret.col1 = _mm_mul_pd(m.col1,temp);
		ret.col2 = _mm_mul_pd(m.col2,temp); ret.col3 = _mm_mul_pd(m.col3,temp);
		ret.col4 = _mm_mul_pd(m.col4,temp); ret.col5 = _mm_mul_pd(m.col5,temp);
	#else
		ret.m00 = m.m00 * v; ret.m01 = m.m01 * v; ret.m02 = m.m02 * v; ret.m03 = m.m03 * v;
		ret.m10 = m.m10 * v; ret.m11 = m.m11 * v; ret.m12 = m.m12 * v; ret.m13 = m.m13 * v;
		ret.m20 = m.m20 * v; ret.m21 = m.m21 * v; ret.m22 = m.m22 * v; ret.m23 = m.m23 * v;
	#endif
	return ret;
}

vec2 &mul(vec2 &ret,const dmat4 &m,const vec2 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m03);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m13);
	return ret;
}

vec2 &mul(vec2 &ret,const vec2 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y);
	return ret;
}

vec3 &mul(vec3 &ret,const dmat4 &m,const vec3 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23);
	return ret;
}

vec3 &mul(vec3 &ret,const vec3 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul(vec4 &ret,const dmat4 &m,const vec4 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w);
	ret.w = v.w;
	return ret;
}

vec4 &mul(vec4 &ret,const vec4 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	ret.w = (float)(m.m03 * v.x + m.m13 * v.y + m.m23 * v.z + v.w);
	return ret;
}

vec2 &mul(vec2 &ret,const dmat4 &m,const dvec2 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m03);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m13);
	return ret;
}

vec2 &mul(vec2 &ret,const dvec2 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y);
	return ret;
}

vec3 &mul(vec3 &ret,const dmat4 &m,const dvec3 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23);
	return ret;
}

vec3 &mul(vec3 &ret,const dvec3 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul(vec4 &ret,const dmat4 &m,const dvec4 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w);
	ret.w = (float)v.w;
	return ret;
}

vec4 &mul(vec4 &ret,const dvec4 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	ret.w = (float)(m.m03 * v.x + m.m13 * v.y + m.m23 * v.z + v.w);
	return ret;
}

dvec2 &mul(dvec2 &ret,const dmat4 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m03;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m13;
	return ret;
}

dvec2 &mul(dvec2 &ret,const dvec2 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

dvec3 &mul(dvec3 &ret,const dmat4 &m,const dvec3 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
	return ret;
}

dvec3 &mul(dvec3 &ret,const dvec3 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dvec4 &mul(dvec4 &ret,const dmat4 &m,const dvec4 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03 * v.w;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13 * v.w;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23 * v.w;
	ret.w = v.w;
	return ret;
}

dvec4 &mul(dvec4 &ret,const dvec4 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	ret.w = m.m03 * v.x + m.m13 * v.y + m.m23 * v.z + v.w;
	return ret;
}

vec2 &mul3(vec2 &ret,const dmat4 &m,const vec2 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y);
	return ret;
}

vec2 &mul3(vec2 &ret,const vec2 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y);
	return ret;
}

vec3 &mul3(vec3 &ret,const dmat4 &m,const vec3 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z);
	return ret;
}

vec3 &mul3(vec3 &ret,const vec3 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul3(vec4 &ret,const dmat4 &m,const vec4 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul3(vec4 &ret,const vec4 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

vec2 &mul3(vec2 &ret,const dmat4 &m,const dvec2 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y);
	return ret;
}

vec2 &mul3(vec2 &ret,const dvec2 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y);
	return ret;
}

vec3 &mul3(vec3 &ret,const dmat4 &m,const dvec3 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z);
	return ret;
}

vec3 &mul3(vec3 &ret,const dvec3 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul3(vec4 &ret,const dmat4 &m,const dvec4 &v) {
	ret.x = (float)(m.m00 * v.x + m.m01 * v.y + m.m02 * v.z);
	ret.y = (float)(m.m10 * v.x + m.m11 * v.y + m.m12 * v.z);
	ret.z = (float)(m.m20 * v.x + m.m21 * v.y + m.m22 * v.z);
	return ret;
}

vec4 &mul3(vec4 &ret,const dvec4 &v,const dmat4 &m) {
	ret.x = (float)(m.m00 * v.x + m.m10 * v.y + m.m20 * v.z);
	ret.y = (float)(m.m01 * v.x + m.m11 * v.y + m.m21 * v.z);
	ret.z = (float)(m.m02 * v.x + m.m12 * v.y + m.m22 * v.z);
	return ret;
}

dvec2 &mul3(dvec2 &ret,const dmat4 &m,const dvec2 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y;
	ret.y = m.m10 * v.x + m.m11 * v.y;
	return ret;
}

dvec2 &mul3(dvec2 &ret,const dvec2 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y;
	ret.y = m.m01 * v.x + m.m11 * v.y;
	return ret;
}

dvec3 &mul3(dvec3 &ret,const dmat4 &m,const dvec3 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	return ret;
}

dvec3 &mul3(dvec3 &ret,const dvec3 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dvec4 &mul3(dvec4 &ret,const dmat4 &m,const dvec4 &v) {
	ret.x = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z;
	ret.y = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z;
	ret.z = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z;
	return ret;
}

dvec4 &mul3(dvec4 &ret,const dvec4 &v,const dmat4 &m) {
	ret.x = m.m00 * v.x + m.m10 * v.y + m.m20 * v.z;
	ret.y = m.m01 * v.x + m.m11 * v.y + m.m21 * v.z;
	ret.z = m.m02 * v.x + m.m12 * v.y + m.m22 * v.z;
	return ret;
}

dmat4 &mul(dmat4 &ret,const dmat4 &m0,const dmat4 &m1) {
	#ifdef USE_SSE2
		__m128d col_01 = _mm_shuffle_pd(m0.col0,m0.col1,_MM_PERM22(Y,X));
		__m128d col_12 = _mm_shuffle_pd(m0.col1,m0.col2,_MM_PERM22(Y,X));
		__m128d col_34 = _mm_shuffle_pd(m0.col3,m0.col4,_MM_PERM22(Y,X));
		__m128d res_0 = _mm_mul_pd(m0.col0,_MM_SWIZZLE2(m1.col0,X,X));
		__m128d res_1 = _mm_mul_pd(col_12,_MM_SWIZZLE2(m1.col0,Y,Y));
		__m128d res_2 = _mm_mul_pd(m0.col3,_MM_SWIZZLE2(m1.col1,X,X));
		ret.col0 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		res_0 = _mm_mul_pd(col_01,_MM_SWIZZLE2(m1.col1,Y,Y));
		res_1 = _mm_mul_pd(m0.col2,_MM_SWIZZLE2(m1.col2,X,X));
		res_2 = _mm_mul_pd(col_34,_MM_SWIZZLE2(m1.col2,Y,Y));
		ret.col2 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		res_0 = _mm_mul_pd(m0.col0,_MM_SWIZZLE2(m1.col3,X,X));
		res_1 = _mm_mul_pd(col_12,_MM_SWIZZLE2(m1.col3,Y,Y));
		res_2 = _mm_mul_pd(m0.col3,_MM_SWIZZLE2(m1.col4,X,X));
		ret.col3 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		res_0 = _mm_mul_pd(col_01,_MM_SWIZZLE2(m1.col4,Y,Y));
		res_1 = _mm_mul_pd(m0.col2,_MM_SWIZZLE2(m1.col5,X,X));
		res_2 = _mm_mul_pd(col_34,_MM_SWIZZLE2(m1.col5,Y,Y));
		ret.col5 = _mm_add_pd(_mm_add_pd(res_0,res_1),_mm_add_pd(res_2,m0.col5));
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
		ret.m03 = m0.m00 * m1.m03 + m0.m01 * m1.m13 + m0.m02 * m1.m23 + m0.m03;
	#else
		ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20;
		ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20;
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21;
		ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21;
		ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22;
		ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
		ret.m03 = m0.m00 * m1.m03 + m0.m01 * m1.m13 + m0.m02 * m1.m23 + m0.m03;
		ret.m13 = m0.m10 * m1.m03 + m0.m11 * m1.m13 + m0.m12 * m1.m23 + m0.m13;
		ret.m23 = m0.m20 * m1.m03 + m0.m21 * m1.m13 + m0.m22 * m1.m23 + m0.m23;
	#endif
	return ret;
}

dmat4 &mul4(dmat4 &ret,const dmat4 &m0,const dmat4 &m1) {
	return mul(ret,m0,m1);
}

dmat4 &mul3(dmat4 &ret,const dmat4 &m0,const dmat4 &m1) {
	#ifdef USE_SSE2
		__m128d col_01 = _mm_shuffle_pd(m0.col0,m0.col1,_MM_PERM22(Y,X));
		__m128d col_12 = _mm_shuffle_pd(m0.col1,m0.col2,_MM_PERM22(Y,X));
		__m128d col_34 = _mm_shuffle_pd(m0.col3,m0.col4,_MM_PERM22(Y,X));
		__m128d res_0 = _mm_mul_pd(m0.col0,_MM_SWIZZLE2(m1.col0,X,X));
		__m128d res_1 = _mm_mul_pd(col_12,_MM_SWIZZLE2(m1.col0,Y,Y));
		__m128d res_2 = _mm_mul_pd(m0.col3,_MM_SWIZZLE2(m1.col1,X,X));
		ret.col0 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		res_0 = _mm_mul_pd(col_01,_MM_SWIZZLE2(m1.col1,Y,Y));
		res_1 = _mm_mul_pd(m0.col2,_MM_SWIZZLE2(m1.col2,X,X));
		res_2 = _mm_mul_pd(col_34,_MM_SWIZZLE2(m1.col2,Y,Y));
		ret.col2 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		res_0 = _mm_mul_pd(m0.col0,_MM_SWIZZLE2(m1.col3,X,X));
		res_1 = _mm_mul_pd(col_12,_MM_SWIZZLE2(m1.col3,Y,Y));
		res_2 = _mm_mul_pd(m0.col3,_MM_SWIZZLE2(m1.col4,X,X));
		ret.col3 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
	#else
		ret.m00 = m0.m00 * m1.m00 + m0.m01 * m1.m10 + m0.m02 * m1.m20;
		ret.m10 = m0.m10 * m1.m00 + m0.m11 * m1.m10 + m0.m12 * m1.m20;
		ret.m20 = m0.m20 * m1.m00 + m0.m21 * m1.m10 + m0.m22 * m1.m20;
		ret.m01 = m0.m00 * m1.m01 + m0.m01 * m1.m11 + m0.m02 * m1.m21;
		ret.m11 = m0.m10 * m1.m01 + m0.m11 * m1.m11 + m0.m12 * m1.m21;
		ret.m21 = m0.m20 * m1.m01 + m0.m21 * m1.m11 + m0.m22 * m1.m21;
		ret.m02 = m0.m00 * m1.m02 + m0.m01 * m1.m12 + m0.m02 * m1.m22;
		ret.m12 = m0.m10 * m1.m02 + m0.m11 * m1.m12 + m0.m12 * m1.m22;
		ret.m22 = m0.m20 * m1.m02 + m0.m21 * m1.m12 + m0.m22 * m1.m22;
	#endif
	return ret;
}

dmat4 &mult(dmat4 &ret,const dmat4 &m,const dvec3 &v) {
	ret.m00 = m.m00; ret.m01 = m.m01; ret.m02 = m.m02;
	ret.m10 = m.m10; ret.m11 = m.m11; ret.m12 = m.m12;
	ret.m20 = m.m20; ret.m21 = m.m21; ret.m22 = m.m22;
	ret.m03 = m.m00 * v.x + m.m01 * v.y + m.m02 * v.z + m.m03;
	ret.m13 = m.m10 * v.x + m.m11 * v.y + m.m12 * v.z + m.m13;
	ret.m23 = m.m20 * v.x + m.m21 * v.y + m.m22 * v.z + m.m23;
	return ret;
}

dmat4 &add(dmat4 &ret,const dmat4 &m0,const dmat4 &m1) {
	#ifdef USE_SSE2
		ret.col0 = _mm_add_pd(m0.col0,m1.col0); ret.col1 = _mm_add_pd(m0.col1,m1.col1);
		ret.col2 = _mm_add_pd(m0.col2,m1.col2); ret.col3 = _mm_add_pd(m0.col3,m1.col3);
		ret.col4 = _mm_add_pd(m0.col4,m1.col4); ret.col5 = _mm_add_pd(m0.col5,m1.col5);
	#else
		ret.m00 = m0.m00 + m1.m00; ret.m01 = m0.m01 + m1.m01; ret.m02 = m0.m02 + m1.m02; ret.m03 = m0.m03 + m1.m03;
		ret.m10 = m0.m10 + m1.m10; ret.m11 = m0.m11 + m1.m11; ret.m12 = m0.m12 + m1.m12; ret.m13 = m0.m13 + m1.m13;
		ret.m20 = m0.m20 + m1.m20; ret.m21 = m0.m21 + m1.m21; ret.m22 = m0.m22 + m1.m22; ret.m23 = m0.m23 + m1.m23;
	#endif
	return ret;
}

dmat4 &sub(dmat4 &ret,const dmat4 &m0,const dmat4 &m1) {
	#ifdef USE_SSE2
		ret.col0 = _mm_sub_pd(m0.col0,m1.col0); ret.col1 = _mm_sub_pd(m0.col1,m1.col1);
		ret.col2 = _mm_sub_pd(m0.col2,m1.col2); ret.col3 = _mm_sub_pd(m0.col3,m1.col3);
		ret.col4 = _mm_sub_pd(m0.col4,m1.col4); ret.col5 = _mm_sub_pd(m0.col5,m1.col5);
	#else
		ret.m00 = m0.m00 - m1.m00; ret.m01 = m0.m01 - m1.m01; ret.m02 = m0.m02 - m1.m02; ret.m03 = m0.m03 - m1.m03;
		ret.m10 = m0.m10 - m1.m10; ret.m11 = m0.m11 - m1.m11; ret.m12 = m0.m12 - m1.m12; ret.m13 = m0.m13 - m1.m13;
		ret.m20 = m0.m20 - m1.m20; ret.m21 = m0.m21 - m1.m21; ret.m22 = m0.m22 - m1.m22; ret.m23 = m0.m23 - m1.m23;
	#endif
	return ret;
}

dmat4 &orthonormalize(dmat4 &ret,const dmat4 &m) {
	dvec3 x = dvec3(m.m00,m.m10,m.m20);
	dvec3 y = dvec3(m.m01,m.m11,m.m21);
	dvec3 z = cross(x,y);
	cross(y,z,x);
	x.normalize();
	y.normalize();
	z.normalize();
	ret.m00 = x.x; ret.m01 = y.x; ret.m02 = z.x; ret.m03 = m.m03;
	ret.m10 = x.y; ret.m11 = y.y; ret.m12 = z.y; ret.m13 = m.m13;
	ret.m20 = x.z; ret.m21 = y.z; ret.m22 = z.z; ret.m23 = m.m23;
	return ret;
}

dmat4 &rotation(dmat4 &ret,const dmat4 &m) {
	ret.m00 = m.m00; ret.m01 = m.m01; ret.m02 = m.m02; ret.m03 = 0.0;
	ret.m10 = m.m10; ret.m11 = m.m11; ret.m12 = m.m12; ret.m13 = 0.0;
	ret.m20 = m.m20; ret.m21 = m.m21; ret.m22 = m.m22; ret.m23 = 0.0;
	return ret;
}

dmat4 &inverse(dmat4 &ret,const dmat4 &m) {
	double idet = Math::rcp(determinant(m));
	ret.m00 =  (m.m11 * m.m22 - m.m12 * m.m21) * idet;
	ret.m10 = -(m.m10 * m.m22 - m.m12 * m.m20) * idet;
	ret.m20 =  (m.m10 * m.m21 - m.m11 * m.m20) * idet;
	ret.m01 = -(m.m01 * m.m22 - m.m02 * m.m21) * idet;
	ret.m11 =  (m.m00 * m.m22 - m.m02 * m.m20) * idet;
	ret.m21 = -(m.m00 * m.m21 - m.m01 * m.m20) * idet;
	ret.m02 =  (m.m01 * m.m12 - m.m02 * m.m11) * idet;
	ret.m12 = -(m.m00 * m.m12 - m.m02 * m.m10) * idet;
	ret.m22 =  (m.m00 * m.m11 - m.m01 * m.m10) * idet;
	ret.m03 = -(ret.m00 * m.m03 + ret.m01 * m.m13 + ret.m02 * m.m23);
	ret.m13 = -(ret.m10 * m.m03 + ret.m11 * m.m13 + ret.m12 * m.m23);
	ret.m23 = -(ret.m20 * m.m03 + ret.m21 * m.m13 + ret.m22 * m.m23);
	return ret;
}

dmat4 &lerp(dmat4 &ret,const dmat4 &m0,const dmat4 &m1,double k) {
	dvec3 positions[3];
	quat rotations[3];
	vec3 scales[3];
	decomposeTransform(m0,positions[0],rotations[0],scales[0]);
	decomposeTransform(m1,positions[1],rotations[1],scales[1]);
	lerp(positions[2],positions[0],positions[1],k);
	slerp(rotations[2],rotations[0],rotations[1],(float)k);
	lerp(scales[2],scales[0],scales[1],(float)k);
	return composeTransform(ret,positions[2],rotations[2],scales[2]);
}

/*
 */
dmat4 orthonormalize(const dmat4 &m) {
	dmat4 ret;
	return orthonormalize(ret,m);
}

dmat4 rotation(const dmat4 &m) {
	dmat4 ret;
	return rotation(ret,m);
}

dmat4 inverse(const dmat4 &m) {
	dmat4 ret;
	return inverse(ret,m);
}

dmat4 lerp(const dmat4 &m0,const dmat4 &m1,double k) {
	dmat4 ret;
	return lerp(ret,m0,m1,k);
}

/*
 */
dmat4 translate(const dvec3 &v) {
	dmat4 ret;
	ret.setTranslate(v);
	return ret;
}

dmat4 translate(double x,double y,double z) {
	dmat4 ret;
	ret.setTranslate(dvec3(x,y,z));
	return ret;
}

/*
 */
dmat4 reflect(const dvec4 &plane) {
	dmat4 ret;
	double x = plane.x;
	double y = plane.y;
	double z = plane.z;
	double x2 = x * 2.0;
	double y2 = y * 2.0;
	double z2 = z * 2.0;
	ret.m00 = 1.0 - x * x2; ret.m01 = -y * x2;      ret.m02 = -z * x2;      ret.m03 = -plane.w * x2;
	ret.m10 = -x * y2;      ret.m11 = 1.0 - y * y2; ret.m12 = -z * y2;      ret.m13 = -plane.w * y2;
	ret.m20 = -x * z2;      ret.m21 = -y * z2;      ret.m22 = 1.0 - z * z2; ret.m23 = -plane.w * z2;
	return ret;
}

dmat4 setTo(const dvec3 &position,const dvec3 &direction,const vec3 &up) {
	dmat4 ret;
	vec3 z = normalize(vec3(position - direction));
	vec3 x = normalize(cross(up,z));
	vec3 y = normalize(cross(z,x));
	ret.m00 = x.x; ret.m01 = y.x; ret.m02 = z.x; ret.m03 = position.x;
	ret.m10 = x.y; ret.m11 = y.y; ret.m12 = z.y; ret.m13 = position.y;
	ret.m20 = x.z; ret.m21 = y.z; ret.m22 = z.z; ret.m23 = position.z;
	return ret;
}

dmat4 lookAt(const dvec3 &position,const dvec3 &direction,const vec3 &up) {
	dmat4 ret,m0,m1;
	vec3 z = normalize(vec3(position - direction));
	vec3 x = normalize(cross(up,z));
	vec3 y = normalize(cross(z,x));
	m0.m00 = x.x; m0.m01 = x.y; m0.m02 = x.z; m0.m03 = 0.0;
	m0.m10 = y.x; m0.m11 = y.y; m0.m12 = y.z; m0.m13 = 0.0;
	m0.m20 = z.x; m0.m21 = z.y; m0.m22 = z.z; m0.m23 = 0.0;
	m1.setTranslate(-position);
	return mul(ret,m0,m1);
}

/*
 */
void decomposeTransform(const dmat4 &m,dvec3 &position,quat &rot,vec3 &s) {
	mat3 rotate,scale;
	mat3 rotation = mat3(m);
	orthonormalize(rotate,rotation);
	mul(scale,transpose(rotate),rotation);
	position.x = m.m03;
	position.y = m.m13;
	position.z = m.m23;
	rot = rotate.getQuat();
	s.x = scale.m00;
	s.y = scale.m11;
	s.z = scale.m22;
}

dmat4 &composeTransform(dmat4 &ret,const dvec3 &position,const quat &rot,const vec3 &s) {
	mat3 rotation,scale;
	scale.setDiagonal(s);
	ret.set(mul(rotation,rot.getMat3(),scale),position);
	return ret;
}

/******************************************************************************\
*
* quat
*
\******************************************************************************/

/*
 */
const quat quat_identity(0.0f,0.0f,0.0f,1.0f);

/*
 */
quat::quat(const vec3 &axis,float angle) {
	set(axis,angle);
}

quat::quat(float x,float y,float z,float angle) {
	set(x,y,z,angle);
}

quat::quat(float angle_x,float angle_y,float angle_z) {
	set(angle_x,angle_y,angle_z);
}

quat::quat(const mat3 &m) {
	quat q = m.getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

quat::quat(const mat4 &m) {
	quat q = mat3(m).getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

quat::quat(const dmat4 &m) {
	quat q = mat3(m).getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

/*
 */
quat quat::operator-() const {
	quat ret;
	ret.x = -x;
	ret.y = -y;
	ret.z = -z;
	ret.w = -w;
	return ret;
}

quat &quat::operator*=(float v) {
	return mul(*this,*this,v);
}

quat &quat::operator*=(const quat &q) {
	return mul(*this,quat(*this),q);
}

quat &quat::operator+=(const quat &q) {
	return add(*this,*this,q);
}

quat &quat::operator-=(const quat &q) {
	return sub(*this,*this,q);
}

/*
 */
void quat::set(const vec3 &v) {
	x = v.x;
	y = v.y;
	z = v.z;
	w = 0.0f;
}

void quat::set(const mat3 &m) {
	quat q = m.getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

void quat::set(const mat4 &m) {
	quat q = mat3(m).getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

void quat::set(const dmat4 &m) {
	quat q = mat3(m).getQuat();
	x = q.x;
	y = q.y;
	z = q.z;
	w = q.w;
}

void quat::set(const vec3 &axis,float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD * 0.5f,s,c);
	float ilength = Math::rsqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
	x = axis.x * ilength * s;
	y = axis.y * ilength * s;
	z = axis.z * ilength * s;
	w = c;
}

void quat::set(float axis_x,float axis_y,float axis_z,float angle) {
	float s,c;
	Math::sincos(angle * DEG2RAD * 0.5f,s,c);
	float ilength = Math::rsqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
	x = axis_x * ilength * s;
	y = axis_y * ilength * s;
	z = axis_z * ilength * s;
	w = c;
}

void quat::set(float angle_x,float angle_y,float angle_z) {
	float sx,cx,sy,cy,sz,cz;
	Math::sincos(angle_x * DEG2RAD * 0.5f,sx,cx);
	Math::sincos(angle_y * DEG2RAD * 0.5f,sy,cy);
	Math::sincos(angle_z * DEG2RAD * 0.5f,sz,cz);
	float sxcy = sx * cy;
	float cxcy = cx * cy;
	float sxsy = sx * sy;
	float cxsy = cx * sy;
	x = cxsy * sz + sxcy * cz;
	y = cxsy * cz - sxcy * sz;
	z = sxsy * cz + cxcy * sz;
	w = cxcy * cz - sxsy * sz;
}

void quat::get(vec3 &axis,float &angle) const {
	float ilength = Math::rsqrt(x * x + y * y + z * z);
	axis.x = x * ilength;
	axis.y = y * ilength;
	axis.z = z * ilength;
	angle = Math::acos(clamp(w,-1.0f,1.0f)) * RAD2DEG * 2.0f;
	if(angle > 180.0f) angle -= 360.0f;
}

/*
 */
mat3 quat::getMat3() const {
	mat3 ret;
	float x2 = x + x;
	float y2 = y + y;
	float z2 = z + z;
	float xx2 = x * x2;
	float yy2 = y * y2;
	float zz2 = z * z2;
	float zx2 = z * x2;
	float xy2 = x * y2;
	float yz2 = y * z2;
	float wx2 = w * x2;
	float wy2 = w * y2;
	float wz2 = w * z2;
	ret.m00 = 1.0f - yy2 - zz2; ret.m01 = xy2 - wz2;        ret.m02 = zx2 + wy2;
	ret.m10 = xy2 + wz2;        ret.m11 = 1.0f - xx2 - zz2; ret.m12 = yz2 - wx2;
	ret.m20 = zx2 - wy2;        ret.m21 = yz2 + wx2;        ret.m22 = 1.0f - xx2 - yy2;
	return ret;
}

float quat::getAngle(const vec3 &axis) const {
	vec3 v0,v1;
	mul(v0,*this,axis);
	add(v1,v0,axis);
	v1.normalize();
	float angle = dot(axis,v1);
	if(Math::abs(angle) > EPSILON) {
		quat q0,q1;
		q0.set(-cross(v0,axis,v1));
		q0.w = angle;
		mul(q1,q0,*this);
		q1.normalize();
		q1.get(v0,angle);
	} else {
		get(v0,angle);
	}
	if(dot(v0,axis) < 0.0f) angle = -angle;
	return angle;
}

/*
 */
int operator==(const quat &q0,const quat &q1) {
	return compare(q0,q1);
}

int operator!=(const quat &q0,const quat &q1) {
	return !compare(q0,q1);
}

quat operator*(const quat &q,float v) {
	quat ret;
	return mul(ret,q,v);
}

vec3 operator*(const quat &q,const vec3 &v) {
	vec3 ret;
	return mul(ret,q,v);
}

vec3 operator*(const vec3 &v,const quat &q) {
	vec3 ret;
	return mul(ret,v,q);
}

dvec3 operator*(const quat &q,const dvec3 &v) {
	dvec3 ret;
	return mul(ret,q,v);
}

dvec3 operator*(const dvec3 &v,const quat &q) {
	dvec3 ret;
	return mul(ret,v,q);
}

quat operator*(const quat &q0,const quat &q1) {
	quat ret;
	return mul(ret,q0,q1);
}

quat operator+(const quat &q0,const quat &q1) {
	quat ret;
	return add(ret,q0,q1);
}

quat operator-(const quat &q0,const quat &q1) {
	quat ret;
	return sub(ret,q0,q1);
}

/*
 */
vec3 &mul(vec3 &ret,const quat &q,const vec3 &v) {
	float x2 = q.x + q.x;
	float y2 = q.y + q.y;
	float z2 = q.z + q.z;
	float xx2 = q.x * x2;
	float yy2 = q.y * y2;
	float zz2 = q.z * z2;
	float zx2 = q.z * x2;
	float xy2 = q.x * y2;
	float yz2 = q.y * z2;
	float wx2 = q.w * x2;
	float wy2 = q.w * y2;
	float wz2 = q.w * z2;
	ret.x = (1.0f - yy2 - zz2) * v.x + (xy2 - wz2) * v.y + (zx2 + wy2) * v.z;
	ret.y = (xy2 + wz2) * v.x + (1.0f - xx2 - zz2) * v.y + (yz2 - wx2) * v.z;
	ret.z = (zx2 - wy2) * v.x + (yz2 + wx2) * v.y + (1.0f - xx2 - yy2) * v.z;
	return ret;
}

vec3 &mul(vec3 &ret,const vec3 &v,const quat &q) {
	float x2 = q.x + q.x;
	float y2 = q.y + q.y;
	float z2 = q.z + q.z;
	float xx2 = q.x * x2;
	float yy2 = q.y * y2;
	float zz2 = q.z * z2;
	float zx2 = q.z * x2;
	float xy2 = q.x * y2;
	float yz2 = q.y * z2;
	float wx2 = q.w * x2;
	float wy2 = q.w * y2;
	float wz2 = q.w * z2;
	ret.x = (1.0f - yy2 - zz2) * v.x + (xy2 + wz2) * v.y + (zx2 - wy2) * v.z;
	ret.y = (xy2 - wz2) * v.x + (1.0f - xx2 - zz2) * v.y + (yz2 + wx2) * v.z;
	ret.z = (zx2 + wy2) * v.x + (yz2 - wx2) * v.y + (1.0f - xx2 - yy2) * v.z;
	return ret;
}

dvec3 &mul(dvec3 &ret,const quat &q,const dvec3 &v) {
	double x2 = q.x + q.x;
	double y2 = q.y + q.y;
	double z2 = q.z + q.z;
	double xx2 = q.x * x2;
	double yy2 = q.y * y2;
	double zz2 = q.z * z2;
	double zx2 = q.z * x2;
	double xy2 = q.x * y2;
	double yz2 = q.y * z2;
	double wx2 = q.w * x2;
	double wy2 = q.w * y2;
	double wz2 = q.w * z2;
	ret.x = (1.0 - yy2 - zz2) * v.x + (xy2 - wz2) * v.y + (zx2 + wy2) * v.z;
	ret.y = (xy2 + wz2) * v.x + (1.0 - xx2 - zz2) * v.y + (yz2 - wx2) * v.z;
	ret.z = (zx2 - wy2) * v.x + (yz2 + wx2) * v.y + (1.0 - xx2 - yy2) * v.z;
	return ret;
}

dvec3 &mul(dvec3 &ret,const dvec3 &v,const quat &q) {
	double x2 = q.x + q.x;
	double y2 = q.y + q.y;
	double z2 = q.z + q.z;
	double xx2 = q.x * x2;
	double yy2 = q.y * y2;
	double zz2 = q.z * z2;
	double zx2 = q.z * x2;
	double xy2 = q.x * y2;
	double yz2 = q.y * z2;
	double wx2 = q.w * x2;
	double wy2 = q.w * y2;
	double wz2 = q.w * z2;
	ret.x = (1.0 - yy2 - zz2) * v.x + (xy2 + wz2) * v.y + (zx2 - wy2) * v.z;
	ret.y = (xy2 - wz2) * v.x + (1.0 - xx2 - zz2) * v.y + (yz2 + wx2) * v.z;
	ret.z = (zx2 + wy2) * v.x + (yz2 - wx2) * v.y + (1.0 - xx2 - yy2) * v.z;
	return ret;
}

quat &mul(quat &ret,const quat &q0,const quat &q1) {
	ret.x = q0.w * q1.x + q0.x * q1.w + q0.y * q1.z - q0.z * q1.y;
	ret.y = q0.w * q1.y + q0.y * q1.w + q0.z * q1.x - q0.x * q1.z;
	ret.z = q0.w * q1.z + q0.z * q1.w + q0.x * q1.y - q0.y * q1.x;
	ret.w = q0.w * q1.w - q0.x * q1.x - q0.y * q1.y - q0.z * q1.z;
	return ret;
}

quat &slerp(quat &ret,const quat &q0,const quat &q1,float k) {
	if(k <= 0.0f) {
		ret = q0;
	} else if(k >= 1.0f) {
		ret = q1;
	} else {
		float k0,k1;
		float c = dot(q0,q1);
		float ac = Math::abs(c);
		if(ac < 1.0f - EPSILON) {
			float angle = Math::acos(ac);
			float is = Math::rcp(Math::sin(angle));
			k0 = Math::sin(angle * (1.0f - k)) * is;
			k1 = Math::sin(angle * k) * is;
		} else {
			k0 = 1.0f - k;
			k1 = k;
		}
		if(c < 0.0f) k1 = -k1;
		ret.x = q0.x * k0 + q1.x * k1;
		ret.y = q0.y * k0 + q1.y * k1;
		ret.z = q0.z * k0 + q1.z * k1;
		ret.w = q0.w * k0 + q1.w * k1;
	}
	return ret;
}

/*
 */
quat normalize(const quat &q) {
	quat ret = q;
	return ret.normalize();
}

quat inverse(const quat &q) {
	quat ret;
	return inverse(ret,q);
}

quat slerp(const quat &q0,const quat &q1,float k) {
	quat ret;
	return slerp(ret,q0,q1,k);
}

/*
 */
void decomposeTransform(const mat4 &m,quat &q0,quat &q1) {
	mat3 rotate;
	q0 = orthonormalize(rotate,mat3(m)).getQuat();
	q1.x = ( m.m03 * q0.w + m.m13 * q0.z - m.m23 * q0.y) * 0.5f;
	q1.y = (-m.m03 * q0.z + m.m13 * q0.w + m.m23 * q0.x) * 0.5f;
	q1.z = ( m.m03 * q0.y - m.m13 * q0.x + m.m23 * q0.w) * 0.5f;
	q1.w = ( m.m03 * q0.x + m.m13 * q0.y + m.m23 * q0.z) * 0.5f;
	mul(q0,q0,Math::sqrt(m.m00 * m.m00 + m.m01 * m.m01 + m.m02 * m.m02));
}

mat4 &composeTransform(mat4 &ret,const quat &q0,const quat &q1) {
	float ilength = Math::rcp(dot(q0,q1)) * 2.0f;
	float nx = q0.x * ilength;
	float ny = q0.y * ilength;
	float nz = q0.z * ilength;
	float nw = q0.w * ilength;
	float xx = q0.x * nx;
	float yy = q0.y * ny;
	float zz = q0.z * nz;
	float ww = q0.w * nw;
	float xy = q0.x * ny;
	float xz = q0.x * nz;
	float xw = q0.x * nw;
	float yz = q0.y * nz;
	float yw = q0.y * nw;
	float zw = q0.z * nw;
	ret.m00 = (ww + xx - yy - zz) * 0.5f;
	ret.m10 = xy + zw;
	ret.m20 = xz - yw;
	ret.m30 = 0.0f;
	ret.m01 = xy - zw;
	ret.m11 = (ww + yy - xx - zz) * 0.5f;
	ret.m21 = yz + xw;
	ret.m31 = 0.0f;
	ret.m02 = xz + yw;
	ret.m12 = yz - xw;
	ret.m22 = (ww + zz - xx - yy) * 0.5f;
	ret.m32 = 0.0f;
	ret.m03 = q1.w * nx + q1.x * nw + q1.z * ny - q1.y * nz;
	ret.m13 = q1.w * ny + q1.x * nz + q1.y * nw - q1.z * nx;
	ret.m23 = q1.w * nz + q1.y * nx + q1.z * nw - q1.x * ny;
	ret.m33 = 1.0f;
	return ret;
}
