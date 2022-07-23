/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    MathLib.h
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

#ifndef __MATH_LIB_H__
#define __MATH_LIB_H__

#ifdef USE_SSE2
	#include <emmintrin.h>
#elif USE_SSE
	#include <xmmintrin.h>
#elif USE_ALTIVEC
	#include <altivec.h>
	#include <vec_types.h>
#elif USE_NEON
	#include <arm_neon.h>
#endif

#ifdef __PPU__
	#include <ppu_intrinsics.h>
#endif

#include "Base.h"
#include "Half.h"

/*
 */
#ifdef PI
	#undef PI
#endif

#ifdef DEG2RAD
	#undef DEG2RAD
#endif

#ifdef RAD2DEG
	#undef RAD2DEG
#endif

#ifdef EPSILON
	#undef EPSILON
#endif

#ifdef INFINITY
	#undef INFINITY
#endif

#ifdef min
	#undef min
#endif

#ifdef max
	#undef max
#endif

/*
 */
#define PI			3.141592654f
#define PI2			6.283185308f
#define PI05		1.570796327f
#define LOG2		0.693147181f
#define LOG10		2.302585093f
#define SQRT2		1.414213562f
#define DEG2RAD		(PI / 180.0f)
#define RAD2DEG		(180.0f / PI)
#define EPSILON		1e-6f
#define INFINITY	1e+9f

/*
 */
#ifdef USE_SSE
	#define _MM_PERM2_X		0
	#define _MM_PERM2_Y		1
	#define _MM_PERM2_Z		2
	#define _MM_PERM2_W		3
	#define _MM_PERM2(X,Y,Z,W) _MM_SHUFFLE(_MM_PERM2_ ## W,_MM_PERM2_ ## Z,_MM_PERM2_ ## Y,_MM_PERM2_ ## X)
	#define _MM_SWIZZLE(V,X,Y,Z,W) _mm_shuffle_ps(V,V,_MM_PERM2(X,Y,Z,W))
	#ifdef USE_SSE2
		#define _MM_PERM22(X,Y) _MM_SHUFFLE2(_MM_PERM2_ ## Y,_MM_PERM2_ ## X)
		#define _MM_SWIZZLE2(V,X,Y) _mm_shuffle_pd(V,V,_MM_PERM22(X,Y))
	#endif
#elif USE_ALTIVEC
	#define VEC_PERM4_LX 	0x00010203
	#define VEC_PERM4_LY 	0x04050607
	#define VEC_PERM4_LZ 	0x08090a0b
	#define VEC_PERM4_LW 	0x0c0d0e0f
	#define VEC_PERM4_RX 	0x10111213
	#define VEC_PERM4_RY 	0x14151617
	#define VEC_PERM4_RZ 	0x18191a1b
	#define VEC_PERM4_RW 	0x1c1d1e1f
	#define VEC_PERM4_LB0	0x0004080c
	#define VEC_PERM4_LB1	0x0105090d
	#define VEC_PERM4_LB2	0x02060a0e
	#define VEC_PERM4_LB3	0x03070b0f
	#define VEC_PERM4_RB0	0x1014181c
	#define VEC_PERM4_RB1	0x1115191d
	#define VEC_PERM4_RB2	0x12161a1e
	#define VEC_PERM4_RB3	0x13171b1f
	#define VEC_PERM4(X,Y,Z,W) (vec_uchar16) (vec_uint4) { VEC_PERM4_ ## X, VEC_PERM4_ ## Y, VEC_PERM4_ ## Z, VEC_PERM4_ ## W }
	#define VEC_PERM2(X,Y,Z,W) (vec_uchar16) (vec_uint4) { VEC_PERM4_L ## X, VEC_PERM4_L ## Y, VEC_PERM4_R ## Z, VEC_PERM4_R ## W }
	#define VEC_SWIZZLE(V,X,Y,Z,W) vec_perm(V,V,VEC_PERM2(X,Y,Z,W))
	#define VEC_FLOAT4(X,Y,Z,W) (vec_float4) { X, Y, Z, W }
	#define VEC_UINT4(X,Y,Z,W) (vec_uint4) { X, Y, Z, W }
#elif USE_NEON
	#if defined(_WINRT) || defined(_IOS)
		#define NEON_UINT2(X,Y) vset_lane_u32(Y,vdup_n_u32(X),1)
		#define NEON_UINT4(X,Y,Z,W) vsetq_lane_u32(W,vsetq_lane_u32(Z,vsetq_lane_u32(Y,vdupq_n_u32(X),1),2),3)
	#else
		#define NEON_UINT2(X,Y) { X, Y }
		#define NEON_UINT4(X,Y,Z,W) { X, Y, Z, W }
	#endif
#endif

/*
 */
#ifdef USE_DOUBLE
	#define Scalar			double
	#define Vec2			dvec2
	#define Vec3			dvec3
	#define Vec4			dvec4
	#define Mat4			dmat4
	#define Vec2_zero		dvec2_zero
	#define Vec3_zero		dvec3_zero
	#define Vec4_zero		dvec4_zero
	#define Vec2_one		dvec2_one
	#define Vec3_one		dvec3_one
	#define Vec4_one		dvec4_one
	#define Vec2_epsilon	dvec2_epsilon
	#define Vec3_epsilon	dvec3_epsilon
	#define Vec4_epsilon	dvec4_epsilon
	#define Vec2_infinity	dvec2_infinity
	#define Vec3_infinity	dvec3_infinity
	#define Vec4_infinity	dvec4_infinity
	#define Mat4_identity	dmat4_identity
#else
	#define Scalar			float
	#define Vec2			vec2
	#define Vec3			vec3
	#define Vec4			vec4
	#define Mat4			mat4
	#define Vec2_zero		vec2_zero
	#define Vec3_zero		vec3_zero
	#define Vec4_zero		vec4_zero
	#define Vec2_one		vec2_one
	#define Vec3_one		vec3_one
	#define Vec4_one		vec4_one
	#define Vec2_epsilon	vec2_epsilon
	#define Vec3_epsilon	vec3_epsilon
	#define Vec4_epsilon	vec4_epsilon
	#define Vec2_infinity	vec2_infinity
	#define Vec3_infinity	vec3_infinity
	#define Vec4_infinity	vec4_infinity
	#define Mat4_identity	mat4_identity
#endif

/*
 */
struct vec2;
struct vec3;
struct vec4;
struct dvec2;
struct dvec3;
struct dvec4;
struct hvec2;
struct hvec3;
struct hvec4;
struct ivec2;
struct ivec3;
struct ivec4;
struct bvec4;
struct mat2;
struct mat3;
struct mat4;
struct dmat4;
struct quat;

/******************************************************************************\
*
* Type conversion
*
\******************************************************************************/

/*
 */
union IntFloat {
	IntFloat() { }
	IntFloat(int i) : i(i) { }
	IntFloat(float f) : f(f) { }
	IntFloat(unsigned int ui) : ui(ui) { }
	int i;
	float f;
	unsigned int ui;
};

/*
 */
union IntDouble {
	IntDouble() { }
	IntDouble(double d) : d(d) { }
	int i[2];
	double d;
	unsigned int ui[2];
};

/******************************************************************************\
*
* MathLib
*
\******************************************************************************/

/*
 */
class Math {
		
		Math();
		
	public:
		
		// projection
		static void setGL();
		static void setD3D9();
		static void setD3D10();
		
		// functions
		static int signMask(int v);
		static float sign(float v);
		static double sign(double v);
		
		static int abs(int v);
		static long long abs(long long v);
		static float abs(float v);
		static double abs(double v);
		
		static float ceil(float v);
		static double ceil(double v);
		
		static float floor(float v);
		static double floor(double v);
		
		static float frac(float v);
		static double frac(double v);
		
		static float sqrt(float v);
		static double sqrt(double v);
		static float sqrtFast(float v);
		
		static float rcp(float v);
		static double rcp(double v);
		static float rcpFast(float v);
		
		static float rsqrt(float v);
		static double rsqrt(double v);
		static float rsqrtFast(float v);
		
		static float mod(float x,float y);
		static double mod(double x,double y);
		
		static float pow(float x,float y);
		static double pow(double x,double y);
		static float powFast(float x,float y);
		
		static float exp(float v);
		static double exp(double v);
		static float expFast(float v);
		
		static float exp2(float v);
		static double exp2(double v);
		static float exp2Fast(float v);
		
		static float log(float v);
		static double log(double v);
		static float logFast(float v);
		
		static int log2(int v);
		static float log2(float v);
		static double log2(double v);
		static float log2Fast(float v);
		
		static float log10(float v);
		static double log10(double v);
		
		// trigonometry
		static float sin(float a);
		static double sin(double a);
		static float sinFast(float a);
		
		static float cos(float a);
		static double cos(double a);
		static float cosFast(float a);
		
		static float tan(float a);
		static double tan(double a);
		
		static float asin(float v);
		static double asin(double v);
		static float asinFast(float v);
		
		static float acos(float v);
		static double acos(double v);
		static float acosFast(float v);
		
		static float atan(float v);
		static double atan(double v);
		
		static float atan2(float y,float x);
		static double atan2(double y,double x);
		
		static void sincos(float a,float &s,float &c);
		static void sincos(double a,double &s,double &c);
		static void sincosFast(float a,float &s,float &c);
		
		// branching
		static int select(int c,int v0,int v1);
		static float select(int c,float v0,float v1);
		static float select(float c,float v0,float v1);
		
		// conversion
		static float itof(int v);
		static int ftoi(float v);
		static int round(float v);
		
		static double itod(int v);
		static int dtoi(double v);
		
		static float ltof(long long v);
		static long long ftol(float v);
		
		static double ltod(long long v);
		static long long dtol(double v);
		
		// nearest power of two
		static int npot(int v);
		
		// rounding up division
		static int udiv(int x,int y);
		
		// bezier curve solver
		static float bezier(const float *t,const float *v,float time);
		static double bezier(const float *t,const double *v,float time);
		
		// memory
		static void prefetch(const void *ptr);
		static void memset(void *dest,int c,size_t size);
		static void memcpy(void *dest,const void *src,size_t size);
		static int memcmp(const void *src_0,const void *src_1,size_t size);
};

/*
 */
INLINE int Math::signMask(int v) {
	return (v >> 31);
}

INLINE float Math::sign(float v) {
	return IntFloat((IntFloat(v).ui & 0x80000000) | 0x3f800000).f;
}

INLINE double Math::sign(double v) {
	if(v >= 0.0) return 1.0;
	return -1.0;
}

/*
 */
INLINE int Math::abs(int v) {
	if(v >= 0) return v;
	return -v;
}

INLINE long long Math::abs(long long v) {
	if(v >= 0) return v;
	return -v;
}

INLINE float Math::abs(float v) {
	#ifdef _CELLOS_LV2
		return __fabsf(v);
	#else
		if(v >= 0.0f) return v;
		return -v;
	#endif
}

INLINE double Math::abs(double v) {
	#ifdef _CELLOS_LV2
		return __fabs(v);
	#else
		if(v >= 0.0) return v;
		return -v;
	#endif
}

/*
 */
INLINE float Math::ceil(float v) {
	return ::ceilf(v);
}

INLINE double Math::ceil(double v) {
	return ::ceil(v);
}

/*
 */
INLINE float Math::floor(float v) {
	return ::floorf(v);
}

INLINE double Math::floor(double v) {
	return ::floor(v);
}

/*
 */
INLINE float Math::frac(float v) {
	return v - ::floorf(v);
}

INLINE double Math::frac(double v) {
	return v - ::floor(v);
}

/*
 */
INLINE float Math::sqrt(float v) {
	#if defined(_CELLOS_LV2) && !defined(__SNC__)
		return __fsqrts(v);
	#elif USE_SSE
		_mm_store_ss(&v,_mm_sqrt_ss(_mm_set_ss(v)));
		return v;
	#else
		return ::sqrtf(v);
	#endif
}

INLINE double Math::sqrt(double v) {
	#if defined(_CELLOS_LV2) && !defined(__SNC__)
		return __fsqrt(v);
	#else
		return ::sqrt(v);
	#endif
}

INLINE float Math::sqrtFast(float v) {
	#ifdef ARCH_ARM
		return ::sqrtf(v);
	#else
		IntFloat i = v;
		i.i = 0x5f3759df - (i.i >> 1);
		v = i.f * v;
		return v * (1.5f - (i.f * v * 0.5f));
	#endif
}

/*
 */
INLINE float Math::rcp(float v) {
	return 1.0f / v;
}

INLINE double Math::rcp(double v) {
	return 1.0 / v;
}

INLINE float Math::rcpFast(float v) {
	#if defined(_CELLOS_LV2) && !defined(__SNC__)
		return __fres(v);
	#elif USE_SSE
		_mm_store_ss(&v,_mm_rcp_ss(_mm_set_ss(v)));
		return v;
	#elif ARCH_ARM
		return 1.0f / v;
	#else
		IntFloat i = v;
		i.i = 0x7f000000 - i.i;
		return i.f * (2.0f - v * i.f);
	#endif
}

/*
 */
INLINE float Math::rsqrt(float v) {
	if(v < 1e-18f) return INFINITY;
	return 1.0f / ::sqrtf(v);
}

INLINE double Math::rsqrt(double v) {
	if(v < 1e-18) return INFINITY;
	return 1.0 / ::sqrt(v);
}

INLINE float Math::rsqrtFast(float v) {
	#ifdef _CELLOS_LV2
		return (float)__frsqrte(v);
	#elif USE_SSE
		_mm_store_ss(&v,_mm_rsqrt_ss(_mm_set_ss(v)));
		return v;
	#else
		IntFloat i = v;
		i.i = 0x5f3759df - (i.i >> 1);
		return i.f * (1.5f - (i.f * i.f * v * 0.5f));
	#endif
}

/*
 */
INLINE float Math::mod(float x,float y) {
	return ::fmodf(x,y);
}

INLINE double Math::mod(double x,double y) {
	return ::fmod(x,y);
}

/*
 */
INLINE float Math::pow(float x,float y) {
	return ::powf(x,y);
}

INLINE double Math::pow(double x,double y) {
	return ::pow(x,y);
}

INLINE float Math::powFast(float x,float y) {
	return Math::exp2Fast(Math::log2Fast(x) * y);
}

/*
 */
INLINE float Math::exp(float v) {
	return ::expf(v);
}

INLINE double Math::exp(double v) {
	return ::exp(v);
}

INLINE float Math::expFast(float v) {
	return exp2Fast(v * (1.0f / LOG2));
}

/*
 */
INLINE float Math::exp2(float v) {
	return ::expf(v * LOG2);
}

INLINE double Math::exp2(double v) {
	return ::exp(v * LOG2);
}

/*
 */
INLINE float Math::log(float v) {
	return ::logf(v);
}

INLINE double Math::log(double v) {
	return ::log(v);
}

INLINE float Math::logFast(float v) {
	return log2Fast(v) * LOG2;
}

/*
 */
INLINE int Math::log2(int v) {
	int ret = 0;
	if(v >= 1 << 16) { v >>= 16; ret |= 16; }
	if(v >= 1 << 8) { v >>= 8; ret |= 8; }
	if(v >= 1 << 4) { v >>= 4; ret |= 4; }
	if(v >= 1 << 2) { v >>= 2; ret |= 2; }
	if(v >= 1 << 1) { ret |= 1; }
	return ret;
}

INLINE float Math::log2(float v) {
	return ::logf(v) * (1.0f / LOG2);
}

INLINE double Math::log2(double v) {
	return ::log(v) * (1.0 / LOG2);
}

/*
 */
INLINE float Math::log10(float v) {
	return ::logf(v) * (1.0f / LOG10);
}

INLINE double Math::log10(double v) {
	return ::log(v) * (1.0 / LOG10);
}

/*
 */
INLINE float Math::sin(float a) {
	return ::sinf(a);
}

INLINE double Math::sin(double a) {
	return ::sin(a);
}

/*
 */
INLINE float Math::cos(float a) {
	return ::cosf(a);
}

INLINE double Math::cos(double a) {
	return ::cos(a);
}

/*
 */
INLINE float Math::tan(float a) {
	return ::tanf(a);
}

INLINE double Math::tan(double a) {
	return ::tan(a);
}

/*
 */
INLINE float Math::asin(float v) {
	return ::asinf(v);
}

INLINE double Math::asin(double v) {
	return ::asin(v);
}

/*
 */
INLINE float Math::acos(float v) {
	return ::acosf(v);
}

INLINE double Math::acos(double v) {
	return ::acos(v);
}

/*
 */
INLINE float Math::atan(float v) {
	return ::atanf(v);
}

INLINE double Math::atan(double v) {
	return ::atan(v);
}

/*
 */
INLINE float Math::atan2(float y,float x) {
	return ::atan2f(y,x);
}

INLINE double Math::atan2(double y,double x) {
	return ::atan2(y,x);
}

/*
 */
INLINE void Math::sincosFast(float a,float &s,float &c) {
	if(a < 0.0f) a -= Math::ftoi(a * (1.0f / PI2)) * PI2 - PI2;
	else if(a >= PI2) a -= Math::ftoi(a * (1.0f / PI2)) * PI2;
	c = 1.0f;
	s = PI - a;
	if(s < -PI05) s = -PI - s;
	else if(s > PI05) s = PI - s;
	else c = -1.0f;
	float a2 = s * s;
	s *= ((0.00761f * a2 - 0.16605f) * a2 + 1.0f);
	c *= ((0.03705f * a2 - 0.49670f) * a2 + 1.0f);
}

/*
 */
INLINE int Math::select(int c,int v0,int v1) {
	int mask = Math::signMask(c | -c);
	return (v0 & mask) | (v1 & ~mask);
}

INLINE float Math::select(int c,float v0,float v1) {
	int mask = Math::signMask(c | -c);
	return IntFloat((IntFloat(v0).i & mask) | (IntFloat(v1).i & ~mask)).f;
}

INLINE float Math::select(float c,float v0,float v1) {
	#ifdef _CELLOS_LV2
		return __fsels(c,v1,v0);
	#else
		int mask = Math::signMask(IntFloat(c).i);
		return IntFloat((IntFloat(v0).i & mask) | (IntFloat(v1).i & ~mask)).f;
	#endif
}

/*
 */
INLINE float Math::itof(int v) {
	return static_cast<float>(v);
}

INLINE int Math::ftoi(float v) {
	#ifdef _CELLOS_LV2
		return __fctiwz(v);
	#else
		return static_cast<int>(v);
	#endif
}

INLINE int Math::round(float v) {
	#ifdef _CELLOS_LV2
		return __fctiw(v);
	#elif defined(_WIN32) && defined(ARCH_X86)
		int i;
		__asm {
			fld v
			fistp i
		}
		return i;
	#elif defined(_LINUX) && !defined(ARCH_ARM)
		int i;
		asm volatile("fistpl %0" : "=m"(i) : "t"(v) : "st");
		return i;
	#elif USE_SSE
		return _mm_cvt_ss2si(_mm_load_ss(&v));
	#else
		return static_cast<int>(v + 0.5f);
	#endif
}

/*
 */
INLINE double Math::itod(int v) {
	return static_cast<double>(v);
}

INLINE int Math::dtoi(double v) {
	return static_cast<int>(v);
}

/*
 */
INLINE float Math::ltof(long long v) {
	return static_cast<float>(v);
}

INLINE long long Math::ftol(float v) {
	return static_cast<long long>(v);
}

/*
 */
INLINE double Math::ltod(long long v) {
	return static_cast<double>(v);
}

INLINE long long Math::dtol(double v) {
	return static_cast<long long>(v);
}

/*
 */
INLINE int Math::npot(int v) {
	int i = 1;
	while(i < v) i += i;
	return i;
}

INLINE int Math::udiv(int x,int y) {
	return x / y + (x % y != 0);
}

/*
 */
INLINE void Math::prefetch(const void *ptr) {
	#ifdef USE_SSE
		_mm_prefetch((const char*)ptr,_MM_HINT_NTA);
	#elif USE_ALTIVEC
		__dcbt(ptr);
	#elif defined(ARCH_ARM) && (defined(_LINUX) || defined(_ANDROID))
		asm volatile("pld [%r0,#0]" : "=r"(ptr));
	#endif
}

/*
 */
#ifdef _WEBGL
	
	/*
	 */
	INLINE void Math::memset(void *dest,int c,size_t size) {
		::memset(dest,c,size);
	}
	
	INLINE void Math::memcpy(void *dest,const void *src,size_t size) {
		::memcpy(dest,src,size);
	}
	
	INLINE int Math::memcmp(const void *src_0,const void *src_1,size_t size) {
		return ::memcmp(src_0,src_1,size);
	}
	
#endif

/******************************************************************************\
*
* Scalars
*
\******************************************************************************/

/*
 */
INLINE int compare(float v0,float v1) {
	float v = Math::abs(v0 - v1);
	return (v < EPSILON);
}

INLINE int compare(float v0,float v1,float epsilon) {
	float v = Math::abs(v0 - v1);
	return (v < (Math::abs(v0) + Math::abs(v1) + 1.0f) * epsilon);
}

/*
 */
INLINE float min(float v0,float v1) {
	#ifdef _CELLOS_LV2
		return __fsels(v1 - v0,v0,v1);
	#else
		return (v0 < v1) ? v0 : v1;
	#endif
}

INLINE float max(float v0,float v1) {
	#ifdef _CELLOS_LV2
		return __fsels(v0 - v1,v0,v1);
	#else
		return (v0 > v1) ? v0 : v1;
	#endif
}

INLINE float clamp(float v,float v0,float v1) {
	#ifdef _CELLOS_LV2
		v = __fsels(v - v0,v,v0);
		return __fsels(v1 - v,v,v1);
	#else
		if(v < v0) return v0;
		if(v > v1) return v1;
		return v;
	#endif
}

INLINE float saturate(float v) {
	#ifdef _CELLOS_LV2
		v = __fsels(v,v,0.0f);
		return __fsels(1.0f - v,v,1.0f);
	#else
		if(v < 0.0f) return 0.0f;
		if(v > 1.0f) return 1.0f;
		return v;
	#endif
}

INLINE float lerp(float v0,float v1,float k) {
	return v0 + (v1 - v0) * k;
}

/*
 */
INLINE int compare(double v0,double v1) {
	double v = Math::abs(v0 - v1);
	return (v < EPSILON);
}

INLINE int compare(double v0,double v1,double epsilon) {
	double v = Math::abs(v0 - v1);
	return (v < (Math::abs(v0) + Math::abs(v1) + 1.0) * epsilon);
}

/*
 */
INLINE double min(double v0,double v1) {
	#ifdef _CELLOS_LV2
		return __fsel(v1 - v0,v0,v1);
	#else
		return (v0 < v1) ? v0 : v1;
	#endif
}

INLINE double max(double v0,double v1) {
	#ifdef _CELLOS_LV2
		return __fsel(v0 - v1,v0,v1);
	#else
		return (v0 > v1) ? v0 : v1;
	#endif
}

INLINE double clamp(double v,double v0,double v1) {
	#ifdef _CELLOS_LV2
		v = __fsel(v - v0,v,v0);
		return __fsel(v1 - v,v,v1);
	#else
		if(v < v0) return v0;
		if(v > v1) return v1;
		return v;
	#endif
}

INLINE double saturate(double v) {
	#ifdef _CELLOS_LV2
		v = __fsel(v,v,0.0);
		return __fsel(1.0 - v,v,1.0);
	#else
		if(v < 0.0) return 0.0;
		if(v > 1.0) return 1.0;
		return v;
	#endif
}

INLINE double lerp(double v0,double v1,double k) {
	return v0 + (v1 - v0) * k;
}

/*
 */
INLINE int min(int v0,int v1) {
	return (v0 < v1) ? v0 : v1;
}

INLINE int max(int v0,int v1) {
	return (v0 > v1) ? v0 : v1;
}

INLINE int clamp(int v,int v0,int v1) {
	if(v < v0) return v0;
	if(v > v1) return v1;
	return v;
}

INLINE int lerp(int v0,int v1,int k) {
	return v0 + (((v1 - v0) * k) >> 16);
}

/******************************************************************************\
*
* Vectors
*
\******************************************************************************/

/*
 */
INLINE float dot33(const float * RESTRICT v0,const float * RESTRICT v1) {
	return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

INLINE float dot34(const float * RESTRICT v0,const float * RESTRICT v1) {
	return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] + v1[3];
}

INLINE float dot43(const float * RESTRICT v0,const float * RESTRICT v1) {
	return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] + v0[3];
}

INLINE float dot44(const float * RESTRICT v0,const float * RESTRICT v1) {
	return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] + v0[3] * v1[3];
}

INLINE float *mul3(float * RESTRICT ret,const float * RESTRICT v0,float v1) {
	ret[0] = v0[0] * v1;
	ret[1] = v0[1] * v1;
	ret[2] = v0[2] * v1;
	return ret;
}

INLINE float *add3(float * RESTRICT ret,const float * RESTRICT v0,const float * RESTRICT v1) {
	ret[0] = v0[0] + v1[0];
	ret[1] = v0[1] + v1[1];
	ret[2] = v0[2] + v1[2];
	return ret;
}

INLINE float *sub3(float * RESTRICT ret,const float * RESTRICT v0,const float * RESTRICT v1) {
	ret[0] = v0[0] - v1[0];
	ret[1] = v0[1] - v1[1];
	ret[2] = v0[2] - v1[2];
	return ret;
}

INLINE float *cross3(float * RESTRICT ret,const float * RESTRICT v0,const float * RESTRICT v1) {
	ret[0] = v0[1] * v1[2] - v0[2] * v1[1];
	ret[1] = v0[2] * v1[0] - v0[0] * v1[2];
	ret[2] = v0[0] * v1[1] - v0[1] * v1[0];
	return ret;
}

/*
 */
#ifdef USE_SSE
	
	/*
	 */
	INLINE __m128 _mm_rcp_ss_nr(__m128 v) {
		__m128 iv = _mm_rcp_ss(v);
		iv = _mm_sub_ss(_mm_add_ss(iv,iv),_mm_mul_ss(v,_mm_mul_ss(iv,iv)));
		return _mm_sub_ss(_mm_add_ss(iv,iv),_mm_mul_ss(v,_mm_mul_ss(iv,iv)));
	}
	
	INLINE __m128 _mm_rcp_ps_nr(__m128 v) {
		__m128 iv = _mm_rcp_ps(v);
		iv = _mm_sub_ps(_mm_add_ps(iv,iv),_mm_mul_ps(v,_mm_mul_ps(iv,iv)));
		return _mm_sub_ps(_mm_add_ps(iv,iv),_mm_mul_ps(v,_mm_mul_ps(iv,iv)));
	}
	
	/*
	 */
	INLINE __m128 _mm_rsqrt_ss_nr(__m128 v) {
		__m128 iv = _mm_rsqrt_ss(v);
		__m128 nr = _mm_mul_ss(_mm_mul_ss(v,iv),iv);
		return _mm_mul_ss(_mm_mul_ss(_mm_set1_ps(0.5f),iv),_mm_sub_ss(_mm_set1_ps(3.0f),nr));
	}
	
	INLINE __m128 _mm_rsqrt_ps_nr(__m128 v) {
		__m128 iv = _mm_rsqrt_ps(v);
		__m128 nr = _mm_mul_ps(_mm_mul_ps(v,iv),iv);
		return _mm_mul_ps(_mm_mul_ps(_mm_set1_ps(0.5f),iv),_mm_sub_ps(_mm_set1_ps(3.0f),nr));
	}
	
	/*
	 */
	INLINE __m128 _mm_dot33_ps(__m128 v0,__m128 v1) {
		__m128 v2 = _mm_mul_ps(v0,v1);
		__m128 v3 = _mm_add_ps(v2,_MM_SWIZZLE(v2,Y,X,Y,W));
		return _mm_add_ps(v3,_MM_SWIZZLE(v2,Z,Z,X,W));
	}
	
	INLINE __m128 _mm_dot44_ps(__m128 v0,__m128 v1) {
		__m128 v2 = _mm_mul_ps(v0,v1);
		v2 = _mm_add_ps(v2,_MM_SWIZZLE(v2,Y,X,W,Z));
		return _mm_add_ps(v2,_MM_SWIZZLE(v2,Z,W,X,Y));
	}
	
	/*
	 */
	INLINE __m128 _mm_normalize3_ps(__m128 v) {
		__m128 length2 = _mm_dot33_ps(v,v);
		return _mm_mul_ps(v,_mm_rsqrt_ps_nr(length2));
	}
	
	INLINE __m128 _mm_normalize4_ps(__m128 v) {
		__m128 length2 = _mm_dot44_ps(v,v);
		return _mm_mul_ps(v,_mm_rsqrt_ps_nr(length2));
	}
	
	/*
	 */
	INLINE __m128 _mm_cross_ps(__m128 v0,__m128 v1) {
		__m128 v0_yzxw = _MM_SWIZZLE(v0,Y,Z,X,W);
		__m128 v1_yzxw = _MM_SWIZZLE(v1,Y,Z,X,W);
		__m128 v2 = _mm_sub_ps(_mm_mul_ps(v0,v1_yzxw),_mm_mul_ps(v1,v0_yzxw));
		return _MM_SWIZZLE(v2,Y,Z,X,W);
	}
	
/*
 */
#elif USE_ALTIVEC
	
	/*
	 */
	INLINE vec_float4 vec_rcp_nr(vec_float4 v) {
		vec_float4 iv = vec_re(v);
		vec_float4 one = vec_splats(1.0f);
		iv = vec_madd(vec_nmsub(iv,v,one),iv,iv);
		return vec_madd(vec_nmsub(iv,v,one),iv,iv);
	}
	
	/*
	 */
	INLINE vec_float4 vec_rsqrt_nr(vec_float4 v) {
		vec_float4 iv = vec_rsqrte(v);
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 nr = vec_madd(vec_madd(v,iv,zero),iv,zero);
		return vec_madd(vec_madd(vec_splats(0.5f),iv,zero),vec_sub(vec_splats(3.0f),nr),zero);
	}
	
	/*
	 */
	INLINE vec_float4 vec_dot33(vec_float4 v0,vec_float4 v1) {
		vec_float4 v2 = vec_madd(v0,v1,vec_splats(0.0f));
		vec_float4 v3 = vec_add(v2,VEC_SWIZZLE(v2,Y,X,Y,W));
		return vec_add(v3,VEC_SWIZZLE(v2,Z,Z,X,W));
	}
	
	INLINE vec_float4 vec_dot44(vec_float4 v0,vec_float4 v1) {
		vec_float4 v2 = vec_madd(v0,v1,vec_splats(0.0f));
		v2 = vec_add(v2,vec_sld(v2,v2,8));
		return vec_add(v2,vec_sld(v2,v2,4));
	}
	
	/*
	 */
	INLINE vec_float4 vec_normalize3(vec_float4 v) {
		vec_float4 length2 = vec_dot33(v,v);
		return vec_madd(v,vec_rsqrt_nr(length2),vec_splats(0.0f));
	}
	
	INLINE vec_float4 vec_normalize4(vec_float4 v) {
		vec_float4 length2 = vec_dot44(v,v);
		return vec_madd(v,vec_rsqrt_nr(length2),vec_splats(0.0f));
	}
	
	/*
	 */
	INLINE vec_float4 vec_cross(vec_float4 v0,vec_float4 v1) {
		vec_uchar16 yzxw = VEC_PERM2(Y,Z,X,W);
		vec_float4 v0_yzxw = vec_perm(v0,v0,yzxw);
		vec_float4 v1_yzxw = vec_perm(v1,v1,yzxw);
		vec_float4 v2 = vec_nmsub(v1,v0_yzxw,vec_madd(v0,v1_yzxw,vec_splats(0.0f)));
		return vec_perm(v2,v2,yzxw);
	}
	
/*
 */
#elif USE_NEON
	
	/*
	 */
	INLINE int vmask_u32(uint32x2_t v) {
		const uint32x2_t mask = NEON_UINT2(1,2);
		#if defined(_WINRT) || defined(_IOS)
			uint32x2_t res_0 = vand_u32(v,mask);
			res_0 = vpadd_u32(res_0,res_0);
			return vget_lane_u32(res_0,0);
		#else
			int ret;
			asm volatile(
				"vand.u32  %P1, %P1, %P2	\n"
				"vpadd.u32 %P1, %P1, %P1	\n"
				"vmov.u32  %r0, %P1[0]		\n"
				: "=r"(ret) : "w"(v), "w"(mask)
			);
			return ret;
		#endif
	}
	
	INLINE int vmaskq_u32(uint32x4_t v) {
		const uint32x4_t mask = NEON_UINT4(1,2,4,8);
		#if defined(_WINRT) || defined(_IOS)
			uint32x4_t res_0 = vandq_u32(v,mask);
			uint32x2_t res_1 = vget_low_u32(res_0);
			uint32x2_t res_2 = vget_high_u32(res_0);
			res_1 = vorr_u32(res_1,res_2);
			res_1 = vpadd_u32(res_1,res_1);
			return vget_lane_u32(res_1,0);
		#else
			int ret;
			asm volatile(
				"vand.u32  %q1, %q1, %q2	\n"
				"vorr.u32  %e1, %e1, %f1	\n"
				"vpadd.u32 %e1, %e1, %e1	\n"
				"vmov.u32  %r0, %e1[0]		\n"
				: "=r"(ret) : "w"(v), "w"(mask)
			);
			return ret;
		#endif
	}
	
	/*
	 */
	INLINE float32x2_t vrcp_nr_f32(float32x2_t v) {
		#if defined(_WINRT) || defined(_IOS)
			float32x2_t res_0 = vrecpe_f32(v);
			float32x2_t res_1 = vrecps_f32(res_0,v);
			res_0 = vmul_f32(res_0,res_1);
			res_1 = vrecps_f32(res_0,v);
			return vmul_f32(res_0,res_1);
		#else
			float32x2_t ret;
			asm volatile(
				"vrecpe.f32  d0, %P1		\n"
				"vrecps.f32  d1,  d0, %P1	\n"
				"vmul.f32    d0,  d0,  d1	\n"
				"vrecps.f32  d1,  d0, %P1	\n"
				"vmul.f32   %P0,  d0,  d1	\n"
				: "=w"(ret) : "w"(v) : "d0", "d1"
			);
			return ret;
		#endif
	}
	
	INLINE float32x4_t vrcpq_nr_f32(float32x4_t v) {
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t res_0 = vrecpeq_f32(v);
			float32x4_t res_1 = vrecpsq_f32(res_0,v);
			res_0 = vmulq_f32(res_0,res_1);
			res_1 = vrecpsq_f32(res_0,v);
			return vmulq_f32(res_0,res_1);
		#else
			float32x4_t ret;
			asm volatile(
				"vrecpe.f32  q0, %q1		\n"
				"vrecps.f32  q1,  q0, %q1	\n"
				"vmul.f32    q0,  q0,  q1	\n"
				"vrecps.f32  q1,  q0, %q1	\n"
				"vmul.f32   %q0,  q0,  q1	\n"
				: "=w"(ret) : "w"(v) : "q0", "q1"
			);
			return ret;
		#endif
	}
	
	/*
	 */
	INLINE float32x2_t vrsqrt_nr_f32(float32x2_t v) {
		#if defined(_WINRT) || defined(_IOS)
			float32x2_t res_0 = vrsqrte_f32(v);
			float32x2_t res_1 = vmul_f32(res_0,v);
			float32x2_t res_2 = vrsqrts_f32(res_0,res_1);
			return vmul_f32(res_0,res_2);
		#else
			float32x2_t ret;
			asm volatile(
				"vrsqrte.f32  d0, %P1		\n"
				"vmul.f32     d1,  d0, %P1	\n"
				"vrsqrts.f32  d2,  d0,  d1	\n"
				"vmul.f32    %P0,  d0,  d2	\n"
				: "=w"(ret) : "w"(v) : "d0", "d1", "d2"
			);
			return ret;
		#endif
	}
	
	INLINE float32x4_t vrsqrtq_nr_f32(float32x4_t v) {
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t res_0 = vrsqrteq_f32(v);
			float32x4_t res_1 = vmulq_f32(res_0,v);
			float32x4_t res_2 = vrsqrtsq_f32(res_0,res_1);
			return vmulq_f32(res_0,res_2);
		#else
			float32x4_t ret;
			asm volatile(
				"vrsqrte.f32  q0, %q1		\n"
				"vmul.f32     q1,  q0, %q1	\n"
				"vrsqrts.f32  q2,  q0,  q1	\n"
				"vmul.f32    %q0,  q0,  q2	\n"
				: "=w"(ret) : "w"(v) : "q0", "q1", "q2"
			);
			return ret;
		#endif
	}
	
	/*
	 */
	INLINE float32x4_t vdot33q_f32(float32x4_t v0,float32x4_t v1) {
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t res_0 = vmulq_f32(v0,v1);
			float32x2_t res_1 = vget_low_f32(res_0);
			float32x2_t res_2 = vget_high_f32(res_0);
			float32x2_t res_3 = vpadd_f32(res_1,res_1);
			res_3 = vadd_f32(res_3,res_2);
			return vdupq_lane_f32(res_3,0);
		#else
			float32x4_t ret;
			asm volatile(
				"vmul.f32   q0, %q1,   %q2	\n"
				"vpadd.f32  d0,  d0,    d0	\n"
				"vadd.f32   d0,  d0,    d1	\n"
				"vdup.32   %q0,  d0[0]		\n"
				: "=w"(ret) : "w"(v0), "w"(v1) : "q0"
			);
			return ret;
		#endif
	}
	
	INLINE float32x4_t vdot44q_f32(float32x4_t v0,float32x4_t v1) {
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t res_0 = vmulq_f32(v0,v1);
			float32x4_t res_1 = vextq_f32(res_0,res_0,2);
			res_0 = vaddq_f32(res_0,res_1);
			res_1 = vextq_f32(res_0,res_0,1);
			return vaddq_f32(res_0,res_1);
		#else
			float32x4_t ret;
			asm volatile(
				"vmul.f32  q0, %q1, %q2		\n"
				"vext.32   q1,  q0,  q0, #2	\n"
				"vadd.f32  q0,  q0,  q1		\n"
				"vext.32   q1,  q0,  q0, #1	\n"
				"vadd.f32 %q0,  q0,  q1		\n"
				: "=w"(ret) : "w"(v0), "w"(v1) : "q0", "q1"
			);
			return ret;
		#endif
	}
	
	/*
	 */
	INLINE float32x4_t vnormalize3q_f32(float32x4_t v) {
		float32x4_t length2 = vdot33q_f32(v,v);
		return vmulq_f32(v,vrsqrtq_nr_f32(length2));
	}
	
	INLINE float32x4_t vnormalize4q_f32(float32x4_t v) {
		float32x4_t length2 = vdot44q_f32(v,v);
		return vmulq_f32(v,vrsqrtq_nr_f32(length2));
	}
	
	/*
	 */
	INLINE float32x4_t vcrossq_f32(float32x4_t v0,float32x4_t v1) {
		#if defined(_WINRT) || defined(_IOS)
			float32x2_t low_0 = vget_low_f32(v0);
			float32x2_t low_1 = vget_low_f32(v1);
			float32x2_t high_0 = vget_high_f32(v0);
			float32x2_t high_1 = vget_high_f32(v1);
			float32x2_t res_0 = vext_f32(low_0,high_0,1);
			float32x2_t res_1 = vext_f32(low_1,high_1,1);
			float32x2_t res_2 = vmul_f32(high_0,low_1);
			float32x2_t res_3 = vmul_f32(low_0,res_1);
			res_2 = vmls_f32(res_2,high_1,low_0);
			res_3 = vmls_f32(res_3,low_1,res_0);
			res_0 = vext_f32(res_3,res_2,1);
			return vcombine_f32(res_0,res_3);
		#else
			float32x4_t ret;
			asm volatile(
				"vext.32   d0, %e1, %f1, #1	\n"
				"vext.32   d1, %e2, %f2, #1	\n"
				"vmul.f32  d2, %f1, %e2		\n"
				"vmul.f32  d3, %e1,  d1		\n"
				"vmls.f32  d2, %f2, %e1		\n"
				"vmls.f32  d3, %e2,  d0		\n"
				"vext.32  %e0,  d3,  d2, #1	\n"
				"vmov     %f0,  d3			\n"
				: "=w"(ret) : "w"(v0), "w"(v1) : "q0", "q1"
			);
			return ret;
		#endif
	}
	
#endif

/******************************************************************************\
*
* vec2
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED8(struct) vec2 {
	
	INLINE vec2() { }
	INLINE vec2(const vec2 &v) : x(v.x), y(v.y) { }
	INLINE vec2(float x,float y) : x(x), y(y) { }
	explicit INLINE vec2(float v) : x(v), y(v) { }
	explicit INLINE vec2(const vec3 &v);
	explicit INLINE vec2(const vec4 &v);
	explicit INLINE vec2(const float *v) : x(v[0]), y(v[1]) { }
	explicit INLINE vec2(const dvec2 &v);
	explicit INLINE vec2(const hvec2 &v);
	explicit INLINE vec2(const ivec2 &v);
	#ifdef USE_SSE
		explicit INLINE vec2(__m64 v) : vec(v) { }
	#elif USE_NEON
		explicit INLINE vec2(float32x2_t v) : vec(v) { }
	#endif
	
	INLINE vec2 &operator=(const vec2 &v) {
		x = v.x; y = v.y;
		return *this;
	}
	INLINE vec2 operator-() const {
		return vec2(-x,-y);
	}
	INLINE vec2 &operator*=(float v) {
		x *= v; y *= v;
		return *this;
	}
	INLINE vec2 &operator*=(const vec2 &v) {
		x *= v.x; y *= v.y;
		return *this;
	}
	INLINE vec2 &operator/=(float v) {
		float iv = Math::rcp(v);
		x *= iv; y *= iv;
		return *this;
	}
	INLINE vec2 &operator/=(const vec2 &v) {
		x /= v.x; y /= v.y;
		return *this;
	}
	INLINE vec2 &operator+=(const vec2 &v) {
		x += v.x; y += v.y;
		return *this;
	}
	INLINE vec2 &operator-=(const vec2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}
	
	INLINE operator float*() { return v; }
	INLINE operator const float*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 2 && "vec2::operator[](): bad index");
		return v[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 2 && "vec2::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(float v) {
		x = v; y = v;
	}
	INLINE void set(float x_,float y_) {
		x = x_; y = y_;
	}
	INLINE void set(const vec2 &v) {
		x = v.x; y = v.y;
	}
	INLINE void set(const float *v) {
		x = v[0]; y = v[1];
	}
	INLINE void get(float *v) const {
		v[0] = x; v[1] = y;
	}
	INLINE float *get() { return v; }
	INLINE const float *get() const { return v; }
	
	INLINE float length2() const {
		return x * x + y * y;
	}
	INLINE float length() const {
		return Math::sqrt(x * x + y * y);
	}
	INLINE vec2 &normalize() {
		float ilength = Math::rsqrt(x * x + y * y);
		x *= ilength; y *= ilength;
		return *this;
	}
	INLINE vec2 &normalizeFast() {
		float ilength = Math::rsqrtFast(x * x + y * y);
		x *= ilength; y *= ilength;
		return *this;
	}
	
	union {
		struct {
			float x,y;
		};
		float v[2];
		#ifdef USE_SSE
			__m64 vec;
		#elif USE_NEON
			float32x2_t vec;
		#endif
	};
};

/*
 */
extern const vec2 vec2_zero;
extern const vec2 vec2_one;
extern const vec2 vec2_epsilon;
extern const vec2 vec2_infinity;

/*
 */
INLINE int operator==(const vec2 &v0,const vec2 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE int operator!=(const vec2 &v0,const vec2 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE vec2 operator*(const vec2 &v0,float v1) {
	return vec2(v0.x * v1,v0.y * v1);
}

INLINE vec2 operator*(const vec2 &v0,const vec2 &v1) {
	return vec2(v0.x * v1.x,v0.y * v1.y);
}

INLINE vec2 operator/(const vec2 &v0,float v1) {
	float iv1 = Math::rcp(v1);
	return vec2(v0.x * iv1,v0.y * iv1);
}

INLINE vec2 operator/(const vec2 &v0,const vec2 &v1) {
	return vec2(v0.x / v1.x,v0.y / v1.y);
}

INLINE vec2 operator+(const vec2 &v0,const vec2 &v1) {
	return vec2(v0.x + v1.x,v0.y + v1.y);
}

INLINE vec2 operator-(const vec2 &v0,const vec2 &v1) {
	return vec2(v0.x - v1.x,v0.y - v1.y);
}

/*
 */
INLINE int compare(const vec2 &v0,const vec2 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE int compare(const vec2 &v0,const vec2 &v1,float epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon));
}

INLINE float dot(const vec2 &v0,const vec2 &v1) {
	return v0.x * v1.x + v0.y * v1.y;
}

INLINE vec2 &mul(vec2 &ret,const vec2 &v0,float v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	return ret;
}

INLINE vec2 &mul(vec2 &ret,const vec2 &v0,const vec2 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	return ret;
}

INLINE vec2 &mad(vec2 &ret,const vec2 &v0,float v1,const vec2 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	return ret;
}

INLINE vec2 &mad(vec2 &ret,const vec2 &v0,const vec2 &v1,const vec2 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	return ret;
}

INLINE vec2 &add(vec2 &ret,const vec2 &v0,const vec2 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	return ret;
}

INLINE vec2 &sub(vec2 &ret,const vec2 &v0,const vec2 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	return ret;
}

INLINE vec2 &lerp(vec2 &ret,const vec2 &v0,const vec2 &v1,float k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	return ret;
}

/*
 */
INLINE float length(const vec2 &v) {
	return v.length();
}

INLINE float length2(const vec2 &v) {
	return v.length2();
}

INLINE vec2 normalize(const vec2 &v) {
	vec2 ret = v;
	return ret.normalize();
}

INLINE float cross(const vec2 &v0,const vec2 &v1) {
	return v0.x * v1.y - v0.y * v1.x;
}

/*
 */
vec2 min(const vec2 &v0,const vec2 &v1);
vec2 max(const vec2 &v0,const vec2 &v1);
vec2 clamp(const vec2 &v,const vec2 &v0,const vec2 &v1);
vec2 saturate(const vec2 &v);
vec2 lerp(const vec2 &v0,const vec2 &v1,float k);

/******************************************************************************\
*
* vec3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) vec3 {
	
	INLINE vec3() { }
	INLINE vec3(const vec2 &v,float z) : x(v.x), y(v.y), z(z), w(0.0f) { }
	INLINE vec3(float x,float y,float z) : x(x), y(y), z(z), w(0.0f) { }
	explicit INLINE vec3(float v) : x(v), y(v), z(v), w(0.0f) { }
	explicit INLINE vec3(const vec2 &v) : x(v.x), y(v.y), z(0.0f), w(0.0f) { }
	explicit INLINE vec3(const vec4 &v);
	explicit INLINE vec3(const float *v) : x(v[0]), y(v[1]), z(v[2]), w(0.0f) { }
	explicit INLINE vec3(const dvec3 &v);
	explicit INLINE vec3(const hvec3 &v);
	explicit INLINE vec3(const ivec3 &v);
	#ifdef USE_SSE
		INLINE vec3(const vec3 &v) : vec(v.vec) { }
		explicit INLINE vec3(__m128 v) : vec(v) { }
	#elif USE_ALTIVEC
		INLINE vec3(const vec3 &v) : vec(v.vec) { }
		explicit INLINE vec3(vec_float4 v) : vec(v) { }
	#elif USE_NEON
		INLINE vec3(const vec3 &v) : vec(v.vec) { }
		explicit INLINE vec3(float32x4_t v) : vec(v) { }
	#else
		INLINE vec3(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(0.0f) { }
	#endif
	
	INLINE vec3 &operator=(const vec3 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE vec3 operator-() const {
		return vec3(-x,-y,-z);
	}
	INLINE vec3 &operator*=(float v) {
		x *= v; y *= v; z *= v;
		return *this;
	}
	INLINE vec3 &operator*=(const vec3 &v) {
		x *= v.x; y *= v.y; z *= v.z;
		return *this;
	}
	INLINE vec3 &operator/=(float v) {
		float iv = Math::rcp(v);
		x *= iv; y *= iv; z *= iv;
		return *this;
	}
	INLINE vec3 &operator/=(const vec3 &v) {
		x /= v.x; y /= v.y; z /= v.z;
		return *this;
	}
	INLINE vec3 &operator+=(const vec3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	INLINE vec3 &operator-=(const vec3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	
	INLINE operator float*() { return v; }
	INLINE operator const float*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 3 && "vec3::operator[](): bad index");
		return v[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 3 && "vec3::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(float v) {
		x = v; y = v; z = v;
	}
	INLINE void set(float x_,float y_,float z_) {
		x = x_; y = y_; z = z_;
	}
	INLINE void set(const vec2 &v,float z_ = 0.0f) {
		x = v.x; y = v.y; z = z_;
	}
	INLINE void set(const vec3 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const float *v) {
		x = v[0]; y = v[1]; z = v[2];
	}
	INLINE void get(float *v) const {
		v[0] = x; v[1] = y; v[2] = z;
	}
	INLINE float *get() { return v; }
	INLINE const float *get() const { return v; }
	
	INLINE float length2() const {
		return x * x + y * y + z * z;
	}
	INLINE float length() const {
		return Math::sqrt(x * x + y * y + z * z);
	}
	INLINE vec3 &normalize() {
		float ilength = Math::rsqrt(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	INLINE vec3 &normalizeFast() {
		float ilength = Math::rsqrtFast(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float v[4];
		#ifdef USE_SSE
			__m128 vec;
		#elif USE_ALTIVEC
			vec_float4 vec;
		#elif USE_NEON
			float32x4_t vec;
		#endif
	};
};

/*
 */
extern const vec3 vec3_zero;
extern const vec3 vec3_one;
extern const vec3 vec3_epsilon;
extern const vec3 vec3_infinity;

/*
 */
INLINE vec2::vec2(const vec3 &v) : x(v.x), y(v.y) { }

/*
 */
INLINE int operator==(const vec3 &v0,const vec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int operator!=(const vec3 &v0,const vec3 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE vec3 operator*(const vec3 &v0,float v1) {
	return vec3(v0.x * v1,v0.y * v1,v0.z * v1);
}

INLINE vec3 operator*(const vec3 &v0,const vec3 &v1) {
	return vec3(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z);
}

INLINE vec3 operator/(const vec3 &v0,float v1) {
	float iv1 = Math::rcp(v1);
	return vec3(v0.x * iv1,v0.y * iv1,v0.z * iv1);
}

INLINE vec3 operator/(const vec3 &v0,const vec3 &v1) {
	return vec3(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z);
}

INLINE vec3 operator+(const vec3 &v0,const vec3 &v1) {
	return vec3(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z);
}

INLINE vec3 operator-(const vec3 &v0,const vec3 &v1) {
	return vec3(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z);
}

/*
 */
INLINE int compare(const vec3 &v0,const vec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int compare(const vec3 &v0,const vec3 &v1,float epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon) && compare(v0.z,v1.z,epsilon));
}

INLINE float dot(const vec2 &v0,const vec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v1.z;
}

INLINE float dot(const vec3 &v0,const vec2 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z;
}

INLINE float dot(const vec3 &v0,const vec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE vec3 &mul(vec3 &ret,const vec3 &v0,float v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	ret.z = v0.z * v1;
	return ret;
}

INLINE vec3 &mul(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	ret.z = v0.z * v1.z;
	return ret;
}

INLINE vec3 &mad(vec3 &ret,const vec3 &v0,float v1,const vec3 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	ret.z = v0.z * v1 + v2.z;
	return ret;
}

INLINE vec3 &mad(vec3 &ret,const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	ret.z = v0.z * v1.z + v2.z;
	return ret;
}

INLINE vec3 &add(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	ret.z = v0.z + v1.z;
	return ret;
}

INLINE vec3 &sub(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	ret.z = v0.z - v1.z;
	return ret;
}

INLINE vec3 &lerp(vec3 &ret,const vec3 &v0,const vec3 &v1,float k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	return ret;
}

INLINE vec3 &cross(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.x = v0.y * v1.z - v0.z * v1.y;
	ret.y = v0.z * v1.x - v0.x * v1.z;
	ret.z = v0.x * v1.y - v0.y * v1.x;
	return ret;
}

INLINE vec3 &reflect(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	float k = dot(v0,v1) * 2.0f;
	ret.x = v0.x - v1.x * k;
	ret.y = v0.y - v1.y * k;
	ret.z = v0.z - v1.z * k;
	return ret;
}

/*
 */
INLINE float length(const vec3 &v) {
	return v.length();
}

INLINE float length2(const vec3 &v) {
	return v.length2();
}

INLINE vec3 normalize(const vec3 &v) {
	vec3 ret = v;
	return ret.normalize();
}

INLINE vec3 cross(const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	return cross(ret,v0,v1);
}

INLINE vec3 reflect(const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	return reflect(ret,v0,v1);
}

/*
 */
vec3 min(const vec3 &v0,const vec3 &v1);
vec3 max(const vec3 &v0,const vec3 &v1);
vec3 clamp(const vec3 &v,const vec3 &v0,const vec3 &v1);
vec3 saturate(const vec3 &v);
vec3 lerp(const vec3 &v0,const vec3 &v1,float k);

/******************************************************************************\
*
* vec4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) vec4 {
	
	INLINE vec4() { }
	INLINE vec4(const vec3 &v,float w) : x(v.x), y(v.y), z(v.z), w(w) { }
	INLINE vec4(const vec2 &v,float z,float w) : x(v.x), y(v.y), z(z), w(w) { }
	INLINE vec4(float x,float y,float z,float w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE vec4(float v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE vec4(const vec2 &v) : x(v.x), y(v.y), z(0.0f), w(1.0f) { }
	explicit INLINE vec4(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(1.0f) { }
	explicit INLINE vec4(const float *v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
	explicit INLINE vec4(const dvec4 &v);
	explicit INLINE vec4(const hvec4 &v);
	explicit INLINE vec4(const ivec4 &v);
	explicit INLINE vec4(const bvec4 &v);
	#ifdef USE_SSE
		INLINE vec4(const vec4 &v) : vec(v.vec) { }
		explicit INLINE vec4(__m128 v) : vec(v) { }
	#elif USE_ALTIVEC
		INLINE vec4(const vec4 &v) : vec(v.vec) { }
		explicit INLINE vec4(vec_float4 v) : vec(v) { }
	#elif USE_NEON
		INLINE vec4(const vec4 &v) : vec(v.vec) { }
		explicit INLINE vec4(float32x4_t v) : vec(v) { }
	#else
		INLINE vec4(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	#endif
	
	INLINE vec4 &operator=(const vec4 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE vec4 operator-() const {
		return vec4(-x,-y,-z,-w);
	}
	INLINE vec4 &operator*=(float v) {
		x *= v; y *= v; z *= v; w *= v;
		return *this;
	}
	INLINE vec4 &operator*=(const vec4 &v) {
		x *= v.x; y *= v.y; z *= v.z; w *= v.w;
		return *this;
	}
	INLINE vec4 &operator/=(float v) {
		float iv = Math::rcp(v);
		x *= iv; y *= iv; z *= iv; w *= iv;
		return *this;
	}
	INLINE vec4 &operator/=(const vec4 &v) {
		x /= v.x; y /= v.y; z /= v.z; w /= v.w;
		return *this;
	}
	INLINE vec4 &operator+=(const vec4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}
	INLINE vec4 &operator-=(const vec4 &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}
	
	INLINE operator float*() { return v; }
	INLINE operator const float*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 4 && "vec4::operator[](): bad index");
		return v[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 4 && "vec4::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(float v) {
		x = v; y = v; z = v; w = v;
	}
	INLINE void set(float x_,float y_,float z_,float w_) {
		x = x_; y = y_; z = z_; w = w_;
	}
	INLINE void set(const vec2 &v,float z_ = 0.0f,float w_ = 1.0f) {
		x = v.x; y = v.y; z = z_; w = w_;
	}
	INLINE void set(const vec3 &v,float w_ = 1.0f) {
		x = v.x; y = v.y; z = v.z; w = w_;
	}
	INLINE void set(const vec4 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const float *v) {
		x = v[0]; y = v[1]; z = v[2]; w = v[3];
	}
	INLINE void get(float *v) const {
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}
	INLINE float *get() { return v; }
	INLINE const float *get() const { return v; }
	
	INLINE float length2() const {
		return x * x + y * y + z * z + w * w;
	}
	INLINE float length() const {
		return Math::sqrt(x * x + y * y + z * z + w * w);
	}
	INLINE vec4 &normalize() {
		float ilength = Math::rsqrt(x * x + y * y + z * z + w * w);
		x *= ilength; y *= ilength; z *= ilength; w *= ilength;
		return *this;
	}
	INLINE vec4 &normalizeFast() {
		float ilength = Math::rsqrtFast(x * x + y * y + z * z + w * w);
		x *= ilength; y *= ilength; z *= ilength; w *= ilength;
		return *this;
	}
	INLINE vec4 &normalize3() {
		float ilength = Math::rsqrt(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	INLINE vec4 &normalizeFast3() {
		float ilength = Math::rsqrtFast(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float v[4];
		#ifdef USE_SSE
			__m128 vec;
		#elif USE_ALTIVEC
			vec_float4 vec;
		#elif USE_NEON
			float32x4_t vec;
		#endif
	};
};

/*
 */
extern const vec4 vec4_zero;
extern const vec4 vec4_one;
extern const vec4 vec4_epsilon;
extern const vec4 vec4_infinity;

/*
 */
INLINE vec2::vec2(const vec4 &v) : x(v.x), y(v.y) { }
INLINE vec3::vec3(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(0.0f) { }

/*
 */
INLINE int operator==(const vec4 &v0,const vec4 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE int operator!=(const vec4 &v0,const vec4 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE vec4 operator*(const vec4 &v0,float v1) {
	return vec4(v0.x * v1,v0.y * v1,v0.z * v1,v0.w * v1);
}

INLINE vec4 operator*(const vec4 &v0,const vec4 &v1) {
	return vec4(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z,v0.w * v1.w);
}

INLINE vec4 operator/(const vec4 &v0,float v1) {
	float iv1 = Math::rcp(v1);
	return vec4(v0.x * iv1,v0.y * iv1,v0.z * iv1,v0.w * iv1);
}

INLINE vec4 operator/(const vec4 &v0,const vec4 &v1) {
	return vec4(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z,v0.w / v1.w);
}

INLINE vec4 operator+(const vec4 &v0,const vec4 &v1) {
	return vec4(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z,v0.w + v1.w);
}

INLINE vec4 operator-(const vec4 &v0,const vec4 &v1) {
	return vec4(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z,v0.w - v1.w);
}

/*
 */
INLINE int compare(const vec4 &v0,const vec4 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE int compare(const vec4 &v0,const vec4 &v1,float epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon) && compare(v0.z,v1.z,epsilon) && compare(v0.w,v1.w,epsilon));
}

INLINE float dot(const vec3 &v0,const vec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v1.w;
}

INLINE float dot(const vec4 &v0,const vec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w;
}

INLINE float dot(const vec4 &v0,const vec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
}

INLINE float dot3(const vec3 &v0,const vec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE float dot3(const vec4 &v0,const vec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE float dot3(const vec4 &v0,const vec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE vec4 &mul(vec4 &ret,const vec4 &v0,float v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	ret.z = v0.z * v1;
	ret.w = v0.w * v1;
	return ret;
}

INLINE vec4 &mul(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	ret.z = v0.z * v1.z;
	ret.w = v0.w * v1.w;
	return ret;
}

INLINE vec4 &mad(vec4 &ret,const vec4 &v0,float v1,const vec4 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	ret.z = v0.z * v1 + v2.z;
	ret.w = v0.w * v1 + v2.w;
	return ret;
}

INLINE vec4 &mad(vec4 &ret,const vec4 &v0,const vec4 &v1,const vec4 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	ret.z = v0.z * v1.z + v2.z;
	ret.w = v0.w * v1.w + v2.w;
	return ret;
}

INLINE vec4 &add(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	ret.z = v0.z + v1.z;
	ret.w = v0.w + v1.w;
	return ret;
}

INLINE vec4 &sub(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	ret.z = v0.z - v1.z;
	ret.w = v0.w - v1.w;
	return ret;
}

INLINE vec4 &lerp(vec4 &ret,const vec4 &v0,const vec4 &v1,float k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	ret.w = lerp(v0.w,v1.w,k);
	return ret;
}

INLINE vec4 &cross(vec4 &ret,const vec3 &v0,const vec3 &v1) {
	ret.x = v0.y * v1.z - v0.z * v1.y;
	ret.y = v0.z * v1.x - v0.x * v1.z;
	ret.z = v0.x * v1.y - v0.y * v1.x;
	return ret;
}

/*
 */
INLINE float length(const vec4 &v) {
	return v.length();
}

INLINE float length2(const vec4 &v) {
	return v.length2();
}

INLINE vec4 normalize(const vec4 &v) {
	vec4 ret = v;
	return ret.normalize();
}

/*
 */
vec4 min(const vec4 &v0,const vec4 &v1);
vec4 max(const vec4 &v0,const vec4 &v1);
vec4 clamp(const vec4 &v,const vec4 &v0,const vec4 &v1);
vec4 saturate(const vec4 &v);
vec4 lerp(const vec4 &v0,const vec4 &v1,float k);

/******************************************************************************\
*
* dvec2
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) dvec2 {
	
	INLINE dvec2() { }
	INLINE dvec2(double x,double y) : x(x), y(y) { }
	explicit INLINE dvec2(double v) : x(v), y(v) { }
	explicit INLINE dvec2(const dvec3 &v);
	explicit INLINE dvec2(const dvec4 &v);
	explicit INLINE dvec2(const double *v) : x(v[0]), y(v[1]) { }
	explicit INLINE dvec2(const vec2 &v) : x(v.x), y(v.y) { }
	explicit INLINE dvec2(const hvec2 &v);
	explicit INLINE dvec2(const ivec2 &v);
	#ifdef USE_SSE2
		INLINE dvec2(const dvec2 &v) : vec(v.vec) { }
		explicit INLINE dvec2(__m128d v) : vec(v) { }
	#else
		INLINE dvec2(const dvec2 &v) : x(v.x), y(v.y) { }
	#endif
	
	INLINE dvec2 &operator=(const dvec2 &v) {
		#ifdef USE_SSE2
			vec = v.vec;
		#else
			x = v.x; y = v.y;
		#endif
		return *this;
	}
	INLINE dvec2 operator-() const {
		return dvec2(-x,-y);
	}
	INLINE dvec2 &operator*=(double v) {
		x *= v; y *= v;
		return *this;
	}
	INLINE dvec2 &operator*=(const dvec2 &v) {
		x *= v.x; y *= v.y;
		return *this;
	}
	INLINE dvec2 &operator/=(double v) {
		double iv = Math::rcp(v);
		x *= iv; y *= iv;
		return *this;
	}
	INLINE dvec2 &operator/=(const dvec2 &v) {
		x /= v.x; y /= v.y;
		return *this;
	}
	INLINE dvec2 &operator+=(const dvec2 &v) {
		x += v.x; y += v.y;
		return *this;
	}
	INLINE dvec2 &operator-=(const dvec2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}
	
	INLINE operator double*() { return v; }
	INLINE operator const double*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE double &operator[](int i) {
		assert((unsigned int)i < 2 && "dvec2::operator[](): bad index");
		return v[i];
	}
	INLINE double operator[](int i) const {
		assert((unsigned int)i < 2 && "dvec2::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(double v) {
		x = v; y = v;
	}
	INLINE void set(double x_,double y_) {
		x = x_; y = y_;
	}
	INLINE void set(const dvec2 &v) {
		#ifdef USE_SSE2
			vec = v.vec;
		#else
			x = v.x; y = v.y;
		#endif
	}
	INLINE void set(const double *v) {
		x = v[0]; y = v[1];
	}
	INLINE void get(double *v) const {
		v[0] = x; v[1] = y;
	}
	INLINE double *get() { return v; }
	INLINE const double *get() const { return v; }
	
	INLINE double length() const {
		return Math::sqrt(x * x + y * y);
	}
	INLINE double length2() const {
		return x * x + y * y;
	}
	INLINE dvec2 &normalize() {
		double ilength = Math::rsqrt(x * x + y * y);
		x *= ilength; y *= ilength;
		return *this;
	}
	
	union {
		struct {
			double x,y;
		};
		double v[2];
		#ifdef USE_SSE2
			__m128d vec;
		#endif
	};
};

/*
 */
extern const dvec2 dvec2_zero;
extern const dvec2 dvec2_one;
extern const dvec2 dvec2_epsilon;
extern const dvec2 dvec2_infinity;

/*
 */
INLINE vec2::vec2(const dvec2 &v) : x((float)v.x), y((float)v.y) { }

/*
 */
INLINE int operator==(const dvec2 &v0,const dvec2 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE int operator!=(const dvec2 &v0,const dvec2 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE dvec2 operator*(const dvec2 &v0,double v1) {
	return dvec2(v0.x * v1,v0.y * v1);
}

INLINE dvec2 operator*(const dvec2 &v0,const dvec2 &v1) {
	return dvec2(v0.x * v1.x,v0.y * v1.y);
}

INLINE dvec2 operator/(const dvec2 &v0,double v1) {
	double iv1 = Math::rcp(v1);
	return dvec2(v0.x * iv1,v0.y * iv1);
}

INLINE dvec2 operator/(const dvec2 &v0,const dvec2 &v1) {
	return dvec2(v0.x / v1.x,v0.y / v1.y);
}

INLINE dvec2 operator+(const dvec2 &v0,const dvec2 &v1) {
	return dvec2(v0.x + v1.x,v0.y + v1.y);
}

INLINE dvec2 operator-(const dvec2 &v0,const dvec2 &v1) {
	return dvec2(v0.x - v1.x,v0.y - v1.y);
}

/*
 */
INLINE int compare(const dvec2 &v0,const dvec2 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y));
}

INLINE int compare(const dvec2 &v0,const dvec2 &v1,double epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon));
}

INLINE double dot(const dvec2 &v0,const dvec2 &v1) {
	return v0.x * v1.x + v0.y * v1.y;
}

INLINE dvec2 &mul(dvec2 &ret,const dvec2 &v0,double v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	return ret;
}

INLINE dvec2 &mul(dvec2 &ret,const dvec2 &v0,const dvec2 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	return ret;
}

INLINE dvec2 &mad(dvec2 &ret,const dvec2 &v0,double v1,const dvec2 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	return ret;
}

INLINE dvec2 &mad(dvec2 &ret,const dvec2 &v0,const dvec2 &v1,const dvec2 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	return ret;
}

INLINE dvec2 &add(dvec2 &ret,const dvec2 &v0,const dvec2 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	return ret;
}

INLINE dvec2 &sub(dvec2 &ret,const dvec2 &v0,const dvec2 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	return ret;
}

INLINE dvec2 &lerp(dvec2 &ret,const dvec2 &v0,const dvec2 &v1,double k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	return ret;
}

/*
 */
INLINE double length(const dvec2 &v) {
	return v.length();
}

INLINE double length2(const dvec2 &v) {
	return v.length2();
}

INLINE dvec2 normalize(const dvec2 &v) {
	dvec2 ret = v;
	return ret.normalize();
}

/*
 */
dvec2 min(const dvec2 &v0,const dvec2 &v1);
dvec2 max(const dvec2 &v0,const dvec2 &v1);
dvec2 clamp(const dvec2 &v,const dvec2 &v0,const dvec2 &v1);
dvec2 saturate(const dvec2 &v);
dvec2 lerp(const dvec2 &v0,const dvec2 &v1,double k);

/******************************************************************************\
*
* dvec3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) dvec3 {
	
	INLINE dvec3() { }
	INLINE dvec3(const dvec3 &v) : x(v.x), y(v.y), z(v.z), w(0.0) { }
	INLINE dvec3(const dvec2 &v,double z) : x(v.x), y(v.y), z(z), w(0.0) { }
	INLINE dvec3(double x,double y,double z) : x(x), y(y), z(z), w(0.0) { }
	explicit INLINE dvec3(double v) : x(v), y(v), z(v), w(0.0) { }
	explicit INLINE dvec3(const dvec2 &v) : x(v.x), y(v.y), z(0.0), w(0.0) { }
	explicit INLINE dvec3(const dvec4 &v);
	explicit INLINE dvec3(const double *v) : x(v[0]), y(v[1]), z(v[2]), w(0.0) { }
	explicit INLINE dvec3(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(0.0) { }
	explicit INLINE dvec3(const hvec3 &v);
	explicit INLINE dvec3(const ivec3 &v);
	
	INLINE dvec3 &operator=(const dvec3 &v) {
		#ifdef USE_SSE2
			vec0 = v.vec0; vec1 = v.vec1;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE dvec3 operator-() const {
		return dvec3(-x,-y,-z);
	}
	INLINE dvec3 &operator*=(double v) {
		x *= v; y *= v; z *= v;
		return *this;
	}
	INLINE dvec3 &operator*=(const dvec3 &v) {
		x *= v.x; y *= v.y; z *= v.z;
		return *this;
	}
	INLINE dvec3 &operator/=(double v) {
		double iv = Math::rcp(v);
		x *= iv; y *= iv; z *= iv;
		return *this;
	}
	INLINE dvec3 &operator/=(const dvec3 &v) {
		x /= v.x; y /= v.y; z /= v.z;
		return *this;
	}
	INLINE dvec3 &operator+=(const dvec3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	INLINE dvec3 &operator-=(const dvec3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	
	INLINE operator double*() { return v; }
	INLINE operator const double*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE double &operator[](int i) {
		assert((unsigned int)i < 3 && "dvec3::operator[](): bad index");
		return v[i];
	}
	INLINE double operator[](int i) const {
		assert((unsigned int)i < 3 && "dvec3::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(double v) {
		x = v; y = v; z = v;
	}
	INLINE void set(double x_,double y_,double z_) {
		x = x_; y = y_; z = z_;
	}
	INLINE void set(const dvec2 &v,double z_ = 0.0) {
		x = v.x; y = v.y; z = z_;
	}
	INLINE void set(const dvec3 &v) {
		#ifdef USE_SSE2
			vec0 = v.vec0; vec1 = v.vec1;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const double *v) {
		x = v[0]; y = v[1]; z = v[2];
	}
	INLINE void get(double *v) const {
		v[0] = x; v[1] = y; v[2] = z;
	}
	INLINE double *get() { return v; }
	INLINE const double *get() const { return v; }
	
	INLINE double length() const {
		return Math::sqrt(x * x + y * y + z * z);
	}
	INLINE double length2() const {
		return x * x + y * y + z * z;
	}
	INLINE dvec3 &normalize() {
		double ilength = Math::rsqrt(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	
	union {
		struct {
			double x,y,z,w;
		};
		double v[4];
		#ifdef USE_SSE2
			struct {
				__m128d vec0;
				__m128d vec1;
			};
		#endif
	};
};

/*
 */
extern const dvec3 dvec3_zero;
extern const dvec3 dvec3_one;
extern const dvec3 dvec3_epsilon;
extern const dvec3 dvec3_infinity;

/*
 */
INLINE dvec2::dvec2(const dvec3 &v) : x(v.x), y(v.y) { }
INLINE vec3::vec3(const dvec3 &v) : x((float)v.x), y((float)v.y), z((float)v.z), w(0.0f) { }

/*
 */
INLINE int operator==(const dvec3 &v0,const dvec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int operator!=(const dvec3 &v0,const dvec3 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE dvec3 operator*(const dvec3 &v0,double v1) {
	return dvec3(v0.x * v1,v0.y * v1,v0.z * v1);
}

INLINE dvec3 operator*(const dvec3 &v0,const dvec3 &v1) {
	return dvec3(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z);
}

INLINE dvec3 operator/(const dvec3 &v0,double v1) {
	double iv1 = Math::rcp(v1);
	return dvec3(v0.x * iv1,v0.y * iv1,v0.z * iv1);
}

INLINE dvec3 operator/(const dvec3 &v0,const dvec3 &v1) {
	return dvec3(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z);
}

INLINE dvec3 operator+(const dvec3 &v0,const dvec3 &v1) {
	return dvec3(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z);
}

INLINE dvec3 operator-(const dvec3 &v0,const dvec3 &v1) {
	return dvec3(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z);
}

/*
 */
INLINE int compare(const dvec3 &v0,const dvec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int compare(const dvec3 &v0,const dvec3 &v1,double epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon) && compare(v0.z,v1.z,epsilon));
}

INLINE double dot(const dvec2 &v0,const dvec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v1.z;
}

INLINE double dot(const dvec3 &v0,const dvec2 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z;
}

INLINE double dot(const dvec3 &v0,const dvec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE dvec3 &mul(dvec3 &ret,const dvec3 &v0,double v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	ret.z = v0.z * v1;
	return ret;
}

INLINE dvec3 &mul(dvec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	ret.z = v0.z * v1.z;
	return ret;
}

INLINE dvec3 &mad(dvec3 &ret,const dvec3 &v0,double v1,const dvec3 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	ret.z = v0.z * v1 + v2.z;
	return ret;
}

INLINE dvec3 &mad(dvec3 &ret,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	ret.z = v0.z * v1.z + v2.z;
	return ret;
}

INLINE vec3 &add(vec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = (float)(v0.x + v1.x);
	ret.y = (float)(v0.y + v1.y);
	ret.z = (float)(v0.z + v1.z);
	return ret;
}

INLINE dvec3 &add(dvec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	ret.z = v0.z + v1.z;
	return ret;
}

INLINE vec3 &sub(vec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = (float)(v0.x - v1.x);
	ret.y = (float)(v0.y - v1.y);
	ret.z = (float)(v0.z - v1.z);
	return ret;
}

INLINE dvec3 &sub(dvec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	ret.z = v0.z - v1.z;
	return ret;
}

INLINE dvec3 &lerp(dvec3 &ret,const dvec3 &v0,const dvec3 &v1,double k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	return ret;
}

INLINE dvec3 &cross(dvec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = v0.y * v1.z - v0.z * v1.y;
	ret.y = v0.z * v1.x - v0.x * v1.z;
	ret.z = v0.x * v1.y - v0.y * v1.x;
	return ret;
}

INLINE dvec3 &reflect(dvec3 &ret,const dvec3 &v0,const dvec3 &v1) {
	double k = 2.0 * (v0.x * v1.x + v0.y * v1.y + v0.z * v1.z);
	ret.x = v0.x - v1.x * k;
	ret.y = v0.y - v1.y * k;
	ret.z = v0.z - v1.z * k;
	return ret;
}

/*
 */
INLINE double length(const dvec3 &v) {
	return v.length();
}

INLINE double length2(const dvec3 &v) {
	return v.length2();
}

INLINE dvec3 normalize(const dvec3 &v) {
	dvec3 ret = v;
	return ret.normalize();
}

INLINE dvec3 cross(const dvec3 &v0,const dvec3 &v1) {
	dvec3 ret;
	return cross(ret,v0,v1);
}

INLINE dvec3 reflect(const dvec3 &v0,const dvec3 &v1) {
	dvec3 ret;
	return reflect(ret,v0,v1);
}

/*
 */
dvec3 min(const dvec3 &v0,const dvec3 &v1);
dvec3 max(const dvec3 &v0,const dvec3 &v1);
dvec3 clamp(const dvec3 &v,const dvec3 &v0,const dvec3 &v1);
dvec3 saturate(const dvec3 &v);
dvec3 lerp(const dvec3 &v0,const dvec3 &v1,double k);

/******************************************************************************\
*
* dvec4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) dvec4 {

	INLINE dvec4() { }
	INLINE dvec4(const dvec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	INLINE dvec4(const dvec3 &v,double w) : x(v.x), y(v.y), z(v.z), w(w) { }
	INLINE dvec4(const dvec2 &v,double z,double w) : x(v.x), y(v.y), z(z), w(w) { }
	INLINE dvec4(double x,double y,double z,double w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE dvec4(double v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE dvec4(const dvec2 &v) : x(v.x), y(v.y), z(0.0), w(1.0) { }
	explicit INLINE dvec4(const dvec3 &v) : x(v.x), y(v.y), z(v.z), w(1.0) { }
	explicit INLINE dvec4(const double *v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
	explicit INLINE dvec4(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	explicit INLINE dvec4(const hvec4 &v);
	explicit INLINE dvec4(const ivec4 &v);
	explicit INLINE dvec4(const bvec4 &v);
	
	INLINE dvec4 &operator=(const dvec4 &v) {
		#ifdef USE_SSE2
			vec0 = v.vec0; vec1 = v.vec1;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE dvec4 operator-() const {
		return dvec4(-x,-y,-z,-w);
	}
	INLINE dvec4 &operator*=(double v) {
		x *= v; y *= v; z *= v; w *= v;
		return *this;
	}
	INLINE dvec4 &operator*=(const dvec4 &v) {
		x *= v.x; y *= v.y; z *= v.z; w *= v.w;
		return *this;
	}
	INLINE dvec4 &operator/=(double v) {
		double iv = Math::rcp(v);
		x *= iv; y *= iv; z *= iv; w *= iv;
		return *this;
	}
	INLINE dvec4 &operator/=(const dvec4 &v) {
		x /= v.x; y /= v.y; z /= v.z; w /= v.w;
		return *this;
	}
	INLINE dvec4 &operator+=(const dvec4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}
	INLINE dvec4 &operator-=(const dvec4 &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}
	
	INLINE operator double*() { return v; }
	INLINE operator const double*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE double &operator[](int i) {
		assert((unsigned int)i < 4 && "dvec4::operator[](): bad index");
		return v[i];
	}
	INLINE double operator[](int i) const {
		assert((unsigned int)i < 4 && "dvec4::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(double v) {
		x = v; y = v; z = v; w = v;
	}
	INLINE void set(double x_,double y_,double z_,double w_) {
		x = x_; y = y_; z = z_; w = w_;
	}
	INLINE void set(const dvec2 &v,double z_ = 0.0,double w_ = 1.0) {
		x = v.x; y = v.y; z = z_; w = w_;
	}
	INLINE void set(const dvec3 &v,double w_ = 1.0) {
		x = v.x; y = v.y; z = v.z; w = w_;
	}
	INLINE void set(const dvec4 &v) {
		#ifdef USE_SSE2
			vec0 = v.vec0; vec1 = v.vec1;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const double *v) {
		x = v[0]; y = v[1]; z = v[2]; w = v[3];
	}
	INLINE void get(double *v) const {
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}
	INLINE double *get() { return v; }
	INLINE const double *get() const { return v; }
	
	INLINE double length() const {
		return Math::sqrt(x * x + y * y + z * z + w * w);
	}
	INLINE double length2() const {
		return x * x + y * y + z * z + w * w;
	}
	INLINE dvec4 &normalize() {
		double ilength = Math::rsqrt(x * x + y * y + z * z + w * w);
		x *= ilength; y *= ilength; z *= ilength; w *= ilength;
		return *this;
	}
	INLINE dvec4 &normalize3() {
		double ilength = Math::rsqrt(x * x + y * y + z * z);
		x *= ilength; y *= ilength; z *= ilength;
		return *this;
	}
	
	union {
		struct {
			double x,y,z,w;
		};
		double v[4];
		#ifdef USE_SSE2
			struct {
				__m128d vec0;
				__m128d vec1;
			};
		#endif
	};
};

/*
 */
extern const dvec4 dvec4_zero;
extern const dvec4 dvec4_one;
extern const dvec4 dvec4_epsilon;
extern const dvec4 dvec4_infinity;

/*
 */
INLINE dvec2::dvec2(const dvec4 &v) : x(v.x), y(v.y) { }
INLINE dvec3::dvec3(const dvec4 &v) : x(v.x), y(v.y), z(v.z), w(0.0) { }
INLINE vec4::vec4(const dvec4 &v) : x((float)v.x), y((float)v.y), z((float)v.z), w((float)v.w) { }

/*
 */
INLINE int operator==(const dvec4 &v0,const dvec4 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE int operator!=(const dvec4 &v0,const dvec4 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE dvec4 operator*(const dvec4 &v0,double v1) {
	return dvec4(v0.x * v1,v0.y * v1,v0.z * v1,v0.w * v1);
}

INLINE dvec4 operator*(const dvec4 &v0,const dvec4 &v1) {
	return dvec4(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z,v0.w * v1.w);
}

INLINE dvec4 operator/(const dvec4 &v0,double v1) {
	double iv1 = Math::rcp(v1);
	return dvec4(v0.x * iv1,v0.y * iv1,v0.z * iv1,v0.w * iv1);
}

INLINE dvec4 operator/(const dvec4 &v0,const dvec4 &v1) {
	return dvec4(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z,v0.w / v1.w);
}

INLINE dvec4 operator+(const dvec4 &v0,const dvec4 &v1) {
	return dvec4(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z,v0.w + v1.w);
}

INLINE dvec4 operator-(const dvec4 &v0,const dvec4 &v1) {
	return dvec4(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z,v0.w - v1.w);
}

/*
 */
INLINE int compare(const dvec4 &v0,const dvec4 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z) && compare(v0.w,v1.w));
}

INLINE int compare(const dvec4 &v0,const dvec4 &v1,double epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon) && compare(v0.z,v1.z,epsilon) && compare(v0.w,v1.w,epsilon));
}

INLINE double dot(const dvec3 &v0,const dvec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v1.w;
}

INLINE double dot(const dvec4 &v0,const dvec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w;
}

INLINE double dot(const dvec4 &v0,const dvec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
}

INLINE double dot3(const dvec3 &v0,const dvec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE double dot3(const dvec4 &v0,const dvec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE double dot3(const dvec4 &v0,const dvec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE dvec4 &mul(dvec4 &ret,const dvec4 &v0,double v1) {
	ret.x = v0.x * v1;
	ret.y = v0.y * v1;
	ret.z = v0.z * v1;
	ret.w = v0.w * v1;
	return ret;
}

INLINE dvec4 &mul(dvec4 &ret,const dvec4 &v0,const dvec4 &v1) {
	ret.x = v0.x * v1.x;
	ret.y = v0.y * v1.y;
	ret.z = v0.z * v1.z;
	ret.w = v0.w * v1.w;
	return ret;
}

INLINE dvec4 &mad(dvec4 &ret,const dvec4 &v0,double v1,const dvec4 &v2) {
	ret.x = v0.x * v1 + v2.x;
	ret.y = v0.y * v1 + v2.y;
	ret.z = v0.z * v1 + v2.z;
	ret.w = v0.w * v1 + v2.w;
	return ret;
}

INLINE dvec4 &mad(dvec4 &ret,const dvec4 &v0,const dvec4 &v1,const dvec4 &v2) {
	ret.x = v0.x * v1.x + v2.x;
	ret.y = v0.y * v1.y + v2.y;
	ret.z = v0.z * v1.z + v2.z;
	ret.w = v0.w * v1.w + v2.w;
	return ret;
}

INLINE dvec4 &add(dvec4 &ret,const dvec4 &v0,const dvec4 &v1) {
	ret.x = v0.x + v1.x;
	ret.y = v0.y + v1.y;
	ret.z = v0.z + v1.z;
	ret.w = v0.w + v1.w;
	return ret;
}

INLINE dvec4 &sub(dvec4 &ret,const dvec4 &v0,const dvec4 &v1) {
	ret.x = v0.x - v1.x;
	ret.y = v0.y - v1.y;
	ret.z = v0.z - v1.z;
	ret.w = v0.w - v1.w;
	return ret;
}

INLINE dvec4 &lerp(dvec4 &ret,const dvec4 &v0,const dvec4 &v1,double k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	ret.w = lerp(v0.w,v1.w,k);
	return ret;
}

INLINE dvec4 &cross(dvec4 &ret,const dvec3 &v0,const dvec3 &v1) {
	ret.x = v0.y * v1.z - v0.z * v1.y;
	ret.y = v0.z * v1.x - v0.x * v1.z;
	ret.z = v0.x * v1.y - v0.y * v1.x;
	return ret;
}

/*
 */
INLINE double length(const dvec4 &v) {
	return v.length();
}

INLINE double length2(const dvec4 &v) {
	return v.length2();
}

INLINE dvec4 normalize(const dvec4 &v) {
	dvec4 ret = v;
	return ret.normalize();
}

/*
 */
dvec4 min(const dvec4 &v0,const dvec4 &v1);
dvec4 max(const dvec4 &v0,const dvec4 &v1);
dvec4 clamp(const dvec4 &v,const dvec4 &v0,const dvec4 &v1);
dvec4 saturate(const dvec4 &v);
dvec4 lerp(const dvec4 &v0,const dvec4 &v1,double k);

/******************************************************************************\
*
* hvec2
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED4(struct) hvec2 {
	
	INLINE hvec2() { }
	INLINE hvec2(const hvec2 &v) : x(v.x), y(v.y) { }
	INLINE hvec2(half x,half y) : x(x), y(y) { }
	explicit INLINE hvec2(half v) : x(v), y(v) { }
	explicit INLINE hvec2(float v) : x(v), y(v) { }
	explicit INLINE hvec2(const vec2 &v) : x(v.x), y(v.y) { }
	explicit INLINE hvec2(const dvec2 &v) : x((float)v.x), y((float)v.y) { }
	
	INLINE hvec2 &operator=(const hvec2 &v) {
		x = v.x; y = v.y;
		return *this;
	}
	
	INLINE operator half*() { return &x; }
	INLINE operator const half*() const { return &x; }
	INLINE operator void*() { return &x; }
	INLINE operator const void*() const { return &x; }
	
	INLINE half &operator[](int i) {
		assert((unsigned int)i < 2 && "hvec2::operator[](): bad index");
		return (&x)[i];
	}
	INLINE half operator[](int i) const {
		assert((unsigned int)i < 2 && "hvec2::operator[](): bad index");
		return (&x)[i];
	}
	
	half x,y;
};

/*
 */
extern const hvec2 hvec2_zero;
extern const hvec2 hvec2_one;

/*
 */
INLINE vec2::vec2(const hvec2 &v) : x(v.x.getFloat()), y(v.y.getFloat()) { }
INLINE dvec2::dvec2(const hvec2 &v) : x(v.x.getFloat()), y(v.y.getFloat()) { }

/******************************************************************************\
*
* hvec3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED8(struct) hvec3 {
	
	INLINE hvec3() { }
	INLINE hvec3(const hvec3 &v) : x(v.x), y(v.y), z(v.z) { }
	INLINE hvec3(half x,half y,half z) : x(x), y(y), z(z) { }
	explicit INLINE hvec3(half v) : x(v), y(v), z(v) { }
	explicit INLINE hvec3(float v) : x(v), y(v), z(v) { }
	explicit INLINE hvec3(const vec3 &v) : x(v.x), y(v.y), z(v.z) { }
	explicit INLINE hvec3(const dvec3 &v) : x((float)v.x), y((float)v.y), z((float)v.z) { }
	
	INLINE hvec3 &operator=(const hvec3 &v) {
		x = v.x; y = v.y; z = v.z; w = v.w;
		return *this;
	}
	
	INLINE operator half*() { return &x; }
	INLINE operator const half*() const { return &x; }
	INLINE operator void*() { return &x; }
	INLINE operator const void*() const { return &x; }
	
	INLINE half &operator[](int i) {
		assert((unsigned int)i < 3 && "hvec3::operator[](): bad index");
		return (&x)[i];
	}
	INLINE half operator[](int i) const {
		assert((unsigned int)i < 3 && "hvec3::operator[](): bad index");
		return (&x)[i];
	}
	
	half x,y,z,w;
};

/*
 */
extern const hvec3 hvec3_zero;
extern const hvec3 hvec3_one;

/*
 */
INLINE vec3::vec3(const hvec3 &v) : x(v.x.getFloat()), y(v.y.getFloat()), z(v.z.getFloat()), w(0.0f) { }
INLINE dvec3::dvec3(const hvec3 &v) : x(v.x.getFloat()), y(v.y.getFloat()), z(v.z.getFloat()), w(0.0) { }

/******************************************************************************\
*
* hvec4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED8(struct) hvec4 {
	
	INLINE hvec4() { }
	INLINE hvec4(const hvec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	INLINE hvec4(half x,half y,half z,half w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE hvec4(half v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE hvec4(float v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE hvec4(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	explicit INLINE hvec4(const dvec4 &v) : x((float)v.x), y((float)v.y), z((float)v.z), w((float)v.w) { }
	
	INLINE hvec4 &operator=(const hvec4 &v) {
		x = v.x; y = v.y; z = v.z; w = v.w;
		return *this;
	}
	
	INLINE operator half*() { return &x; }
	INLINE operator const half*() const { return &x; }
	INLINE operator void*() { return &x; }
	INLINE operator const void*() const { return &x; }
	
	INLINE half &operator[](int i) {
		assert((unsigned int)i < 4 && "hvec4::operator[](): bad index");
		return (&x)[i];
	}
	INLINE half operator[](int i) const {
		assert((unsigned int)i < 4 && "hvec4::operator[](): bad index");
		return (&x)[i];
	}
	
	half x,y,z,w;
};

/*
 */
extern const hvec4 hvec4_zero;
extern const hvec4 hvec4_one;

/*
 */
INLINE vec4::vec4(const hvec4 &v) : x(v.x.getFloat()), y(v.y.getFloat()), z(v.z.getFloat()), w(v.w.getFloat()) { }
INLINE dvec4::dvec4(const hvec4 &v) : x(v.x.getFloat()), y(v.y.getFloat()), z(v.z.getFloat()), w(v.w.getFloat()) { }

/******************************************************************************\
*
* ivec2
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED8(struct) ivec2 {
	
	INLINE ivec2() { }
	INLINE ivec2(const ivec2 &v) : x(v.x), y(v.y) { }
	INLINE ivec2(int x,int y) : x(x), y(y) { }
	explicit INLINE ivec2(int v) : x(v), y(v) { }
	explicit INLINE ivec2(const int *v) : x(v[0]), y(v[1]) { }
	explicit INLINE ivec2(const vec2 &v) : x(Math::ftoi(v.x)), y(Math::ftoi(v.y)) { }
	explicit INLINE ivec2(const dvec2 &v) : x(Math::dtoi(v.x)), y(Math::dtoi(v.y)) { }
	#ifdef USE_SSE
		explicit INLINE ivec2(__m64 v) : vec(v) { }
	#elif USE_NEON
		explicit INLINE ivec2(int32x2_t v) : vec(v) { }
	#endif
	
	INLINE ivec2 &operator=(const ivec2 &v) {
		x = v.x; y = v.y;
		return *this;
	}
	
	INLINE ivec2 operator-() const {
		return ivec2(-x,-y);
	}
	INLINE ivec2 &operator*=(int v) {
		x *= v; y *= v;
		return *this;
	}
	INLINE ivec2 &operator*=(const ivec2 &v) {
		x *= v.x; y *= v.y;
		return *this;
	}
	INLINE ivec2 &operator/=(int v) {
		x /= v; y /= v;
		return *this;
	}
	INLINE ivec2 &operator/=(const ivec2 &v) {
		x /= v.x; y /= v.y;
		return *this;
	}
	INLINE ivec2 &operator+=(const ivec2 &v) {
		x += v.x; y += v.y;
		return *this;
	}
	INLINE ivec2 &operator-=(const ivec2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}
	INLINE ivec2 &operator<<=(int v) {
		x <<= v; y <<= v;
		return *this;
	}
	INLINE ivec2 &operator>>=(int v) {
		x >>= v; y >>= v;
		return *this;
	}
	
	INLINE operator int*() { return v; }
	INLINE operator const int*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE int &operator[](int i) {
		assert((unsigned int)i < 2 && "ivec2::operator[](): bad index");
		return v[i];
	}
	INLINE int operator[](int i) const {
		assert((unsigned int)i < 2 && "ivec2::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(int v) {
		x = v; y = v;
	}
	INLINE void set(int x_,int y_) {
		x = x_; y = y_;
	}
	INLINE void set(const ivec2 &v) {
		x = v.x; y = v.y;
	}
	INLINE void set(const int *v) {
		x = v[0]; y = v[1];
	}
	INLINE void get(int *v) const {
		v[0] = x; v[1] = y;
	}
	INLINE int *get() { return v; }
	INLINE const int *get() const { return v; }
	
	INLINE int length2() const {
		return x * x + y * y;
	}
	
	union {
		struct {
			int x,y;
		};
		int v[2];
		#ifdef USE_SSE
			__m64 vec;
		#elif USE_NEON
			int32x2_t vec;
		#endif
	};
};

/*
 */
extern const ivec2 ivec2_zero;
extern const ivec2 ivec2_one;

/*
 */
INLINE vec2::vec2(const ivec2 &v) : x(Math::itof(v.x)), y(Math::itof(v.y)) { }
INLINE dvec2::dvec2(const ivec2 &v) : x(Math::itod(v.x)), y(Math::itod(v.y)) { }

/*
 */
INLINE int operator==(const ivec2 &v0,const ivec2 &v1) {
	return (v0.x == v1.x && v0.y == v1.y);
}

INLINE int operator!=(const ivec2 &v0,const ivec2 &v1) {
	return (v0.x != v1.x || v0.y != v1.y);
}

INLINE ivec2 operator*(const ivec2 &v0,int v1) {
	return ivec2(v0.x * v1,v0.y * v1);
}

INLINE ivec2 operator*(const ivec2 &v0,const ivec2 &v1) {
	return ivec2(v0.x * v1.x,v0.y * v1.y);
}

INLINE ivec2 operator/(const ivec2 &v0,int v1) {
	return ivec2(v0.x / v1,v0.y / v1);
}

INLINE ivec2 operator/(const ivec2 &v0,const ivec2 &v1) {
	return ivec2(v0.x / v1.x,v0.y / v1.y);
}

INLINE ivec2 operator+(const ivec2 &v0,const ivec2 &v1) {
	return ivec2(v0.x + v1.x,v0.y + v1.y);
}

INLINE ivec2 operator-(const ivec2 &v0,const ivec2 &v1) {
	return ivec2(v0.x - v1.x,v0.y - v1.y);
}

INLINE ivec2 operator<<(const ivec2 &v0,int v1) {
	return ivec2(v0.x << v1,v0.y << v1);
}

INLINE ivec2 operator>>(const ivec2 &v0,int v1) {
	return ivec2(v0.x >> v1,v0.y >> v1);
}

/*
 */
INLINE int dot(const ivec2 &v0,const ivec2 &v1) {
	return v0.x * v1.x + v0.y * v1.y;
}

INLINE ivec2 &lerp(ivec2 &ret,const ivec2 &v0,const ivec2 &v1,int k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	return ret;
}

/*
 */
INLINE int length2(const ivec2 &v) {
	return v.length2();
}

/*
 */
ivec2 min(const ivec2 &v0,const ivec2 &v1);
ivec2 max(const ivec2 &v0,const ivec2 &v1);
ivec2 clamp(const ivec2 &v,const ivec2 &v0,const ivec2 &v1);
ivec2 lerp(const ivec2 &v0,const ivec2 &v1,int k);

/******************************************************************************\
*
* ivec3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) ivec3 {
	
	INLINE ivec3() { }
	INLINE ivec3(int x,int y,int z) : x(x), y(y), z(z), w(0) { }
	explicit INLINE ivec3(int v) : x(v), y(v), z(v), w(0) { }
	explicit INLINE ivec3(const int *v) : x(v[0]), y(v[1]), z(v[2]) { }
	explicit INLINE ivec3(const vec3 &v) : x(Math::ftoi(v.x)), y(Math::ftoi(v.y)), z(Math::ftoi(v.z)) { }
	explicit INLINE ivec3(const dvec3 &v) : x(Math::dtoi(v.x)), y(Math::dtoi(v.y)), z(Math::dtoi(v.z)) { }
	#ifdef USE_SSE2
		INLINE ivec3(const ivec3 &v) : vec(v.vec) { }
		explicit INLINE ivec3(__m128i v) : vec(v) { }
	#elif USE_SSE
		INLINE ivec3(const ivec3 &v) : vec(v.vec) { }
		explicit INLINE ivec3(__m128 v) : vec(v) { }
	#elif USE_ALTIVEC
		INLINE ivec3(const ivec3 &v) : vec(v.vec) { }
		explicit INLINE ivec3(vec_int4 v) : vec(v) { }
	#elif USE_NEON
		INLINE ivec3(const ivec3 &v) : vec(v.vec) { }
		explicit INLINE ivec3(int32x4_t v) : vec(v) { }
	#else
		INLINE ivec3(const ivec3 &v) : x(v.x), y(v.y), z(v.z), w(0) { }
	#endif
	
	INLINE ivec3 &operator=(const ivec3 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE ivec3 operator-() const {
		return ivec3(-x,-y,-z);
	}
	INLINE ivec3 &operator*=(int v) {
		x *= v; y *= v; z *= v;
		return *this;
	}
	INLINE ivec3 &operator*=(const ivec3 &v) {
		x *= v.x; y *= v.y; z *= v.z;
		return *this;
	}
	INLINE ivec3 &operator/=(int v) {
		x /= v; y /= v; z /= v;
		return *this;
	}
	INLINE ivec3 &operator/=(const ivec3 &v) {
		x /= v.x; y /= v.y; z /= v.z;
		return *this;
	}
	INLINE ivec3 &operator+=(const ivec3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	INLINE ivec3 &operator-=(const ivec3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	INLINE ivec3 &operator<<=(int v) {
		x <<= v; y <<= v; z <<= v;
		return *this;
	}
	INLINE ivec3 &operator>>=(int v) {
		x >>= v; y >>= v; z >>= v;
		return *this;
	}
	
	INLINE operator int*() { return v; }
	INLINE operator const int*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE int &operator[](int i) {
		assert((unsigned int)i < 3 && "ivec3::operator[](): bad index");
		return v[i];
	}
	INLINE int operator[](int i) const {
		assert((unsigned int)i < 3 && "ivec3::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(int v) {
		x = v; y = v; z = v;
	}
	INLINE void set(int x_,int y_,int z_) {
		x = x_; y = y_; z = z_;
	}
	INLINE void set(const ivec3 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const int *v) {
		x = v[0]; y = v[1]; z = v[2];
	}
	INLINE void get(int *v) const {
		v[0] = x; v[1] = y; v[2] = z;
	}
	INLINE int *get() { return v; }
	INLINE const int *get() const { return v; }
	
	INLINE int length2() const {
		return x * x + y * y + z * z;
	}
	
	union {
		struct {
			int x,y,z,w;
		};
		int v[4];
		#ifdef USE_SSE2
			__m128i vec;
		#elif USE_SSE
			__m128 vec;
		#elif USE_ALTIVEC
			vec_int4 vec;
		#elif USE_NEON
			int32x4_t vec;
		#endif
	};
};

/*
 */
extern const ivec3 ivec3_zero;
extern const ivec3 ivec3_one;

/*
 */
INLINE vec3::vec3(const ivec3 &v) : x(Math::itof(v.x)), y(Math::itof(v.y)), z(Math::itof(v.z)), w(0.0f) { }
INLINE dvec3::dvec3(const ivec3 &v) : x(Math::itod(v.x)), y(Math::itod(v.y)), z(Math::itod(v.z)), w(0.0) { }

/*
 */
INLINE int operator==(const ivec3 &v0,const ivec3 &v1) {
	return (v0.x == v1.x && v0.y == v1.y && v0.z == v1.z);
}

INLINE int operator!=(const ivec3 &v0,const ivec3 &v1) {
	return (v0.x != v1.x || v0.y != v1.y || v0.z != v1.z);
}

INLINE ivec3 operator*(const ivec3 &v0,int v1) {
	return ivec3(v0.x * v1,v0.y * v1,v0.z * v1);
}

INLINE ivec3 operator*(const ivec3 &v0,const ivec3 &v1) {
	return ivec3(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z);
}

INLINE ivec3 operator/(const ivec3 &v0,int v1) {
	return ivec3(v0.x / v1,v0.y / v1,v0.z / v1);
}

INLINE ivec3 operator/(const ivec3 &v0,const ivec3 &v1) {
	return ivec3(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z);
}

INLINE ivec3 operator+(const ivec3 &v0,const ivec3 &v1) {
	return ivec3(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z);
}

INLINE ivec3 operator-(const ivec3 &v0,const ivec3 &v1) {
	return ivec3(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z);
}

INLINE ivec3 operator<<(const ivec3 &v0,int v1) {
	return ivec3(v0.x << v1,v0.y << v1,v0.z << v1);
}

INLINE ivec3 operator>>(const ivec3 &v0,int v1) {
	return ivec3(v0.x >> v1,v0.y >> v1,v0.z >> v1);
}

/*
 */
INLINE int dot(const ivec3 &v0,const ivec3 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

INLINE ivec3 &lerp(ivec3 &ret,const ivec3 &v0,const ivec3 &v1,int k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	return ret;
}

INLINE ivec3 &cross(ivec3 &ret,const ivec3 &v0,const ivec3 &v1) {
	ret.x = v0.y * v1.z - v0.z * v1.y;
	ret.y = v0.z * v1.x - v0.x * v1.z;
	ret.z = v0.x * v1.y - v0.y * v1.x;
	return ret;
}

/*
 */
INLINE int length2(const ivec3 &v) {
	return v.length2();
}

INLINE ivec3 cross(const ivec3 &v0,const ivec3 &v1) {
	ivec3 ret;
	return cross(ret,v0,v1);
}

/*
 */
ivec3 min(const ivec3 &v0,const ivec3 &v1);
ivec3 max(const ivec3 &v0,const ivec3 &v1);
ivec3 clamp(const ivec3 &v,const ivec3 &v0,const ivec3 &v1);
ivec3 lerp(const ivec3 &v0,const ivec3 &v1,int k);

/******************************************************************************\
*
* ivec4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) ivec4 {
	
	INLINE ivec4() { }
	INLINE ivec4(int x,int y,int z,int w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE ivec4(int v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE ivec4(const int *v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
	explicit INLINE ivec4(const vec4 &v) : x(Math::ftoi(v.x)), y(Math::ftoi(v.y)), z(Math::ftoi(v.z)), w(Math::ftoi(v.w)) { }
	explicit INLINE ivec4(const dvec4 &v) : x(Math::dtoi(v.x)), y(Math::dtoi(v.y)), z(Math::dtoi(v.z)), w(Math::dtoi(v.w)) { }
	explicit INLINE ivec4(const bvec4 &v);
	#ifdef USE_SSE2
		INLINE ivec4(const ivec4 &v) : vec(v.vec) { }
		explicit INLINE ivec4(__m128i v) : vec(v) { }
	#elif USE_SSE
		INLINE ivec4(const ivec4 &v) : vec(v.vec) { }
		explicit INLINE ivec4(__m128 v) : vec(v) { }
	#elif USE_ALTIVEC
		INLINE ivec4(const ivec4 &v) : vec(v.vec) { }
		explicit INLINE ivec4(vec_int4 v) : vec(v) { }
	#elif USE_NEON
		INLINE ivec4(const ivec4 &v) : vec(v.vec) { }
		explicit INLINE ivec4(int32x4_t v) : vec(v) { }
	#else
		INLINE ivec4(const ivec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	#endif
	
	INLINE ivec4 &operator=(const ivec4 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
		return *this;
	}
	INLINE ivec4 operator-() const {
		return ivec4(-x,-y,-z,-w);
	}
	INLINE ivec4 &operator*=(int v) {
		x *= v; y *= v; z *= v; w *= v;
		return *this;
	}
	INLINE ivec4 &operator*=(const ivec4 &v) {
		x *= v.x; y *= v.y; z *= v.z; w *= v.w;
		return *this;
	}
	INLINE ivec4 &operator/=(int v) {
		x /= v; y /= v; z /= v; w /= v;
		return *this;
	}
	INLINE ivec4 &operator/=(const ivec4 &v) {
		x /= v.x; y /= v.y; z /= v.z; w /= v.w;
		return *this;
	}
	INLINE ivec4 &operator+=(const ivec4 &v) {
		x += v.x; y += v.y; z += v.z; w += v.w;
		return *this;
	}
	INLINE ivec4 &operator-=(const ivec4 &v) {
		x -= v.x; y -= v.y; z -= v.z; w -= v.w;
		return *this;
	}
	INLINE ivec4 &operator<<=(int v) {
		x <<= v; y <<= v; z <<= v; w <<= v;
		return *this;
	}
	INLINE ivec4 &operator>>=(int v) {
		x >>= v; y >>= v; z >>= v; w >>= v;
		return *this;
	}
	
	INLINE operator int*() { return v; }
	INLINE operator const int*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE int &operator[](int i) {
		assert((unsigned int)i < 4 && "ivec4::operator[](): bad index");
		return v[i];
	}
	INLINE int operator[](int i) const {
		assert((unsigned int)i < 4 && "ivec4::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(int v) {
		x = v; y = v; z = v; w = v;
	}
	INLINE void set(int x_,int y_,int z_,int w_) {
		x = x_; y = y_; z = z_; w = w_;
	}
	INLINE void set(const ivec4 &v) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = v.vec;
		#else
			x = v.x; y = v.y; z = v.z; w = v.w;
		#endif
	}
	INLINE void set(const int *v) {
		x = v[0]; y = v[1]; z = v[2]; w = v[3];
	}
	INLINE void get(int *v) const {
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}
	INLINE int *get() { return v; }
	INLINE const int *get() const { return v; }
	
	INLINE int length2() const {
		return x * x + y * y + z * z + w * w;
	}
	
	union {
		struct {
			int x,y,z,w;
		};
		int v[4];
		#ifdef USE_SSE2
			__m128i vec;
		#elif USE_SSE
			__m128 vec;
		#elif USE_ALTIVEC
			vec_int4 vec;
		#elif USE_NEON
			int32x4_t vec;
		#endif
	};
};

/*
 */
extern const ivec4 ivec4_zero;
extern const ivec4 ivec4_one;

/*
 */
INLINE vec4::vec4(const ivec4 &v) : x(Math::itof(v.x)), y(Math::itof(v.y)), z(Math::itof(v.z)), w(Math::itof(v.w)) { }
INLINE dvec4::dvec4(const ivec4 &v) : x(Math::itod(v.x)), y(Math::itod(v.y)), z(Math::itod(v.z)), w(Math::itod(v.w)) { }

/*
 */
INLINE int operator==(const ivec4 &v0,const ivec4 &v1) {
	return (v0.x == v1.x && v0.y == v1.y && v0.z == v1.z && v0.w == v1.w);
}

INLINE int operator!=(const ivec4 &v0,const ivec4 &v1) {
	return (v0.x != v1.x || v0.y != v1.y || v0.z != v1.z || v0.w != v1.w);
}

INLINE ivec4 operator*(const ivec4 &v0,int v1) {
	return ivec4(v0.x * v1,v0.y * v1,v0.z * v1,v0.w * v1);
}

INLINE ivec4 operator*(const ivec4 &v0,const ivec4 &v1) {
	return ivec4(v0.x * v1.x,v0.y * v1.y,v0.z * v1.z,v0.w * v1.w);
}

INLINE ivec4 operator/(const ivec4 &v0,int v1) {
	return ivec4(v0.x / v1,v0.y / v1,v0.z / v1,v0.w / v1);
}

INLINE ivec4 operator/(const ivec4 &v0,const ivec4 &v1) {
	return ivec4(v0.x / v1.x,v0.y / v1.y,v0.z / v1.z,v0.w / v1.w);
}

INLINE ivec4 operator+(const ivec4 &v0,const ivec4 &v1) {
	return ivec4(v0.x + v1.x,v0.y + v1.y,v0.z + v1.z,v0.w + v1.w);
}

INLINE ivec4 operator-(const ivec4 &v0,const ivec4 &v1) {
	return ivec4(v0.x - v1.x,v0.y - v1.y,v0.z - v1.z,v0.w - v1.w);
}

INLINE ivec4 operator<<(const ivec4 &v0,int v1) {
	return ivec4(v0.x << v1,v0.y << v1,v0.z << v1,v0.w << v1);
}

INLINE ivec4 operator>>(const ivec4 &v0,int v1) {
	return ivec4(v0.x >> v1,v0.y >> v1,v0.z >> v1,v0.w >> v1);
}

/*
 */
INLINE int dot(const ivec4 &v0,const ivec4 &v1) {
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + v0.w * v1.w;
}

INLINE ivec4 &lerp(ivec4 &ret,const ivec4 &v0,const ivec4 &v1,int k) {
	ret.x = lerp(v0.x,v1.x,k);
	ret.y = lerp(v0.y,v1.y,k);
	ret.z = lerp(v0.z,v1.z,k);
	ret.w = lerp(v0.w,v1.w,k);
	return ret;
}

/*
 */
INLINE int length2(const ivec4 &v) {
	return v.length2();
}

/*
 */
ivec4 min(const ivec4 &v0,const ivec4 &v1);
ivec4 max(const ivec4 &v0,const ivec4 &v1);
ivec4 clamp(const ivec4 &v,const ivec4 &v0,const ivec4 &v1);
ivec4 lerp(const ivec4 &v0,const ivec4 &v1,int k);

/******************************************************************************\
*
* bvec4
*
\******************************************************************************/

/*
 */
struct bvec4 {
	
	INLINE bvec4() { }
	INLINE bvec4(unsigned char x,unsigned char y,unsigned char z,unsigned char w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE bvec4(unsigned char v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE bvec4(const unsigned char *v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) { }
	explicit INLINE bvec4(const vec4 &v) : x((unsigned char)Math::ftoi(v.x)), y((unsigned char)Math::ftoi(v.y)), z((unsigned char)Math::ftoi(v.z)), w((unsigned char)Math::ftoi(v.w)) { }
	explicit INLINE bvec4(const dvec4 &v) : x((unsigned char)Math::dtoi(v.x)), y((unsigned char)Math::dtoi(v.y)), z((unsigned char)Math::dtoi(v.z)), w((unsigned char)Math::dtoi(v.w)) { }
	explicit INLINE bvec4(const ivec4 &v) : x((unsigned char)v.x), y((unsigned char)v.y), z((unsigned char)v.z), w((unsigned char)v.w) { }
	INLINE bvec4(const bvec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	
	INLINE bvec4 &operator=(const bvec4 &v) {
		vec = v.vec;
		return *this;
	}
	
	INLINE operator unsigned char*() { return v; }
	INLINE operator const unsigned char*() const { return v; }
	INLINE operator void*() { return v; }
	INLINE operator const void*() const { return v; }
	
	INLINE unsigned char &operator[](int i) {
		assert((unsigned int)i < 4 && "bvec4::operator[](): bad index");
		return v[i];
	}
	INLINE unsigned char operator[](int i) const {
		assert((unsigned int)i < 4 && "bvec4::operator[](): bad index");
		return v[i];
	}
	
	INLINE void set(unsigned char v) {
		x = v; y = v; z = v; w = v;
	}
	INLINE void set(unsigned char x_,unsigned char y_,unsigned char z_,unsigned char w_) {
		x = x_; y = y_; z = z_; w = w_;
	}
	INLINE void set(const bvec4 &v) {
		x = v.x; y = v.y; z = v.z; w = v.w;
	}
	INLINE void set(const unsigned char *v) {
		x = v[0]; y = v[1]; z = v[2]; w = v[3];
	}
	INLINE void get(unsigned char *v) const {
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}
	INLINE unsigned char *get() { return v; }
	INLINE const unsigned char *get() const { return v; }
	
	union {
		struct {
			unsigned char x,y,z,w;
		};
		unsigned char v[4];
		unsigned int vec;
	};
};

/*
 */
extern const bvec4 bvec4_zero;
extern const bvec4 bvec4_one;

/*
 */
INLINE vec4::vec4(const bvec4 &v) : x(Math::itof(v.x)), y(Math::itof(v.y)), z(Math::itof(v.z)), w(Math::itof(v.w)) { }
INLINE dvec4::dvec4(const bvec4 &v) : x(Math::itod(v.x)), y(Math::itod(v.y)), z(Math::itod(v.z)), w(Math::itod(v.w)) { }
INLINE ivec4::ivec4(const bvec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }

/*
 */
INLINE int operator==(const bvec4 &v0,const bvec4 &v1) {
	return (v0.vec == v1.vec);
}

INLINE int operator!=(const bvec4 &v0,const bvec4 &v1) {
	return (v0.vec != v1.vec);
}

/*
 */
bvec4 min(const bvec4 &v0,const bvec4 &v1);
bvec4 max(const bvec4 &v0,const bvec4 &v1);
bvec4 clamp(const bvec4 &v,const bvec4 &v0,const bvec4 &v1);

/******************************************************************************\
*
* mat2
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) mat2 {
	
	INLINE mat2() { }
	mat2(const mat2 &m);
	explicit mat2(float v);
	explicit mat2(const mat3 &m);
	explicit mat2(const mat4 &m);
	explicit mat2(const dmat4 &m);
	explicit mat2(const float *m);
	
	INLINE mat2 &operator=(const mat2 &m) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			col = m.col;
		#else
			m00 = m.m00; m10 = m.m10;
			m01 = m.m01; m11 = m.m11;
		#endif
		return *this;
	}
	
	mat2 operator-() const;
	mat2 &operator*=(float v);
	mat2 &operator*=(const mat2 &m);
	mat2 &operator+=(const mat2 &m);
	mat2 &operator-=(const mat2 &m);
	
	INLINE operator float*() { return mat; }
	INLINE operator const float*() const { return mat; }
	INLINE operator void*() { return mat; }
	INLINE operator const void*() const { return mat; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 4 && "mat2::operator[](): bad index");
		return mat[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 4 && "mat2::operator[](): bad index");
		return mat[i];
	}
	
	INLINE void set(int row,int column,float v) {
		assert((unsigned int)row < 2 && "mat2::set(): bad row");
		assert((unsigned int)column < 2 && "mat2::set(): bad column");
		mat[column * 2 + row] = v;
	}
	INLINE float &get(int row,int column) {
		assert((unsigned int)row < 2 && "mat2::get(): bad row");
		assert((unsigned int)column < 2 && "mat2::get(): bad column");
		return mat[column * 2 + row];
	}
	INLINE float get(int row,int column) const {
		assert((unsigned int)row < 2 && "mat2::get(): bad row");
		assert((unsigned int)column < 2 && "mat2::get(): bad column");
		return mat[column * 2 + row];
	}
	
	void set(const mat2 &m);
	void set(const mat3 &m);
	void set(const mat4 &m);
	void set(const dmat4 &m);
	void set(const float *m);
	void get(float *m) const;
	INLINE float *get() { return mat; }
	INLINE const float *get() const { return mat; }
	
	void setRow(int row,const vec2 &v);
	vec2 getRow(int row) const;
	
	void setColumn(int column,const vec2 &v);
	vec2 getColumn(int column) const;
	
	void setZero();
	void setIdentity();
	void setRotate(float angle);
	void setScale(const vec2 &v);
	
	union {
		struct {
			float m00,m10;
			float m01,m11;
		};
		float mat[4];
		#ifdef USE_SSE
			__m128 col;
		#elif USE_ALTIVEC
			vec_float4 col;
		#elif USE_NEON
			float32x4_t col;
		#endif
	};
};

/*
 */
extern const mat2 mat2_zero;
extern const mat2 mat2_one;
extern const mat2 mat2_identity;

/*
 */
int operator==(const mat2 &m0,const mat2 &m1);
int operator!=(const mat2 &m0,const mat2 &m1);
mat2 operator*(const mat2 &m,float v);
vec2 operator*(const mat2 &m,const vec2 &v);
vec2 operator*(const vec2 &v,const mat2 &m);
dvec2 operator*(const mat2 &m,const dvec2 &v);
dvec2 operator*(const dvec2 &v,const mat2 &m);
mat2 operator*(const mat2 &m0,const mat2 &m1);
mat2 operator+(const mat2 &m0,const mat2 &m1);
mat2 operator-(const mat2 &m0,const mat2 &m1);

/*
 */
int compare(const mat2 &m0,const mat2 &m1);
int compare(const mat2 &m0,const mat2 &m1,float epsilon);
float trace(const mat2 &m);
float determinant(const mat2 &m);
mat2 &mul(mat2 &ret,const mat2 &m,float v);
vec2 &mul(vec2 &ret,const mat2 &m,const vec2 &v);
vec2 &mul(vec2 &ret,const vec2 &v,const mat2 &m);
dvec2 &mul(dvec2 &ret,const mat2 &m,const dvec2 &v);
dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat2 &m);
mat2 &mul(mat2 &ret,const mat2 &m0,const mat2 &m1);
mat2 &add(mat2 &ret,const mat2 &m0,const mat2 &m1);
mat2 &sub(mat2 &ret,const mat2 &m0,const mat2 &m1);
mat2 &transpose(mat2 &ret,const mat2 &m);
mat2 &inverse(mat2 &ret,const mat2 &m);
mat2 &inverse(mat2 &ret,const mat2 &m,float det);

/*
 */
mat2 transpose(const mat2 &m);
mat2 inverse(const mat2 &m);
mat2 inverse(const mat2 &m,float det);

/******************************************************************************\
*
* mat3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) mat3 {
	
	INLINE mat3() { }
	mat3(const mat3 &m);
	explicit mat3(float v);
	explicit mat3(const mat2 &m);
	explicit mat3(const mat4 &m);
	explicit mat3(const dmat4 &m);
	explicit mat3(const quat &q);
	explicit mat3(const float *m);
	
	INLINE mat3 &operator=(const mat3 &m) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			col0 = m.col0;
			col1 = m.col1;
			col2 = m.col2;
		#else
			m00 = m.m00; m10 = m.m10; m20 = m.m20; m30 = m.m30;
			m01 = m.m01; m11 = m.m11; m21 = m.m21; m31 = m.m31;
			m02 = m.m02; m12 = m.m12; m22 = m.m22; m32 = m.m32;
		#endif
		return *this;
	}
	
	mat3 operator-() const;
	mat3 &operator*=(float v);
	mat3 &operator*=(const mat3 &m);
	mat3 &operator+=(const mat3 &m);
	mat3 &operator-=(const mat3 &m);
	
	INLINE operator float*() { return mat; }
	INLINE operator const float*() const { return mat; }
	INLINE operator void*() { return mat; }
	INLINE operator const void*() const { return mat; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 12 && "mat3::operator[](): bad index");
		return mat[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 12 && "mat3::operator[](): bad index");
		return mat[i];
	}
	
	INLINE void set(int row,int column,float v) {
		assert((unsigned int)row < 3 && "mat3::set(): bad row");
		assert((unsigned int)column < 3 && "mat3::set(): bad column");
		mat[column * 4 + row] = v;
	}
	INLINE float &get(int row,int column) {
		assert((unsigned int)row < 3 && "mat3::get(): bad row");
		assert((unsigned int)column < 3 && "mat3::get(): bad column");
		return mat[column * 4 + row];
	}
	INLINE float get(int row,int column) const {
		assert((unsigned int)row < 3 && "mat3::get(): bad row");
		assert((unsigned int)column < 3 && "mat3::get(): bad column");
		return mat[column * 4 + row];
	}
	
	void set(const mat2 &m);
	void set(const mat3 &m);
	void set(const mat4 &m);
	void set(const dmat4 &m);
	void set(const quat &q);
	void set(const float *m);
	void get(float *m) const;
	INLINE float *get() { return mat; }
	INLINE const float *get() const { return mat; }
	
	void setRow(int row,const vec3 &v);
	vec3 getRow(int row) const;
	
	void setColumn(int column,const vec3 &v);
	vec3 getColumn(int column) const;
	
	void setDiagonal(const vec3 &v);
	vec3 getDiagonal() const;
	
	void setZero();
	void setIdentity();
	void setSkewSymmetric(const vec3 &v);
	void setRotate(const vec3 &axis,float angle);
	void setRotateX(float angle);
	void setRotateY(float angle);
	void setRotateZ(float angle);
	void setScale(const vec3 &v);
	
	quat getQuat() const;
	
	union {
		struct {
			float m00,m10,m20,m30;
			float m01,m11,m21,m31;
			float m02,m12,m22,m32;
		};
		float mat[12];
		#ifdef USE_SSE
			struct {
				__m128 col0;
				__m128 col1;
				__m128 col2;
			};
		#elif USE_ALTIVEC
			struct {
				vec_float4 col0;
				vec_float4 col1;
				vec_float4 col2;
			};
		#elif USE_NEON
			struct {
				float32x4_t col0;
				float32x4_t col1;
				float32x4_t col2;
			};
		#endif
	};
};

/*
 */
extern const mat3 mat3_zero;
extern const mat3 mat3_one;
extern const mat3 mat3_identity;

/*
 */
int operator==(const mat3 &m0,const mat3 &m1);
int operator!=(const mat3 &m0,const mat3 &m1);
mat3 operator*(const mat3 &m,float v);
vec2 operator*(const mat3 &m,const vec2 &v);
vec2 operator*(const vec2 &v,const mat3 &m);
vec3 operator*(const mat3 &m,const vec3 &v);
vec3 operator*(const vec3 &v,const mat3 &m);
dvec2 operator*(const mat3 &m,const dvec2 &v);
dvec2 operator*(const dvec2 &v,const mat3 &m);
dvec3 operator*(const mat3 &m,const dvec3 &v);
dvec3 operator*(const dvec3 &v,const mat3 &m);
mat3 operator*(const mat3 &m0,const mat3 &m1);
mat3 operator+(const mat3 &m0,const mat3 &m1);
mat3 operator-(const mat3 &m0,const mat3 &m1);

/*
 */
int compare(const mat3 &m0,const mat3 &m1);
int compare(const mat3 &m0,const mat3 &m1,float epsilon);
float trace(const mat3 &m);
float determinant(const mat3 &m);
mat3 &mul(mat3 &ret,const mat3 &m,float v);
vec2 &mul(vec2 &ret,const mat3 &m,const vec2 &v);
vec2 &mul(vec2 &ret,const vec2 &v,const mat3 &m);
vec3 &mul(vec3 &ret,const mat3 &m,const vec3 &v);
vec3 &mul(vec3 &ret,const vec3 &v,const mat3 &m);
dvec2 &mul(dvec2 &ret,const mat3 &m,const dvec2 &v);
dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat3 &m);
dvec3 &mul(dvec3 &ret,const mat3 &m,const dvec3 &v);
dvec3 &mul(dvec3 &ret,const dvec3 &v,const mat3 &m);
mat3 &mul(mat3 &ret,const mat3 &m,const vec3 &v);
mat3 &mul(mat3 &ret,const vec3 &v,const mat3 &m);
mat3 &mul(mat3 &ret,const mat3 &m0,const mat3 &m1);
mat3 &add(mat3 &ret,const mat3 &m0,const mat3 &m1);
mat3 &sub(mat3 &ret,const mat3 &m0,const mat3 &m1);
mat3 &orthonormalize(mat3 &ret,const mat3 &m);
mat3 &transpose(mat3 &ret,const mat3 &m);
mat3 &inverse(mat3 &ret,const mat3 &m);
mat3 &inverse(mat3 &ret,const mat3 &m,float det);

/*
 */
mat3 orthonormalize(const mat3 &m);
mat3 transpose(const mat3 &m);
mat3 inverse(const mat3 &m);
mat3 inverse(const mat3 &m,float det);

/*
 */
mat3 rotate3(const vec3 &axis,float angle);
mat3 rotate3(float x,float y,float z,float angle);
mat3 rotate3(const quat &q);
mat3 rotateX3(float angle);
mat3 rotateY3(float angle);
mat3 rotateZ3(float angle);
mat3 scale3(const vec3 &v);
mat3 scale3(float x,float y,float z);

/*
 */
mat3 jacobi(const mat3 &m,mat3 &v);

/******************************************************************************\
*
* mat4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) mat4 {
	
	INLINE mat4() { }
	mat4(const mat4 &m);
	explicit mat4(float v);
	explicit mat4(const mat2 &m);
	explicit mat4(const mat3 &m);
	explicit mat4(const dmat4 &m);
	explicit mat4(const quat &q);
	explicit mat4(const float *m);
	mat4(const mat3 &m,const vec3 &v);
	mat4(const quat &q,const vec3 &v);
	
	INLINE mat4 &operator=(const mat4 &m) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			col0 = m.col0;
			col1 = m.col1;
			col2 = m.col2;
			col3 = m.col3;
		#else
			m00 = m.m00; m10 = m.m10; m20 = m.m20; m30 = m.m30;
			m01 = m.m01; m11 = m.m11; m21 = m.m21; m31 = m.m31;
			m02 = m.m02; m12 = m.m12; m22 = m.m22; m32 = m.m32;
			m03 = m.m03; m13 = m.m13; m23 = m.m23; m33 = m.m33;
		#endif
		return *this;
	}
	
	mat4 operator-() const;
	mat4 &operator*=(float v);
	mat4 &operator*=(const mat4 &m);
	mat4 &operator+=(const mat4 &m);
	mat4 &operator-=(const mat4 &m);
	
	INLINE operator float*() { return mat; }
	INLINE operator const float*() const { return mat; }
	INLINE operator void*() { return mat; }
	INLINE operator const void*() const { return mat; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 16 && "mat4::operator[](): bad index");
		return mat[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 16 && "mat4::operator[](): bad index");
		return mat[i];
	}
	
	INLINE void set(int row,int column,float v) {
		assert((unsigned int)row < 4 && "mat4::set(): bad row");
		assert((unsigned int)column < 4 && "mat4::set(): bad column");
		mat[column * 4 + row] = v;
	}
	INLINE float &get(int row,int column) {
		assert((unsigned int)row < 4 && "mat4::get(): bad row");
		assert((unsigned int)column < 4 && "mat4::get(): bad column");
		return mat[column * 4 + row];
	}
	INLINE float get(int row,int column) const {
		assert((unsigned int)row < 4 && "mat4::get(): bad row");
		assert((unsigned int)column < 4 && "mat4::get(): bad column");
		return mat[column * 4 + row];
	}
	
	void set(const mat2 &m);
	void set(const mat3 &m);
	void set(const mat4 &m);
	void set(const dmat4 &m);
	void set(const quat &q);
	void set(const float *m);
	void set(const mat3 &m,const vec3 &v);
	void set(const quat &q,const vec3 &v);
	void get(float *m) const;
	INLINE float *get() { return mat; }
	INLINE const float *get() const { return mat; }
	
	void setRow(int row,const vec4 &v);
	void setRow3(int row,const vec3 &v);
	vec4 getRow(int row) const;
	vec3 getRow3(int row) const;
	
	void setColumn(int column,const vec4 &v);
	void setColumn3(int column,const vec3 &v);
	vec4 getColumn(int column) const;
	vec3 getColumn3(int column) const;
	
	void setDiagonal(const vec4 &v);
	vec4 getDiagonal() const;
	
	void setZero();
	void setIdentity();
	void setTranslate(const vec3 &v);
	void setRotate(const vec3 &axis,float angle);
	void setRotateX(float angle);
	void setRotateY(float angle);
	void setRotateZ(float angle);
	void setScale(const vec3 &v);
	
	union {
		struct {
			float m00,m10,m20,m30;
			float m01,m11,m21,m31;
			float m02,m12,m22,m32;
			float m03,m13,m23,m33;
		};
		float mat[16];
		#ifdef USE_SSE
			struct {
				__m128 col0;
				__m128 col1;
				__m128 col2;
				__m128 col3;
			};
		#elif USE_ALTIVEC
			struct {
				vec_float4 col0;
				vec_float4 col1;
				vec_float4 col2;
				vec_float4 col3;
			};
		#elif USE_NEON
			struct {
				float32x4_t col0;
				float32x4_t col1;
				float32x4_t col2;
				float32x4_t col3;
			};
		#endif
	};
};

/*
 */
extern const mat4 mat4_zero;
extern const mat4 mat4_one;
extern const mat4 mat4_identity;

/*
 */
int operator==(const mat4 &m0,const mat4 &m1);
int operator!=(const mat4 &m0,const mat4 &m1);
mat4 operator*(const mat4 &m,float v);
vec2 operator*(const mat4 &m,const vec2 &v);
vec2 operator*(const vec2 &v,const mat4 &m);
vec3 operator*(const mat4 &m,const vec3 &v);
vec3 operator*(const vec3 &v,const mat4 &m);
vec4 operator*(const mat4 &m,const vec4 &v);
vec4 operator*(const vec4 &v,const mat4 &m);
dvec2 operator*(const mat4 &m,const dvec2 &v);
dvec2 operator*(const dvec2 &v,const mat4 &m);
dvec3 operator*(const mat4 &m,const dvec3 &v);
dvec3 operator*(const dvec3 &v,const mat4 &m);
dvec4 operator*(const mat4 &m,const dvec4 &v);
dvec4 operator*(const dvec4 &v,const mat4 &m);
mat4 operator*(const mat4 &m0,const mat4 &m1);
mat4 operator+(const mat4 &m0,const mat4 &m1);
mat4 operator-(const mat4 &m0,const mat4 &m1);

/*
 */
int compare(const mat4 &m0,const mat4 &m1);
int compare(const mat4 &m0,const mat4 &m1,float epsilon);
float trace(const mat4 &m);
float determinant(const mat4 &m);
float determinant3(const mat4 &m);
mat4 &mul(mat4 &ret,const mat4 &m,float v);
vec2 &mul(vec2 &ret,const mat4 &m,const vec2 &v);
vec2 &mul(vec2 &ret,const vec2 &v,const mat4 &m);
vec3 &mul(vec3 &ret,const mat4 &m,const vec3 &v);
vec3 &mul(vec3 &ret,const vec3 &v,const mat4 &m);
vec4 &mul(vec4 &ret,const mat4 &m,const vec4 &v);
vec4 &mul(vec4 &ret,const vec4 &v,const mat4 &m);
dvec2 &mul(dvec2 &ret,const mat4 &m,const dvec2 &v);
dvec2 &mul(dvec2 &ret,const dvec2 &v,const mat4 &m);
dvec3 &mul(dvec3 &ret,const mat4 &m,const dvec3 &v);
dvec3 &mul(dvec3 &ret,const dvec3 &v,const mat4 &m);
dvec4 &mul(dvec4 &ret,const mat4 &m,const dvec4 &v);
dvec4 &mul(dvec4 &ret,const dvec4 &v,const mat4 &m);
vec2 &mul3(vec2 &ret,const mat4 &m,const vec2 &v);
vec2 &mul3(vec2 &ret,const vec2 &v,const mat4 &m);
vec3 &mul3(vec3 &ret,const mat4 &m,const vec3 &v);
vec3 &mul3(vec3 &ret,const vec3 &v,const mat4 &m);
vec4 &mul3(vec4 &ret,const mat4 &m,const vec4 &v);
vec4 &mul3(vec4 &ret,const vec4 &v,const mat4 &m);
dvec2 &mul3(dvec2 &ret,const mat4 &m,const dvec2 &v);
dvec2 &mul3(dvec2 &ret,const dvec2 &v,const mat4 &m);
dvec3 &mul3(dvec3 &ret,const mat4 &m,const dvec3 &v);
dvec3 &mul3(dvec3 &ret,const dvec3 &v,const mat4 &m);
dvec4 &mul3(dvec4 &ret,const mat4 &m,const dvec4 &v);
dvec4 &mul3(dvec4 &ret,const dvec4 &v,const mat4 &m);
mat4 &mul(mat4 &ret,const mat4 &m0,const mat4 &m1);
mat4 &mul4(mat4 &ret,const mat4 &m0,const mat4 &m1);
mat4 &mul3(mat4 &ret,const mat4 &m0,const mat4 &m1);
mat4 &mult(mat4 &ret,const mat4 &m,const vec3 &v);
vec3 &proj(vec3 &ret,const mat4 &m,const vec3 &v);
vec4 &proj(vec4 &ret,const mat4 &m,const vec4 &v);
dvec3 &proj(dvec3 &ret,const mat4 &m,const dvec3 &v);
dvec4 &proj(dvec4 &ret,const mat4 &m,const dvec4 &v);
mat4 &add(mat4 &ret,const mat4 &m0,const mat4 &m1);
mat4 &sub(mat4 &ret,const mat4 &m0,const mat4 &m1);
mat4 &orthonormalize(mat4 &ret,const mat4 &m);
mat4 &rotation(mat4 &ret,const mat4 &m);
mat4 &transpose(mat4 &ret,const mat4 &m);
mat4 &transpose3(mat4 &ret,const mat4 &m);
mat4 &inverse(mat4 &ret,const mat4 &m);
mat4 &inverse4(mat4 &ret,const mat4 &m);
mat4 &lerp(mat4 &ret,const mat4 &m0,const mat4 &q1,float k);

/*
 */
mat4 orthonormalize(const mat4 &m);
mat4 rotation(const mat4 &m);
mat4 transpose(const mat4 &m);
mat4 transpose3(const mat4 &m);
mat4 inverse(const mat4 &m);
mat4 inverse4(const mat4 &m);
mat4 lerp(const mat4 &m0,const mat4 &q1,float k);

/*
 */
mat4 translate(const vec3 &v);
mat4 translate(float x,float y,float z);
mat4 rotate(const vec3 &axis,float angle);
mat4 rotate(float x,float y,float z,float angle);
mat4 rotate(const quat &q);
mat4 rotateX(float angle);
mat4 rotateY(float angle);
mat4 rotateZ(float angle);
mat4 scale(const vec3 &v);
mat4 scale(float x,float y,float z);

/*
 */
mat4 reflect(const vec4 &plane);
mat4 ortho(float left,float right,float bottom,float top,float znear,float zfar);
mat4 frustum(float left,float right,float bottom,float top,float znear,float zfar);
mat4 perspective(float fov,float aspect,float znear,float zfar);
mat4 setTo(const vec3 &position,const vec3 &direction,const vec3 &up);
mat4 lookAt(const vec3 &position,const vec3 &direction,const vec3 &up);
mat4 obliqueProjection(const mat4 &projection,const vec4 &plane);
mat4 symmetryProjection(const mat4 &projection);
mat4 cubeTransform(int face);

/*
 */
void decomposeTransform(const mat4 &m,vec4 &position,quat &rot);
mat4 &composeTransform(mat4 &ret,const vec4 &position,const quat &rot);
void decomposeTransform(const mat4 &m,vec3 &position,quat &rot,vec3 &scale);
mat4 &composeTransform(mat4 &ret,const vec3 &position,const quat &rot,const vec3 &scale);
void decomposeProjection(const mat4 &projection,float &znear,float &zfar);

/*
 */
const mat4 &hardwareProjectionGL(const mat4 &projection,int width,int height);
const mat4 &hardwareProjectionD3D9(const mat4 &projection,int width,int height);
const mat4 &hardwareProjectionD3D10(const mat4 &projection,int width,int height);
extern const mat4 &(*hardwareProjection)(const mat4 &projection,int width,int height);

/******************************************************************************\
*
* dmat4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) dmat4 {
	
	INLINE dmat4() { }
	dmat4(const dmat4 &m);
	explicit dmat4(double v);
	explicit dmat4(const mat2 &m);
	explicit dmat4(const mat3 &m);
	explicit dmat4(const mat4 &m);
	explicit dmat4(const quat &q);
	dmat4(const double *m);
	dmat4(const mat3 &m,const dvec3 &v);
	dmat4(const quat &q,const dvec3 &v);
	
	INLINE dmat4 &operator=(const dmat4 &m) {
		#ifdef USE_SSE2
			col0 = m.col0; col1 = m.col1;
			col2 = m.col2; col3 = m.col3;
			col4 = m.col4; col5 = m.col5;
		#else
			m00 = m.m00; m10 = m.m10; m20 = m.m20;
			m01 = m.m01; m11 = m.m11; m21 = m.m21;
			m02 = m.m02; m12 = m.m12; m22 = m.m22;
			m03 = m.m03; m13 = m.m13; m23 = m.m23;
		#endif
		return *this;
	}
	
	dmat4 operator-() const;
	dmat4 &operator*=(double v);
	dmat4 &operator*=(const dmat4 &m);
	dmat4 &operator+=(const dmat4 &m);
	dmat4 &operator-=(const dmat4 &m);
	
	INLINE operator double*() { return mat; }
	INLINE operator const double*() const { return mat; }
	INLINE operator void*() { return mat; }
	INLINE operator const void*() const { return mat; }
	
	INLINE double &operator[](int i) {
		assert((unsigned int)i < 12 && "dmat4::operator[](): bad index");
		return mat[i];
	}
	INLINE double operator[](int i) const {
		assert((unsigned int)i < 12 && "dmat4::operator[](): bad index");
		return mat[i];
	}
	
	INLINE void set(int row,int column,double v) {
		assert((unsigned int)row < 3 && "dmat4::set(): bad row");
		assert((unsigned int)column < 4 && "dmat4::set(): bad column");
		mat[column * 3 + row] = v;
	}
	INLINE double &get(int row,int column) {
		assert((unsigned int)row < 3 && "dmat4::get(): bad row");
		assert((unsigned int)column < 4 && "dmat4::get(): bad column");
		return mat[column * 3 + row];
	}
	INLINE double get(int row,int column) const {
		assert((unsigned int)row < 3 && "dmat4::get(): bad row");
		assert((unsigned int)column < 4 && "dmat4::get(): bad column");
		return mat[column * 3 + row];
	}
	
	void set(const mat2 &m);
	void set(const mat3 &m);
	void set(const mat4 &m);
	void set(const dmat4 &m);
	void set(const quat &q);
	void set(const double *v);
	void set(const mat3 &m,const dvec3 &v);
	void set(const quat &q,const dvec3 &v);
	void get(double *m) const;
	INLINE double *get() { return mat; }
	INLINE const double *get() const { return mat; }
	
	void setRow(int row,const dvec4 &v);
	void setRow3(int row,const dvec3 &v);
	dvec4 getRow(int row) const;
	dvec3 getRow3(int row) const;
	
	void setColumn(int column,const dvec4 &v);
	void setColumn3(int column,const dvec3 &v);
	dvec4 getColumn(int column) const;
	dvec3 getColumn3(int column) const;
	
	void setZero();
	void setIdentity();
	void setTranslate(const dvec3 &v);
	void setRotate(const dvec3 &axis,double angle);
	void setRotateX(double angle);
	void setRotateY(double angle);
	void setRotateZ(double angle);
	void setScale(const dvec3 &v);
	
	union {
		struct {
			double m00,m10,m20;
			double m01,m11,m21;
			double m02,m12,m22;
			double m03,m13,m23;
		};
		double mat[12];
		#ifdef USE_SSE2
			struct {
				__m128d col0;
				__m128d col1;
				__m128d col2;
				__m128d col3;
				__m128d col4;
				__m128d col5;
			};
		#endif
	};
};

/*
 */
extern const dmat4 dmat4_zero;
extern const dmat4 dmat4_one;
extern const dmat4 dmat4_identity;

/*
 */
int operator==(const dmat4 &m0,const dmat4 &m1);
int operator!=(const dmat4 &m0,const dmat4 &m1);
dmat4 operator*(const dmat4 &m,double v);
vec2 operator*(const dmat4 &m,const vec2 &v);
vec2 operator*(const vec2 &v,const dmat4 &m);
vec3 operator*(const dmat4 &m,const vec3 &v);
vec3 operator*(const vec3 &v,const dmat4 &m);
vec4 operator*(const dmat4 &m,const vec4 &v);
vec4 operator*(const vec4 &v,const dmat4 &m);
dvec2 operator*(const dmat4 &m,const dvec2 &v);
dvec2 operator*(const dvec2 &v,const dmat4 &m);
dvec3 operator*(const dmat4 &m,const dvec3 &v);
dvec3 operator*(const dvec3 &v,const dmat4 &m);
dvec4 operator*(const dmat4 &m,const dvec4 &v);
dvec4 operator*(const dvec4 &v,const dmat4 &m);
dmat4 operator*(const dmat4 &m0,const dmat4 &m1);
dmat4 operator+(const dmat4 &m0,const dmat4 &m1);
dmat4 operator-(const dmat4 &m0,const dmat4 &m1);

/*
 */
int compare(const dmat4 &m0,const dmat4 &m1);
int compare(const dmat4 &m0,const dmat4 &m1,double epsilon);
double determinant(const dmat4 &m);
dmat4 &mul(dmat4 &ret,const dmat4 &m,double v);
vec2 &mul(vec2 &ret,const dmat4 &m,const vec2 &v);
vec2 &mul(vec2 &ret,const vec2 &v,const dmat4 &m);
vec3 &mul(vec3 &ret,const dmat4 &m,const vec3 &v);
vec3 &mul(vec3 &ret,const vec3 &v,const dmat4 &m);
vec4 &mul(vec4 &ret,const dmat4 &m,const vec4 &v);
vec4 &mul(vec4 &ret,const vec4 &v,const dmat4 &m);
vec2 &mul(vec2 &ret,const dmat4 &m,const dvec2 &v);
vec2 &mul(vec2 &ret,const dvec2 &v,const dmat4 &m);
vec3 &mul(vec3 &ret,const dmat4 &m,const dvec3 &v);
vec3 &mul(vec3 &ret,const dvec3 &v,const dmat4 &m);
vec4 &mul(vec4 &ret,const dmat4 &m,const dvec4 &v);
vec4 &mul(vec4 &ret,const dvec4 &v,const dmat4 &m);
dvec2 &mul(dvec2 &ret,const dmat4 &m,const dvec2 &v);
dvec2 &mul(dvec2 &ret,const dvec2 &v,const dmat4 &m);
dvec3 &mul(dvec3 &ret,const dmat4 &m,const dvec3 &v);
dvec3 &mul(dvec3 &ret,const dvec3 &v,const dmat4 &m);
dvec4 &mul(dvec4 &ret,const dmat4 &m,const dvec4 &v);
dvec4 &mul(dvec4 &ret,const dvec4 &v,const dmat4 &m);
vec2 &mul3(vec2 &ret,const dmat4 &m,const vec2 &v);
vec2 &mul3(vec2 &ret,const vec2 &v,const dmat4 &m);
vec3 &mul3(vec3 &ret,const dmat4 &m,const vec3 &v);
vec3 &mul3(vec3 &ret,const vec3 &v,const dmat4 &m);
vec4 &mul3(vec4 &ret,const dmat4 &m,const vec4 &v);
vec4 &mul3(vec4 &ret,const vec4 &v,const dmat4 &m);
vec2 &mul3(vec2 &ret,const dmat4 &m,const dvec2 &v);
vec2 &mul3(vec2 &ret,const dvec2 &v,const dmat4 &m);
vec3 &mul3(vec3 &ret,const dmat4 &m,const dvec3 &v);
vec3 &mul3(vec3 &ret,const dvec3 &v,const dmat4 &m);
vec4 &mul3(vec4 &ret,const dmat4 &m,const dvec4 &v);
vec4 &mul3(vec4 &ret,const dvec4 &v,const dmat4 &m);
dvec2 &mul3(dvec2 &ret,const dmat4 &m,const dvec2 &v);
dvec2 &mul3(dvec2 &ret,const dvec2 &v,const dmat4 &m);
dvec3 &mul3(dvec3 &ret,const dmat4 &m,const dvec3 &v);
dvec3 &mul3(dvec3 &ret,const dvec3 &v,const dmat4 &m);
dvec4 &mul3(dvec4 &ret,const dmat4 &m,const dvec4 &v);
dvec4 &mul3(dvec4 &ret,const dvec4 &v,const dmat4 &m);
dmat4 &mul(dmat4 &ret,const dmat4 &m0,const dmat4 &m1);
dmat4 &mul4(dmat4 &ret,const dmat4 &m0,const dmat4 &m1);
dmat4 &mul3(dmat4 &ret,const dmat4 &m0,const dmat4 &m1);
dmat4 &mult(dmat4 &ret,const dmat4 &m,const dvec3 &v);
dmat4 &add(dmat4 &ret,const dmat4 &m0,const dmat4 &m1);
dmat4 &sub(dmat4 &ret,const dmat4 &m0,const dmat4 &m1);
dmat4 &orthonormalize(dmat4 &ret,const dmat4 &m);
dmat4 &rotation(dmat4 &ret,const dmat4 &m);
dmat4 &inverse(dmat4 &ret,const dmat4 &m);
dmat4 &lerp(dmat4 &ret,const dmat4 &m0,const dmat4 &q1,double k);

/*
 */
dmat4 orthonormalize(const dmat4 &m);
dmat4 rotation(const dmat4 &m);
dmat4 inverse(const dmat4 &m);
dmat4 lerp(const dmat4 &m0,const dmat4 &m1,double k);

/*
 */
dmat4 translate(const dvec3 &v);
dmat4 translate(double x,double y,double z);

/*
 */
dmat4 reflect(const dvec4 &plane);
dmat4 setTo(const dvec3 &position,const dvec3 &direction,const vec3 &up);
dmat4 lookAt(const dvec3 &position,const dvec3 &direction,const vec3 &up);

/*
 */
void decomposeTransform(const dmat4 &m,dvec3 &position,quat &rot,vec3 &scale);
dmat4 &composeTransform(dmat4 &ret,const dvec3 &position,const quat &rot,const vec3 &scale);

/******************************************************************************\
*
* quat
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) quat {
	
	INLINE quat() { }
	quat(const vec3 &axis,float angle);
	quat(float x,float y,float z,float angle);
	quat(float angle_x,float angle_y,float angle_z);
	explicit INLINE quat(const float *q) : x(q[0]), y(q[1]), z(q[2]), w(q[3]) { }
	explicit INLINE quat(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(0.0f) { }
	explicit INLINE quat(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	explicit quat(const mat3 &m);
	explicit quat(const mat4 &m);
	explicit quat(const dmat4 &m);
	#ifdef USE_SSE
		INLINE quat(const quat &q) : vec(q.vec) { }
		explicit INLINE quat(__m128 vec) : vec(vec) { }
	#elif USE_ALTIVEC
		INLINE quat(const quat &q) : vec(q.vec) { }
		explicit INLINE quat(vec_float4 vec) : vec(vec) { }
	#elif USE_NEON
		INLINE quat(const quat &q) : vec(q.vec) { }
		explicit INLINE quat(float32x4_t vec) : vec(vec) { }
	#else
		INLINE quat(const quat &q) : x(q.x), y(q.y), z(q.z), w(q.w) { }
	#endif
	
	INLINE quat &operator=(const quat &q) {
		#if defined(USE_SSE) || defined(USE_ALTIVEC) || defined(USE_NEON)
			vec = q.vec;
		#else
			x = q.x; y = q.y; z = q.z; w = q.w;
		#endif
		return *this;
	}

	quat operator-() const;
	quat &operator*=(float v);
	quat &operator*=(const quat &q);
	quat &operator+=(const quat &q);
	quat &operator-=(const quat &q);
	
	INLINE operator float*() { return q; }
	INLINE operator const float*() const { return q; }
	INLINE operator void*() { return q; }
	INLINE operator const void*() const { return q; }
	
	INLINE float &operator[](int i) {
		assert((unsigned int)i < 4 && "quat::operator[](): bad index");
		return q[i];
	}
	INLINE float operator[](int i) const {
		assert((unsigned int)i < 4 && "quat::operator[](): bad index");
		return q[i];
	}
	
	void set(const vec3 &v);
	void set(const mat3 &m);
	void set(const mat4 &m);
	void set(const dmat4 &m);
	void set(const vec3 &axis,float angle);
	void set(float x,float y,float z,float angle);
	void set(float angle_x,float angle_y,float angle_z);
	void get(vec3 &axis,float &angle) const;
	
	INLINE void set(const float *q) {
		x = q[0]; y = q[1]; z = q[2]; w = q[3];
	}
	INLINE void get(float *q) const {
		q[0] = x; q[1] = y; q[2] = z; q[3] = w;
	}
	INLINE float *get() { return q; }
	INLINE const float *get() const { return q; }
	
	mat3 getMat3() const;
	float getAngle(const vec3 &axis) const;
	
	INLINE quat &normalize() {
		float ilength = Math::rsqrt(x * x + y * y + z * z + w * w);
		x *= ilength; y *= ilength; z *= ilength; w *= ilength;
		return *this;
	}
	INLINE quat &normalizeFast() {
		float ilength = Math::rsqrtFast(x * x + y * y + z * z + w * w);
		x *= ilength; y *= ilength; z *= ilength; w *= ilength;
		return *this;
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float q[4];
		#ifdef USE_SSE
			__m128 vec;
		#elif USE_ALTIVEC
			vec_float4 vec;
		#elif USE_NEON
			float32x4_t vec;
		#endif
	};
};

/*
 */
extern const quat quat_identity;

/*
 */
int operator==(const quat &q0,const quat &q1);
int operator!=(const quat &q0,const quat &q1);
quat operator*(const quat &q,float v);
vec3 operator*(const quat &q,const vec3 &v);
vec3 operator*(const vec3 &v,const quat &q);
dvec3 operator*(const quat &q,const dvec3 &v);
dvec3 operator*(const dvec3 &v,const quat &q);
quat operator*(const quat &q0,const quat &q1);
quat operator+(const quat &q0,const quat &q1);
quat operator-(const quat &q0,const quat &q1);

/*
 */
INLINE int compare(const quat &q0,const quat &q1) {
	return (compare(q0.x,q1.x) && compare(q0.y,q1.y) && compare(q0.z,q1.z) && compare(q0.w,q1.w));
}

INLINE int compare(const quat &q0,const quat &q1,float epsilon) {
	return (compare(q0.x,q1.x,epsilon) && compare(q0.y,q1.y,epsilon) && compare(q0.z,q1.z,epsilon) && compare(q0.w,q1.w,epsilon));
}

INLINE float dot(const quat &q0,const quat &q1) {
	return q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
}

INLINE quat &mul(quat &ret,const quat &q,float v) {
	ret.x = q.x * v;
	ret.y = q.y * v;
	ret.z = q.z * v;
	ret.w = q.w * v;
	return ret;
}

INLINE quat &mad(quat &ret,const quat &q0,float v,const quat &q1) {
	ret.x = q0.x * v + q1.x;
	ret.y = q0.y * v + q1.y;
	ret.z = q0.z * v + q1.z;
	ret.w = q0.w * v + q1.w;
	return ret;
}

INLINE quat &add(quat &ret,const quat &q0,const quat &q1) {
	ret.x = q0.x + q1.x;
	ret.y = q0.y + q1.y;
	ret.z = q0.z + q1.z;
	ret.w = q0.w + q1.w;
	return ret;
}

INLINE quat &sub(quat &ret,const quat &q0,const quat &q1) {
	ret.x = q0.x - q1.x;
	ret.y = q0.y - q1.y;
	ret.z = q0.z - q1.z;
	ret.w = q0.w - q1.w;
	return ret;
}

INLINE quat &inverse(quat &ret,const quat &q) {
	ret.x = -q.x;
	ret.y = -q.y;
	ret.z = -q.z;
	ret.w = q.w;
	return ret;
}

/*
 */
vec3 &mul(vec3 &ret,const quat &q,const vec3 &v);
vec3 &mul(vec3 &ret,const vec3 &v,const quat &q);
dvec3 &mul(dvec3 &ret,const quat &q,const dvec3 &v);
dvec3 &mul(dvec3 &ret,const dvec3 &v,const quat &q);
quat &mul(quat &ret,const quat &q0,const quat &q1);
quat &slerp(quat &ret,const quat &q0,const quat &q1,float k);

/*
 */
quat normalize(const quat &q);
quat inverse(const quat &q);
quat slerp(const quat &q0,const quat &q1,float k);

/*
 */
void decomposeTransform(const mat4 &m,quat &q0,quat &q1);
mat4 &composeTransform(mat4 &ret,const quat &q0,const quat &q1);

#endif /* __MATH_LIB_H__ */
