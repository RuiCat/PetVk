/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    MathLib.spu.h
 * Desc:    Math spu library
 * Version: 1.03
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

#ifndef __MATH_LIB_SPU_H__
#define __MATH_LIB_SPU_H__

#include "../utils/Base.spu.h"
#include "../utils/Half.h"

/*
 */
#ifdef INFINITY
	#undef INFINITY
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
#define SPU_PERM4_LX 	0x00010203
#define SPU_PERM4_LY 	0x04050607
#define SPU_PERM4_LZ 	0x08090a0b
#define SPU_PERM4_LW 	0x0c0d0e0f
#define SPU_PERM4_RX 	0x10111213
#define SPU_PERM4_RY 	0x14151617
#define SPU_PERM4_RZ 	0x18191a1b
#define SPU_PERM4_RW 	0x1c1d1e1f
#define SPU_PERM4_LB0	0x0004080c
#define SPU_PERM4_LB1	0x0105090d
#define SPU_PERM4_LB2	0x02060a0e
#define SPU_PERM4_LB3	0x03070b0f
#define SPU_PERM4_RB0	0x1014181c
#define SPU_PERM4_RB1	0x1115191d
#define SPU_PERM4_RB2	0x12161a1e
#define SPU_PERM4_RB3	0x13171b1f
#define SPU_PERM4_LPXY	0x02030607
#define SPU_PERM4_LPZW	0x0a0b0e0f
#define SPU_PERM4_RPXY	0x12131617
#define SPU_PERM4_RPZW	0x1a1b1e1f
#define SPU_PERM4_LULX	0x80800001
#define SPU_PERM4_LULY	0x80800203
#define SPU_PERM4_LULZ	0x80800405
#define SPU_PERM4_LULW	0x80800607
#define SPU_PERM4_RULX	0x80801011
#define SPU_PERM4_RULY	0x80801213
#define SPU_PERM4_RULZ	0x80801415
#define SPU_PERM4_RULW	0x80801617
#define SPU_PERM4_LURX	0x80800809
#define SPU_PERM4_LURY	0x80800a0b
#define SPU_PERM4_LURZ	0x80800c0d
#define SPU_PERM4_LURW	0x80800e0f
#define SPU_PERM4_RURX	0x80801819
#define SPU_PERM4_RURY	0x80801a1b
#define SPU_PERM4_RURZ	0x80801c1d
#define SPU_PERM4_RURW	0x80801e1f
#define SPU_PERM4(X,Y,Z,W) (vec_uchar16) (vec_uint4) { SPU_PERM4_ ## X, SPU_PERM4_ ## Y, SPU_PERM4_ ## Z, SPU_PERM4_ ## W }
#define SPU_PERM2(X,Y,Z,W) (vec_uchar16) (vec_uint4) { SPU_PERM4_L ## X, SPU_PERM4_L ## Y, SPU_PERM4_R ## Z, SPU_PERM4_R ## W }
#define SPU_SWIZZLE(V,X,Y,Z,W) spu_shuffle(V,V,SPU_PERM2(X,Y,Z,W))
#define SPU_FLOAT4(X,Y,Z,W) (vec_float4) { X, Y, Z, W }
#define SPU_UINT4(X,Y,Z,W) (vec_uint4) { X, Y, Z, W }

/*
 */
struct vec3;
struct vec4;
struct mat3;
struct mat4;
struct quat;

/*
 */
INLINE vec_float4 spu_rcp_nr(vec_float4 v);
INLINE vec_float4 spu_rsqrt_nr(vec_float4 v);

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

/******************************************************************************\
*
* MathLib
*
\******************************************************************************/

/*
 */
class Math {
		
		Math() { }
		
	public:
		
		// functions
		static int signMask(int v);
		static float sign(float v);
		static vec_float4 sign4(vec_float4 v);
		
		static int abs(int v);
		static float abs(float v);
		static vec_float4 abs4(vec_float4 v);
		
		static float ceil(float v);
		static vec_float4 ceil4(vec_float4 v);
		
		static float floor(float v);
		static vec_float4 floor4(vec_float4 v);
		
		static float frac(float v);
		static vec_float4 frac4(vec_float4 v);
		
		static float sqrt(float v);
		static vec_float4 sqrt4(vec_float4 v);
		
		static float rcp(float v);
		static vec_float4 rcp4(vec_float4 v);
		
		static float rsqrt(float v);
		static float rsqrtFast(float v);
		static vec_float4 rsqrt4(vec_float4 v);
		static vec_float4 rsqrtFast4(vec_float4 v);
		
		static float mod(float x,float y);
		static vec_float4 mod(vec_float4 x,vec_float4 y);
		
		static float pow(float x,float y);
		static float powFast(float x,float y);
		static vec_float4 pow4(vec_float4 x,vec_float4 y);
		
		static float exp(float v);
		static float expFast(float v);
		static vec_float4 exp4(vec_float4 v);
		
		static float exp2(float v);
		static float exp2Fast(float v);
		static vec_float4 exp24(vec_float4 v);
		
		static float log(float v);
		static float logFast(float v);
		static vec_float4 log4(vec_float4 v);
		
		static float log2(float v);
		static float log2Fast(float v);
		static vec_float4 log24(vec_float4 v);
		
		static float log10(float v);
		static vec_float4 log104(vec_float4 v);
		
		// trigonometry
		static float sin(float a);
		static vec_float4 sin4(vec_float4 a);
		
		static float cos(float a);
		static vec_float4 cos4(vec_float4 a);
		
		static float tan(float a);
		static vec_float4 tan4(vec_float4 a);
		
		static float asin(float v);
		static vec_float4 asin4(vec_float4 v);
		
		static float acos(float v);
		static vec_float4 acos4(vec_float4 v);
		
		static float atan(float v);
		static vec_float4 atan4(vec_float4 v);
		
		static float atan2(float y,float x);
		static vec_float4 atan24(vec_float4 y,vec_float4 x);
		
		static void sincos(float a,float &s,float &c);
		static void sincosFast(float a,float &s,float &c);
		static void sincos4(vec_float4 a,vec_float4 &s,vec_float4 &c);
		static void sincosFast4(vec_float4 a,vec_float4 &s,vec_float4 &c);
		
		// branching
		static int select(int c,int v0,int v1);
		static float select(int c,float v0,float v1);
		static float select(float c,float v0,float v1);
		
		// conversion
		static float itof(int v);
		static int ftoi(float v);
		static int round(float v);
};

/*
 */
INLINE int Math::signMask(int v) {
	return (v >> 31);
}

INLINE float Math::sign(float v) {
	return IntFloat((IntFloat(v).ui & 0x80000000) | 0x3f800000).f;
}

INLINE vec_float4 Math::sign4(vec_float4 v) {
	return (vec_float4)spu_or(spu_and((vec_uint4)v,0x80000000),0x3f800000);
}

/*
 */
INLINE int Math::abs(int v) {
	return ::abs(v);
}

INLINE float Math::abs(float v) {
	return IntFloat(IntFloat(v).ui & 0x7fffffff).f;
}

INLINE vec_float4 Math::abs4(vec_float4 v) {
	return ::fabsf4(v);
}

/*
 */
INLINE float Math::ceil(float v) {
	return ::ceilf(v);
}

INLINE vec_float4 Math::ceil4(vec_float4 v) {
	return ::ceilf4(v);
}

/*
 */
INLINE float Math::floor(float v) {
	return ::floorf(v);
}

INLINE vec_float4 Math::floor4(vec_float4 v) {
	return ::floorf4(v);
}

/*
 */
INLINE float Math::frac(float v) {
	return v - ::floorf(v);
}

INLINE vec_float4 Math::frac4(vec_float4 v) {
	return spu_sub(v,::floorf4(v));
}

/*
 */
INLINE float Math::sqrt(float v) {
	return ::sqrtf(v);
}

INLINE vec_float4 Math::sqrt4(vec_float4 v) {
	return ::sqrtf4(v);
}

/*
 */
INLINE float Math::rcp(float v) {
	return 1.0f / v;
}

INLINE vec_float4 Math::rcp4(vec_float4 v) {
	return spu_rcp_nr(v);
}

/*
 */
INLINE float Math::rsqrt(float v) {
	if(v < 1e-8f) return INFINITY;
	return 1.0f / ::sqrtf(v);
}

INLINE float Math::rsqrtFast(float v) {
	IntFloat i = v;
	i.i = 0x5f3759df - (i.i >> 1);
	return i.f * (1.5f - (i.f * i.f * v * 0.5f));
}

INLINE vec_float4 Math::rsqrt4(vec_float4 v) {
	return spu_rsqrt_nr(v);
}

INLINE vec_float4 Math::rsqrtFast4(vec_float4 v) {
	return spu_rsqrte(v);
}

/*
 */
INLINE float Math::pow(float x,float y) {
	return ::powf(x,y);
}

INLINE float Math::powFast(float x,float y) {
	return Math::exp2Fast(Math::log2Fast(x) * y);
}

INLINE vec_float4 Math::pow4(vec_float4 x,vec_float4 y) {
	return ::powf4(x,y);
}

/*
 */
INLINE float Math::mod(float x,float y) {
	return ::fmodf(x,y);
}

INLINE vec_float4 Math::mod(vec_float4 x,vec_float4 y) {
	return ::fmodf4(x,y);
}

/*
 */
INLINE float Math::exp(float v) {
	return ::expf(v);
}

INLINE float Math::expFast(float v) {
	return exp2Fast(v * (1.0f / LOG2));
}

INLINE vec_float4 Math::exp4(vec_float4 v) {
	return ::expf4(v);
}

/*
 */
INLINE float Math::exp2(float v) {
	return ::exp2f(v);
}

INLINE float Math::exp2Fast(float v) {
	int i = ftoi(v - 0.5f);
	v = v - itof(i);
	return IntFloat((i + 127) << 23).f * (((((0.0018775767f * v + 0.0089893397f) * v + 0.055826318f) * v + 0.24015361f) * v + 0.69315308f) * v + 0.99999994f);
}

INLINE vec_float4 Math::exp24(vec_float4 v) {
	return ::exp2f4(v);
}

/*
 */
INLINE float Math::log(float v) {
	return ::logf(v);
}

INLINE float Math::logFast(float v) {
	return log2Fast(v) * LOG2;
}

INLINE vec_float4 Math::log4(vec_float4 v) {
	return ::logf4(v);
}

/*
 */
INLINE float Math::log2(float v) {
	return ::log2f(v);
}

INLINE float Math::log2Fast(float v) {
	int i = IntFloat(v).i;
	int e = ((i >> 23) & 0xff) - 127;
	v = IntFloat((i & 0x007fffff) | 0x3f800000).f;
	return itof(e) + (v - 1.0f) * (((((-0.034436006f * v + 0.31821337f) * v - 1.2315303f) * v + 2.5988452f) * v - 3.3241990f) * v + 3.1157899f);
}

INLINE vec_float4 Math::log24(vec_float4 v) {
	return ::log2f4(v);
}

/*
 */
INLINE float Math::log10(float v) {
	return ::log10f(v);
}

INLINE vec_float4 Math::log104(vec_float4 v) {
	return ::log10f4(v);
}

/*
 */
INLINE float Math::sin(float a) {
	return ::sinf(a);
}

INLINE vec_float4 Math::sin4(vec_float4 a) {
	return ::sinf4(a);
}

/*
 */
INLINE float Math::cos(float a) {
	return ::cosf(a);
}

INLINE vec_float4 Math::cos4(vec_float4 a) {
	return ::cosf4(a);
}

/*
 */
INLINE float Math::tan(float a) {
	return ::tanf(a);
}

INLINE vec_float4 Math::tan4(vec_float4 a) {
	return ::tanf4(a);
}

/*
 */
INLINE float Math::asin(float v) {
	return ::asinf(v);
}

INLINE vec_float4 Math::asin4(vec_float4 v) {
	return ::asinf4(v);
}

/*
 */
INLINE float Math::acos(float v) {
	return ::acosf(v);
}

INLINE vec_float4 Math::acos4(vec_float4 v) {
	return ::acosf4(v);
}

/*
 */
INLINE float Math::atan(float v) {
	return ::atanf(v);
}

INLINE vec_float4 Math::atan4(vec_float4 v) {
	return ::atanf4(v);
}

/*
 */
INLINE float Math::atan2(float y,float x) {
	return ::atan2f(y,x);
}

INLINE vec_float4 Math::atan24(vec_float4 y,vec_float4 x) {
	return ::atan2f4(y,x);
}

/*
 */
INLINE void Math::sincos(float a,float &s,float &c) {
	s = Math::sin(a);
	c = Math::cos(a);
}

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

INLINE void Math::sincos4(vec_float4 a,vec_float4 &s,vec_float4 &c) {
	::sincosf4(a,&s,&c);
}

INLINE void Math::sincosFast4(vec_float4 a,vec_float4 &s,vec_float4 &c) {
	vec_float4 one = spu_splats(1.0f);
	vec_float4 ia = floorf4(spu_mul(a,spu_splats(1.0f / PI2)));
	a = spu_sub(a,spu_mul(ia,spu_splats(PI2)));
	s = spu_sub(spu_splats(PI),a);
	vec_uint4 sign = spu_and((vec_uint4)s,0x80000000);
	vec_uint4 mask = spu_cmpabsgt(s,spu_splats(PI05));
	vec_float4 is = spu_sub((vec_float4)spu_or((vec_uint4)spu_splats(PI),sign),s);
	s = spu_sel(s,is,mask);
	c = spu_sel(spu_splats(-1.0f),one,mask);
	vec_float4 a2 = spu_mul(s,s);
	s = spu_mul(s,spu_madd(spu_madd(spu_splats(7.610e-03f),a2,spu_splats(-1.6605e-01f)),a2,one));
	c = spu_mul(c,spu_madd(spu_madd(spu_splats(3.705e-02f),a2,spu_splats(-4.9670e-01f)),a2,one));
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
	int mask = Math::signMask(IntFloat(c).i);
	return IntFloat((IntFloat(v0).i & mask) | (IntFloat(v1).i & ~mask)).f;
}

/*
 */
INLINE float Math::itof(int v) {
	return static_cast<float>(v);
}

INLINE int Math::ftoi(float v) {
	return static_cast<int>(v);
}

INLINE int Math::round(float v) {
	return static_cast<int>(v + 0.5f);
}

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
	return (v0 < v1) ? v0 : v1;
}

INLINE float max(float v0,float v1) {
	return (v0 > v1) ? v0 : v1;
}

INLINE float clamp(float v,float v0,float v1) {
	if(v < v0) return v0;
	if(v > v1) return v1;
	return v;
}

INLINE float saturate(float v) {
	if(v < 0.0f) return 0.0f;
	if(v > 1.0f) return 1.0f;
	return v;
}

INLINE float lerp(float v0,float v1,float k) {
	return v0 + (v1 - v0) * k;
}

/*
 */
INLINE int min(int v0,int v1) {
	int i = v1 - v0;
	i &= ~Math::signMask(i);
	return v1 - i;
}

INLINE int max(int v0,int v1) {
	int i = v1 - v0;
	i &= ~Math::signMask(i);
	return v0 + i;
}

INLINE int clamp(int v,int v0,int v1) {
	int i = v0 - v;
	i &= ~Math::signMask(i);
	i = v1 - v - i;
	i &= ~Math::signMask(i);
	return v1 - i;
}

/******************************************************************************\
*
* Vectors
*
\******************************************************************************/

/*
 */
INLINE vec_float4 spu_rcp_nr(vec_float4 v) {
	vec_float4 iv = spu_re(v);
	vec_float4 one = spu_splats(1.0f);
	iv = spu_madd(spu_nmsub(iv,v,one),iv,iv);
	return spu_madd(spu_nmsub(iv,v,one),iv,iv);
}

/*
 */
INLINE vec_float4 spu_rsqrt_nr(vec_float4 v) {
	vec_float4 iv = spu_rsqrte(v);
	vec_float4 nr = spu_mul(spu_mul(v,iv),iv);
	return spu_mul(spu_mul(spu_splats(0.5f),iv),spu_sub(spu_splats(3.0f),nr));
}

/*
 */
INLINE vec_float4 spu_dot33(vec_float4 v0,vec_float4 v1) {
	vec_float4 v2 = spu_mul(v0,v1);
	vec_float4 v3 = spu_add(v2,SPU_SWIZZLE(v2,Y,X,Y,W));
	return spu_add(v3,SPU_SWIZZLE(v2,Z,Z,X,W));
}

INLINE vec_float4 spu_dot44(vec_float4 v0,vec_float4 v1) {
	vec_float4 v2 = spu_mul(v0,v1);
	v2 = spu_add(v2,spu_rlqwbyte(v2,8));
	return spu_add(v2,spu_rlqwbyte(v2,4));
}

/*
 */
INLINE vec_float4 spu_normalize3(vec_float4 v) {
	vec_float4 ilength = spu_rsqrt_nr(spu_dot33(v,v));
	return spu_mul(v,ilength);
}

INLINE vec_float4 spu_normalize4(vec_float4 v) {
	vec_float4 ilength = spu_rsqrt_nr(spu_dot44(v,v));
	return spu_mul(v,ilength);
}

/*
 */
INLINE vec_float4 spu_cross(vec_float4 v0,vec_float4 v1) {
	vec_uchar16 yzxw = SPU_PERM2(Y,Z,X,W);
	vec_float4 v0_yzxw = spu_shuffle(v0,v0,yzxw);
	vec_float4 v1_yzxw = spu_shuffle(v1,v1,yzxw);
	vec_float4 v2 = spu_nmsub(v1,v0_yzxw,spu_mul(v0,v1_yzxw));
	return spu_shuffle(v2,v2,yzxw);
}

/*
 */
INLINE vec_float4 spu_half_to_float(vec_uint4 v) {
	vec_uint4 em = spu_and(v,0x7fff);
	vec_uint4 fn = spu_add(spu_rl(em,13),(127 - 15) << 23);
	vec_uint4 fd = (vec_uint4)spu_convtf(em,24);
	vec_uint4 f = spu_sel(fd,fn,spu_cmpgt(em,spu_splats((unsigned int)0x03ff)));
	return (vec_float4)spu_sel(f,spu_rl(v,16),spu_splats((unsigned int)0x80000000));
}

INLINE vec_uint4 spu_float_to_half(vec_float4 v) {
	vec_uint4 i = (vec_uint4)v;
	vec_uint4 e = spu_and(spu_rl(i,32 - 23),0x00ff);
	vec_uint4 m = spu_and(i,0x007fffff);
	vec_uint4 hn = spu_and(spu_rl(i,32 - 13),0x3fff);
	vec_uint4 hd = spu_rlmask(spu_or(m,0x00800000),spu_sub((vec_int4)e,spu_splats(127 - 1)));
	vec_uint4 h = spu_sel(hn,hd,spu_cmpgt(spu_splats((unsigned int)(127 - 14)),e));
	return spu_sel(h,spu_rl(i,16),spu_splats((unsigned int)0xc000));
}

/******************************************************************************\
*
* vec3
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) vec3 {
	
	INLINE vec3() { }
	INLINE vec3(const vec3 &v) : vec(v.vec) { }
	INLINE vec3(float x,float y,float z) : x(x), y(y), z(z), w(0.0f) { }
	explicit INLINE vec3(float v) : x(v), y(v), z(v), w(0.0f) { }
	explicit INLINE vec3(const vec4 &v);
	explicit INLINE vec3(vec_float4 v) : vec(v) { }
	
	INLINE vec3 &operator=(const vec3 &v) {
		vec = v.vec;
		return *this;
	}
	INLINE vec3 operator-() const {
		return vec3((vec_float4)spu_xor((vec_uint4)vec,0x80000000));
	}
	INLINE vec3 &operator*=(float v) {
		vec = spu_mul(vec,spu_splats(v));
		return *this;
	}
	INLINE vec3 &operator*=(const vec3 &v) {
		vec = spu_mul(vec,v.vec);
		return *this;
	}
	INLINE vec3 &operator/=(float v) {
		vec = spu_mul(vec,spu_splats(Math::rcp(v)));
		return *this;
	}
	INLINE vec3 &operator/=(const vec3 &v) {
		vec = spu_mul(vec,Math::rcp4(v.vec));
		return *this;
	}
	INLINE vec3 &operator+=(const vec3 &v) {
		vec = spu_add(vec,v.vec);
		return *this;
	}
	INLINE vec3 &operator-=(const vec3 &v) {
		vec = spu_sub(vec,v.vec);
		return *this;
	}
	
	INLINE float &operator[](short i) {
		assert((unsigned short)i < 3 && "vec3::operator[](): bad index");
		return v[i];
	}
	INLINE float operator[](short i) const {
		assert((unsigned short)i < 3 && "vec3::operator[](): bad index");
		return v[i];
	}
	
	INLINE float length2() const {
		return spu_extract(spu_dot33(vec,vec),0);
	}
	INLINE float length() const {
		return Math::sqrt(spu_extract(spu_dot33(vec,vec),0));
	}
	INLINE vec3 &normalize() {
		vec = spu_normalize3(vec);
		return *this;
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float v[4];
		vec_float4 vec;
	};
};

/*
 */
INLINE int operator==(const vec3 &v0,const vec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int operator!=(const vec3 &v0,const vec3 &v1) {
	return !(compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE vec3 operator*(const vec3 &v0,float v1) {
	return vec3(spu_mul(v0.vec,spu_splats(v1)));
}

INLINE vec3 operator*(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_mul(v0.vec,v1.vec));
}

INLINE vec3 operator/(const vec3 &v0,float v1) {
	return vec3(spu_mul(v0.vec,spu_splats(Math::rcp(v1))));
}

INLINE vec3 operator/(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_mul(v0.vec,Math::rcp4(v1.vec)));
}

INLINE vec3 operator+(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_add(v0.vec,v1.vec));
}

INLINE vec3 operator-(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_sub(v0.vec,v1.vec));
}

/*
 */
INLINE int compare(const vec3 &v0,const vec3 &v1) {
	return (compare(v0.x,v1.x) && compare(v0.y,v1.y) && compare(v0.z,v1.z));
}

INLINE int compare(const vec3 &v0,const vec3 &v1,float epsilon) {
	return (compare(v0.x,v1.x,epsilon) && compare(v0.y,v1.y,epsilon) && compare(v0.z,v1.z,epsilon));
}

INLINE float dot(const vec3 &v0,const vec3 &v1) {
	return spu_extract(spu_dot33(v0.vec,v1.vec),0);
}

INLINE vec3 &mul(vec3 &ret,const vec3 &v0,float v1) {
	ret.vec = spu_mul(v0.vec,spu_splats(v1));
	return ret;
}

INLINE vec3 &mul(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.vec = spu_mul(v0.vec,v1.vec);
	return ret;
}

INLINE vec3 &mad(vec3 &ret,const vec3 &v0,float v1,const vec3 &v2) {
	ret.vec = spu_madd(v0.vec,spu_splats(v1),v2.vec);
	return ret;
}

INLINE vec3 &mad(vec3 &ret,const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	ret.vec = spu_madd(v0.vec,v1.vec,v2.vec);
	return ret;
}

INLINE vec3 &add(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.vec = spu_add(v0.vec,v1.vec);
	return ret;
}

INLINE vec3 &sub(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.vec = spu_sub(v0.vec,v1.vec);
	return ret;
}

INLINE vec3 &lerp(vec3 &ret,const vec3 &v0,const vec3 &v1,float k) {
	ret.vec = spu_madd(spu_sub(v1.vec,v0.vec),spu_splats(k),v0.vec);
	return ret;
}

INLINE vec3 &cross(vec3 &ret,const vec3 &v0,const vec3 &v1) {
	ret.vec = spu_cross(v0.vec,v1.vec);
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

INLINE vec3 lerp(const vec3 &v0,const vec3 &v1,float k) {
	vec3 ret;
	return lerp(ret,v0,v1,k);
}

INLINE vec3 cross(const vec3 &v0,const vec3 &v1) {
	vec3 ret;
	return cross(ret,v0,v1);
}

/*
 */
INLINE vec3 min(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_sel(v0.vec,v1.vec,spu_cmpgt(v0.vec,v1.vec)));
}

INLINE vec3 max(const vec3 &v0,const vec3 &v1) {
	return vec3(spu_sel(v0.vec,v1.vec,spu_cmpgt(v1.vec,v0.vec)));
}

INLINE vec3 clamp(const vec3 &v,const vec3 &v0,const vec3 &v1) {
	vec_uint4 mask0 = spu_cmpgt(v.vec,v0.vec);
	vec_uint4 mask1 = spu_cmpgt(v.vec,v1.vec);
	return vec3(spu_sel(v0.vec,spu_sel(v.vec,v1.vec,mask1),mask0));
}

INLINE vec3 saturate(const vec3 &v) {
	vec_float4 one = spu_splats(1.0f);
	vec_float4 zero = spu_splats(0.0f);
	vec_uint4 mask_0 = spu_cmpgt(v.vec,one);
	vec_uint4 mask_1 = spu_cmpgt(v.vec,zero);
	return vec3(spu_sel(zero,spu_sel(v.vec,one,mask_0),mask_1));
}

/******************************************************************************\
*
* vec4
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) vec4 {
	
	INLINE vec4() { }
	INLINE vec4(const vec4 &v) : vec(v.vec) { }
	INLINE vec4(const vec3 &v,float w) : x(v.x), y(v.y), z(v.z), w(w) { }
	INLINE vec4(float x,float y,float z,float w) : x(x), y(y), z(z), w(w) { }
	explicit INLINE vec4(float v) : x(v), y(v), z(v), w(v) { }
	explicit INLINE vec4(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(1.0f) { }
	explicit INLINE vec4(vec_float4 v) : vec(v) { }
	
	INLINE vec4 &operator=(const vec4 &v) {
		vec = v.vec;
		return *this;
	}
	INLINE vec4 operator-() const {
		return vec4((vec_float4)spu_xor((vec_uint4)vec,0x80000000));
	}
	INLINE vec4 &operator*=(float v) {
		vec = spu_mul(vec,spu_splats(v));
		return *this;
	}
	INLINE vec4 &operator*=(const vec4 &v) {
		vec = spu_mul(vec,v.vec);
		return *this;
	}
	INLINE vec4 &operator/=(float v) {
		vec = spu_mul(vec,spu_splats(Math::rcp(v)));
		return *this;
	}
	INLINE vec4 &operator/=(const vec4 &v) {
		vec = spu_mul(vec,Math::rcp4(v.vec));
		return *this;
	}
	INLINE vec4 &operator+=(const vec4 &v) {
		vec = spu_add(vec,v.vec);
		return *this;
	}
	INLINE vec4 &operator-=(const vec4 &v) {
		vec = spu_sub(vec,v.vec);
		return *this;
	}
	
	INLINE float &operator[](short i) {
		assert((unsigned short)i < 4 && "vec4::operator[](): bad index");
		return v[i];
	}
	INLINE float operator[](short i) const {
		assert((unsigned short)i < 4 && "vec4::operator[](): bad index");
		return v[i];
	}
	
	INLINE float length2() const {
		return spu_extract(spu_dot44(vec,vec),0);
	}
	INLINE float length() const {
		return Math::sqrt(spu_extract(spu_dot44(vec,vec),0));
	}
	INLINE vec4 &normalize() {
		vec = spu_normalize4(vec);
		return *this;
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float v[4];
		vec_float4 vec;
	};
};

/*
 */
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
	return vec4(spu_mul(v0.vec,spu_splats(v1)));
}

INLINE vec4 operator*(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_mul(v0.vec,v1.vec));
}

INLINE vec4 operator/(const vec4 &v0,float v1) {
	return vec4(spu_mul(v0.vec,spu_splats(Math::rcp(v1))));
}

INLINE vec4 operator/(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_mul(v0.vec,Math::rcp4(v1.vec)));
}

INLINE vec4 operator+(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_add(v0.vec,v1.vec));
}

INLINE vec4 operator-(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_sub(v0.vec,v1.vec));
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
	return spu_extract(spu_dot33(v0.vec,v1.vec),0) + v1.w;
}

INLINE float dot(const vec4 &v0,const vec3 &v1) {
	return spu_extract(spu_dot33(v0.vec,v1.vec),0) + v0.w;
}

INLINE float dot(const vec4 &v0,const vec4 &v1) {
	return spu_extract(spu_dot44(v0.vec,v1.vec),0);
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
	ret.vec = spu_mul(v0.vec,spu_splats(v1));
	return ret;
}

INLINE vec4 &mul(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.vec = spu_mul(v0.vec,v1.vec);
	return ret;
}

INLINE vec4 &mad(vec4 &ret,const vec4 &v0,float v1,const vec4 &v2) {
	ret.vec = spu_madd(v0.vec,spu_splats(v1),v2.vec);
	return ret;
}

INLINE vec4 &mad(vec4 &ret,const vec4 &v0,const vec4 &v1,const vec4 &v2) {
	ret.vec = spu_madd(v0.vec,v1.vec,v2.vec);
	return ret;
}

INLINE vec4 &add(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.vec = spu_add(v0.vec,v1.vec);
	return ret;
}

INLINE vec4 &sub(vec4 &ret,const vec4 &v0,const vec4 &v1) {
	ret.vec = spu_sub(v0.vec,v1.vec);
	return ret;
}

INLINE vec4 &lerp(vec4 &ret,const vec4 &v0,const vec4 &v1,float k) {
	ret.vec = spu_madd(spu_sub(v1.vec,v0.vec),spu_splats(k),v0.vec);
	return ret;
}

INLINE vec4 &cross(vec4 &ret,const vec3 &v0,const vec3 &v1) {
	ret.vec = spu_cross(v0.vec,v1.vec);
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

INLINE vec4 lerp(const vec4 &v0,const vec4 &v1,float k) {
	vec4 ret;
	return lerp(ret,v0,v1,k);
}

/*
 */
INLINE vec4 min(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_sel(v0.vec,v1.vec,spu_cmpgt(v0.vec,v1.vec)));
}

INLINE vec4 max(const vec4 &v0,const vec4 &v1) {
	return vec4(spu_sel(v0.vec,v1.vec,spu_cmpgt(v1.vec,v0.vec)));
}

INLINE vec4 clamp(const vec4 &v,const vec4 &v0,const vec4 &v1) {
	vec_uint4 mask0 = spu_cmpgt(v.vec,v0.vec);
	vec_uint4 mask1 = spu_cmpgt(v.vec,v1.vec);
	return vec4(spu_sel(v0.vec,spu_sel(v.vec,v1.vec,mask1),mask0));
}

INLINE vec4 saturate(const vec4 &v) {
	vec_float4 one = spu_splats(1.0f);
	vec_float4 zero = spu_splats(0.0f);
	vec_uint4 mask0 = spu_cmpgt(v.vec,zero);
	vec_uint4 mask1 = spu_cmpgt(v.vec,one);
	return vec4(spu_sel(zero,spu_sel(v.vec,one,mask1),mask0));
}

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
	explicit mat3(const mat4 &m);
	explicit mat3(const quat &q);
	
	INLINE mat3 &operator=(const mat3 &v) {
		col0 = v.col0;
		col1 = v.col1;
		col2 = v.col2;
		return *this;
	}
	
	mat3 operator-() const;
	mat3 &operator*=(const mat3 &m);
	mat3 &operator+=(const mat3 &m);
	mat3 &operator-=(const mat3 &m);
	
	INLINE float &operator[](short i) {
		assert((unsigned short)i < 12 && "mat3::operator[](): bad index");
		return mat[i];
	}
	INLINE float operator[](short i) const {
		assert((unsigned short)i < 12 && "mat3::operator[](): bad index");
		return mat[i];
	}
	
	void setRow(short row,const vec3 &v);
	vec3 getRow(short row) const;
	
	void setColumn(short column,const vec3 &v);
	vec3 getColumn(short column) const;
	
	union {
		struct {
			float m00,m10,m20,m30;
			float m01,m11,m21,m31;
			float m02,m12,m22,m32;
		};
		float mat[12];
		struct {
			vec_float4 col0;
			vec_float4 col1;
			vec_float4 col2;
		};
	};
};

/*
 */
INLINE mat3::mat3(float v) {
	m00 = v;    m01 = v;    m02 = v;
	m10 = v;    m11 = v;    m12 = v;
	m20 = v;    m21 = v;    m22 = v;
	m30 = 0.0f; m31 = 0.0f; m32 = 0.0f;
}

INLINE mat3::mat3(const mat3 &m) {
	col0 = m.col0;
	col1 = m.col1;
	col2 = m.col2;
}

/*
 */
INLINE vec3 &mul(vec3 &ret,const mat3 &m,const vec3 &v) {
	vec_float4 res_0 = spu_mul(m.col0,SPU_SWIZZLE(v.vec,X,X,X,W));
	vec_float4 res_1 = spu_madd(m.col1,SPU_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
	ret.vec = spu_madd(m.col2,SPU_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	return ret;
}

INLINE mat3 &mul(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	vec_uchar16 xxxx = SPU_PERM2(X,X,X,X);
	vec_uchar16 yyyy = SPU_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = SPU_PERM2(Z,Z,Z,Z);
	ret.col0 = spu_mul(m0.col0,spu_shuffle(m1.col0,m1.col0,xxxx));
	ret.col1 = spu_mul(m0.col0,spu_shuffle(m1.col1,m1.col1,xxxx));
	ret.col2 = spu_mul(m0.col0,spu_shuffle(m1.col2,m1.col2,xxxx));
	ret.col0 = spu_madd(m0.col1,spu_shuffle(m1.col0,m1.col0,yyyy),ret.col0);
	ret.col1 = spu_madd(m0.col1,spu_shuffle(m1.col1,m1.col1,yyyy),ret.col1);
	ret.col2 = spu_madd(m0.col1,spu_shuffle(m1.col2,m1.col2,yyyy),ret.col2);
	ret.col0 = spu_madd(m0.col2,spu_shuffle(m1.col0,m1.col0,zzzz),ret.col0);
	ret.col1 = spu_madd(m0.col2,spu_shuffle(m1.col1,m1.col1,zzzz),ret.col1);
	ret.col2 = spu_madd(m0.col2,spu_shuffle(m1.col2,m1.col2,zzzz),ret.col2);
	return ret;
}

INLINE mat3 &add(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	ret.col0 = spu_add(m0.col0,m1.col0);
	ret.col1 = spu_add(m0.col1,m1.col1);
	ret.col2 = spu_add(m0.col2,m1.col2);
	return ret;
}

INLINE mat3 &sub(mat3 &ret,const mat3 &m0,const mat3 &m1) {
	ret.col0 = spu_sub(m0.col0,m1.col0);
	ret.col1 = spu_sub(m0.col1,m1.col1);
	ret.col2 = spu_sub(m0.col2,m1.col2);
	return ret;
}

/*
 */
INLINE mat3 mat3::operator-() const {
	mat3 ret;
	ret.col0 = (vec_float4)spu_xor((vec_uint4)col0,spu_splats(0x80000000));
	ret.col1 = (vec_float4)spu_xor((vec_uint4)col1,spu_splats(0x80000000));
	ret.col2 = (vec_float4)spu_xor((vec_uint4)col2,spu_splats(0x80000000));
	return ret;
}

INLINE mat3 &mat3::operator*=(const mat3 &m) {
	return mul(*this,mat3(*this),m);
}

INLINE mat3 &mat3::operator+=(const mat3 &m) {
	return add(*this,*this,m);
}

INLINE mat3 &mat3::operator-=(const mat3 &m) {
	return sub(*this,*this,m);
}

/*
 */
INLINE void mat3::setRow(short row,const vec3 &v) {
	assert((unsigned short)row < 3 && "mat3::setRow(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
}

INLINE vec3 mat3::getRow(short row) const {
	assert((unsigned short)row < 3 && "mat3::getRow(): bad row");
	return vec3(mat[row + 0],mat[row + 4],mat[row + 8]);
}

/*
 */
INLINE void mat3::setColumn(short column,const vec3 &v) {
	assert((unsigned short)column < 3 && "mat3::setColumn(): bad column");
	mat[column * 4 + 0] = v.x;
	mat[column * 4 + 1] = v.y;
	mat[column * 4 + 2] = v.z;
}

INLINE vec3 mat3::getColumn(short column) const {
	assert((unsigned short)column < 3 && "mat3::getColumn(): bad column");
	return vec3(mat[column * 4 + 0],mat[column * 4 + 1],mat[column * 4 + 2]);
}

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
	explicit mat4(const mat3 &m);
	explicit mat4(const quat &q);
	explicit mat4(const mat3 &m,const vec3 &v);
	explicit mat4(const quat &q,const vec3 &v);
	
	INLINE mat4 &operator=(const mat4 &v) {
		col0 = v.col0;
		col1 = v.col1;
		col2 = v.col2;
		col3 = v.col3;
		return *this;
	}
	
	mat4 operator-() const;
	mat4 &operator*=(const mat4 &m);
	mat4 &operator+=(const mat4 &m);
	mat4 &operator-=(const mat4 &m);
	
	INLINE float &operator[](short i) {
		assert((unsigned short)i < 16 && "mat4::operator[](): bad index");
		return mat[i];
	}
	INLINE float operator[](short i) const {
		assert((unsigned short)i < 16 && "mat4::operator[](): bad index");
		return mat[i];
	}
	
	void setRow(short row,const vec4 &v);
	vec4 getRow(short row) const;
	void setRow3(short row,const vec3 &v);
	vec3 getRow3(short row) const;
	
	void setColumn(short column,const vec4 &v);
	vec4 getColumn(short column) const;
	void setColumn3(short column,const vec3 &v);
	vec3 getColumn3(short column) const;
	
	void setZero();
	void setIdentity();
	void setTranslate(const vec3 &v);
	void setScale(const vec3 &v);
	
	union {
		struct {
			float m00,m10,m20,m30;
			float m01,m11,m21,m31;
			float m02,m12,m22,m32;
			float m03,m13,m23,m33;
		};
		float mat[16];
		struct {
			vec_float4 col0;
			vec_float4 col1;
			vec_float4 col2;
			vec_float4 col3;
		};
	};
};

/*
 */
INLINE mat3::mat3(const mat4 &m) {
	col0 = m.col0;
	col1 = m.col1;
	col2 = m.col2;
}

/*
 */
INLINE mat4::mat4(float v) {
	vec_float4 temp = spu_splats(v);
	col0 = temp;
	col1 = temp;
	col2 = temp;
	col3 = temp;
}

INLINE mat4::mat4(const mat3 &m) {
	m00 = m.m00; m01 = m.m01; m02 = m.m02; m03 = 0.0f;
	m10 = m.m10; m11 = m.m11; m12 = m.m12; m13 = 0.0f;
	m20 = m.m20; m21 = m.m21; m22 = m.m22; m23 = 0.0f;
	m30 = 0.0f;  m31 = 0.0f;  m32 = 0.0f;  m33 = 1.0f;
}

INLINE mat4::mat4(const mat4 &m) {
	col0 = m.col0;
	col1 = m.col1;
	col2 = m.col2;
	col3 = m.col3;
}

/*
 */
INLINE vec3 &mul(vec3 &ret,const mat4 &m,const vec3 &v) {
	vec_float4 res_0 = spu_madd(m.col0,SPU_SWIZZLE(v.vec,X,X,X,W),m.col3);
	vec_float4 res_1 = spu_madd(m.col1,SPU_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
	ret.vec = spu_madd(m.col2,SPU_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	return ret;
}

INLINE vec4 &mul(vec4 &ret,const mat4 &m,const vec4 &v) {
	vec_float4 res_0 = spu_mul(m.col0,SPU_SWIZZLE(v.vec,X,X,X,X));
	vec_float4 res_1 = spu_madd(m.col1,SPU_SWIZZLE(v.vec,Y,Y,Y,Y),res_0);
	vec_float4 res_2 = spu_madd(m.col2,SPU_SWIZZLE(v.vec,Z,Z,Z,Z),res_1);
	ret.vec = spu_madd(m.col3,SPU_SWIZZLE(v.vec,W,W,W,W),res_2);
	return ret;
}

INLINE vec3 &mul3(vec3 &ret,const mat4 &m,const vec3 &v) {
	vec_float4 res_0 = spu_mul(m.col0,SPU_SWIZZLE(v.vec,X,X,X,W));
	vec_float4 res_1 = spu_madd(m.col1,SPU_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
	ret.vec = spu_madd(m.col2,SPU_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	return ret;
}

INLINE vec4 &mul3(vec4 &ret,const mat4 &m,const vec4 &v) {
	vec_float4 res_0 = spu_mul(m.col0,SPU_SWIZZLE(v.vec,X,X,X,W));
	vec_float4 res_1 = spu_madd(m.col1,SPU_SWIZZLE(v.vec,Y,Y,Y,W),res_0);
	ret.vec = spu_madd(m.col2,SPU_SWIZZLE(v.vec,Z,Z,Z,W),res_1);
	ret.vec = spu_insert(spu_extract(v.vec,3),ret.vec,3);
	return ret;
}

INLINE mat4 &mul(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	vec_uchar16 xxxx = SPU_PERM2(X,X,X,X);
	vec_uchar16 yyyy = SPU_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = SPU_PERM2(Z,Z,Z,Z);
	vec_uchar16 wwww = SPU_PERM2(W,W,W,W);
	ret.col0 = spu_mul(m0.col0,spu_shuffle(m1.col0,m1.col0,xxxx));
	ret.col1 = spu_mul(m0.col0,spu_shuffle(m1.col1,m1.col1,xxxx));
	ret.col2 = spu_mul(m0.col0,spu_shuffle(m1.col2,m1.col2,xxxx));
	ret.col3 = spu_mul(m0.col0,spu_shuffle(m1.col3,m1.col3,xxxx));
	ret.col0 = spu_madd(m0.col1,spu_shuffle(m1.col0,m1.col0,yyyy),ret.col0);
	ret.col1 = spu_madd(m0.col1,spu_shuffle(m1.col1,m1.col1,yyyy),ret.col1);
	ret.col2 = spu_madd(m0.col1,spu_shuffle(m1.col2,m1.col2,yyyy),ret.col2);
	ret.col3 = spu_madd(m0.col1,spu_shuffle(m1.col3,m1.col3,yyyy),ret.col3);
	ret.col0 = spu_madd(m0.col2,spu_shuffle(m1.col0,m1.col0,zzzz),ret.col0);
	ret.col1 = spu_madd(m0.col2,spu_shuffle(m1.col1,m1.col1,zzzz),ret.col1);
	ret.col2 = spu_madd(m0.col2,spu_shuffle(m1.col2,m1.col2,zzzz),ret.col2);
	ret.col3 = spu_madd(m0.col2,spu_shuffle(m1.col3,m1.col3,zzzz),ret.col3);
	ret.col0 = spu_madd(m0.col3,spu_shuffle(m1.col0,m1.col0,wwww),ret.col0);
	ret.col1 = spu_madd(m0.col3,spu_shuffle(m1.col1,m1.col1,wwww),ret.col1);
	ret.col2 = spu_madd(m0.col3,spu_shuffle(m1.col2,m1.col2,wwww),ret.col2);
	ret.col3 = spu_madd(m0.col3,spu_shuffle(m1.col3,m1.col3,wwww),ret.col3);
	return ret;
}

INLINE mat4 &add(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	ret.col0 = spu_add(m0.col0,m1.col0);
	ret.col1 = spu_add(m0.col1,m1.col1);
	ret.col2 = spu_add(m0.col2,m1.col2);
	ret.col3 = spu_add(m0.col3,m1.col3);
	return ret;
}

INLINE mat4 &sub(mat4 &ret,const mat4 &m0,const mat4 &m1) {
	ret.col0 = spu_sub(m0.col0,m1.col0);
	ret.col1 = spu_sub(m0.col1,m1.col1);
	ret.col2 = spu_sub(m0.col2,m1.col2);
	ret.col3 = spu_sub(m0.col3,m1.col3);
	return ret;
}

/*
 */
INLINE mat4 mat4::operator-() const {
	mat4 ret;
	ret.col0 = (vec_float4)spu_xor((vec_uint4)col0,spu_splats(0x80000000));
	ret.col1 = (vec_float4)spu_xor((vec_uint4)col1,spu_splats(0x80000000));
	ret.col2 = (vec_float4)spu_xor((vec_uint4)col2,spu_splats(0x80000000));
	ret.col3 = (vec_float4)spu_xor((vec_uint4)col3,spu_splats(0x80000000));
	return ret;
}

INLINE mat4 &mat4::operator*=(const mat4 &m) {
	return mul(*this,mat4(*this),m);
}

INLINE mat4 &mat4::operator+=(const mat4 &m) {
	return add(*this,*this,m);
}

INLINE mat4 &mat4::operator-=(const mat4 &m) {
	return sub(*this,*this,m);
}

/*
 */
INLINE void mat4::setRow(short row,const vec4 &v) {
	assert((unsigned short)row < 4 && "mat4::setRow(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
	mat[row + 12] = v.w;
}

INLINE vec4 mat4::getRow(short row) const {
	assert((unsigned short)row < 4 && "mat4::getRow(): bad row");
	return vec4(mat[row + 0],mat[row + 4],mat[row + 8],mat[row + 12]);
}

INLINE void mat4::setRow3(short row,const vec3 &v) {
	assert((unsigned short)row < 4 && "mat4::setRow3(): bad row");
	mat[row + 0] = v.x;
	mat[row + 4] = v.y;
	mat[row + 8] = v.z;
}

INLINE vec3 mat4::getRow3(short row) const {
	assert((unsigned short)row < 4 && "mat4::getRow(): bad row");
	return vec3(mat[row + 0],mat[row + 4],mat[row + 8]);
}

/*
 */
INLINE void mat4::setColumn(short column,const vec4 &v) {
	assert((unsigned short)column < 4 && "mat4::setColumn(): bad column");
	mat[column * 4 + 0] = v.x;
	mat[column * 4 + 1] = v.y;
	mat[column * 4 + 2] = v.z;
	mat[column * 4 + 3] = v.w;
}

INLINE vec4 mat4::getColumn(short column) const {
	assert((unsigned short)column < 4 && "mat4::getColumn(): bad column");
	return vec4(mat[column * 4 + 0],mat[column * 4 + 1],mat[column * 4 + 2],mat[column * 4 + 3]);
}

INLINE void mat4::setColumn3(short column,const vec3 &v) {
	assert((unsigned short)column < 4 && "mat4::setColumn3(): bad column");
	mat[column * 4 + 0] = v.x;
	mat[column * 4 + 1] = v.y;
	mat[column * 4 + 2] = v.z;
}

INLINE vec3 mat4::getColumn3(short column) const {
	assert((unsigned short)column < 4 && "mat4::getColumn3(): bad column");
	return vec3(mat[column * 4 + 0],mat[column * 4 + 1],mat[column * 4 + 2]);
}

/*
 */
INLINE void mat4::setZero() {
	vec_float4 zero = spu_splats(0.0f);
	col0 = zero;
	col1 = zero;
	col2 = zero;
	col3 = zero;
}

INLINE void mat4::setIdentity() {
	col0 = SPU_FLOAT4(1.0f,0.0f,0.0f,0.0f);
	col1 = SPU_FLOAT4(0.0f,1.0f,0.0f,0.0f);
	col2 = SPU_FLOAT4(0.0f,0.0f,1.0f,0.0f);
	col3 = SPU_FLOAT4(0.0f,0.0f,0.0f,1.0f);
}

INLINE void mat4::setTranslate(const vec3 &v) {
	col0 = SPU_FLOAT4(1.0f,0.0f,0.0f,0.0f);
	col1 = SPU_FLOAT4(0.0f,1.0f,0.0f,0.0f);
	col2 = SPU_FLOAT4(0.0f,0.0f,1.0f,0.0f);
	col3 = spu_insert(1.0f,v.vec,3);
}

INLINE void mat4::setScale(const vec3 &v) {
	vec_float4 zero = spu_splats(0.0f);
	col0 = spu_insert(spu_extract(v.vec,0),zero,0);
	col1 = spu_insert(spu_extract(v.vec,1),zero,1);
	col2 = spu_insert(spu_extract(v.vec,2),zero,2);
	col3 = SPU_FLOAT4(0.0f,0.0f,0.0f,1.0f);
}

/*
 */
INLINE vec3 operator*(const mat4 &m,const vec3 &v) {
	vec3 ret;
	return mul(ret,m,v);
}

INLINE vec4 operator*(const mat4 &m,const vec4 &v) {
	vec4 ret;
	return mul(ret,m,v);
}

INLINE mat4 operator*(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return mul(ret,m0,m1);
}

INLINE mat4 operator+(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return add(ret,m0,m1);
}

INLINE mat4 operator-(const mat4 &m0,const mat4 &m1) {
	mat4 ret;
	return sub(ret,m0,m1);
}

/*
 */
INLINE mat4 &transpose(mat4 &ret,const mat4 &m) {
	vec_float4 res_0 = spu_shuffle(m.col0,m.col1,SPU_PERM2(X,Y,X,Y));
	vec_float4 res_1 = spu_shuffle(m.col0,m.col1,SPU_PERM2(Z,W,Z,W));
	vec_float4 res_2 = spu_shuffle(m.col2,m.col3,SPU_PERM2(X,Y,X,Y));
	vec_float4 res_3 = spu_shuffle(m.col2,m.col3,SPU_PERM2(Z,W,Z,W));
	ret.col0 = spu_shuffle(res_0,res_2,SPU_PERM2(X,Z,X,Z));
	ret.col1 = spu_shuffle(res_0,res_2,SPU_PERM2(Y,W,Y,W));
	ret.col2 = spu_shuffle(res_1,res_3,SPU_PERM2(X,Z,X,Z));
	ret.col3 = spu_shuffle(res_1,res_3,SPU_PERM2(Y,W,Y,W));
	return ret;
}

INLINE mat4 &inverse(mat4 &ret,const mat4 &m) {
	vec_uchar16 yxwz = SPU_PERM2(Y,X,W,Z);
	vec_uchar16 zwxy = SPU_PERM2(Z,W,X,Y);
	vec_float4 res_0 = spu_shuffle(m.col0,m.col1,SPU_PERM2(X,Y,X,Y));
	vec_float4 res_1 = spu_shuffle(m.col0,m.col1,SPU_PERM2(Z,W,Z,W));
	vec_float4 res_2 = spu_shuffle(m.col2,m.col3,SPU_PERM2(X,Y,X,Y));
	vec_float4 res_3 = spu_shuffle(m.col2,m.col3,SPU_PERM2(Z,W,Z,W));
	vec_float4 row_0 = spu_shuffle(res_0,res_2,SPU_PERM2(X,Z,X,Z));
	vec_float4 row_1 = spu_shuffle(res_2,res_0,SPU_PERM2(Y,W,Y,W));
	vec_float4 row_2 = spu_shuffle(res_1,res_3,SPU_PERM2(X,Z,X,Z));
	vec_float4 row_3 = spu_shuffle(res_3,res_1,SPU_PERM2(Y,W,Y,W));
	vec_float4 temp = spu_mul(row_2,row_3);
	temp = spu_shuffle(temp,temp,yxwz);
	res_0 = spu_mul(row_1,temp);
	res_1 = spu_mul(row_0,temp);
	temp = spu_shuffle(temp,temp,zwxy);
	res_0 = spu_sub(spu_mul(row_1,temp),res_0);
	res_1 = spu_sub(spu_mul(row_0,temp),res_1);
	res_1 = spu_shuffle(res_1,res_1,zwxy);
	temp = spu_mul(row_1,row_2);
	temp = spu_shuffle(temp,temp,yxwz);
	res_0 = spu_madd(row_3,temp,res_0);
	res_3 = spu_mul(row_0,temp);
	temp = spu_shuffle(temp,temp,zwxy);
	res_0 = spu_nmsub(row_3,temp,res_0);
	res_3 = spu_sub(spu_mul(row_0,temp),res_3);
	res_3 = spu_shuffle(res_3,res_3,zwxy);
	temp = spu_mul(spu_shuffle(row_1,row_1,zwxy),row_3);
	temp = spu_shuffle(temp,temp,yxwz);
	row_2 = spu_shuffle(row_2,row_2,zwxy);
	res_0 = spu_madd(row_2,temp,res_0);
	res_2 = spu_mul(row_0,temp);
	temp = spu_shuffle(temp,temp,zwxy);
	res_0 = spu_nmsub(row_2,temp,res_0);
	res_2 = spu_sub(spu_mul(row_0,temp),res_2);
	res_2 = spu_shuffle(res_2,res_2,zwxy);
	temp = spu_mul(row_0,row_1);
	temp = spu_shuffle(temp,temp,yxwz);
	res_2 = spu_madd(row_3,temp,res_2);
	res_3 = spu_sub(spu_mul(row_2,temp),res_3);
	temp = spu_shuffle(temp,temp,zwxy);
	res_2 = spu_sub(spu_mul(row_3,temp),res_2);
	res_3 = spu_nmsub(row_2,temp,res_3);
	temp = spu_mul(row_0,row_3);
	temp = spu_shuffle(temp,temp,yxwz);
	res_1 = spu_nmsub(row_2,temp,res_1);
	res_2 = spu_madd(row_1,temp,res_2);
	temp = spu_shuffle(temp,temp,zwxy);
	res_1 = spu_madd(row_2,temp,res_1);
	res_2 = spu_nmsub(row_1,temp,res_2);
	temp = spu_mul(row_0,row_2);
	temp = spu_shuffle(temp,temp,yxwz);
	res_1 = spu_madd(row_3,temp,res_1);
	res_3 = spu_nmsub(row_1,temp,res_3);
	temp = spu_shuffle(temp,temp,zwxy);
	res_1 = spu_nmsub(row_3,temp,res_1);
	res_3 = spu_madd(row_1,temp,res_3);
	vec_float4 det = spu_mul(row_0,res_0);
	det = spu_add(det,spu_rlqwbyte(det,8));
	det = spu_add(det,spu_rlqwbyte(det,4));
	temp = spu_rcp_nr(det);
	ret.col0 = spu_mul(res_0,temp);
	ret.col1 = spu_mul(res_1,temp);
	ret.col2 = spu_mul(res_2,temp);
	ret.col3 = spu_mul(res_3,temp);
	return ret;
}

/*
 */
INLINE mat4 transpose(const mat4 &m) {
	mat4 ret;
	return transpose(ret,m);
}

INLINE mat4 inverse(const mat4 &m) {
	mat4 ret;
	return inverse(ret,m);
}

/*
 */
INLINE mat4 translate(const vec3 &v) {
	mat4 ret;
	ret.setTranslate(v);
	return ret;
}

INLINE mat4 translate(float x,float y,float z) {
	return translate(vec3(x,y,z));
}

INLINE mat4 scale(const vec3 &v) {
	mat4 ret;
	ret.setScale(v);
	return ret;
}

INLINE mat4 scale(float x,float y,float z) {
	return scale(vec3(x,y,z));
}

/******************************************************************************\
*
* quat
*
\******************************************************************************/

/*
 */
ATTRIBUTE_ALIGNED16(struct) quat {
	
	INLINE quat() : x(0.0f), y(0.0f), z(0.0f), w(1.0f) { }
	INLINE quat(const quat &q) : x(q.x), y(q.y), z(q.z), w(q.w) { }
	explicit INLINE quat(const vec3 &v) : x(v.x), y(v.y), z(v.z), w(0.0f) { }
	explicit INLINE quat(const vec4 &v) : x(v.x), y(v.y), z(v.z), w(v.w) { }
	explicit INLINE quat(vec_float4 vec) : vec(vec) { }
	explicit quat(const mat3 &m);
	explicit quat(const mat4 &m);
	
	INLINE float &operator[](short i) {
		assert((unsigned short)i < 4 && "quat::operator[](): bad index");
		return q[i];
	}
	INLINE float operator[](short i) const {
		assert((unsigned short)i < 4 && "quat::operator[](): bad index");
		return q[i];
	}
	
	union {
		struct {
			float x,y,z,w;
		};
		float q[4];
		vec_float4 vec;
	};
};

#endif /* __MATH_LIB_SPU_H__ */
