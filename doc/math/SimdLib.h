/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    SimdLib.h
 * Desc:    Simd library
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

#ifndef __SIMD_LIB_H__
#define __SIMD_LIB_H__

#include "MathLib.h"

/******************************************************************************\
*
* SimdLib
*
\******************************************************************************/

/*
 */
class Simd {
		
		Simd();
		
	public:
		
		// shutdown
		static void shutdown();
		
		// codepath
		static void setGeneric();
		static void setSSE();
		static void setSSE2();
		static void setAltiVec();
		static void setNeon();
		
		// math
		static void dot(float &ret,const float *v0,const float *v1,int num);
		static void mul(float *ret,const float *v0,float v1,int num);
		static void mul(float *ret,const float *v0,const float *v1,int num);
		static void mad(float *ret,const float *v0,float v1,const float *v2,int num);
		static void add(float *ret,const float *v0,const float *v1,int num);
		static void sub(float *ret,const float *v0,const float *v1,int num);
		
		// bounds
		static void minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num);
		static void minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num);
		static void minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num);
		static void minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num);
		
		// transformed bounds
		static void minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num);
		static void minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num);
		
		// vector dot products
		static void dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num);
		static void dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num);
		
		// vector normalizations
		static void normalizeVec3(void *ret,int ret_stride,int num);
		static void normalizeVec4(void *ret,int ret_stride,int num);
		
		// vector multiplications
		static void mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num);
		static void mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num);
		static void mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num);
		static void mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num);
		static void mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num);
		static void mulMat4Vec4(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num);
		
		// vector projections
		static void projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num);
		static void projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num);
		
		// matrix multiplications
		static void mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num);
		static void mulMat4Mat4(vec4 *ret,const dmat4 &m,const dmat4 **matrices,int num);
		
		// matrix palette skinning
		static void skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src);
		static void skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src);
		
		// matrix decomposition
		static void eliminate(float *ret,const float *column,const float *factor,int rows,int num);
};

#endif /* __SIMD_LIB_H__ */
