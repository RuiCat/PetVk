/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    SimdLib.cpp
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

#include "Log.h"
#include "SimdLib.h"

/******************************************************************************\
*
* Generic simd processor
*
\******************************************************************************/

/*
 */
class SimdGeneric {
		
	public:
		
		virtual ~SimdGeneric() { }
		
		// codepath
		virtual const char *name() const { return "Generic"; }
		
		// math
		virtual void FASTCALL dot(float &ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,float v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mad(float *ret,const float *v0,float v1,const float *v2,int num) const;
		virtual void FASTCALL add(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL sub(float *ret,const float *v0,const float *v1,int num) const;
		
		// bounds
		virtual void FASTCALL minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num) const;
		
		// transformed bounds
		virtual void FASTCALL minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		
		// vector dot products
		virtual void FASTCALL dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const;
		virtual void FASTCALL dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const;
		
		// vector normalizations
		virtual void FASTCALL normalizeVec3(void *ret,int ret_stride,int num) const;
		virtual void FASTCALL normalizeVec4(void *ret,int ret_stride,int num) const;
		
		// vector multiplications
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec4(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const;
		
		// vector projections
		virtual void FASTCALL projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		// matrix multiplications
		virtual void FASTCALL mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const;
		virtual void FASTCALL mulMat4Mat4(vec4 *ret,const dmat4 &m,const dmat4 **matrices,int num) const;
		
		// matrix palette skinning
		virtual void FASTCALL skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		virtual void FASTCALL skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		
		// matrix decomposition
		virtual void FASTCALL eliminate(float *ret,const float *column,const float *factor,int rows,int num) const;
};

/*
 */
void SimdGeneric::dot(float &ret,const float *v0,const float *v1,int num) const {
	
	ret = 0.0f;
	
	if(num & ~3) {
		
		float res_0 = 0.0f;
		float res_1 = 0.0f;
		float res_2 = 0.0f;
		float res_3 = 0.0f;
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			res_0 += v0[0] * v1[0];
			res_1 += v0[1] * v1[1];
			res_2 += v0[2] * v1[2];
			res_3 += v0[3] * v1[3];
			
			v0 += 4;
			v1 += 4;
		}
		
		ret = res_0 + res_1 + res_2 + res_3;
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		ret += *v0++ * *v1++;
	}
}

void SimdGeneric::mul(float *ret,const float *v0,float v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			ret[0] = v0[0] * v1;
			ret[1] = v0[1] * v1;
			ret[2] = v0[2] * v1;
			ret[3] = v0[3] * v1;
			
			ret += 4;
			v0 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1;
	}
}

void SimdGeneric::mul(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			ret[0] = v0[0] * v1[0];
			ret[1] = v0[1] * v1[1];
			ret[2] = v0[2] * v1[2];
			ret[3] = v0[3] * v1[3];
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * *v1++;
	}
}

void SimdGeneric::mad(float *ret,const float *v0,float v1,const float *v2,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			ret[0] = v0[0] * v1 + v2[0];
			ret[1] = v0[1] * v1 + v2[1];
			ret[2] = v0[2] * v1 + v2[2];
			ret[3] = v0[3] * v1 + v2[3];
			
			ret += 4;
			v0 += 4;
			v2 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1 + *v2++;
	}
}

void SimdGeneric::add(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			ret[0] = v0[0] + v1[0];
			ret[1] = v0[1] + v1[1];
			ret[2] = v0[2] + v1[2];
			ret[3] = v0[3] + v1[3];
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ + *v1++;
	}
}

void SimdGeneric::sub(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			ret[0] = v0[0] - v1[0];
			ret[1] = v0[1] - v1[1];
			ret[2] = v0[2] - v1[2];
			ret[3] = v0[3] - v1[3];
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ - *v1++;
	}
}

/*
 */
void SimdGeneric::minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		const vec3 &v = *(const vec3*)s;
		if(min.x > v.x) min.x = v.x;
		if(max.x < v.x) max.x = v.x;
		if(min.y > v.y) min.y = v.y;
		if(max.y < v.y) max.y = v.y;
		if(min.z > v.z) min.z = v.z;
		if(max.z < v.z) max.z = v.z;
		
		s += src_stride;
	}
}

void SimdGeneric::minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		const vec4 &v = *(const vec4*)s;
		if(min.x > v.x) min.x = v.x;
		if(max.x < v.x) max.x = v.x;
		if(min.y > v.y) min.y = v.y;
		if(max.y < v.y) max.y = v.y;
		if(min.z > v.z) min.z = v.z;
		if(max.z < v.z) max.z = v.z;
		if(min.w > v.w) min.w = v.w;
		if(max.w < v.w) max.w = v.w;
		
		s += src_stride;
	}
}

void SimdGeneric::minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		const dvec3 &v = *(const dvec3*)s;
		if(min.x > v.x) min.x = v.x;
		if(max.x < v.x) max.x = v.x;
		if(min.y > v.y) min.y = v.y;
		if(max.y < v.y) max.y = v.y;
		if(min.z > v.z) min.z = v.z;
		if(max.z < v.z) max.z = v.z;
		
		s += src_stride;
	}
}

void SimdGeneric::minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		const dvec4 &v = *(const dvec4*)s;
		if(min.x > v.x) min.x = v.x;
		if(max.x < v.x) max.x = v.x;
		if(min.y > v.y) min.y = v.y;
		if(max.y < v.y) max.y = v.y;
		if(min.z > v.z) min.z = v.z;
		if(max.z < v.z) max.z = v.z;
		if(min.w > v.w) min.w = v.w;
		if(max.w < v.w) max.w = v.w;
		
		s += src_stride;
	}
}

/*
 */
void SimdGeneric::minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	vec3 res;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const vec3*)s);
		if(min.x > res.x) min.x = res.x;
		if(max.x < res.x) max.x = res.x;
		if(min.y > res.y) min.y = res.y;
		if(max.y < res.y) max.y = res.y;
		if(min.z > res.z) min.z = res.z;
		if(max.z < res.z) max.z = res.z;
		
		s += src_stride;
	}
}

void SimdGeneric::minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	vec4 res;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const vec4*)s);
		if(min.x > res.x) min.x = res.x;
		if(max.x < res.x) max.x = res.x;
		if(min.y > res.y) min.y = res.y;
		if(max.y < res.y) max.y = res.y;
		if(min.z > res.z) min.z = res.z;
		if(max.z < res.z) max.z = res.z;
		if(min.w > res.w) min.w = res.w;
		if(max.w < res.w) max.w = res.w;
		
		s += src_stride;
	}
}

/*
 */
void SimdGeneric::dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		*(float*)d = ::dot(v,*(const vec3*)s);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		*(float*)d = ::dot(v,*(const vec4*)s);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdGeneric::normalizeVec3(void *ret,int ret_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		((vec3*)d)->normalize();
		
		d += ret_stride;
	}
}

void SimdGeneric::normalizeVec4(void *ret,int ret_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		((vec4*)d)->normalize();
		
		d += ret_stride;
	}
}

/*
 */
void SimdGeneric::mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	vec3 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		mul3(res,m,*(const vec3*)s);
		*(vec3*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	vec3 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const vec3*)s);
		*(vec3*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	vec4 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const vec4*)s);
		*(vec4*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const {
	
	dvec3 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		mul3(res,m,*(const dvec3*)s);
		*(dvec3*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const {
	
	dvec3 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const dvec3*)s);
		*(dvec3*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::mulMat4Vec4(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const {
	
	dvec4 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::mul(res,m,*(const dvec4*)s);
		*(dvec4*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdGeneric::projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	vec3 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::proj(res,m,*(const vec3*)s);
		*(vec3*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdGeneric::projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	vec4 res;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		::proj(res,m,*(const vec4*)s);
		*(vec4*)d = res;
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdGeneric::mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const {
	
	float m00 = m.m00;
	float m01 = m.m01;
	float m02 = m.m02;
	float m03 = m.m03;
	float m10 = m.m10;
	float m11 = m.m11;
	float m12 = m.m12;
	float m13 = m.m13;
	float m20 = m.m20;
	float m21 = m.m21;
	float m22 = m.m22;
	float m23 = m.m23;
	
	for(int i = 0; i < num; i++) {
		const mat4 *m = matrices[i];
		
		ret[0].x = m00 * m->m00 + m01 * m->m10 + m02 * m->m20;
		ret[0].y = m00 * m->m01 + m01 * m->m11 + m02 * m->m21;
		ret[0].z = m00 * m->m02 + m01 * m->m12 + m02 * m->m22;
		ret[0].w = m00 * m->m03 + m01 * m->m13 + m02 * m->m23 + m03;
		
		ret[1].x = m10 * m->m00 + m11 * m->m10 + m12 * m->m20;
		ret[1].y = m10 * m->m01 + m11 * m->m11 + m12 * m->m21;
		ret[1].z = m10 * m->m02 + m11 * m->m12 + m12 * m->m22;
		ret[1].w = m10 * m->m03 + m11 * m->m13 + m12 * m->m23 + m13;
		
		ret[2].x = m20 * m->m00 + m21 * m->m10 + m22 * m->m20;
		ret[2].y = m20 * m->m01 + m21 * m->m11 + m22 * m->m21;
		ret[2].z = m20 * m->m02 + m21 * m->m12 + m22 * m->m22;
		ret[2].w = m20 * m->m03 + m21 * m->m13 + m22 * m->m23 + m23;
		
		ret += 3;
	}
}

void SimdGeneric::mulMat4Mat4(vec4 *ret,const dmat4 &m,const dmat4 **matrices,int num) const {
	
	double m00 = m.m00;
	double m01 = m.m01;
	double m02 = m.m02;
	double m03 = m.m03;
	double m10 = m.m10;
	double m11 = m.m11;
	double m12 = m.m12;
	double m13 = m.m13;
	double m20 = m.m20;
	double m21 = m.m21;
	double m22 = m.m22;
	double m23 = m.m23;
	
	for(int i = 0; i < num; i++) {
		const dmat4 *m = matrices[i];
		
		ret[0].x = (float)(m00 * m->m00 + m01 * m->m10 + m02 * m->m20);
		ret[0].y = (float)(m00 * m->m01 + m01 * m->m11 + m02 * m->m21);
		ret[0].z = (float)(m00 * m->m02 + m01 * m->m12 + m02 * m->m22);
		ret[0].w = (float)(m00 * m->m03 + m01 * m->m13 + m02 * m->m23 + m03);
		
		ret[1].x = (float)(m10 * m->m00 + m11 * m->m10 + m12 * m->m20);
		ret[1].y = (float)(m10 * m->m01 + m11 * m->m11 + m12 * m->m21);
		ret[1].z = (float)(m10 * m->m02 + m11 * m->m12 + m12 * m->m22);
		ret[1].w = (float)(m10 * m->m03 + m11 * m->m13 + m12 * m->m23 + m13);
		
		ret[2].x = (float)(m20 * m->m00 + m21 * m->m10 + m22 * m->m20);
		ret[2].y = (float)(m20 * m->m01 + m21 * m->m11 + m22 * m->m21);
		ret[2].z = (float)(m20 * m->m02 + m21 * m->m12 + m22 * m->m22);
		ret[2].w = (float)(m20 * m->m03 + m21 * m->m13 + m22 * m->m23 + m23);
		
		ret += 3;
	}
}

/*
 */
void SimdGeneric::skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	vec3 row_0,row_1,row_2;
	
	float weight = weights[0];
	const mat4 *matrix = matrices[0];
	row_0.x = matrix->m00 * weight;
	row_1.x = matrix->m10 * weight;
	row_2.x = matrix->m20 * weight;
	row_0.y = matrix->m01 * weight;
	row_1.y = matrix->m11 * weight;
	row_2.y = matrix->m21 * weight;
	row_0.z = matrix->m02 * weight;
	row_1.z = matrix->m12 * weight;
	row_2.z = matrix->m22 * weight;
	
	for(int i = 1; i < num; i++) {
		float weight = weights[i];
		const mat4 *matrix = matrices[i];
		row_0.x += matrix->m00 * weight;
		row_1.x += matrix->m10 * weight;
		row_2.x += matrix->m20 * weight;
		row_0.y += matrix->m01 * weight;
		row_1.y += matrix->m11 * weight;
		row_2.y += matrix->m21 * weight;
		row_0.z += matrix->m02 * weight;
		row_1.z += matrix->m12 * weight;
		row_2.z += matrix->m22 * weight;
	}
	
	ret.x = ::dot(src,row_0);
	ret.y = ::dot(src,row_1);
	ret.z = ::dot(src,row_2);
}

void SimdGeneric::skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	vec4 row_0,row_1,row_2;
	
	float weight = weights[0];
	const mat4 *matrix = matrices[0];
	row_0.x = matrix->m00 * weight;
	row_1.x = matrix->m10 * weight;
	row_2.x = matrix->m20 * weight;
	row_0.y = matrix->m01 * weight;
	row_1.y = matrix->m11 * weight;
	row_2.y = matrix->m21 * weight;
	row_0.z = matrix->m02 * weight;
	row_1.z = matrix->m12 * weight;
	row_2.z = matrix->m22 * weight;
	row_0.w = matrix->m03 * weight;
	row_1.w = matrix->m13 * weight;
	row_2.w = matrix->m23 * weight;
	
	for(int i = 1; i < num; i++) {
		float weight = weights[i];
		const mat4 *matrix = matrices[i];
		row_0.x += matrix->m00 * weight;
		row_1.x += matrix->m10 * weight;
		row_2.x += matrix->m20 * weight;
		row_0.y += matrix->m01 * weight;
		row_1.y += matrix->m11 * weight;
		row_2.y += matrix->m21 * weight;
		row_0.z += matrix->m02 * weight;
		row_1.z += matrix->m12 * weight;
		row_2.z += matrix->m22 * weight;
		row_0.w += matrix->m03 * weight;
		row_1.w += matrix->m13 * weight;
		row_2.w += matrix->m23 * weight;
	}
	
	ret.x = ::dot(src,row_0);
	ret.y = ::dot(src,row_1);
	ret.z = ::dot(src,row_2);
}

/*
 */
void SimdGeneric::eliminate(float *ret,const float *column,const float *factor,int rows,int num) const {
	
	if(num > 3) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		float *d2 = d1 + rows;
		float *d3 = d2 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		float f2 = factor[rows + rows];
		float f3 = factor[rows + rows + rows];
		
		for(size_t i = num; i > 0; i--) {
			float c = *column++;
			*d0++ -= c * f0;
			*d1++ -= c * f1;
			*d2++ -= c * f2;
			*d3++ -= c * f3;
		}
	}
	else if(num == 3) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		float *d2 = d1 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		float f2 = factor[rows + rows];
		
		float c0 = *column++;
		float c1 = *column++;
		float c2 = *column++;
		
		*d0++ -= c0 * f0;
		*d1++ -= c0 * f1;
		*d2++ -= c0 * f2;
		
		*d0++ -= c1 * f0;
		*d1++ -= c1 * f1;
		*d2++ -= c1 * f2;
		
		*d0++ -= c2 * f0;
		*d1++ -= c2 * f1;
		*d2++ -= c2 * f2;
	}
	else if(num == 2) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		
		float c0 = *column++;
		float c1 = *column++;
		
		*d0++ -= c0 * f0;
		*d1++ -= c0 * f1;
		
		*d0++ -= c1 * f0;
		*d1++ -= c1 * f1;
	}
	else if(num == 1) {
		
		*ret -= *column * *factor;
	}
}

/******************************************************************************\
*
* SSE simd processor
*
\******************************************************************************/

/*
 */
#ifdef USE_SSE

/*
 */
class SimdSSE : public SimdGeneric {
		
	public:
		
		virtual const char *name() const { return "SSE"; }
		
		virtual void FASTCALL dot(float &ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,float v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mad(float *ret,const float *v0,float v1,const float *v2,int num) const;
		virtual void FASTCALL add(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL sub(float *ret,const float *v0,const float *v1,int num) const;
		
		virtual void FASTCALL minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const;
		virtual void FASTCALL dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL normalizeVec3(void *ret,int ret_stride,int num) const;
		virtual void FASTCALL normalizeVec4(void *ret,int ret_stride,int num) const;
		
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const;
		
		virtual void FASTCALL skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		virtual void FASTCALL skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		
		virtual void FASTCALL eliminate(float *ret,const float *column,const float *factor,int rows,int num) const;
};

/*
 */
void SimdSSE::dot(float &ret,const float *v0,const float *v1,int num) const {
	
	ret = 0.0f;
	
	if(num & ~3) {
		
		__m128 res_0,res_1,res_2;
		
		if(IS_ALIGNED16(v0) && IS_ALIGNED16(v1)) {
			
			res_0 = _mm_load_ps(v0);
			res_1 = _mm_load_ps(v1);
			res_2 = _mm_mul_ps(res_0,res_1);
			
			v0 += 4;
			v1 += 4;
			
			for(size_t i = (num >> 2) - 1; i > 0; i--) {
				
				res_0 = _mm_load_ps(v0);
				res_1 = _mm_load_ps(v1);
				res_2 = _mm_add_ps(res_2,_mm_mul_ps(res_0,res_1));
				
				v0 += 4;
				v1 += 4;
			}
		}
		else {
			
			res_0 = _mm_loadu_ps(v0);
			res_1 = _mm_loadu_ps(v1);
			res_2 = _mm_mul_ps(res_0,res_1);
			
			v0 += 4;
			v1 += 4;
			
			for(size_t i = (num >> 2) - 1; i > 0; i--) {
				
				res_0 = _mm_loadu_ps(v0);
				res_1 = _mm_loadu_ps(v1);
				res_2 = _mm_add_ps(res_2,_mm_mul_ps(res_0,res_1));
				
				v0 += 4;
				v1 += 4;
			}
		}
		
		res_2 = _mm_add_ps(res_2,_MM_SWIZZLE(res_2,Y,X,W,Z));
		res_2 = _mm_add_ss(res_2,_MM_SWIZZLE(res_2,Z,W,X,Y));
		_mm_store_ss(&ret,res_2);
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		ret += *v0++ * *v1++;
	}
}

void SimdSSE::mul(float *ret,const float *v0,float v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		
		__m128 res_1 = _mm_set1_ps(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__m128 res_0 = _mm_load_ps(v0);
			_mm_stream_ps(ret,_mm_mul_ps(res_0,res_1));
			
			ret += 4;
			v0 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1;
	}
}

void SimdSSE::mul(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__m128 res_0 = _mm_load_ps(v0);
			__m128 res_1 = _mm_load_ps(v1);
			_mm_stream_ps(ret,_mm_mul_ps(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * *v1++;
	}
}

void SimdSSE::mad(float *ret,const float *v0,float v1,const float *v2,int num) const {
	
	if(num & ~3) {
		
		__m128 res_1 = _mm_set1_ps(v1);
		
		if(IS_ALIGNED16(ret) && IS_ALIGNED16(v0) && IS_ALIGNED16(v2)) {
			
			for(size_t i = num >> 2; i > 0; i--) {
				
				__m128 res_0 = _mm_load_ps(v0);
				__m128 res_2 = _mm_load_ps(v2);
				_mm_stream_ps(ret,_mm_add_ps(_mm_mul_ps(res_0,res_1),res_2));
				
				ret += 4;
				v0 += 4;
				v2 += 4;
			}
		}
		else {
			
			for(size_t i = num >> 2; i > 0; i--) {
				
				__m128 res_0 = _mm_loadu_ps(v0);
				__m128 res_2 = _mm_loadu_ps(v2);
				_mm_storeu_ps(ret,_mm_add_ps(_mm_mul_ps(res_0,res_1),res_2));
				
				ret += 4;
				v0 += 4;
				v2 += 4;
			}
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1 + *v2++;
	}
}

void SimdSSE::add(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__m128 res_0 = _mm_load_ps(v0);
			__m128 res_1 = _mm_load_ps(v1);
			_mm_stream_ps(ret,_mm_add_ps(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ + *v1++;
	}
}

void SimdSSE::sub(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__m128 res_0 = _mm_load_ps(v0);
			__m128 res_1 = _mm_load_ps(v1);
			_mm_stream_ps(ret,_mm_sub_ps(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ - *v1++;
	}
}

/*
 */
static INLINE void simd_sse_min_max_vec(__m128 &min,__m128 &max,const void *src,int src_stride,int num) {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
			
			__m128 res_0 = _mm_load_ps((const float*)s);
			s += src_stride;
			__m128 res_1 = _mm_load_ps((const float*)s);
			s += src_stride;
			__m128 res_2 = _mm_load_ps((const float*)s);
			s += src_stride;
			__m128 res_3 = _mm_load_ps((const float*)s);
			s += src_stride;
			
			min = _mm_min_ps(min,res_0);
			min = _mm_min_ps(min,res_1);
			min = _mm_min_ps(min,res_2);
			min = _mm_min_ps(min,res_3);
			
			max = _mm_max_ps(max,res_0);
			max = _mm_max_ps(max,res_1);
			max = _mm_max_ps(max,res_2);
			max = _mm_max_ps(max,res_3);
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		__m128 res = _mm_load_ps((const float*)s);
		
		min = _mm_min_ps(min,res);
		max = _mm_max_ps(max,res);
		
		s += src_stride;
	}
}

void SimdSSE::minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_sse_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

void SimdSSE::minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_sse_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

/*
 */
void SimdSSE::minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	__m128 min_vec = min.vec;
	__m128 max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,W));
		temp = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3));
		
		min_vec = _mm_min_ps(min_vec,temp);
		max_vec = _mm_max_ps(max_vec,temp);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

void SimdSSE::minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	__m128 min_vec = min.vec;
	__m128 max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(col_3,_MM_SWIZZLE(temp,W,W,W,W));
		temp = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		
		min_vec = _mm_min_ps(min_vec,temp);
		max_vec = _mm_max_ps(max_vec,temp);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

/*
 */
void SimdSSE::dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 vec = v.vec;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
			
			__m128 res_0 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_1 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_2 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_3 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			
			__m128 res_4 = _mm_add_ps(res_0,_MM_SWIZZLE(res_0,Y,X,Y,W));
			__m128 res_5 = _mm_add_ps(res_1,_MM_SWIZZLE(res_1,Y,X,Y,W));
			__m128 res_6 = _mm_add_ps(res_2,_MM_SWIZZLE(res_2,Y,X,Y,W));
			__m128 res_7 = _mm_add_ps(res_3,_MM_SWIZZLE(res_3,Y,X,Y,W));
			
			res_0 = _mm_add_ss(res_4,_MM_SWIZZLE(res_0,Z,Z,X,W));
			res_1 = _mm_add_ss(res_5,_MM_SWIZZLE(res_1,Z,Z,X,W));
			res_2 = _mm_add_ss(res_6,_MM_SWIZZLE(res_2,Z,Z,X,W));
			res_3 = _mm_add_ss(res_7,_MM_SWIZZLE(res_3,Z,Z,X,W));
			
			_mm_store_ss((float*)d,res_0);
			d += ret_stride;
			_mm_store_ss((float*)d,res_1);
			d += ret_stride;
			_mm_store_ss((float*)d,res_2);
			d += ret_stride;
			_mm_store_ss((float*)d,res_3);
			d += ret_stride;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		__m128 res_0 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
		__m128 res_1 = _mm_add_ps(res_0,_MM_SWIZZLE(res_0,Y,X,Y,W));
		res_0 = _mm_add_ss(res_1,_MM_SWIZZLE(res_0,Z,Z,X,W));
		
		_mm_store_ss((float*)d,res_0);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdSSE::dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 vec = v.vec;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
			
			__m128 res_0 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_1 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_2 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			__m128 res_3 = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
			s += src_stride;
			
			res_0 = _mm_add_ps(res_0,_MM_SWIZZLE(res_0,Y,X,W,Z));
			res_1 = _mm_add_ps(res_1,_MM_SWIZZLE(res_1,Y,X,W,Z));
			res_2 = _mm_add_ps(res_2,_MM_SWIZZLE(res_2,Y,X,W,Z));
			res_3 = _mm_add_ps(res_3,_MM_SWIZZLE(res_3,Y,X,W,Z));
			
			res_0 = _mm_add_ss(res_0,_MM_SWIZZLE(res_0,Z,W,X,Y));
			res_1 = _mm_add_ss(res_1,_MM_SWIZZLE(res_1,Z,W,X,Y));
			res_2 = _mm_add_ss(res_2,_MM_SWIZZLE(res_2,Z,W,X,Y));
			res_3 = _mm_add_ss(res_3,_MM_SWIZZLE(res_3,Z,W,X,Y));
			
			_mm_store_ss((float*)d,res_0);
			d += ret_stride;
			_mm_store_ss((float*)d,res_1);
			d += ret_stride;
			_mm_store_ss((float*)d,res_2);
			d += ret_stride;
			_mm_store_ss((float*)d,res_3);
			d += ret_stride;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		__m128 res = _mm_mul_ps(vec,_mm_load_ps((const float*)s));
		res = _mm_add_ps(res,_MM_SWIZZLE(res,Y,X,W,Z));
		res = _mm_add_ss(res,_MM_SWIZZLE(res,Z,W,X,Y));
		
		_mm_store_ss((float*)d,res);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdSSE::normalizeVec3(void *ret,int ret_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)d + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)d);
		
		temp = _mm_normalize3_ps(temp);
		
		_mm_store_ps((float*)d,temp);
		
		d += ret_stride;
	}
}

void SimdSSE::normalizeVec4(void *ret,int ret_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)d + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)d);
		
		temp = _mm_normalize4_ps(temp);
		
		_mm_store_ps((float*)d,temp);
		
		d += ret_stride;
	}
}

/*
 */
void SimdSSE::mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,W));
		temp = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
		
		_mm_stream_ps((float*)d,temp);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdSSE::mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,W));
		temp = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3));
		
		_mm_stream_ps((float*)d,temp);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdSSE::mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(col_3,_MM_SWIZZLE(temp,W,W,W,W));
		temp = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		
		_mm_stream_ps((float*)d,temp);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdSSE::projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,Z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3));
		temp = _mm_div_ps(res_3,_MM_SWIZZLE(res_3,W,W,W,W));
		
		_mm_stream_ps((float*)d,temp);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdSSE::projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128 col_0 = m.col0;
	__m128 col_1 = m.col1;
	__m128 col_2 = m.col2;
	__m128 col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128 temp = _mm_load_ps((const float*)s);
		
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(temp,X,X,X,X));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(temp,Y,Y,Y,Y));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(temp,Z,Z,Z,Z));
		__m128 res_3 = _mm_mul_ps(col_3,_MM_SWIZZLE(temp,W,W,W,W));
		__m128 res_4 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,res_3));
		temp = _mm_div_ps(res_4,_MM_SWIZZLE(res_4,W,W,W,W));
		
		_mm_stream_ps((float*)d,temp);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdSSE::mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const {
	
	__m128 col_0_x = _MM_SWIZZLE(m.col0,X,X,X,X);
	__m128 col_0_y = _MM_SWIZZLE(m.col0,Y,Y,Y,Y);
	__m128 col_0_z = _MM_SWIZZLE(m.col0,Z,Z,Z,Z);
	__m128 col_1_x = _MM_SWIZZLE(m.col1,X,X,X,X);
	__m128 col_1_y = _MM_SWIZZLE(m.col1,Y,Y,Y,Y);
	__m128 col_1_z = _MM_SWIZZLE(m.col1,Z,Z,Z,Z);
	__m128 col_2_x = _MM_SWIZZLE(m.col2,X,X,X,X);
	__m128 col_2_y = _MM_SWIZZLE(m.col2,Y,Y,Y,Y);
	__m128 col_2_z = _MM_SWIZZLE(m.col2,Z,Z,Z,Z);
	__m128 col_3_x = vec4(0.0f,0.0f,0.0f,m.m03).vec;
	__m128 col_3_y = vec4(0.0f,0.0f,0.0f,m.m13).vec;
	__m128 col_3_z = vec4(0.0f,0.0f,0.0f,m.m23).vec;
	
	for(int i = 0; i < num; i++) {
		const mat4 *m = matrices[i];
		
		if(i + 4 < num) _mm_prefetch((const char*)matrices[i + 4],_MM_HINT_NTA);
		
		__m128 res_0 = _mm_shuffle_ps(m->col0,m->col1,_MM_PERM2(X,Y,X,Y));
		__m128 res_1 = _mm_shuffle_ps(m->col0,m->col1,_MM_PERM2(Z,W,Z,W));
		__m128 res_2 = _mm_shuffle_ps(m->col2,m->col3,_MM_PERM2(X,Y,X,Y));
		__m128 res_3 = _mm_shuffle_ps(m->col2,m->col3,_MM_PERM2(Z,W,Z,W));
		__m128 col_0 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(X,Z,X,Z));
		__m128 col_1 = _mm_shuffle_ps(res_0,res_2,_MM_PERM2(Y,W,Y,W));
		__m128 col_2 = _mm_shuffle_ps(res_1,res_3,_MM_PERM2(X,Z,X,Z));
		
		res_0 = _mm_mul_ps(col_0,col_0_x);
		res_1 = _mm_mul_ps(col_1,col_1_x);
		res_2 = _mm_mul_ps(col_2,col_2_x);
		_mm_store_ps(ret[0].v,_mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3_x)));
		
		res_0 = _mm_mul_ps(col_0,col_0_y);
		res_1 = _mm_mul_ps(col_1,col_1_y);
		res_2 = _mm_mul_ps(col_2,col_2_y);
		_mm_store_ps(ret[1].v,_mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3_y)));
		
		res_0 = _mm_mul_ps(col_0,col_0_z);
		res_1 = _mm_mul_ps(col_1,col_1_z);
		res_2 = _mm_mul_ps(col_2,col_2_z);
		_mm_store_ps(ret[2].v,_mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3_z)));
		
		ret += 3;
	}
}

/*
 */
void SimdSSE::skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	__m128 weight = _mm_set1_ps(weights[0]);
	__m128 col_0 = _mm_mul_ps(matrices[0]->col0,weight);
	__m128 col_1 = _mm_mul_ps(matrices[0]->col1,weight);
	__m128 col_2 = _mm_mul_ps(matrices[0]->col2,weight);
	
	for(int i = 1; i < num; i++) {
		weight = _mm_set1_ps(weights[i]);
		col_0 = _mm_add_ps(col_0,_mm_mul_ps(matrices[i]->col0,weight));
		col_1 = _mm_add_ps(col_1,_mm_mul_ps(matrices[i]->col1,weight));
		col_2 = _mm_add_ps(col_2,_mm_mul_ps(matrices[i]->col2,weight));
	}
	
	__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(src.vec,X,X,X,W));
	__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(src.vec,Y,Y,Y,W));
	__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(src.vec,Z,Z,Z,W));
	ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),res_2);
}

void SimdSSE::skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	__m128 weight = _mm_set1_ps(weights[0]);
	__m128 col_0 = _mm_mul_ps(matrices[0]->col0,weight);
	__m128 col_1 = _mm_mul_ps(matrices[0]->col1,weight);
	__m128 col_2 = _mm_mul_ps(matrices[0]->col2,weight);
	__m128 col_3 = _mm_mul_ps(matrices[0]->col3,weight);
	
	for(int i = 1; i < num; i++) {
		weight = _mm_set1_ps(weights[i]);
		col_0 = _mm_add_ps(col_0,_mm_mul_ps(matrices[i]->col0,weight));
		col_1 = _mm_add_ps(col_1,_mm_mul_ps(matrices[i]->col1,weight));
		col_2 = _mm_add_ps(col_2,_mm_mul_ps(matrices[i]->col2,weight));
		col_3 = _mm_add_ps(col_3,_mm_mul_ps(matrices[i]->col3,weight));
	}
	
	__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(src.vec,X,X,X,W));
	__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(src.vec,Y,Y,Y,W));
	__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(src.vec,Z,Z,Z,W));
	ret.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3));
}

/*
 */
void SimdSSE::eliminate(float *ret,const float *column,const float *factor,int rows,int num) const {
	
	if(num > 3) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		float *d2 = d1 + rows;
		float *d3 = d2 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		float f2 = factor[rows + rows];
		float f3 = factor[rows + rows + rows];
		
		if(num & ~7) {
			
			__m128 factor_0 = _mm_set1_ps(f0);
			__m128 factor_1 = _mm_set1_ps(f1);
			__m128 factor_2 = _mm_set1_ps(f2);
			__m128 factor_3 = _mm_set1_ps(f3);
			
			for(size_t i = num >> 3; i > 0; i--) {
				
				__m128 col_0 = _mm_loadu_ps(column + 0);
				__m128 col_1 = _mm_loadu_ps(column + 4);
				
				__m128 res_0 = _mm_loadu_ps(d0 + 0);
				__m128 res_1 = _mm_loadu_ps(d0 + 4);
				res_0 = _mm_sub_ps(res_0,_mm_mul_ps(col_0,factor_0));
				res_1 = _mm_sub_ps(res_1,_mm_mul_ps(col_1,factor_0));
				_mm_storeu_ps(d0 + 0,res_0);
				_mm_storeu_ps(d0 + 4,res_1);
				
				res_0 = _mm_loadu_ps(d1 + 0);
				res_1 = _mm_loadu_ps(d1 + 4);
				res_0 = _mm_sub_ps(res_0,_mm_mul_ps(col_0,factor_1));
				res_1 = _mm_sub_ps(res_1,_mm_mul_ps(col_1,factor_1));
				_mm_storeu_ps(d1 + 0,res_0);
				_mm_storeu_ps(d1 + 4,res_1);
				
				res_0 = _mm_loadu_ps(d2 + 0);
				res_1 = _mm_loadu_ps(d2 + 4);
				res_0 = _mm_sub_ps(res_0,_mm_mul_ps(col_0,factor_2));
				res_1 = _mm_sub_ps(res_1,_mm_mul_ps(col_1,factor_2));
				_mm_storeu_ps(d2 + 0,res_0);
				_mm_storeu_ps(d2 + 4,res_1);
				
				res_0 = _mm_loadu_ps(d3 + 0);
				res_1 = _mm_loadu_ps(d3 + 4);
				res_0 = _mm_sub_ps(res_0,_mm_mul_ps(col_0,factor_3));
				res_1 = _mm_sub_ps(res_1,_mm_mul_ps(col_1,factor_3));
				_mm_storeu_ps(d3 + 0,res_0);
				_mm_storeu_ps(d3 + 4,res_1);
				
				column += 8;
				d0 += 8;
				d1 += 8;
				d2 += 8;
				d3 += 8;
			}
			
			num &= 7;
		}
		
		for(size_t i = num; i > 0; i--) {
			float c = *column++;
			*d0++ -= c * f0;
			*d1++ -= c * f1;
			*d2++ -= c * f2;
			*d3++ -= c * f3;
		}
	}
	else if(num == 3) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		float *d2 = d1 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		float f2 = factor[rows + rows];
		
		float c0 = *column++;
		float c1 = *column++;
		float c2 = *column++;
		
		*d0++ -= c0 * f0;
		*d1++ -= c0 * f1;
		*d2++ -= c0 * f2;
		
		*d0++ -= c1 * f0;
		*d1++ -= c1 * f1;
		*d2++ -= c1 * f2;
		
		*d0++ -= c2 * f0;
		*d1++ -= c2 * f1;
		*d2++ -= c2 * f2;
	}
	else if(num == 2) {
		
		float *d0 = ret;
		float *d1 = d0 + rows;
		
		float f0 = factor[0];
		float f1 = factor[rows];
		
		float c0 = *column++;
		float c1 = *column++;
		
		*d0++ -= c0 * f0;
		*d1++ -= c0 * f1;
		
		*d0++ -= c1 * f0;
		*d1++ -= c1 * f1;
	}
	else if(num == 1) {
		
		*ret -= *column * *factor;
	}
}

#endif /* USE_SSE */

/******************************************************************************\
*
* SSE2 simd processor
*
\******************************************************************************/

/*
 */
#ifdef USE_SSE2

/*
 */
class SimdSSE2 : public SimdSSE {
		
	public:
		
		virtual const char *name() const { return "SSE2"; }
		
		virtual void FASTCALL minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const;
};

/*
 */
static INLINE void simd_sse2_min_max_vec(__m128d &min0,__m128d &min1,__m128d &max0,__m128d &max1,const void *src,int src_stride,int num) {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
			
			__m128d res_00 = _mm_load_pd((const double*)s);
			__m128d res_01 = _mm_load_pd((const double*)s + 2);
			s += src_stride;
			__m128d res_10 = _mm_load_pd((const double*)s);
			__m128d res_11 = _mm_load_pd((const double*)s + 2);
			s += src_stride;
			__m128d res_20 = _mm_load_pd((const double*)s);
			__m128d res_21 = _mm_load_pd((const double*)s + 2);
			s += src_stride;
			__m128d res_30 = _mm_load_pd((const double*)s);
			__m128d res_31 = _mm_load_pd((const double*)s + 2);
			s += src_stride;
			
			min0 = _mm_min_pd(min0,res_00);
			min1 = _mm_min_pd(min1,res_01);
			min0 = _mm_min_pd(min0,res_10);
			min1 = _mm_min_pd(min1,res_11);
			min0 = _mm_min_pd(min0,res_20);
			min1 = _mm_min_pd(min1,res_21);
			min0 = _mm_min_pd(min0,res_30);
			min1 = _mm_min_pd(min1,res_31);
			
			max0 = _mm_max_pd(max0,res_00);
			max1 = _mm_max_pd(max1,res_01);
			max0 = _mm_max_pd(max0,res_10);
			max1 = _mm_max_pd(max1,res_11);
			max0 = _mm_max_pd(max0,res_20);
			max1 = _mm_max_pd(max1,res_21);
			max0 = _mm_max_pd(max0,res_30);
			max1 = _mm_max_pd(max1,res_31);
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		__m128d res0 = _mm_load_pd((const double*)s);
		__m128d res1 = _mm_load_pd((const double*)s + 2);
		
		min0 = _mm_min_pd(min0,res0);
		min1 = _mm_min_pd(min1,res1);
		max0 = _mm_max_pd(max0,res0);
		max1 = _mm_max_pd(max1,res1);
		
		s += src_stride;
	}
}

void SimdSSE2::minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_sse2_min_max_vec(min.vec0,min.vec1,max.vec0,max.vec1,src,src_stride,num);
}

void SimdSSE2::minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_sse2_min_max_vec(min.vec0,min.vec1,max.vec0,max.vec1,src,src_stride,num);
}

/*
 */
void SimdSSE2::mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128d col_0 = m.col0;
	__m128d col_1 = m.col1;
	__m128d col_12 = _mm_shuffle_pd(m.col1,m.col2,_MM_PERM22(Y,X));
	__m128d col_22 = _mm_shuffle_pd(m.col2,m.col2,_MM_PERM22(Y,Y));
	__m128d col_3 = m.col3;
	__m128d col_4 = m.col4;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128d temp0 = _mm_load_pd((const double*)s);
		__m128d temp1 = _mm_load_sd((const double*)s + 2);
		
		__m128d x = _MM_SWIZZLE2(temp0,X,X);
		__m128d y = _MM_SWIZZLE2(temp0,Y,Y);
		__m128d z = _MM_SWIZZLE2(temp1,X,X);
		
		__m128d res_0 = _mm_mul_pd(col_0,x);
		__m128d res_1 = _mm_mul_pd(col_12,y);
		__m128d res_2 = _mm_mul_pd(col_3,z);
		
		temp0 = _mm_add_pd(_mm_add_pd(res_0,res_1),res_2);
		
		res_0 = _mm_mul_sd(col_1,x);
		res_1 = _mm_mul_sd(col_22,y);
		res_2 = _mm_mul_sd(col_4,z);
		
		temp1 = _mm_add_sd(_mm_add_sd(res_0,res_1),res_2);
		
		_mm_stream_pd((double*)d,temp0);
		_mm_stream_pd((double*)d + 2,temp1);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdSSE2::mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	__m128d col_0 = m.col0;
	__m128d col_1 = m.col1;
	__m128d col_12 = _mm_shuffle_pd(m.col1,m.col2,_MM_PERM22(Y,X));
	__m128d col_22 = _mm_shuffle_pd(m.col2,m.col2,_MM_PERM22(Y,Y));
	__m128d col_3 = m.col3;
	__m128d col_4 = m.col4;
	__m128d col_45 = _mm_shuffle_pd(m.col4,m.col5,_MM_PERM22(Y,X));
	__m128d col_55 = _mm_shuffle_pd(m.col5,m.col5,_MM_PERM22(Y,Y));
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		_mm_prefetch((const char*)s + 1024,_MM_HINT_NTA);
		
		__m128d temp0 = _mm_load_pd((const double*)s);
		__m128d temp1 = _mm_load_sd((const double*)s + 2);
		
		__m128d x = _MM_SWIZZLE2(temp0,X,X);
		__m128d y = _MM_SWIZZLE2(temp0,Y,Y);
		__m128d z = _MM_SWIZZLE2(temp1,X,X);
		
		__m128d res_0 = _mm_mul_pd(col_0,x);
		__m128d res_1 = _mm_mul_pd(col_12,y);
		__m128d res_2 = _mm_mul_pd(col_3,z);
		
		temp0 = _mm_add_pd(_mm_add_pd(res_0,res_1),_mm_add_pd(res_2,col_45));
		
		res_0 = _mm_mul_sd(col_1,x);
		res_1 = _mm_mul_sd(col_22,y);
		res_2 = _mm_mul_sd(col_4,z);
		
		temp1 = _mm_add_sd(_mm_add_sd(res_0,res_1),_mm_add_sd(res_2,col_55));
		
		_mm_stream_pd((double*)d,temp0);
		_mm_stream_pd((double*)d + 2,temp1);
		
		s += src_stride;
		d += ret_stride;
	}
}

#endif /* USE_SSE2 */

/******************************************************************************\
*
* AltiVec simd processor
*
\******************************************************************************/

/*
 */
#ifdef USE_ALTIVEC

/*
 */
class SimdAltiVec : public SimdGeneric {
		
	public:
		
		virtual const char *name() const { return "AltiVec"; }
		
		virtual void FASTCALL dot(float &ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,float v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mad(float *ret,const float *v0,float v1,const float *v2,int num) const;
		virtual void FASTCALL add(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL sub(float *ret,const float *v0,const float *v1,int num) const;
		
		virtual void FASTCALL minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const;
		virtual void FASTCALL dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL normalizeVec3(void *ret,int ret_stride,int num) const;
		virtual void FASTCALL normalizeVec4(void *ret,int ret_stride,int num) const;
		
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const;
		
		virtual void FASTCALL skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		virtual void FASTCALL skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
};

/*
 */
void SimdAltiVec::dot(float &ret,const float *v0,const float *v1,int num) const {
	
	ret = 0.0f;
	
	if(num & ~3) {
		
		if(IS_ALIGNED16(v0) && IS_ALIGNED16(v1)) {
			
			vec_float4 res_0 = vec_ld(0,v0);
			vec_float4 res_1 = vec_ld(0,v1);
			vec_float4 res_2 = vec_madd(res_0,res_1,vec_splats(0.0f));
			
			v0 += 4;
			v1 += 4;
			
			for(size_t i = (num >> 2) - 1; i > 0; i--) {
				
				res_0 = vec_ld(0,v0);
				res_1 = vec_ld(0,v1);
				res_2 = vec_madd(res_0,res_1,res_2);
				
				v0 += 4;
				v1 += 4;
			}
			
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,8));
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,4));
			
			vec4 temp;
			vec_ste(res_2,0,&temp.x);
			ret = temp.x;
		}
		else {
			
			vec_uchar16 mask_0 = vec_lvsl(0,(const unsigned char*)v0);
			vec_uchar16 mask_1 = vec_lvsl(0,(const unsigned char*)v1);
			
			vec_float4 res_0 = vec_perm(vec_ld(0,v0),vec_ld(15,v0),mask_0);
			vec_float4 res_1 = vec_perm(vec_ld(0,v1),vec_ld(15,v1),mask_1);
			vec_float4 res_2 = vec_madd(res_0,res_1,vec_splats(0.0f));
			
			v0 += 4;
			v1 += 4;
			
			for(size_t i = (num >> 2) - 1; i > 0; i--) {
				
				vec_float4 res_0 = vec_perm(vec_ld(0,v0),vec_ld(15,v0),mask_0);
				vec_float4 res_1 = vec_perm(vec_ld(0,v1),vec_ld(15,v1),mask_1);
				res_2 = vec_madd(res_0,res_1,res_2);
				
				v0 += 4;
				v1 += 4;
			}
			
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,8));
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,4));
			
			vec4 temp;
			vec_ste(res_2,0,&temp.x);
			ret = temp.x;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		ret += *v0++ * *v1++;
	}
}

void SimdAltiVec::mul(float *ret,const float *v0,float v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		
		vec_float4 res_1 = vec_splats(v1);
		vec_float4 zero = vec_splats(0.0f);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			vec_float4 res_0 = vec_ld(0,v0);
			vec_st(vec_madd(res_0,res_1,zero),0,ret);
			
			ret += 4;
			v0 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1;
	}
}

void SimdAltiVec::mul(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		vec_float4 zero = vec_splats(0.0f);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			vec_float4 res_0 = vec_ld(0,v0);
			vec_float4 res_1 = vec_ld(0,v1);
			vec_st(vec_madd(res_0,res_1,zero),0,ret);
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * *v1++;
	}
}

void SimdAltiVec::mad(float *ret,const float *v0,float v1,const float *v2,int num) const {
	
	if(num & ~3) {
		
		if(IS_ALIGNED16(ret) && IS_ALIGNED16(v0) && IS_ALIGNED16(v2)) {
			
			vec_float4 res_1 = vec_splats(v1);
			
			for(size_t i = num >> 2; i > 0; i--) {
				
				vec_float4 res_0 = vec_ld(0,v0);
				vec_float4 res_2 = vec_ld(0,v2);
				vec_st(vec_madd(res_0,res_1,res_2),0,ret);
				
				ret += 4;
				v0 += 4;
				v2 += 4;
			}
		}
		else {
			
			for(size_t i = num >> 2; i > 0; i--) {
				
				ret[0] = v0[0] * v1 + v2[0];
				ret[1] = v0[1] * v1 + v2[1];
				ret[2] = v0[2] * v1 + v2[2];
				ret[3] = v0[3] * v1 + v2[3];
				
				ret += 4;
				v0 += 4;
				v2 += 4;
			}
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1 + *v2++;
	}
}

void SimdAltiVec::add(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			vec_float4 res_0 = vec_ld(0,v0);
			vec_float4 res_1 = vec_ld(0,v1);
			vec_st(vec_add(res_0,res_1),0,ret);
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ + *v1++;
	}
}

void SimdAltiVec::sub(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		ASSERT_ALIGNED16(ret);
		ASSERT_ALIGNED16(v0);
		ASSERT_ALIGNED16(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			vec_float4 res_0 = vec_ld(0,v0);
			vec_float4 res_1 = vec_ld(0,v1);
			vec_st(vec_sub(res_0,res_1),0,ret);
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ - *v1++;
	}
}

/*
 */
static INLINE void simd_altivec_min_max_vec(vec_float4 &min,vec_float4 &max,const void *src,int src_stride,int num) {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__dcbt((const void*)(s + 1024));
			
			vec_float4 res_0 = vec_ld(0,(const float*)s);
			s += src_stride;
			vec_float4 res_1 = vec_ld(0,(const float*)s);
			s += src_stride;
			vec_float4 res_2 = vec_ld(0,(const float*)s);
			s += src_stride;
			vec_float4 res_3 = vec_ld(0,(const float*)s);
			s += src_stride;
			
			min = vec_min(min,res_0);
			min = vec_min(min,res_1);
			min = vec_min(min,res_2);
			min = vec_min(min,res_3);
			
			max = vec_max(max,res_0);
			max = vec_max(max,res_1);
			max = vec_max(max,res_2);
			max = vec_max(max,res_3);
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		vec_float4 res = vec_ld(0,(const float*)s);
		
		min = vec_min(min,res);
		max = vec_max(max,res);
		
		s += src_stride;
	}
}

void SimdAltiVec::minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_altivec_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

void SimdAltiVec::minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_altivec_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

/*
 */
void SimdAltiVec::minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_uchar16 xxxw = VEC_PERM2(X,X,X,W);
	vec_uchar16 yyyw = VEC_PERM2(Y,Y,Y,W);
	vec_uchar16 zzzw = VEC_PERM2(Z,Z,Z,W);
	
	vec_float4 min_vec = min.vec;
	vec_float4 max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxw),col_3);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyw),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzw),res_1);
		
		min_vec = vec_min(min_vec,res_2);
		max_vec = vec_max(max_vec,res_2);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

void SimdAltiVec::minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_float4 zero = vec_splats(0.0f);
	vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
	vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
	vec_uchar16 wwww = VEC_PERM2(W,W,W,W);
	
	vec_float4 min_vec = min.vec;
	vec_float4 max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxx),zero);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyy),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzz),res_1);
		vec_float4 res_3 = vec_madd(col_3,vec_perm(temp,temp,wwww),res_2);
		
		min_vec = vec_min(min_vec,res_3);
		max_vec = vec_max(max_vec,res_3);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

/*
 */
void SimdAltiVec::dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 vec = v.vec;
	vec_float4 zero = vec_splats(0.0f);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__dcbt((const void*)(s + 1024));
			
			vec_float4 res_0 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_1 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_2 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_3 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			
			vec_float4 res_4 = vec_add(res_0,VEC_SWIZZLE(res_0,Y,X,Y,W));
			vec_float4 res_5 = vec_add(res_1,VEC_SWIZZLE(res_1,Y,X,Y,W));
			vec_float4 res_6 = vec_add(res_2,VEC_SWIZZLE(res_2,Y,X,Y,W));
			vec_float4 res_7 = vec_add(res_3,VEC_SWIZZLE(res_3,Y,X,Y,W));
			
			res_0 = vec_add(res_4,VEC_SWIZZLE(res_0,Z,Z,X,W));
			res_1 = vec_add(res_5,VEC_SWIZZLE(res_1,Z,Z,X,W));
			res_2 = vec_add(res_6,VEC_SWIZZLE(res_2,Z,Z,X,W));
			res_3 = vec_add(res_7,VEC_SWIZZLE(res_3,Z,Z,X,W));
			
			vec_ste(res_0,0,(float*)d);
			d += ret_stride;
			vec_ste(res_1,0,(float*)d);
			d += ret_stride;
			vec_ste(res_2,0,(float*)d);
			d += ret_stride;
			vec_ste(res_3,0,(float*)d);
			d += ret_stride;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		vec_float4 res_0 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
		vec_float4 res_1 = vec_add(res_0,VEC_SWIZZLE(res_0,Y,X,Y,W));
		res_0 = vec_add(res_1,VEC_SWIZZLE(res_0,Z,Z,X,W));
		
		vec_ste(res_0,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdAltiVec::dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 vec = v.vec;
	vec_float4 zero = vec_splats(0.0f);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			__dcbt((const void*)(s + 1024));
			
			vec_float4 res_0 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_1 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_2 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			vec_float4 res_3 = vec_madd(vec,vec_ld(0,(const float*)s),zero);
			s += src_stride;
			
			res_0 = vec_add(res_0,vec_sld(res_0,res_0,8));
			res_1 = vec_add(res_1,vec_sld(res_1,res_1,8));
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,8));
			res_3 = vec_add(res_3,vec_sld(res_3,res_3,8));
			
			res_0 = vec_add(res_0,vec_sld(res_0,res_0,4));
			res_1 = vec_add(res_1,vec_sld(res_1,res_1,4));
			res_2 = vec_add(res_2,vec_sld(res_2,res_2,4));
			res_3 = vec_add(res_3,vec_sld(res_3,res_3,4));
			
			vec_ste(res_0,0,(float*)d);
			d += ret_stride;
			vec_ste(res_1,0,(float*)d);
			d += ret_stride;
			vec_ste(res_2,0,(float*)d);
			d += ret_stride;
			vec_ste(res_3,0,(float*)d);
			d += ret_stride;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		vec_float4 res = vec_madd(vec,vec_ld(0,(const float*)s),zero);
		res = vec_add(res,vec_sld(res,res,8));
		res = vec_add(res,vec_sld(res,res,4));
		
		vec_ste(res,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdAltiVec::normalizeVec3(void *ret,int ret_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		__dcbtst((void*)(d + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)d);
		
		temp = vec_normalize3(temp);
		
		vec_st(temp,0,(float*)d);
		
		d += ret_stride;
	}
}

void SimdAltiVec::normalizeVec4(void *ret,int ret_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		__dcbtst((void*)(d + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)d);
		
		temp = vec_normalize4(temp);
		
		vec_st(temp,0,(float*)d);
		
		d += ret_stride;
	}
}

/*
 */
void SimdAltiVec::mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	
	vec_float4 zero = vec_splats(0.0f);
	vec_uchar16 xxxw = VEC_PERM2(X,X,X,W);
	vec_uchar16 yyyw = VEC_PERM2(Y,Y,Y,W);
	vec_uchar16 zzzw = VEC_PERM2(Z,Z,Z,W);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxw),zero);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyw),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzw),res_1);
		
		vec_st(res_2,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdAltiVec::mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_uchar16 xxxw = VEC_PERM2(X,X,X,W);
	vec_uchar16 yyyw = VEC_PERM2(Y,Y,Y,W);
	vec_uchar16 zzzw = VEC_PERM2(Z,Z,Z,W);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxw),col_3);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyw),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzw),res_1);
		
		vec_st(res_2,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdAltiVec::mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_float4 zero = vec_splats(0.0f);
	vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
	vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
	vec_uchar16 wwww = VEC_PERM2(W,W,W,W);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxx),zero);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyy),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzz),res_1);
		vec_float4 res_3 = vec_madd(col_3,vec_perm(temp,temp,wwww),res_2);
		
		vec_st(res_3,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdAltiVec::projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_float4 zero = vec_splats(0.0f);
	vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
	vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
	vec_uchar16 wwww = VEC_PERM2(W,W,W,W);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxx),col_3);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyy),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzz),res_1);
		vec_float4 res_3 = vec_madd(res_2,vec_rcp_nr(vec_perm(res_2,res_2,wwww)),zero);
		
		vec_st(res_3,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdAltiVec::projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	ASSERT_ALIGNED16(ret);
	ASSERT_ALIGNED16(ret_stride);
	ASSERT_ALIGNED16(src);
	ASSERT_ALIGNED16(src_stride);
	
	vec_float4 col_0 = m.col0;
	vec_float4 col_1 = m.col1;
	vec_float4 col_2 = m.col2;
	vec_float4 col_3 = m.col3;
	
	vec_float4 zero = vec_splats(0.0f);
	vec_uchar16 xxxx = VEC_PERM2(X,X,X,X);
	vec_uchar16 yyyy = VEC_PERM2(Y,Y,Y,Y);
	vec_uchar16 zzzz = VEC_PERM2(Z,Z,Z,Z);
	vec_uchar16 wwww = VEC_PERM2(W,W,W,W);
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		__dcbt((const void*)(s + 1024));
		
		vec_float4 temp = vec_ld(0,(const float*)s);
		
		vec_float4 res_0 = vec_madd(col_0,vec_perm(temp,temp,xxxx),zero);
		vec_float4 res_1 = vec_madd(col_1,vec_perm(temp,temp,yyyy),res_0);
		vec_float4 res_2 = vec_madd(col_2,vec_perm(temp,temp,zzzz),res_1);
		vec_float4 res_3 = vec_madd(col_3,vec_perm(temp,temp,wwww),res_2);
		vec_float4 res_4 = vec_madd(res_3,vec_rcp_nr(vec_perm(res_3,res_3,wwww)),zero);
		
		vec_st(res_4,0,(float*)d);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdAltiVec::mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const {
	
	vec_float4 col_0_x = VEC_SWIZZLE(m.col0,X,X,X,X);
	vec_float4 col_0_y = VEC_SWIZZLE(m.col0,Y,Y,Y,Y);
	vec_float4 col_0_z = VEC_SWIZZLE(m.col0,Z,Z,Z,Z);
	vec_float4 col_1_x = VEC_SWIZZLE(m.col1,X,X,X,X);
	vec_float4 col_1_y = VEC_SWIZZLE(m.col1,Y,Y,Y,Y);
	vec_float4 col_1_z = VEC_SWIZZLE(m.col1,Z,Z,Z,Z);
	vec_float4 col_2_x = VEC_SWIZZLE(m.col2,X,X,X,X);
	vec_float4 col_2_y = VEC_SWIZZLE(m.col2,Y,Y,Y,Y);
	vec_float4 col_2_z = VEC_SWIZZLE(m.col2,Z,Z,Z,Z);
	vec_float4 col_3_x = vec4(0.0f,0.0f,0.0f,m.m03).vec;
	vec_float4 col_3_y = vec4(0.0f,0.0f,0.0f,m.m13).vec;
	vec_float4 col_3_z = vec4(0.0f,0.0f,0.0f,m.m23).vec;
	
	for(int i = 0; i < num; i++) {
		const mat4 *m = matrices[i];
		
		if(i + 4 < num) __dcbt((const void*)matrices[i + 4]);
		
		vec_float4 res_0 = vec_perm(m->col0,m->col1,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_1 = vec_perm(m->col0,m->col1,VEC_PERM2(Z,W,Z,W));
		vec_float4 res_2 = vec_perm(m->col2,m->col3,VEC_PERM2(X,Y,X,Y));
		vec_float4 res_3 = vec_perm(m->col2,m->col3,VEC_PERM2(Z,W,Z,W));
		vec_float4 col_0 = vec_perm(res_0,res_2,VEC_PERM2(X,Z,X,Z));
		vec_float4 col_1 = vec_perm(res_0,res_2,VEC_PERM2(Y,W,Y,W));
		vec_float4 col_2 = vec_perm(res_1,res_3,VEC_PERM2(X,Z,X,Z));
		
		res_0 = vec_madd(col_0,col_0_x,col_3_x);
		res_1 = vec_madd(col_0,col_0_y,col_3_y);
		res_2 = vec_madd(col_0,col_0_z,col_3_z);
		res_0 = vec_madd(col_1,col_1_x,res_0);
		res_1 = vec_madd(col_1,col_1_y,res_1);
		res_2 = vec_madd(col_1,col_1_z,res_2);
		res_0 = vec_madd(col_2,col_2_x,res_0);
		res_1 = vec_madd(col_2,col_2_y,res_1);
		res_2 = vec_madd(col_2,col_2_z,res_2);
		
		ret[0].vec = res_0;
		ret[1].vec = res_1;
		ret[2].vec = res_2;
		
		ret += 3;
	}
}

/*
 */
void SimdAltiVec::skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	vec_float4 zero = vec_splats(0.0f);
	vec_float4 weight = vec_splats(weights[0]);
	vec_float4 col_0 = vec_madd(matrices[0]->col0,weight,zero);
	vec_float4 col_1 = vec_madd(matrices[0]->col1,weight,zero);
	vec_float4 col_2 = vec_madd(matrices[0]->col2,weight,zero);
	
	for(int i = 1; i < num; i++) {
		weight = vec_splats(weights[i]);
		col_0 = vec_madd(matrices[i]->col0,weight,col_0);
		col_1 = vec_madd(matrices[i]->col1,weight,col_1);
		col_2 = vec_madd(matrices[i]->col2,weight,col_2);
	}
	
	vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(src.vec,X,X,X,W),vec_splats(0.0f));
	vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(src.vec,Y,Y,Y,W),res_0);
	ret.vec = vec_madd(col_2,VEC_SWIZZLE(src.vec,Z,Z,Z,W),res_1);
}

void SimdAltiVec::skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	vec_float4 zero = vec_splats(0.0f);
	vec_float4 weight = vec_splats(weights[0]);
	vec_float4 col_0 = vec_madd(matrices[0]->col0,weight,zero);
	vec_float4 col_1 = vec_madd(matrices[0]->col1,weight,zero);
	vec_float4 col_2 = vec_madd(matrices[0]->col2,weight,zero);
	vec_float4 col_3 = vec_madd(matrices[0]->col3,weight,zero);
	
	for(int i = 1; i < num; i++) {
		weight = vec_splats(weights[i]);
		col_0 = vec_madd(matrices[i]->col0,weight,col_0);
		col_1 = vec_madd(matrices[i]->col1,weight,col_1);
		col_2 = vec_madd(matrices[i]->col2,weight,col_2);
		col_3 = vec_madd(matrices[i]->col3,weight,col_3);
	}
	
	vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(src.vec,X,X,X,W),col_3);
	vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(src.vec,Y,Y,Y,W),res_0);
	ret.vec = vec_madd(col_2,VEC_SWIZZLE(src.vec,Z,Z,Z,W),res_1);
}

#endif /* USE_ALTIVEC */

/******************************************************************************\
*
* Neon simd processor
*
\******************************************************************************/

/*
 */
#ifdef USE_NEON

/*
 */
class SimdNeon : public SimdGeneric {
		
	public:
		
		virtual const char *name() const { return "Neon"; }
		
		virtual void FASTCALL dot(float &ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,float v1,int num) const;
		virtual void FASTCALL mul(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL mad(float *ret,const float *v0,float v1,const float *v2,int num) const;
		virtual void FASTCALL add(float *ret,const float *v0,const float *v1,int num) const;
		virtual void FASTCALL sub(float *ret,const float *v0,const float *v1,int num) const;
		
		virtual void FASTCALL minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const;
		virtual void FASTCALL dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL normalizeVec3(void *ret,int ret_stride,int num) const;
		virtual void FASTCALL normalizeVec4(void *ret,int ret_stride,int num) const;
		
		virtual void FASTCALL mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		virtual void FASTCALL projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const;
		
		virtual void FASTCALL mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const;
		
		virtual void FASTCALL skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
		virtual void FASTCALL skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const;
};

/*
 */
void SimdNeon::dot(float &ret,const float *v0,const float *v1,int num) const {
	
	ret = 0.0f;
	
	if(num & ~3) {
		
		float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
		float32x4_t res_1 = vld1q_f32((const float32_t*)v1);
		float32x4_t res_2 = vmulq_f32(res_0,res_1);
		
		v0 += 4;
		v1 += 4;
		
		for(size_t i = (num >> 2) - 1; i > 0; i--) {
			
			res_0 = vld1q_f32((const float32_t*)v0);
			res_1 = vld1q_f32((const float32_t*)v1);
			res_2 = vmlaq_f32(res_2,res_0,res_1);
			
			v0 += 4;
			v1 += 4;
		}
		
		res_0 = vextq_f32(res_2,res_2,2);
		res_2 = vaddq_f32(res_2,res_0);
		res_0 = vextq_f32(res_2,res_2,1);
		res_2 = vaddq_f32(res_2,res_0);
		
		ret = vgetq_lane_f32(res_2,0);
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		ret += *v0++ * *v1++;
	}
}

void SimdNeon::mul(float *ret,const float *v0,float v1,int num) const {
	
	if(num & ~3) {
		
		float32x4_t res_1 = vdupq_n_f32(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
			vst1q_f32((float32_t*)ret,vmulq_f32(res_0,res_1));
			
			ret += 4;
			v0 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1;
	}
}

void SimdNeon::mul(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
			float32x4_t res_1 = vld1q_f32((const float32_t*)v1);
			vst1q_f32((float32_t*)ret,vmulq_f32(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * *v1++;
	}
}

void SimdNeon::mad(float *ret,const float *v0,float v1,const float *v2,int num) const {
	
	if(num & ~3) {
		
		float32x4_t res_1 = vdupq_n_f32(v1);
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
			float32x4_t res_2 = vld1q_f32((const float32_t*)v2);
			vst1q_f32((float32_t*)ret,vmlaq_f32(res_2,res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v2 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ * v1 + *v2++;
	}
}

void SimdNeon::add(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
			float32x4_t res_1 = vld1q_f32((const float32_t*)v1);
			vst1q_f32((float32_t*)ret,vaddq_f32(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ + *v1++;
	}
}

void SimdNeon::sub(float *ret,const float *v0,const float *v1,int num) const {
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)v0);
			float32x4_t res_1 = vld1q_f32((const float32_t*)v1);
			vst1q_f32((float32_t*)ret,vsubq_f32(res_0,res_1));
			
			ret += 4;
			v0 += 4;
			v1 += 4;
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		*ret++ = *v0++ - *v1++;
	}
}

/*
 */
static INLINE void simd_neon_min_max_vec(float32x4_t &min,float32x4_t &max,const void *src,int src_stride,int num) {
	
	const unsigned char *s = (const unsigned char*)src;
	
	if(num & ~3) {
		
		for(size_t i = num >> 2; i > 0; i--) {
			
			float32x4_t res_0 = vld1q_f32((const float32_t*)s);
			s += src_stride;
			float32x4_t res_1 = vld1q_f32((const float32_t*)s);
			s += src_stride;
			float32x4_t res_2 = vld1q_f32((const float32_t*)s);
			s += src_stride;
			float32x4_t res_3 = vld1q_f32((const float32_t*)s);
			s += src_stride;
			
			min = vminq_f32(min,res_0);
			min = vminq_f32(min,res_1);
			min = vminq_f32(min,res_2);
			min = vminq_f32(min,res_3);
			
			max = vmaxq_f32(max,res_0);
			max = vmaxq_f32(max,res_1);
			max = vmaxq_f32(max,res_2);
			max = vmaxq_f32(max,res_3);
		}
		
		num &= 3;
	}
	
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t res = vld1q_f32((const float32_t*)s);
		
		min = vminq_f32(min,res);
		max = vmaxq_f32(max,res);
		
		s += src_stride;
	}
}

void SimdNeon::minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_neon_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

void SimdNeon::minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	simd_neon_min_max_vec(min.vec,max.vec,src,src_stride,num);
}

/*
 */
void SimdNeon::minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	float32x4_t min_vec = min.vec;
	float32x4_t max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmlaq_lane_f32(col_3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		
		min_vec = vminq_f32(min_vec,res_2);
		max_vec = vmaxq_f32(max_vec,res_2);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

void SimdNeon::minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) const {
	
	min.set(INFINITY);
	max.set(-INFINITY);
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	float32x4_t min_vec = min.vec;
	float32x4_t max_vec = max.vec;
	
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmulq_lane_f32(col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		float32x4_t res_3 = vmlaq_lane_f32(res_2,col_3,high,1);
		
		min_vec = vminq_f32(min_vec,res_3);
		max_vec = vmaxq_f32(max_vec,res_3);
		
		s += src_stride;
	}
	
	min.vec = min_vec;
	max.vec = max_vec;
}

/*
 */
void SimdNeon::dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) const {
	
	float32x4_t vec = v.vec;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	for(size_t i = num; i > 0; i--) {
		
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t temp = vld1q_f32((const float32_t*)s);
			float32x4_t res_0 = vmulq_f32(temp,vec);
			float32x2_t res_1 = vget_low_f32(res_0);
			float32x2_t res_2 = vget_high_f32(res_0);
			float32x2_t res_3 = vpadd_f32(res_1,res_1);
			res_3 = vadd_f32(res_3,res_2);
			*(float*)d = vget_lane_f32(res_3,0);
		#else
			asm volatile(
				"vld1.32  { d0, d1 }, [%r0]	\n"
				"vmul.f32   q0, q0, %q2		\n"
				"vpadd.f32  d0, d0,  d0		\n"
				"vadd.f32   d0, d0,  d1		\n"
				"vst1.32  { d0[0] }, [%r1]	\n"
				: : "r"(s), "r"(d), "w"(vec) : "q0"
			);
		#endif
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdNeon::dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) const {
	
	float32x4_t vec = v.vec;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	
	for(size_t i = num; i > 0; i--) {
		
		#if defined(_WINRT) || defined(_IOS)
			float32x4_t temp = vld1q_f32((const float32_t*)s);
			float32x4_t res_0 = vmulq_f32(temp,vec);
			float32x4_t res_1 = vextq_f32(res_0,res_0,2);
			res_0 = vaddq_f32(res_0,res_1);
			res_1 = vextq_f32(res_0,res_0,1);
			res_0 = vaddq_f32(res_0,res_1);
			*(float*)d = vgetq_lane_f32(res_0,0);
		#else
			asm volatile(
				"vld1.32  { d0, d1 }, [%r0]	\n"
				"vmul.f32   q0, q0, %q2		\n"
				"vext.32    q1, q0,  q0, #2	\n"
				"vadd.f32   d0, d0,  d2		\n"
				"vext.32    d2, d0,  d0, #1	\n"
				"vadd.f32   d0, d0,  d2		\n"
				"vst1.32  { d0[0] }, [%r1]	\n"
				: : "r"(s), "r"(d), "w"(vec) : "q0", "q1"
			);
		#endif
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdNeon::normalizeVec3(void *ret,int ret_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)d);
		
		temp = vnormalize3q_f32(temp);
		
		vst1q_f32((float32_t*)d,temp);
		
		d += ret_stride;
	}
}

void SimdNeon::normalizeVec4(void *ret,int ret_stride,int num) const {
	
	unsigned char *d = (unsigned char*)ret;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)d);
		
		temp = vnormalize4q_f32(temp);
		
		vst1q_f32((float32_t*)d,temp);
		
		d += ret_stride;
	}
}

/*
 */
void SimdNeon::mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmulq_lane_f32(col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		
		vst1q_f32((float32_t*)d,res_2);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdNeon::mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmlaq_lane_f32(col_3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		
		vst1q_f32((float32_t*)d,res_2);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdNeon::mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmulq_lane_f32(col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		float32x4_t res_3 = vmlaq_lane_f32(res_2,col_3,high,1);
		
		vst1q_f32((float32_t*)d,res_3);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdNeon::projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmlaq_lane_f32(col_3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		float32x4_t res_3 = vmulq_lane_f32(res_2,vrcp_nr_f32(vget_high_f32(res_2)),1);
		
		vst1q_f32((float32_t*)d,res_3);
		
		s += src_stride;
		d += ret_stride;
	}
}

void SimdNeon::projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) const {
	
	float32x4_t col_0 = m.col0;
	float32x4_t col_1 = m.col1;
	float32x4_t col_2 = m.col2;
	float32x4_t col_3 = m.col3;
	
	unsigned char *d = (unsigned char*)ret;
	const unsigned char *s = (const unsigned char*)src;
	for(size_t i = num; i > 0; i--) {
		
		float32x4_t temp = vld1q_f32((const float32_t*)s);
		
		float32x2_t low = vget_low_f32(temp);
		float32x2_t high = vget_high_f32(temp);
		float32x4_t res_0 = vmulq_lane_f32(col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0);
		float32x4_t res_3 = vmlaq_lane_f32(res_2,col_3,high,1);
		float32x4_t res_4 = vmulq_lane_f32(res_3,vrcp_nr_f32(vget_high_f32(res_3)),1);
		
		vst1q_f32((float32_t*)d,res_4);
		
		s += src_stride;
		d += ret_stride;
	}
}

/*
 */
void SimdNeon::mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) const {
	
	float32x4_t col_0_x = vdupq_n_f32(m.m00);
	float32x4_t col_0_y = vdupq_n_f32(m.m10);
	float32x4_t col_0_z = vdupq_n_f32(m.m20);
	float32x4_t col_1_x = vdupq_n_f32(m.m01);
	float32x4_t col_1_y = vdupq_n_f32(m.m11);
	float32x4_t col_1_z = vdupq_n_f32(m.m21);
	float32x4_t col_2_x = vdupq_n_f32(m.m02);
	float32x4_t col_2_y = vdupq_n_f32(m.m12);
	float32x4_t col_2_z = vdupq_n_f32(m.m22);
	float32x4_t col_3_x = vec4(0.0f,0.0f,0.0f,m.m03).vec;
	float32x4_t col_3_y = vec4(0.0f,0.0f,0.0f,m.m13).vec;
	float32x4_t col_3_z = vec4(0.0f,0.0f,0.0f,m.m23).vec;
	
	for(int i = 0; i < num; i++) {
		
		float32x4x4_t m = vld4q_f32((const float32_t*)matrices[i]->mat);
		
		float32x4_t res_0 = vmlaq_f32(col_3_x,m.val[0],col_0_x);
		float32x4_t res_1 = vmlaq_f32(col_3_y,m.val[0],col_0_y);
		float32x4_t res_2 = vmlaq_f32(col_3_z,m.val[0],col_0_z);
		res_0 = vmlaq_f32(res_0,m.val[1],col_1_x);
		res_1 = vmlaq_f32(res_1,m.val[1],col_1_y);
		res_2 = vmlaq_f32(res_2,m.val[1],col_1_z);
		res_0 = vmlaq_f32(res_0,m.val[2],col_2_x);
		res_1 = vmlaq_f32(res_1,m.val[2],col_2_y);
		res_2 = vmlaq_f32(res_2,m.val[2],col_2_z);
		
		ret[0].vec = res_0;
		ret[1].vec = res_1;
		ret[2].vec = res_2;
		
		ret += 3;
	}
}

/*
 */
void SimdNeon::skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	float weight = weights[0];
	float32x4_t col_0 = vmulq_n_f32(matrices[0]->col0,weight);
	float32x4_t col_1 = vmulq_n_f32(matrices[0]->col1,weight);
	float32x4_t col_2 = vmulq_n_f32(matrices[0]->col2,weight);
	
	for(int i = 1; i < num; i++) {
		float weight = weights[i];
		col_0 = vmlaq_n_f32(col_0,matrices[i]->col0,weight);
		col_1 = vmlaq_n_f32(col_1,matrices[i]->col1,weight);
		col_2 = vmlaq_n_f32(col_2,matrices[i]->col2,weight);
	}
	
	float32x2_t low = vget_low_f32(src.vec);
	float32x2_t high = vget_high_f32(src.vec);
	float32x4_t res_0 = vmulq_lane_f32(col_0,low,0);
	float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
	ret.vec = vmlaq_lane_f32(res_1,col_2,high,0);
}

void SimdNeon::skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) const {
	
	float weight = weights[0];
	float32x4_t col_0 = vmulq_n_f32(matrices[0]->col0,weight);
	float32x4_t col_1 = vmulq_n_f32(matrices[0]->col1,weight);
	float32x4_t col_2 = vmulq_n_f32(matrices[0]->col2,weight);
	float32x4_t col_3 = vmulq_n_f32(matrices[0]->col3,weight);
	
	for(int i = 1; i < num; i++) {
		float weight = weights[i];
		col_0 = vmlaq_n_f32(col_0,matrices[i]->col0,weight);
		col_1 = vmlaq_n_f32(col_1,matrices[i]->col1,weight);
		col_2 = vmlaq_n_f32(col_2,matrices[i]->col2,weight);
		col_3 = vmlaq_n_f32(col_3,matrices[i]->col3,weight);
	}
	
	float32x2_t low = vget_low_f32(src.vec);
	float32x2_t high = vget_high_f32(src.vec);
	float32x4_t res_0 = vmlaq_lane_f32(col_3,col_0,low,0);
	float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
	ret.vec = vmlaq_lane_f32(res_1,col_2,high,0);
}

#endif /* USE_NEON */

/******************************************************************************\
*
* SimdLib
*
\******************************************************************************/

/*
 */
static SimdGeneric simd_generic;
static SimdGeneric *simd = &simd_generic;

/*
 */
Simd::Simd() {
	
}

/*
 */
void Simd::shutdown() {
	if(simd != &simd_generic) delete simd;
	simd = &simd_generic;
}

/*
 */
void Simd::setGeneric() {
	if(simd != &simd_generic) delete simd;
	simd = &simd_generic;
	Log::message("Set %s simd processor\n",simd->name());
}

void Simd::setSSE() {
	if(simd != &simd_generic) delete simd;
	#ifdef USE_SSE
		simd = new SimdSSE();
	#else
		simd = &simd_generic;
		Log::error("Simd::setSSE(): this binary is built without SSE support\n");
	#endif
	Log::message("Set %s simd processor\n",simd->name());
}

void Simd::setSSE2() {
	if(simd != &simd_generic) delete simd;
	#ifdef USE_SSE2
		simd = new SimdSSE2();
	#else
		simd = &simd_generic;
		Log::error("Simd::setSSE2(): this binary is built without SSE2 support\n");
	#endif
	Log::message("Set %s simd processor\n",simd->name());
}

void Simd::setAltiVec() {
	if(simd != &simd_generic) delete simd;
	#ifdef USE_ALTIVEC
		simd = new SimdAltiVec();
	#else
		simd = &simd_generic;
		Log::error("Simd::setAltiVec(): this binary is built without AltiVec support\n");
	#endif
	Log::message("Set %s simd processor\n",simd->name());
}

void Simd::setNeon() {
	if(simd != &simd_generic) delete simd;
	#ifdef USE_NEON
		simd = new SimdNeon();
	#else
		simd = &simd_generic;
		Log::error("Simd::setNeon(): this binary is built without Neon support\n");
	#endif
	Log::message("Set %s simd processor\n",simd->name());
}

/*
 */
void Simd::dot(float &ret,const float *v0,const float *v1,int num) {
	simd->dot(ret,v0,v1,num);
}

void Simd::mul(float *ret,const float *v0,float v1,int num) {
	simd->mul(ret,v0,v1,num);
}

void Simd::mul(float *ret,const float *v0,const float *v1,int num) {
	simd->mul(ret,v0,v1,num);
}

void Simd::mad(float *ret,const float *v0,float v1,const float *v2,int num) {
	simd->mad(ret,v0,v1,v2,num);
}

void Simd::add(float *ret,const float *v0,const float *v1,int num) {
	simd->add(ret,v0,v1,num);
}

void Simd::sub(float *ret,const float *v0,const float *v1,int num) {
	simd->sub(ret,v0,v1,num);
}

/*
 */
void Simd::minMaxVec3(vec3 &min,vec3 &max,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(vec3));
	simd->minMaxVec3(min,max,src,src_stride,num);
}

void Simd::minMaxVec4(vec4 &min,vec4 &max,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(vec4));
	simd->minMaxVec4(min,max,src,src_stride,num);
}

void Simd::minMaxVec3(dvec3 &min,dvec3 &max,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(dvec3));
	simd->minMaxVec3(min,max,src,src_stride,num);
}

void Simd::minMaxVec4(dvec4 &min,dvec4 &max,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(dvec4));
	simd->minMaxVec4(min,max,src,src_stride,num);
}

/*
 */
void Simd::minMaxMat4Vec3(vec3 &min,vec3 &max,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(vec3));
	simd->minMaxMat4Vec3(min,max,m,src,src_stride,num);
}

void Simd::minMaxMat4Vec4(vec4 &min,vec4 &max,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)src_stride >= sizeof(vec4));
	simd->minMaxMat4Vec4(min,max,m,src,src_stride,num);
}

/*
 */
void Simd::dotVec3Vec3(void *ret,int ret_stride,const vec3 &v,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(float));
	assert((size_t)src_stride >= sizeof(vec3));
	simd->dotVec3Vec3(ret,ret_stride,v,src,src_stride,num);
}

void Simd::dotVec4Vec4(void *ret,int ret_stride,const vec4 &v,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(float));
	assert((size_t)src_stride >= sizeof(vec4));
	simd->dotVec4Vec4(ret,ret_stride,v,src,src_stride,num);
}

/*
 */
void Simd::normalizeVec3(void *ret,int ret_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec3));
	simd->normalizeVec3(ret,ret_stride,num);
}

void Simd::normalizeVec4(void *ret,int ret_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec4));
	simd->normalizeVec4(ret,ret_stride,num);
}

/*
 */
void Simd::mulMat3Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec3));
	assert((size_t)src_stride >= sizeof(vec3));
	simd->mulMat3Vec3(ret,ret_stride,m,src,src_stride,num);
}

void Simd::mulMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec3));
	assert((size_t)src_stride >= sizeof(vec3));
	simd->mulMat4Vec3(ret,ret_stride,m,src,src_stride,num);
}

void Simd::mulMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec4));
	assert((size_t)src_stride >= sizeof(vec4));
	simd->mulMat4Vec4(ret,ret_stride,m,src,src_stride,num);
}

void Simd::mulMat3Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(dvec3));
	assert((size_t)src_stride >= sizeof(dvec3));
	simd->mulMat3Vec3(ret,ret_stride,m,src,src_stride,num);
}

void Simd::mulMat4Vec3(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(dvec3));
	assert((size_t)src_stride >= sizeof(dvec3));
	simd->mulMat4Vec3(ret,ret_stride,m,src,src_stride,num);
}

void Simd::mulMat4Vec4(void *ret,int ret_stride,const dmat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(dvec4));
	assert((size_t)src_stride >= sizeof(dvec4));
	simd->mulMat4Vec4(ret,ret_stride,m,src,src_stride,num);
}

/*
 */
void Simd::projMat4Vec3(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec3));
	assert((size_t)src_stride >= sizeof(vec3));
	simd->projMat4Vec3(ret,ret_stride,m,src,src_stride,num);
}

void Simd::projMat4Vec4(void *ret,int ret_stride,const mat4 &m,const void *src,int src_stride,int num) {
	assert((size_t)ret_stride >= sizeof(vec4));
	assert((size_t)src_stride >= sizeof(vec4));
	simd->projMat4Vec4(ret,ret_stride,m,src,src_stride,num);
}

/*
 */
void Simd::mulMat4Mat4(vec4 *ret,const mat4 &m,const mat4 **matrices,int num) {
	simd->mulMat4Mat4(ret,m,matrices,num);
}

void Simd::mulMat4Mat4(vec4 *ret,const dmat4 &m,const dmat4 **matrices,int num) {
	simd->mulMat4Mat4(ret,m,matrices,num);
}

/*
 */
void Simd::skinningMat3(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) {
	simd->skinningMat3(ret,matrices,weights,num,src);
}

void Simd::skinningMat4(vec3 &ret,const mat4 **matrices,const float *weights,int num,const vec3 &src) {
	simd->skinningMat4(ret,matrices,weights,num,src);
}

/*
 */
void Simd::eliminate(float *ret,const float *column,const float *factor,int rows,int num) {
	simd->eliminate(ret,column,factor,rows,num);
}