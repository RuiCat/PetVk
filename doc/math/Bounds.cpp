/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Bounds.cpp
 * Desc:    Bounding volumes
 * Version: 1.46
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
#include "Bounds.h"

/*
 */
#define BOUNDS_EPSILON 1e-4f

/******************************************************************************\
*
* BoundSphere
*
\******************************************************************************/

/*
 */
BoundSphere::BoundSphere() {
	clear();
}

BoundSphere::BoundSphere(const vec3 &center,float radius) {
	set(center,radius);
}

BoundSphere::BoundSphere(const vec3 *points,int num_points,int optimal) {
	set(points,num_points,optimal);
}

BoundSphere::BoundSphere(const BoundSphere &bs) {
	set(bs);
}

BoundSphere::BoundSphere(const BoundSphere &bs,const mat4 &transform) {
	set(bs,transform);
}

BoundSphere::BoundSphere(const BoundBox &bb) {
	set(bb);
}

BoundSphere::~BoundSphere() {
	
}

/*
 */
void BoundSphere::clear() {
	center = vec3_zero;
	center.w = -1.0f;
}

/*
 */
void BoundSphere::set(const vec3 &c,float r) {
	center = c;
	center.w = r;
}

void BoundSphere::set(const vec3 *points,int num_points,int optimal) {
	clear();
	if(num_points > 2 && optimal) {
		center.w = INFINITY;
		for(int i = 0; i < num_points; i++) {
			for(int j = i + 1; j < num_points; j++) {
				float radius2 = -INFINITY;
				vec3 point = (points[i] + points[j]) * 0.5f;
				for(int k = 0; k < num_points; k++) {
					float length2 = (points[k] - point).length2();
					if(radius2 < length2) radius2 = length2;
				}
				if(center.w > radius2) {
					center = point;
					center.w = radius2;
				}
			}
		}
		center.w = (center.w > 0.0f) ? Math::sqrt(center.w) : -1.0f;
	} else {
		expand(points,num_points);
	}
}

void BoundSphere::set(const BoundSphere &bs) {
	center = bs.center;
	center.w = bs.center.w;
}

void BoundSphere::set(const BoundSphere &bs,const mat4 &transform) {
	center = bs.center;
	center.w = bs.center.w;
	setTransform(transform);
}

void BoundSphere::set(const BoundBox &bb) {
	clear();
	expand(bb);
}

/*
 */
void BoundSphere::setTransform(const mat4 &transform) {
	float radius = center.w;
	#ifdef USE_SSE
		__m128 col_0 = transform.col0;
		__m128 col_1 = transform.col1;
		__m128 col_2 = transform.col2;
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(center.vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(center.vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(center.vec,Z,Z,Z,W));
		center.vec = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,transform.col3));
		col_0 = _mm_mul_ps(col_0,col_0);
		col_1 = _mm_mul_ps(col_1,col_1);
		col_2 = _mm_mul_ps(col_2,col_2);
		col_0 = _mm_add_ps(col_0,_MM_SWIZZLE(col_0,Y,X,W,Z));
		col_1 = _mm_add_ps(col_1,_MM_SWIZZLE(col_1,Y,X,W,Z));
		col_2 = _mm_add_ps(col_2,_MM_SWIZZLE(col_2,Y,X,W,Z));
		col_0 = _mm_add_ss(col_0,_MM_SWIZZLE(col_0,Z,W,X,Y));
		col_1 = _mm_add_ss(col_1,_MM_SWIZZLE(col_1,Z,W,X,Y));
		col_2 = _mm_add_ss(col_2,_MM_SWIZZLE(col_2,Z,W,X,Y));
		col_0 = _mm_max_ss(_mm_max_ss(col_0,col_1),col_2);
		col_0 = _mm_mul_ss(col_0,_mm_rsqrt_ss(col_0));
		center.w = radius * _mm_cvtss_f32(col_0);
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 col_0 = transform.col0;
		vec_float4 col_1 = transform.col1;
		vec_float4 col_2 = transform.col2;
		vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(center.vec,X,X,X,W),transform.col3);
		vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(center.vec,Y,Y,Y,W),res_0);
		center.vec = vec_madd(col_2,VEC_SWIZZLE(center.vec,Z,Z,Z,W),res_1);
		col_0 = vec_madd(col_0,col_0,zero);
		col_1 = vec_madd(col_1,col_1,zero);
		col_2 = vec_madd(col_2,col_2,zero);
		col_0 = vec_add(col_0,vec_sld(col_0,col_0,8));
		col_1 = vec_add(col_1,vec_sld(col_1,col_1,8));
		col_2 = vec_add(col_2,vec_sld(col_2,col_2,8));
		col_0 = vec_add(col_0,vec_sld(col_0,col_0,4));
		col_1 = vec_add(col_1,vec_sld(col_1,col_1,4));
		col_2 = vec_add(col_2,vec_sld(col_2,col_2,4));
		col_0 = vec_max(vec_max(col_0,col_1),col_2);
		col_0 = vec_madd(col_0,vec_rsqrte(col_0),zero);
		center.w = radius * vec_extract(col_0,0);
	#elif USE_NEON
		float32x4_t col_0 = transform.col0;
		float32x4_t col_1 = transform.col1;
		float32x4_t col_2 = transform.col2;
		float32x2_t low = vget_low_f32(center.vec);
		float32x2_t high = vget_high_f32(center.vec);
		float32x4_t res_0 = vmlaq_lane_f32(transform.col3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		center.vec = vmlaq_lane_f32(res_1,col_2,high,0);
		col_0 = vdot33q_f32(col_0,col_0);
		col_1 = vdot33q_f32(col_1,col_1);
		col_2 = vdot33q_f32(col_2,col_2);
		col_0 = vmaxq_f32(col_0,vmaxq_f32(col_1,col_2));
		col_0 = vmulq_f32(col_0,vrsqrteq_f32(col_0));
		center.w = radius * vgetq_lane_f32(col_0,0);
	#else
		center = transform * center;
		float x = transform.m00 * transform.m00 + transform.m10 * transform.m10 + transform.m20 * transform.m20;
		float y = transform.m01 * transform.m01 + transform.m11 * transform.m11 + transform.m21 * transform.m21;
		float z = transform.m02 * transform.m02 + transform.m12 * transform.m12 + transform.m22 * transform.m22;
		float scale = Math::sqrtFast(::max(::max(x,y),z));
		center.w = radius * scale;
	#endif
}

void BoundSphere::setTransform(const dmat4 &transform) {
	setTransform(mat4(transform));
}

/*
 */
int BoundSphere::compare(const BoundSphere &bs) const {
	return (center == bs.center && ::compare(center.w,bs.center.w));
}

/*
 */
void BoundSphere::expand(const vec3 &point) {
	if(isValid()) {
		vec3 direction = point - center;
		float length = direction.length();
		if(length > center.w) {
			float delta = (length - center.w) * 0.5f;
			center += direction * (delta / length);
			center.w += delta;
		}
	} else {
		center = point;
		center.w = BOUNDS_EPSILON;
	}
}

void BoundSphere::expand(const vec3 *points,int num_points) {
	if(isValid()) {
		for(int i = 0; i < num_points; i++) {
			vec3 direction = points[i] - center;
			float length = direction.length();
			if(length > center.w) {
				float delta = (length - center.w) * 0.5f;
				center += direction * (delta / length);
				center.w += delta;
			}
		}
	} else {
		vec3 min,max;
		Simd::minMaxVec3(min,max,points,sizeof(vec3),num_points);
		center = (min + max) * 0.5f;
		float radius2 = -INFINITY;
		for(int i = 0; i < num_points; i++) {
			float length2 = (points[i] - center).length2();
			if(radius2 < length2) radius2 = length2;
		}
		center.w = (radius2 > 0.0f) ? Math::sqrt(radius2) : -1.0f;
	}
}

void BoundSphere::expand(const BoundSphere &bs) {
	if(bs.isValid()) {
		if(isValid()) {
			vec3 direction = bs.center - center;
			float length = direction.length();
			if(length > EPSILON) {
				if(length + center.w < bs.center.w) {
					center = bs.center;
					center.w = bs.center.w;
				} else if(length + bs.center.w > center.w) {
					vec3 p0 = center - direction * (center.w / length);
					vec3 p1 = bs.center + direction * (bs.center.w / length);
					center = (p0 + p1) * 0.5f;
					center.w = (p1 - center).length();
				}
			} else {
				if(center.w < bs.center.w) center.w = bs.center.w;
			}
		} else {
			center = bs.center;
			center.w = bs.center.w;
		}
	}
}

void BoundSphere::expand(const BoundBox &bb) {
	if(bb.isValid()) {
		const vec3 &min = bb.getMin();
		const vec3 &max = bb.getMax();
		if(isValid()) {
			expand(vec3(min.x,min.y,min.z));
			expand(vec3(max.x,min.y,min.z));
			expand(vec3(min.x,max.y,min.z));
			expand(vec3(max.x,max.y,min.z));
			expand(vec3(min.x,min.y,max.z));
			expand(vec3(max.x,min.y,max.z));
			expand(vec3(min.x,max.y,max.z));
			expand(vec3(max.x,max.y,max.z));
		} else {
			center = (min + max) * 0.5f;
			center.w = length(max - center);
		}
	}
}

/*
 */
void BoundSphere::expandRadius(const vec3 &point) {
	if(isValid()) {
		float radius = length(center - point);
		if(center.w < radius) center.w = radius;
	} else {
		center = point;
		center.w = BOUNDS_EPSILON;
	}
}

void BoundSphere::expandRadius(const vec3 *points,int num_points) {
	if(isValid()) {
		for(int i = 0; i < num_points; i++) {
			float radius = length(center - points[i]);
			if(center.w < radius) center.w = radius;
		}
	} else {
		vec3 min,max;
		Simd::minMaxVec3(min,max,points,sizeof(vec3),num_points);
		center = (min + max) * 0.5f;
		float radius2 = -INFINITY;
		for(int i = 0; i < num_points; i++) {
			float length2 = (points[i] - center).length2();
			if(radius2 < length2) radius2 = length2;
		}
		center.w = (radius2 > 0.0f) ? Math::sqrt(radius2) : -1.0f;
	}
}

void BoundSphere::expandRadius(const BoundSphere &bs) {
	if(bs.isValid()) {
		if(isValid()) {
			float radius = length(bs.center - center) + bs.center.w;
			if(center.w < radius) center.w = radius;
		} else {
			center = bs.center;
			center.w = bs.center.w;
		}
	}
}

void BoundSphere::expandRadius(const BoundBox &bb) {
	if(bb.isValid()) {
		const vec3 &min = bb.getMin();
		const vec3 &max = bb.getMax();
		if(isValid()) {
			expandRadius(vec3(min.x,min.y,min.z));
			expandRadius(vec3(max.x,min.y,min.z));
			expandRadius(vec3(min.x,max.y,min.z));
			expandRadius(vec3(max.x,max.y,min.z));
			expandRadius(vec3(min.x,min.y,max.z));
			expandRadius(vec3(max.x,min.y,max.z));
			expandRadius(vec3(min.x,max.y,max.z));
			expandRadius(vec3(max.x,max.y,max.z));
		} else {
			center = (min + max) * 0.5f;
			center.w = length(max - center);
		}
	}
}

/*
 */
int BoundSphere::inside(const vec3 &point) const {
	if(isValid()) return insideValid(point);
	return 0;
}

int BoundSphere::inside(const vec3 &point,float radius) const {
	if(isValid()) return insideValid(point,radius);
	return 0;
}

int BoundSphere::inside(const vec3 &min,const vec3 &max) const {
	if(isValid()) return insideValid(min,max);
	return 0;
}

/*
 */
int BoundSphere::inside(const BoundSphere &bs) const {
	if(isValid() && bs.isValid()) return insideValid(bs.center,bs.center.w);
	return 0;
}

int BoundSphere::inside(const BoundBox &bb) const {
	if(isValid() && bb.isValid()) return insideValid(bb.getMin(),bb.getMax());
	return 0;
}

/*
 */
int BoundSphere::insideAll(const BoundSphere &bs) const {
	if(isValid() && bs.isValid()) return insideAllValid(bs);
	return 0;
}

int BoundSphere::insideAll(const BoundBox &bb) const {
	if(isValid() && bb.isValid()) return insideAllValid(bb);
	return 0;
}

/*
 */
int BoundSphere::rayIntersection(const vec3 &p,const vec3 &direction) const {
	if(isValid()) return rayIntersectionValid(p,direction);
	return 0;
}

int BoundSphere::getIntersection(const vec3 &p0,const vec3 &p1) const {
	if(isValid()) return getIntersectionValid(p0,p1);
	return 0;
}

/*
 */
float BoundSphere::distance() const {
	if(isValid()) return distanceValid();
	return INFINITY;
}

float BoundSphere::distance(const vec3 &p) const {
	if(isValid()) return distanceValid(p);
	return INFINITY;
}

/*
 */
BoundSphere operator*(const mat4 &m,const BoundSphere &bs) {
	BoundSphere ret = bs;
	ret.setTransform(m);
	return ret;
}

BoundSphere operator*(const dmat4 &m,const BoundSphere &bs) {
	BoundSphere ret = bs;
	ret.setTransform(m);
	return ret;
}

/******************************************************************************\
*
* BoundBox
*
\******************************************************************************/

/*
 */
BoundBox::BoundBox() {
	clear();
}

BoundBox::BoundBox(const vec3 &min,const vec3 &max) {
	set(min,max);
}

BoundBox::BoundBox(const vec3 *points,int num_points) {
	set(points,num_points);
}

BoundBox::BoundBox(const BoundBox &bb) {
	set(bb);
}

BoundBox::BoundBox(const BoundBox &bb,const mat4 &transform) {
	set(bb,transform);
}

BoundBox::BoundBox(const BoundSphere &bs) {
	set(bs);
}

BoundBox::~BoundBox() {
	
}

/*
 */
void BoundBox::clear() {
	min = vec3_infinity;
	max = -vec3_infinity;
}

/*
 */
void BoundBox::set(const vec3 &min_,const vec3 &max_) {
	min = min_;
	max = max_;
}

void BoundBox::set(const vec3 *points,int num_points) {
	clear();
	expand(points,num_points);
}

void BoundBox::set(const BoundSphere &bs) {
	clear();
	expand(bs);
}

void BoundBox::set(const BoundBox &bb) {
	min = bb.min;
	max = bb.max;
}

void BoundBox::set(const BoundBox &bb,const mat4 &transform) {
	min = bb.min;
	max = bb.max;
	setTransform(transform);
}

/*
 */
void BoundBox::setTransform(const mat4 &transform) {
	#ifdef USE_SSE
		__m128 sign = _mm_set1_ps(IntFloat(0x7fffffff).f);
		__m128 center = _mm_mul_ps(_mm_add_ps(min.vec,max.vec),_mm_set1_ps(0.5f));
		__m128 axis = _mm_sub_ps(max.vec,center);
		__m128 col_0 = transform.col0;
		__m128 col_1 = transform.col1;
		__m128 col_2 = transform.col2;
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(center,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(center,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(center,Z,Z,Z,W));
		center = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,transform.col3));
		col_0 = _mm_mul_ps(_mm_and_ps(col_0,sign),_MM_SWIZZLE(axis,X,X,X,W));
		col_1 = _mm_mul_ps(_mm_and_ps(col_1,sign),_MM_SWIZZLE(axis,Y,Y,Y,W));
		col_2 = _mm_mul_ps(_mm_and_ps(col_2,sign),_MM_SWIZZLE(axis,Z,Z,Z,W));
		axis = _mm_add_ps(_mm_add_ps(col_0,col_1),col_2);
		min.vec = _mm_sub_ps(center,axis);
		max.vec = _mm_add_ps(center,axis);
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 center = vec_madd(vec_add(min.vec,max.vec),vec_splats(0.5f),zero);
		vec_float4 axis = vec_sub(max.vec,center);
		vec_float4 col_0 = transform.col0;
		vec_float4 col_1 = transform.col1;
		vec_float4 col_2 = transform.col2;
		vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(center,X,X,X,W),transform.col3);
		vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(center,Y,Y,Y,W),res_0);
		center = vec_madd(col_2,VEC_SWIZZLE(center,Z,Z,Z,W),res_1);
		col_0 = vec_madd(vec_abs(col_0),VEC_SWIZZLE(axis,X,X,X,W),zero);
		col_1 = vec_madd(vec_abs(col_1),VEC_SWIZZLE(axis,Y,Y,Y,W),col_0);
		axis = vec_madd(vec_abs(col_2),VEC_SWIZZLE(axis,Z,Z,Z,W),col_1);
		min.vec = vec_sub(center,axis);
		max.vec = vec_add(center,axis);
	#elif USE_NEON
		float32x4_t center = vmulq_n_f32(vaddq_f32(min.vec,max.vec),0.5f);
		float32x4_t axis = vsubq_f32(max.vec,center);
		float32x4_t col_0 = transform.col0;
		float32x4_t col_1 = transform.col1;
		float32x4_t col_2 = transform.col2;
		float32x2_t low = vget_low_f32(center);
		float32x2_t high = vget_high_f32(center);
		float32x4_t res_0 = vmlaq_lane_f32(transform.col3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		center = vmlaq_lane_f32(res_1,col_2,high,0);
		low = vget_low_f32(axis);
		high = vget_high_f32(axis);
		col_0 = vmulq_lane_f32(vabsq_f32(col_0),low,0);
		col_1 = vmlaq_lane_f32(col_0,vabsq_f32(col_1),low,1);
		axis = vmlaq_lane_f32(col_1,vabsq_f32(col_2),high,0);
		min.vec = vsubq_f32(center,axis);
		max.vec = vaddq_f32(center,axis);
	#else
		vec3 center = (min + max) * 0.5f;
		vec3 axis = max - center;
		float x = Math::abs(transform.m00) * axis.x + Math::abs(transform.m01) * axis.y + Math::abs(transform.m02) * axis.z;
		float y = Math::abs(transform.m10) * axis.x + Math::abs(transform.m11) * axis.y + Math::abs(transform.m12) * axis.z;
		float z = Math::abs(transform.m20) * axis.x + Math::abs(transform.m21) * axis.y + Math::abs(transform.m22) * axis.z;
		center = transform * center;
		min = center - vec3(x,y,z);
		max = center + vec3(x,y,z);
	#endif
}

void BoundBox::setTransform(const dmat4 &transform) {
	setTransform(mat4(transform));
}

void BoundBox::setTransform(const BoundSphere &bs,const mat4 &transform) {
	setTransform(transform);
	#ifdef USE_SSE
		__m128 col_0 = transform.col0;
		__m128 col_1 = transform.col1;
		__m128 col_2 = transform.col2;
		__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(bs.getCenter().vec,X,X,X,W));
		__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(bs.getCenter().vec,Y,Y,Y,W));
		__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(bs.getCenter().vec,Z,Z,Z,W));
		__m128 center = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,transform.col3));
		col_0 = _mm_mul_ps(col_0,col_0);
		col_1 = _mm_mul_ps(col_1,col_1);
		col_2 = _mm_mul_ps(col_2,col_2);
		col_0 = _mm_add_ps(col_0,_MM_SWIZZLE(col_0,Y,X,W,Z));
		col_1 = _mm_add_ps(col_1,_MM_SWIZZLE(col_1,Y,X,W,Z));
		col_2 = _mm_add_ps(col_2,_MM_SWIZZLE(col_2,Y,X,W,Z));
		col_0 = _mm_add_ps(col_0,_MM_SWIZZLE(col_0,Z,W,X,Y));
		col_1 = _mm_add_ps(col_1,_MM_SWIZZLE(col_1,Z,W,X,Y));
		col_2 = _mm_add_ps(col_2,_MM_SWIZZLE(col_2,Z,W,X,Y));
		col_0 = _mm_max_ps(_mm_max_ps(col_0,col_1),col_2);
		col_0 = _mm_mul_ps(col_0,_mm_rsqrt_ps(col_0));
		__m128 radius = _mm_mul_ps(col_0,vec3(bs.getRadius()).vec);
		min.vec = _mm_max_ps(min.vec,_mm_sub_ps(center,radius));
		max.vec = _mm_min_ps(max.vec,_mm_add_ps(center,radius));
	#elif USE_ALTIVEC
		vec_float4 zero = vec_splats(0.0f);
		vec_float4 col_0 = transform.col0;
		vec_float4 col_1 = transform.col1;
		vec_float4 col_2 = transform.col2;
		vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(bs.getCenter().vec,X,X,X,W),transform.col3);
		vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(bs.getCenter().vec,Y,Y,Y,W),res_0);
		vec_float4 center = vec_madd(col_2,VEC_SWIZZLE(bs.getCenter().vec,Z,Z,Z,W),res_1);
		col_0 = vec_madd(col_0,col_0,zero);
		col_1 = vec_madd(col_1,col_1,zero);
		col_2 = vec_madd(col_2,col_2,zero);
		col_0 = vec_add(col_0,vec_sld(col_0,col_0,8));
		col_1 = vec_add(col_1,vec_sld(col_1,col_1,8));
		col_2 = vec_add(col_2,vec_sld(col_2,col_2,8));
		col_0 = vec_add(col_0,vec_sld(col_0,col_0,4));
		col_1 = vec_add(col_1,vec_sld(col_1,col_1,4));
		col_2 = vec_add(col_2,vec_sld(col_2,col_2,4));
		col_0 = vec_max(vec_max(col_0,col_1),col_2);
		col_0 = vec_madd(col_0,vec_rsqrte(col_0),zero);
		vec_float4 radius = vec_madd(col_0,vec_splats(bs.getRadius()),zero);
		min.vec = vec_max(min.vec,vec_sub(center,radius));
		max.vec = vec_min(max.vec,vec_add(center,radius));
	#elif USE_NEON
		float32x4_t col_0 = transform.col0;
		float32x4_t col_1 = transform.col1;
		float32x4_t col_2 = transform.col2;
		float32x2_t low = vget_low_f32(bs.getCenter().vec);
		float32x2_t high = vget_high_f32(bs.getCenter().vec);
		float32x4_t res_0 = vmlaq_lane_f32(transform.col3,col_0,low,0);
		float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1);
		float32x4_t center = vmlaq_lane_f32(res_1,col_2,high,0);
		col_0 = vdot33q_f32(col_0,col_0);
		col_1 = vdot33q_f32(col_1,col_1);
		col_2 = vdot33q_f32(col_2,col_2);
		col_0 = vmaxq_f32(col_0,vmaxq_f32(col_1,col_2));
		col_0 = vmulq_f32(col_0,vrsqrteq_f32(col_0));
		float32x4_t radius = vmulq_n_f32(col_0,bs.getRadius());
		min.vec = vmaxq_f32(min.vec,vsubq_f32(center,radius));
		max.vec = vminq_f32(max.vec,vaddq_f32(center,radius));
	#else
		vec3 center = transform * bs.getCenter();
		float x = transform.m00 * transform.m00 + transform.m10 * transform.m10 + transform.m20 * transform.m20;
		float y = transform.m01 * transform.m01 + transform.m11 * transform.m11 + transform.m21 * transform.m21;
		float z = transform.m02 * transform.m02 + transform.m12 * transform.m12 + transform.m22 * transform.m22;
		float radius = Math::sqrtFast(::max(::max(x,y),z)) * bs.getRadius();
		min = ::max(min,center - vec3(radius));
		max = ::min(max,center + vec3(radius));
	#endif
}

void BoundBox::setTransform(const BoundSphere &bs,const dmat4 &transform) {
	setTransform(bs,mat4(transform));
}

/*
 */
int BoundBox::compare(const BoundBox &bb) const {
	return (min == bb.min && max == bb.max);
}

/*
 */
void BoundBox::expand(const vec3 &point) {
	if(isValid()) {
		#ifdef USE_SSE
			min.vec = _mm_min_ps(min.vec,point.vec);
			max.vec = _mm_max_ps(max.vec,point.vec);
		#elif USE_ALTIVEC
			min.vec = vec_min(min.vec,point.vec);
			max.vec = vec_max(max.vec,point.vec);
		#elif USE_NEON
			min.vec = vminq_f32(min.vec,point.vec);
			max.vec = vmaxq_f32(max.vec,point.vec);
		#else
			min = ::min(min,point);
			max = ::max(max,point);
		#endif
	} else {
		min = point - vec3(BOUNDS_EPSILON);
		max = point + vec3(BOUNDS_EPSILON);
	}
}

void BoundBox::expand(const vec3 *points,int num_points) {
	if(isValid()) {
		vec3 min_,max_;
		Simd::minMaxVec3(min_,max_,points,sizeof(vec3),num_points);
		#ifdef USE_SSE
			min.vec = _mm_min_ps(min.vec,min_.vec);
			max.vec = _mm_max_ps(max.vec,max_.vec);
		#elif USE_ALTIVEC
			min.vec = vec_min(min.vec,min_.vec);
			max.vec = vec_max(max.vec,max_.vec);
		#elif USE_NEON
			min.vec = vminq_f32(min.vec,min_.vec);
			max.vec = vmaxq_f32(max.vec,max_.vec);
		#else
			min = ::min(min,min_);
			max = ::max(max,max_);
		#endif
	} else {
		Simd::minMaxVec3(min,max,points,sizeof(vec3),num_points);
	}
}

void BoundBox::expand(const BoundSphere &bs) {
	if(bs.isValid()) {
		const vec3 &center = bs.getCenter();
		float radius = bs.getRadius();
		if(isValid()) {
			#ifdef USE_SSE
				__m128 r = vec3(radius).vec;
				min.vec = _mm_min_ps(min.vec,_mm_sub_ps(center.vec,r));
				max.vec = _mm_max_ps(max.vec,_mm_add_ps(center.vec,r));
			#elif USE_ALTIVEC
				vec_float4 r = vec_splats(radius);
				min.vec = vec_min(min.vec,vec_sub(center.vec,r));
				max.vec = vec_max(max.vec,vec_add(center.vec,r));
			#elif USE_NEON
				float32x4_t r = vdupq_n_f32(radius);
				min.vec = vminq_f32(min.vec,vsubq_f32(center.vec,r));
				max.vec = vmaxq_f32(max.vec,vaddq_f32(center.vec,r));
			#else
				min = ::min(min,center - vec3(radius));
				max = ::max(max,center + vec3(radius));
			#endif
		} else {
			min = center - vec3(radius);
			max = center + vec3(radius);
		}
	}
}

void BoundBox::expand(const BoundBox &bb) {
	if(bb.isValid()) {
		if(isValid()) {
			#ifdef USE_SSE
				min.vec = _mm_min_ps(min.vec,bb.min.vec);
				max.vec = _mm_max_ps(max.vec,bb.max.vec);
			#elif USE_ALTIVEC
				min.vec = vec_min(min.vec,bb.min.vec);
				max.vec = vec_max(max.vec,bb.max.vec);
			#elif USE_NEON
				min.vec = vminq_f32(min.vec,bb.min.vec);
				max.vec = vmaxq_f32(max.vec,bb.max.vec);
			#else
				min = ::min(min,bb.min);
				max = ::max(max,bb.max);
			#endif
		} else {
			min = bb.min;
			max = bb.max;
		}
	}
}

/*
 */
int BoundBox::inside(const vec3 &point) const {
	if(isValid()) return insideValid(point);
	return 0;
}

int BoundBox::inside(const vec3 &point,float radius) const {
	if(isValid()) return insideValid(point,radius);
	return 0;
}

int BoundBox::inside(const vec3 &min,const vec3 &max) const {
	if(isValid()) return insideValid(min,max);
	return 0;
}

/*
 */
int BoundBox::inside(const BoundSphere &bs) const {
	if(isValid() && bs.isValid()) return insideValid(bs.getCenter(),bs.getRadius());
	return 0;
}

int BoundBox::inside(const BoundBox &bb) const {
	if(isValid() && bb.isValid()) return insideValid(bb.min,bb.max);
	return 0;
}

/*
 */
int BoundBox::insideAll(const BoundSphere &bs) const {
	if(isValid() && bs.isValid()) return insideAllValid(bs);
	return 0;
}

int BoundBox::insideAll(const BoundBox &bb) const {
	if(isValid() && bb.isValid()) return insideAllValid(bb);
	return 0;
}

/*
 */
int BoundBox::insideCube(int face,const vec3 &offset) const {
	if(isValid()) {
		vec3 min = getMin() - offset;
		vec3 max = getMax() - offset;
		switch(face) {
			case 0: return (max.x >= 0.0f && min.y <=  max.x && max.y >= -max.x && min.z <=  max.x && max.z >= -max.x);
			case 1: return (min.x <= 0.0f && min.y <= -min.x && max.y >=  min.x && min.z <= -min.x && max.z >=  min.x);
			case 2: return (max.y >= 0.0f && min.x <=  max.y && max.x >= -max.y && min.z <=  max.y && max.z >= -max.y);
			case 3: return (min.y <= 0.0f && min.x <= -min.y && max.x >=  min.y && min.z <= -min.y && max.z >=  min.y);
			case 4: return (max.z >= 0.0f && min.x <=  max.z && max.x >= -max.z && min.y <=  max.z && max.y >= -max.z);
			case 5: return (min.z <= 0.0f && min.x <= -min.z && max.x >=  min.z && min.y <= -min.z && max.y >=  min.z);
		};
		assert(0 && "BoundBox::insideCube(): bad face number");
	}
	return 0;
}

/*
 */
int BoundBox::rayIntersection(const vec3 &p,const vec3 &direction) const {
	if(isValid()) return rayIntersectionValid(p,direction);
	return 0;
}

int BoundBox::irayIntersection(const vec3 &p,const vec3 &idirection) const {
	if(isValid()) return irayIntersectionValid(p,idirection);
	return 0;
}

int BoundBox::getIntersection(const vec3 &p0,const vec3 &p1) const {
	if(isValid()) return getIntersectionValid(p0,p1);
	return 0;
}

/*
 */
float BoundBox::distance() const {
	if(isValid()) return distanceValid();
	return INFINITY;
}

float BoundBox::distance(const vec3 &p) const {
	if(isValid()) return distanceValid(p);
	return INFINITY;
}

/*
 */
void BoundBox::getPoints(vec3 *points,int num_points) const {
	assert(num_points == 8 && "BoundBox::getPoints(): bad points number");
	points[0].set(min.x,min.y,min.z);
	points[1].set(max.x,min.y,min.z);
	points[2].set(min.x,max.y,min.z);
	points[3].set(max.x,max.y,min.z);
	points[4].set(min.x,min.y,max.z);
	points[5].set(max.x,min.y,max.z);
	points[6].set(min.x,max.y,max.z);
	points[7].set(max.x,max.y,max.z);
}

void BoundBox::getPlanes(vec4 *planes,int num_planes) const {
	assert(num_planes == 6 && "BoundBox::getPlanes(): bad planes number");
	planes[0].set( 1.0f, 0.0f, 0.0f,-max.x);
	planes[1].set(-1.0f, 0.0f, 0.0f, min.x);
	planes[2].set( 0.0f, 1.0f, 0.0f,-max.y);
	planes[3].set( 0.0f,-1.0f, 0.0f, min.y);
	planes[4].set( 0.0f, 0.0f, 1.0f,-max.z);
	planes[5].set( 0.0f, 0.0f,-1.0f, min.z);
}

/*
 */
BoundBox operator*(const mat4 &m,const BoundBox &bb) {
	BoundBox ret = bb;
	ret.setTransform(m);
	return ret;
}

BoundBox operator*(const dmat4 &m,const BoundBox &bb) {
	BoundBox ret = bb;
	ret.setTransform(m);
	return ret;
}

/******************************************************************************\
*
* BoundFrustum
*
\******************************************************************************/

/*
 */
BoundFrustum::BoundFrustum() : valid(0) {
	
}

BoundFrustum::BoundFrustum(const mat4 &projection,const mat4 &modelview) {
	set(projection,modelview);
}

BoundFrustum::BoundFrustum(const BoundFrustum &bf) {
	set(bf);
}

BoundFrustum::BoundFrustum(const BoundFrustum &bf,const mat4 &itransform) {
	set(bf,itransform);
}

BoundFrustum::~BoundFrustum() {
	
}

/*
 */
void BoundFrustum::clear() {
	
	valid = 0;
	
	// camera
	camera = vec3_zero;
	
	// clipping planes
	for(int i = 0; i < 6; i++) {
		planes[i] = vec4_zero;
	}
	
	// clipping planes and points
	for(int i = 0; i < 8; i++) {
		tplanes[i] = vec4_zero;
		points[i] = vec3_zero;
	}
	
	// portals
	portals.clear();
}

/*
 */
void BoundFrustum::set(const mat4 &projection,const mat4 &modelview) {
	
	valid = 1;
	
	// camera
	camera = inverse(modelview).getColumn3(3);
	
	// modelview projection matrix
	mat4 mvp = projection * modelview;
	
	// points
	points[0].set(-1.0f,-1.0f,-1.0f);
	points[1].set( 1.0f,-1.0f,-1.0f);
	points[2].set(-1.0f, 1.0f,-1.0f);
	points[3].set( 1.0f, 1.0f,-1.0f);
	points[4].set(-1.0f,-1.0f, 1.0f);
	points[5].set( 1.0f,-1.0f, 1.0f);
	points[6].set(-1.0f, 1.0f, 1.0f);
	points[7].set( 1.0f, 1.0f, 1.0f);
	Simd::projMat4Vec3(points,sizeof(vec3),inverse(mvp),points,sizeof(vec3),8);
	
	// clipping planes
	planes[0].set(mvp.m30 + mvp.m00,mvp.m31 + mvp.m01,mvp.m32 + mvp.m02,mvp.m33 + mvp.m03);
	planes[1].set(mvp.m30 - mvp.m00,mvp.m31 - mvp.m01,mvp.m32 - mvp.m02,mvp.m33 - mvp.m03);
	planes[2].set(mvp.m30 + mvp.m10,mvp.m31 + mvp.m11,mvp.m32 + mvp.m12,mvp.m33 + mvp.m13);
	planes[3].set(mvp.m30 - mvp.m10,mvp.m31 - mvp.m11,mvp.m32 - mvp.m12,mvp.m33 - mvp.m13);
	planes[4].set(mvp.m30 + mvp.m20,mvp.m31 + mvp.m21,mvp.m32 + mvp.m22,mvp.m33 + mvp.m23);
	planes[5].set(mvp.m30 - mvp.m20,mvp.m31 - mvp.m21,mvp.m32 - mvp.m22,mvp.m33 - mvp.m23);
	for(int i = 0; i < 6; i++) {
		planes[i] /= length(vec3(planes[i]));
	}
	for(int i = 0; i < 4; i++) {
		tplanes[0][i] = planes[i].x;
		tplanes[1][i] = planes[i].y;
		tplanes[2][i] = planes[i].z;
		tplanes[3][i] = planes[i].w;
	}
	for(int i = 0, j = 4; i < 2; i++, j++) {
		tplanes[4][i] = planes[j].x;
		tplanes[5][i] = planes[j].y;
		tplanes[6][i] = planes[j].z;
		tplanes[7][i] = planes[j].w;
	}
	
	// portals
	portals.clear();
}

void BoundFrustum::set(const BoundFrustum &bf) {
	
	valid = bf.valid;
	
	// camera
	camera = bf.camera;
	
	// clipping planes and points
	for(int i = 0; i < 6; i++) {
		planes[i] = bf.planes[i];
	}
	for(int i = 0; i < 8; i++) {
		tplanes[i] = bf.tplanes[i];
		points[i] = bf.points[i];
	}
	
	// portals
	portals = bf.portals;
}

void BoundFrustum::set(const BoundFrustum &bf,const mat4 &itransform) {
	set(bf);
	setITransform(itransform);
}

/*
 */
void BoundFrustum::setITransform(const mat4 &itransform) {
	
	mat4 transform = inverse(itransform);
	
	// camera
	camera = transform * camera;
	
	// points
	Simd::mulMat4Vec3(points,sizeof(vec3),transform,points,sizeof(vec3),8);
	
	// clipping planes
	for(int i = 0; i < 6; i++) {
		planes[i] = planes[i] * itransform;
		planes[i] /= length(vec3(planes[i]));
	}
	for(int i = 0; i < 4; i++) {
		tplanes[0][i] = planes[i].x;
		tplanes[1][i] = planes[i].y;
		tplanes[2][i] = planes[i].z;
		tplanes[3][i] = planes[i].w;
	}
	for(int i = 0, j = 4; i < 2; i++, j++) {
		tplanes[4][i] = planes[j].x;
		tplanes[5][i] = planes[j].y;
		tplanes[6][i] = planes[j].z;
		tplanes[7][i] = planes[j].w;
	}
	
	// portals
	for(int i = 0; i < portals.size(); i++) {
		Portal &portal = portals[i];
		portal.plane = portal.plane * itransform;
		portal.plane /= length(vec3(portal.plane));
		Simd::mulMat4Vec3(portal.points,sizeof(vec3),transform,portal.points,sizeof(vec3),4);
		for(int j = 0; j < 4; j++) {
			portal.planes[j] = portal.planes[j] * itransform;
			portal.planes[j] /= length(vec3(portal.planes[j]));
			portal.tplanes[0][j] = portal.planes[j].x;
			portal.tplanes[1][j] = portal.planes[j].y;
			portal.tplanes[2][j] = portal.planes[j].z;
			portal.tplanes[3][j] = portal.planes[j].w;
		}
	}
}

void BoundFrustum::setITransform(const dmat4 &itransform) {
	setITransform(mat4(itransform));
}

/*
 */
int BoundFrustum::compare(const BoundFrustum &bf) const {
	
	if(valid != bf.valid) return 0;
	if(camera != bf.camera) return 0;
	if(portals.size() != bf.portals.size()) return 0;
	
	// clipping planes
	if(planes[0] != bf.planes[0]) return 0;
	if(planes[1] != bf.planes[1]) return 0;
	if(planes[2] != bf.planes[2]) return 0;
	if(planes[3] != bf.planes[3]) return 0;
	if(planes[4] != bf.planes[4]) return 0;
	if(planes[5] != bf.planes[5]) return 0;
	
	// portals
	for(int i = 0; i < portals.size(); i++) {
		if(portals[i].plane != bf.portals[i].plane) return 0;
		if(portals[i].planes[0] != bf.portals[i].planes[0]) return 0;
		if(portals[i].planes[1] != bf.portals[i].planes[1]) return 0;
		if(portals[i].planes[2] != bf.portals[i].planes[2]) return 0;
		if(portals[i].planes[3] != bf.portals[i].planes[3]) return 0;
	}
	
	return 1;
}

/*
 */
int BoundFrustum::addPortal(const vec3 *points,int num_points,const mat4 &transform) {
	
	assert(num_points == 4 && "BoundFrustum::addPortal(): bad points number");
	
	Portal portal;
	
	// portal points
	for(int i = 0; i < num_points; i++) {
		portal.points[i] = transform * points[i];
	}
	
	// check visibility
	if(inside(portal.points,num_points) == 0) return 0;
	
	// portal normal
	vec3 normal = cross(portal.points[1] - portal.points[0],portal.points[2] - portal.points[0]);
	float length = normal.length();
	if(length < EPSILON) return 0;
	normal /= length;
	
	// portal plane
	float angle = dot(camera - portal.points[0],normal);
	if(angle > 0.0f) normal = -normal;
	portal.plane = vec4(normal,-dot(normal,portal.points[0]));
	
	// clipping planes
	int j = (angle > 0.0f) ? num_points - 1 : 1;
	for(int i = 0; i < num_points; i++) {
		normal = normalize(cross(portal.points[i] - camera,portal.points[j] - camera));
		portal.planes[i] = vec4(normal,-dot(normal,camera));
		if(++j == num_points) j = 0;
	}
	
	// copy planes
	for(int i = 0; i < 4; i++) {
		portal.tplanes[0][i] = portal.planes[i].x;
		portal.tplanes[1][i] = portal.planes[i].y;
		portal.tplanes[2][i] = portal.planes[i].z;
		portal.tplanes[3][i] = portal.planes[i].w;
	}
	
	// add portal
	portals.append(portal);
	
	return 1;
}

void BoundFrustum::removePortal() {
	assert(portals.size() > 0 && "BoundFrustum::removePortal(): portals underflow");
	portals.remove();
}

int BoundFrustum::getNumPortals() const {
	return portals.size();
}

const vec3 *BoundFrustum::getPortalPoints(int num) const {
	assert(num >= 0 && num < portals.size() && "BoundFrustum::getPortalPoints(): bad portal number");
	return portals[num].points;
}

/*
 */
void BoundFrustum::expand(float radius) {
	
	// clipping planes
	for(int i = 0; i < 6; i++) {
		planes[i].w += radius;
	}
	for(int i = 0; i < 4; i++) {
		tplanes[0][i] = planes[i].x;
		tplanes[1][i] = planes[i].y;
		tplanes[2][i] = planes[i].z;
		tplanes[3][i] = planes[i].w;
	}
	for(int i = 0, j = 4; i < 2; i++, j++) {
		tplanes[4][i] = planes[j].x;
		tplanes[5][i] = planes[j].y;
		tplanes[6][i] = planes[j].z;
		tplanes[7][i] = planes[j].w;
	}
}

/*
 */
static INLINE int inside_plane(const vec4 &plane,const vec3 &min,const vec3 &max) {
	#ifdef USE_SSE
		__m128 min_xyz = _mm_mul_ps(min.vec,plane.vec);
		__m128 max_xyz = _mm_mul_ps(max.vec,plane.vec);
		__m128 min_max_x = _mm_shuffle_ps(min_xyz,max_xyz,_MM_PERM2(X,X,X,X));
		__m128 min_max_y = _mm_shuffle_ps(min_xyz,max_xyz,_MM_PERM2(Y,Y,Y,Y));
		min_max_x = _MM_SWIZZLE(min_max_x,X,Z,X,Z);
		min_max_y = _mm_add_ps(min_max_y,_MM_SWIZZLE(plane.vec,W,W,W,W));
		__m128 min_max_xy = _mm_add_ps(min_max_x,min_max_y);
		__m128 res_0 = _mm_add_ps(min_max_xy,_MM_SWIZZLE(min_xyz,Z,Z,Z,Z));
		__m128 res_1 = _mm_add_ps(min_max_xy,_MM_SWIZZLE(max_xyz,Z,Z,Z,Z));
		if(_mm_movemask_ps(_mm_and_ps(res_0,res_1)) != 0x0f) return 1;
	#elif USE_ALTIVEC
		vec_float4 plane_w = vec_perm(vec_splats(0.0f),plane.vec,VEC_PERM4(LX,LX,RW,LX));
		vec_float4 min_xyz = vec_madd(min.vec,plane.vec,plane_w);
		vec_float4 max_xyz = vec_madd(max.vec,plane.vec,plane_w);
		vec_float4 min_max_x = vec_perm(min_xyz,max_xyz,VEC_PERM4(LX,RX,LX,RX));
		vec_float4 min_max_y = vec_perm(min_xyz,max_xyz,VEC_PERM4(LY,LY,RY,RY));
		vec_float4 min_max_xy = vec_add(min_max_x,min_max_y);
		vec_float4 res_0 = vec_add(min_max_xy,VEC_SWIZZLE(min_xyz,Z,Z,Z,Z));
		vec_float4 res_1 = vec_add(min_max_xy,VEC_SWIZZLE(max_xyz,Z,Z,Z,Z));
		vec_uint4 res_2 = (vec_uint4)vec_and(res_0,res_1);
		if((VEC_SWIZZLE(res_2,B0,B0,B0,B0)[0] & 0x80808080) != 0x80808080) return 1;
	#else
		float min_x = min.x * plane.x;
		float min_y = min.y * plane.y;
		float min_zw = min.z * plane.z + plane.w;
		float min_min_xy = min_x + min_y;
		if(min_min_xy + min_zw > 0.0f) return 1;
		float max_x = max.x * plane.x;
		float max_min_xy = max_x + min_y;
		if(max_min_xy + min_zw > 0.0f) return 1;
		float max_y = max.y * plane.y;
		float min_max_xy = min_x + max_y;
		if(min_max_xy + min_zw > 0.0f) return 1;
		float max_max_xy = max_x + max_y;
		if(max_max_xy + min_zw > 0.0f) return 1;
		float max_zw = max.z * plane.z + plane.w;
		if(min_min_xy + max_zw > 0.0f) return 1;
		if(max_min_xy + max_zw > 0.0f) return 1;
		if(min_max_xy + max_zw > 0.0f) return 1;
		if(max_max_xy + max_zw > 0.0f) return 1;
	#endif
	return 0;
}

static INLINE int inside_plane(const vec4 &plane,const vec3 *points,int num_points) {
	for(int i = 0; i < num_points; i++) {
		if(dot(plane,points[i]) > 0.0f) return 1;
	}
	return 0;
}

/*
 */
int BoundFrustum::inside_planes(const vec3 &point) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[3].vec));
		if(_mm_movemask_ps(res_3) & 0x0f) return 0;
		res_0 = _mm_mul_ps(tplanes[4].vec,_mm_set1_ps(point.x));
		res_1 = _mm_mul_ps(tplanes[5].vec,_mm_set1_ps(point.y));
		res_2 = _mm_mul_ps(tplanes[6].vec,_mm_set1_ps(point.z));
		res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[7].vec));
		if(_mm_movemask_ps(res_3) & 0x03) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(tplanes[0].vec,vec_splats(point.x),tplanes[3].vec);
		vec_float4 res_1 = vec_madd(tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(tplanes[2].vec,vec_splats(point.z),res_1);
		if(VEC_SWIZZLE((vec_uint4)res_2,B0,B0,B0,B0)[0] & 0x80808080) return 0;
		res_0 = vec_madd(tplanes[4].vec,vec_splats(point.x),tplanes[7].vec);
		res_1 = vec_madd(tplanes[5].vec,vec_splats(point.y),res_0);
		res_2 = vec_madd(tplanes[6].vec,vec_splats(point.z),res_1);
		if(VEC_SWIZZLE((vec_uint4)res_2,B0,B0,B0,B0)[0] & 0x80800000) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(tplanes[3].vec,tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(0.0f)))) return 0;
		float32x2_t res_3 = vmla_n_f32(vget_low_f32(tplanes[7].vec),vget_low_f32(tplanes[4].vec),point.x);
		float32x2_t res_4 = vmla_n_f32(res_3,vget_low_f32(tplanes[5].vec),point.y);
		float32x2_t res_5 = vmla_n_f32(res_4,vget_low_f32(tplanes[6].vec),point.z);
		if(vmask_u32(vclt_f32(res_5,vdup_n_f32(0.0f)))) return 0;
	#else
		if(dot(planes[0],point) < 0.0f) return 0;
		if(dot(planes[1],point) < 0.0f) return 0;
		if(dot(planes[2],point) < 0.0f) return 0;
		if(dot(planes[3],point) < 0.0f) return 0;
		if(dot(planes[4],point) < 0.0f) return 0;
		if(dot(planes[5],point) < 0.0f) return 0;
	#endif
	return 1;
}

int BoundFrustum::inside_planes(const vec3 &point,float radius) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[3].vec));
		if(_mm_movemask_ps(_mm_add_ps(res_3,_mm_set1_ps(radius))) & 0x0f) return 0;
		res_0 = _mm_mul_ps(tplanes[4].vec,_mm_set1_ps(point.x));
		res_1 = _mm_mul_ps(tplanes[5].vec,_mm_set1_ps(point.y));
		res_2 = _mm_mul_ps(tplanes[6].vec,_mm_set1_ps(point.z));
		res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[7].vec));
		if(_mm_movemask_ps(_mm_add_ps(res_3,_mm_set1_ps(radius))) & 0x03) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(tplanes[0].vec,vec_splats(point.x),tplanes[3].vec);
		vec_float4 res_1 = vec_madd(tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(tplanes[2].vec,vec_splats(point.z),res_1);
		vec_uint4 res_3 = (vec_uint4)vec_add(res_2,vec_splats(radius));
		if(VEC_SWIZZLE(res_3,B0,B0,B0,B0)[0] & 0x80808080) return 0;
		res_0 = vec_madd(tplanes[4].vec,vec_splats(point.x),tplanes[7].vec);
		res_1 = vec_madd(tplanes[5].vec,vec_splats(point.y),res_0);
		res_2 = vec_madd(tplanes[6].vec,vec_splats(point.z),res_1);
		res_3 = (vec_uint4)vec_add(res_2,vec_splats(radius));
		if(VEC_SWIZZLE(res_3,B0,B0,B0,B0)[0] & 0x80800000) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(tplanes[3].vec,tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(-radius)))) return 0;
		float32x2_t res_3 = vmla_n_f32(vget_low_f32(tplanes[7].vec),vget_low_f32(tplanes[4].vec),point.x);
		float32x2_t res_4 = vmla_n_f32(res_3,vget_low_f32(tplanes[5].vec),point.y);
		float32x2_t res_5 = vmla_n_f32(res_4,vget_low_f32(tplanes[6].vec),point.z);
		if(vmask_u32(vclt_f32(res_5,vdup_n_f32(-radius)))) return 0;
	#else
		if(dot(planes[0],point) < -radius) return 0;
		if(dot(planes[1],point) < -radius) return 0;
		if(dot(planes[2],point) < -radius) return 0;
		if(dot(planes[3],point) < -radius) return 0;
		if(dot(planes[4],point) < -radius) return 0;
		if(dot(planes[5],point) < -radius) return 0;
	#endif
	return 1;
}

int BoundFrustum::inside_planes(const vec3 &min,const vec3 &max) const {
	if(inside_plane(planes[0],min,max) == 0) return 0;
	if(inside_plane(planes[1],min,max) == 0) return 0;
	if(inside_plane(planes[2],min,max) == 0) return 0;
	if(inside_plane(planes[3],min,max) == 0) return 0;
	if(inside_plane(planes[4],min,max) == 0) return 0;
	if(inside_plane(planes[5],min,max) == 0) return 0;
	return 1;
}

int BoundFrustum::inside_planes(const vec3 *points,int num_points) const {
	if(inside_plane(planes[0],points,num_points) == 0) return 0;
	if(inside_plane(planes[1],points,num_points) == 0) return 0;
	if(inside_plane(planes[2],points,num_points) == 0) return 0;
	if(inside_plane(planes[3],points,num_points) == 0) return 0;
	if(inside_plane(planes[4],points,num_points) == 0) return 0;
	if(inside_plane(planes[5],points,num_points) == 0) return 0;
	return 1;
}

/*
 */
int BoundFrustum::inside_planes_fast(const vec3 &point) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[3].vec));
		if(_mm_movemask_ps(res_3) & 0x0f) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(tplanes[0].vec,vec_splats(point.x),tplanes[3].vec);
		vec_float4 res_1 = vec_madd(tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(tplanes[2].vec,vec_splats(point.z),res_1);
		if(VEC_SWIZZLE((vec_uint4)res_2,B0,B0,B0,B0)[0] & 0x80808080) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(tplanes[3].vec,tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(0.0f)))) return 0;
	#else
		if(dot(planes[0],point) < 0.0f) return 0;
		if(dot(planes[1],point) < 0.0f) return 0;
		if(dot(planes[2],point) < 0.0f) return 0;
		if(dot(planes[3],point) < 0.0f) return 0;
	#endif
	return 1;
}

int BoundFrustum::inside_planes_fast(const vec3 &point,float radius) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,tplanes[3].vec));
		if(_mm_movemask_ps(_mm_add_ps(res_3,_mm_set1_ps(radius))) & 0x0f) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(tplanes[0].vec,vec_splats(point.x),tplanes[3].vec);
		vec_float4 res_1 = vec_madd(tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(tplanes[2].vec,vec_splats(point.z),res_1);
		vec_uint4 res_3 = (vec_uint4)vec_add(res_2,vec_splats(radius));
		if(VEC_SWIZZLE(res_3,B0,B0,B0,B0)[0] & 0x80808080) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(tplanes[3].vec,tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(-radius)))) return 0;
	#else
		if(dot(planes[0],point) < -radius) return 0;
		if(dot(planes[1],point) < -radius) return 0;
		if(dot(planes[2],point) < -radius) return 0;
		if(dot(planes[3],point) < -radius) return 0;
	#endif
	return 1;
}

int BoundFrustum::inside_planes_fast(const vec3 &min,const vec3 &max) const {
	if(inside_plane(planes[0],min,max) == 0) return 0;
	if(inside_plane(planes[1],min,max) == 0) return 0;
	if(inside_plane(planes[2],min,max) == 0) return 0;
	if(inside_plane(planes[3],min,max) == 0) return 0;
	return 1;
}

int BoundFrustum::inside_planes_fast(const vec3 *points,int num_points) const {
	if(inside_plane(planes[0],points,num_points) == 0) return 0;
	if(inside_plane(planes[1],points,num_points) == 0) return 0;
	if(inside_plane(planes[2],points,num_points) == 0) return 0;
	if(inside_plane(planes[3],points,num_points) == 0) return 0;
	return 1;
}

/*
 */
INLINE int BoundFrustum::inside_portal(const Portal &portal,const vec3 &point) const {
	if(dot(portal.plane,point) < 0.0f) return 0;
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(portal.tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(portal.tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(portal.tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,portal.tplanes[3].vec));
		if(_mm_movemask_ps(res_3)) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(portal.tplanes[0].vec,vec_splats(point.x),portal.tplanes[3].vec);
		vec_float4 res_1 = vec_madd(portal.tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(portal.tplanes[2].vec,vec_splats(point.z),res_1);
		if(VEC_SWIZZLE((vec_uint4)res_2,B0,B0,B0,B0)[0] & 0x80808080) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(portal.tplanes[3].vec,portal.tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,portal.tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,portal.tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(0.0f)))) return 0;
	#else
		if(dot(portal.planes[0],point) < 0.0f) return 0;
		if(dot(portal.planes[1],point) < 0.0f) return 0;
		if(dot(portal.planes[2],point) < 0.0f) return 0;
		if(dot(portal.planes[3],point) < 0.0f) return 0;
	#endif
	return 1;
}

INLINE int BoundFrustum::inside_portal(const Portal &portal,const vec3 &point,float radius) const {
	if(dot(portal.plane,point) < -radius) return 0;
	#ifdef USE_SSE
		__m128 res_0 = _mm_mul_ps(portal.tplanes[0].vec,_mm_set1_ps(point.x));
		__m128 res_1 = _mm_mul_ps(portal.tplanes[1].vec,_mm_set1_ps(point.y));
		__m128 res_2 = _mm_mul_ps(portal.tplanes[2].vec,_mm_set1_ps(point.z));
		__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,portal.tplanes[3].vec));
		if(_mm_movemask_ps(_mm_add_ps(res_3,_mm_set1_ps(radius)))) return 0;
	#elif USE_ALTIVEC
		vec_float4 res_0 = vec_madd(portal.tplanes[0].vec,vec_splats(point.x),portal.tplanes[3].vec);
		vec_float4 res_1 = vec_madd(portal.tplanes[1].vec,vec_splats(point.y),res_0);
		vec_float4 res_2 = vec_madd(portal.tplanes[2].vec,vec_splats(point.z),res_1);
		vec_uint4 res_3 = (vec_uint4)vec_add(res_2,vec_splats(radius));
		if(VEC_SWIZZLE(res_3,B0,B0,B0,B0)[0] & 0x80808080) return 0;
	#elif USE_NEON
		float32x4_t res_0 = vmlaq_n_f32(portal.tplanes[3].vec,portal.tplanes[0].vec,point.x);
		float32x4_t res_1 = vmlaq_n_f32(res_0,portal.tplanes[1].vec,point.y);
		float32x4_t res_2 = vmlaq_n_f32(res_1,portal.tplanes[2].vec,point.z);
		if(vmaskq_u32(vcltq_f32(res_2,vdupq_n_f32(-radius)))) return 0;
	#else
		if(dot(portal.planes[0],point) < -radius) return 0;
		if(dot(portal.planes[1],point) < -radius) return 0;
		if(dot(portal.planes[2],point) < -radius) return 0;
		if(dot(portal.planes[3],point) < -radius) return 0;
	#endif
	return 1;
}

INLINE int BoundFrustum::inside_portal(const Portal &portal,const vec3 &min,const vec3 &max) const {
	if(inside_plane(portal.plane,min,max) == 0) return 0;
	if(inside_plane(portal.planes[0],min,max) == 0) return 0;
	if(inside_plane(portal.planes[1],min,max) == 0) return 0;
	if(inside_plane(portal.planes[2],min,max) == 0) return 0;
	if(inside_plane(portal.planes[3],min,max) == 0) return 0;
	return 1;
}

INLINE int BoundFrustum::inside_portal(const Portal &portal,const vec3 *points,int num_points) const {
	if(inside_plane(portal.plane,points,num_points) == 0) return 0;
	if(inside_plane(portal.planes[0],points,num_points) == 0) return 0;
	if(inside_plane(portal.planes[1],points,num_points) == 0) return 0;
	if(inside_plane(portal.planes[2],points,num_points) == 0) return 0;
	if(inside_plane(portal.planes[3],points,num_points) == 0) return 0;
	return 1;
}

/*
 */
int BoundFrustum::inside_portals(const vec3 &point) const {
	for(int i = 0; i < portals.size(); i++) {
		if(inside_portal(portals[i],point) == 0) return 0;
	}
	return 1;
}

int BoundFrustum::inside_portals(const vec3 &point,float radius) const {
	for(int i = 0; i < portals.size(); i++) {
		if(inside_portal(portals[i],point,radius) == 0) return 0;
	}
	return 1;
}

int BoundFrustum::inside_portals(const vec3 &min,const vec3 &max) const {
	for(int i = 0; i < portals.size(); i++) {
		if(inside_portal(portals[i],min,max) == 0) return 0;
	}
	return 1;
}

int BoundFrustum::inside_portals(const vec3 *points,int num_points) const {
	for(int i = 0; i < portals.size(); i++) {
		if(inside_portal(portals[i],points,num_points) == 0) return 0;
	}
	return 1;
}

/*
 */
int BoundFrustum::inside(const BoundSphere &bs) const {
	if(bs.isValid()) return insideValid(bs);
	return 0;
}

int BoundFrustum::inside(const BoundBox &bb) const {
	if(bb.isValid()) return insideValid(bb);
	return 0;
}

int BoundFrustum::inside(const BoundFrustum &bf) const {
	if(bf.isValid()) return insideValid(bf);
	return 0;
}

/*
 */
int BoundFrustum::insideAll(const BoundSphere &bs) const {
	if(bs.isValid()) return insideAllValid(bs);
	return 0;
}

int BoundFrustum::insideAll(const BoundBox &bb) const {
	if(bb.isValid()) return insideAllValid(bb);
	return 0;
}

int BoundFrustum::insideAll(const BoundFrustum &bf) const {
	if(bf.isValid()) return insideAllValid(bf);
	return 0;
}

/*
 */
int BoundFrustum::insidePlanes(const BoundSphere &bs) const {
	if(bs.isValid()) return insidePlanesValid(bs);
	return 0;
}

int BoundFrustum::insidePlanes(const BoundBox &bb) const {
	if(bb.isValid()) return insidePlanesValid(bb);
	return 0;
}

int BoundFrustum::insidePlanes(const BoundFrustum &bf) const {
	if(bf.isValid()) return insidePlanesValid(bf);
	return 0;
}

/*
 */
int BoundFrustum::insidePortals(const BoundSphere &bs) const {
	if(bs.isValid()) return insidePortalsValid(bs);
	return 0;
}

int BoundFrustum::insidePortals(const BoundBox &bb) const {
	if(bb.isValid()) return insidePortalsValid(bb);
	return 0;
}

int BoundFrustum::insidePortals(const BoundFrustum &bf) const {
	if(bf.isValid()) return insidePortalsValid(bf);
	return 0;
}

/*
 */
int BoundFrustum::insideShadow(const vec3 &direction,const BoundSphere &object) const {
	
	if(object.isValid()) {
		
		// object is inside the bound frustum
		if(inside_planes(object.getCenter(),object.getRadius())) return 1;
		
		// shadow volume is inside the bound frustum
		for(int i = 0; i < 6; i++) {
			float k = dot3(planes[i],direction);
			if(Math::abs(k) < EPSILON) continue;
			k = -dot(planes[i],object.getCenter()) / k;
			if(k > object.getRadius()) continue;
			if(inside_planes(object.getCenter() + direction * k,object.getRadius())) return 1;
		}
		
		return 0;
	}
	return 0;
}

int BoundFrustum::insideShadow(const BoundSphere &light,const vec3 &offset,const BoundSphere &object) const {
	
	if(light.isValid() && object.isValid()) {
		
		// object is outside the light bounds
		if(light.inside(object) == 0) return 0;
		
		// object is inside the bound frustum
		if(inside_planes(object.getCenter(),object.getRadius())) return 1;
		
		// direction from light center to object center
		vec3 direction = object.getCenter() - offset;
		float distance = length(direction);
		if(distance < object.getRadius() + EPSILON) return 1;
		direction /= distance;
		
		// basis
		vec3 x,y;
		orthoBasis(direction,x,y);
		
		// near points
		vec3 x0 = x * object.getRadius();
		vec3 y0 = y * object.getRadius();
		vec3 z0 = offset + direction * (distance - object.getRadius());
		
		// far points
		float radius = light.getRadius() + dot(direction,light.getCenter() - offset);
		float k = object.getRadius() * radius / (distance - object.getRadius());
		vec3 x1 = x * k;
		vec3 y1 = y * k;
		vec3 z1 = offset + direction * radius;
		
		// check visibility
		vec3 points[8] = {
			x0 + y0 + z0,x0 - y0 + z0,-x0 - y0 + z0,-x0 + y0 + z0,
			x1 + y1 + z1,x1 - y1 + z1,-x1 - y1 + z1,-x1 + y1 + z1,
		};
		
		return inside_planes(points,8);
	}
	return 0;
}
