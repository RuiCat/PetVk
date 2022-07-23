/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Bounds.h
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

#ifndef __BOUNDS_H__
#define __BOUNDS_H__

#include "Vector.h"
#include "Geometry.h"

/*
 */
class BoundSphere;
class BoundBox;
class BoundFrustum;

/******************************************************************************\
*
* BoundSphere
*
\******************************************************************************/

/*
 */
class BoundSphere {
		
	public:
		
		BoundSphere();
		BoundSphere(const vec3 &center,float radius);
		BoundSphere(const vec3 *points,int num_points,int optimal);
		BoundSphere(const BoundSphere &bs);
		BoundSphere(const BoundSphere &bs,const mat4 &transform);
		explicit BoundSphere(const BoundBox &bb);
		~BoundSphere();
		
		void clear();
		
		void set(const vec3 &center,float radius);
		void set(const vec3 *points,int num_points,int optimal);
		void set(const BoundSphere &bs);
		void set(const BoundSphere &bs,const mat4 &transform);
		void set(const BoundBox &bb);
		
		// transformation
		void setTransform(const mat4 &transform);
		void setTransform(const dmat4 &transform);
		
		// compare
		int compare(const BoundSphere &bs) const;
		INLINE int operator==(const BoundSphere &bs) const { return compare(bs); }
		INLINE int operator!=(const BoundSphere &bs) const { return !compare(bs); }
		
		// expand
		void expand(const vec3 &point);
		void expand(const vec3 *points,int num_points);
		void expand(const BoundSphere &bs);
		void expand(const BoundBox &bb);
		
		// radius expand
		void expandRadius(const vec3 &point);
		void expandRadius(const vec3 *points,int num_points);
		void expandRadius(const BoundSphere &bs);
		void expandRadius(const BoundBox &bb);
		
		// inside points
		int inside(const vec3 &point) const;
		int inside(const vec3 &point,float radius) const;
		int inside(const vec3 &min,const vec3 &max) const;
		
		int insideValid(const vec3 &point) const;
		int insideValid(const vec3 &point,float radius) const;
		int insideValid(const vec3 &min,const vec3 &max) const;
		
		// inside bounds
		int inside(const BoundSphere &bs) const;
		int inside(const BoundBox &bb) const;
		
		int insideValid(const BoundSphere &bs) const;
		int insideValid(const BoundBox &bb) const;
		
		int insideAll(const BoundSphere &bs) const;
		int insideAll(const BoundBox &bb) const;
		
		int insideAllValid(const BoundSphere &bs) const;
		int insideAllValid(const BoundBox &bb) const;
		
		// intersections
		int rayIntersection(const vec3 &point,const vec3 &direction) const;
		int getIntersection(const vec3 &p0,const vec3 &p1) const;
		
		int rayIntersectionValid(const vec3 &point,const vec3 &direction) const;
		int getIntersectionValid(const vec3 &p0,const vec3 &p1) const;
		
		// distance
		float distance() const;
		float distance(const vec3 &point) const;
		
		float distanceValid() const;
		float distanceValid(const vec3 &point) const;
		
		// parameters
		INLINE int isValid() const { return (center.w > 0.0f); }
		INLINE const vec3 &getCenter() const { return center; }
		INLINE float getRadius() const { return center.w; }
		
	private:
		
		vec3 center;	// bounding sphere center and radius
};

/*
 */
BoundSphere operator*(const mat4 &m,const BoundSphere &bs);
BoundSphere operator*(const dmat4 &m,const BoundSphere &bs);

/*
 */
INLINE int BoundSphere::insideValid(const vec3 &point) const {
	return (length2(center - point) <= center.w * center.w);
}

INLINE int BoundSphere::insideValid(const vec3 &point,float radius) const {
	radius += center.w;
	return (length2(center - point) <= radius * radius);
}

INLINE int BoundSphere::insideValid(const vec3 &min,const vec3 &max) const {
	#ifdef USE_SSE
		__m128 radius = _MM_SWIZZLE(center.vec,W,W,W,W);
		__m128 res_0 = _mm_cmplt_ps(_mm_add_ps(center.vec,radius),min.vec);
		__m128 res_1 = _mm_cmpgt_ps(_mm_sub_ps(center.vec,radius),max.vec);
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_float4 radius = VEC_SWIZZLE(center.vec,W,W,W,W);
		vec_uint4 res_0 = (vec_uint4)vec_cmplt(vec_add(center.vec,radius),min.vec);
		vec_uint4 res_1 = (vec_uint4)vec_cmpgt(vec_sub(center.vec,radius),max.vec);
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		float32x4_t radius = vdupq_lane_f32(vget_high_f32(center.vec),1);
		uint32x4_t res_0 = vcltq_f32(vaddq_f32(center.vec,radius),min.vec);
		uint32x4_t res_1 = vcgtq_f32(vsubq_f32(center.vec,radius),max.vec);
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		if(center.x + center.w < min.x || center.x - center.w > max.x) return 0;
		if(center.y + center.w < min.y || center.y - center.w > max.y) return 0;
		if(center.z + center.w < min.z || center.z - center.w > max.z) return 0;
		return 1;
	#endif
}

/*
 */
INLINE int BoundSphere::insideValid(const BoundSphere &bs) const {
	return insideValid(bs.center,bs.center.w);
}

/*
 */
INLINE int BoundSphere::insideAllValid(const BoundSphere &bs) const {
	float radius = center.w - bs.center.w;
	if(radius > 0.0f) return (length2(center - bs.center) <= radius * radius);
	return 0;
}

/*
 */
INLINE int BoundSphere::rayIntersectionValid(const vec3 &point,const vec3 &direction) const {
	float k = saturate(dot(direction,center - point) / length2(direction));
	return (length2(center - point - direction * k) <= center.w * center.w);
}

INLINE int BoundSphere::getIntersectionValid(const vec3 &p0,const vec3 &p1) const {
	return rayIntersectionValid(p0,p1 - p0);
}

/*
 */
INLINE float BoundSphere::distanceValid() const {
	#ifdef USE_SSE
		float ret;
		__m128 direction = _mm_rcp_ss(_mm_rsqrt_ss(_mm_dot33_ps(center.vec,center.vec)));
		direction = _mm_sub_ss(direction,_MM_SWIZZLE(center.vec,W,W,W,W));
		_mm_store_ss(&ret,direction);
		return ret;
	#elif USE_ALTIVEC
		vec4 ret;
		vec_float4 direction = vec_re(vec_rsqrte(vec_dot33(center.vec,center.vec)));
		direction = vec_sub(direction,VEC_SWIZZLE(center.vec,W,W,W,W));
		vec_ste(direction,0,&ret.x);
		return ret.x;
	#elif USE_NEON
		float32x4_t direction = vrecpeq_f32(vrsqrteq_f32(vdot33q_f32(center.vec,center.vec)));
		direction = vsubq_f32(direction,center.vec);
		return vgetq_lane_f32(direction,3);
	#else
		return Math::sqrtFast(center.length2()) - center.w;
	#endif
}

INLINE float BoundSphere::distanceValid(const vec3 &point) const {
	#ifdef USE_SSE
		float ret;
		__m128 direction = _mm_sub_ps(center.vec,point.vec);
		direction = _mm_rcp_ss(_mm_rsqrt_ss(_mm_dot33_ps(direction,direction)));
		direction = _mm_sub_ss(direction,_MM_SWIZZLE(center.vec,W,W,W,W));
		_mm_store_ss(&ret,direction);
		return ret;
	#elif USE_ALTIVEC
		vec4 ret;
		vec_float4 direction = vec_sub(center.vec,point.vec);
		direction = vec_re(vec_rsqrte(vec_dot33(direction,direction)));
		direction = vec_sub(direction,VEC_SWIZZLE(center.vec,W,W,W,W));
		vec_ste(direction,0,&ret.x);
		return ret.x;
	#elif USE_NEON
		float32x4_t direction = vsubq_f32(center.vec,point.vec);
		direction = vrecpeq_f32(vrsqrteq_f32(vdot33q_f32(direction,direction)));
		direction = vsubq_f32(direction,center.vec);
		return vgetq_lane_f32(direction,3);
	#else
		vec3 direction;
		sub(direction,center,point);
		return Math::sqrtFast(direction.length2()) - center.w;
	#endif
}

/******************************************************************************\
*
* BoundBox
*
\******************************************************************************/

/*
 */
class BoundBox {
		
	public:
		
		BoundBox();
		BoundBox(const vec3 &min,const vec3 &max);
		BoundBox(const vec3 *points,int num_points);
		BoundBox(const BoundBox &bb);
		BoundBox(const BoundBox &bb,const mat4 &transform);
		explicit BoundBox(const BoundSphere &bs);
		~BoundBox();
		
		void clear();
		
		void set(const vec3 &min,const vec3 &max);
		void set(const vec3 *points,int num_points);
		void set(const BoundSphere &bs);
		void set(const BoundBox &bb);
		void set(const BoundBox &bb,const mat4 &transform);
		
		// transformation
		void setTransform(const mat4 &transform);
		void setTransform(const dmat4 &transform);
		void setTransform(const BoundSphere &bs,const mat4 &transform);
		void setTransform(const BoundSphere &bs,const dmat4 &transform);
		
		// compare
		int compare(const BoundBox &bb) const;
		INLINE int operator==(const BoundBox &bb) const { return compare(bb); }
		INLINE int operator!=(const BoundBox &bb) const { return !compare(bb); }
		
		// expand
		void expand(const vec3 &point);
		void expand(const vec3 *points,int num_points);
		void expand(const BoundSphere &bs);
		void expand(const BoundBox &bb);
		
		// inside points
		int inside(const vec3 &point) const;
		int inside(const vec3 &point,float radius) const;
		int inside(const vec3 &min,const vec3 &max) const;
		
		int insideValid(const vec3 &point) const;
		int insideValid(const vec3 &point,float radius) const;
		int insideValid(const vec3 &min,const vec3 &max) const;
		
		// inside bounds
		int inside(const BoundSphere &bs) const;
		int inside(const BoundBox &bb) const;
		
		int insideValid(const BoundSphere &bs) const;
		int insideValid(const BoundBox &bb) const;
		
		int insideAll(const BoundSphere &bs) const;
		int insideAll(const BoundBox &bb) const;
		
		int insideAllValid(const BoundSphere &bs) const;
		int insideAllValid(const BoundBox &bb) const;
		
		// inside cube
		int insideCube(int face,const vec3 &offset) const;
		
		// intersections
		int rayIntersection(const vec3 &point,const vec3 &direction) const;
		int irayIntersection(const vec3 &point,const vec3 &idirection) const;
		int getIntersection(const vec3 &p0,const vec3 &p1) const;
		
		int rayIntersectionValid(const vec3 &point,const vec3 &direction) const;
		int irayIntersectionValid(const vec3 &point,const vec3 &idirection) const;
		int getIntersectionValid(const vec3 &p0,const vec3 &p1) const;
		
		// distance
		float distance() const;
		float distance(const vec3 &point) const;
		
		float distanceValid() const;
		float distanceValid(const vec3 &point) const;
		
		// parameters
		INLINE int isValid() const { return (min.x <= max.x); }
		INLINE const vec3 &getMin() const { return min; }
		INLINE const vec3 &getMax() const { return max; }
		void getPoints(vec3 *points,int num_points) const;
		void getPlanes(vec4 *planes,int num_planes) const;
		
	private:
		
		vec3 min;		// bounding box minimum
		vec3 max;		// bounding box maximum
};

/*
 */
BoundBox operator*(const mat4 &m,const BoundBox &bb);
BoundBox operator*(const dmat4 &m,const BoundBox &bb);

/*
 */
INLINE int BoundSphere::insideValid(const BoundBox &bb) const {
	return insideValid(bb.getMin(),bb.getMax());
}

/*
 */
INLINE int BoundSphere::insideAllValid(const BoundBox &bb) const {
	const vec3 &min = bb.getMin();
	const vec3 &max = bb.getMax();
	if(insideValid(vec3(min.x,min.y,min.z)) == 0) return 0;
	if(insideValid(vec3(max.x,min.y,min.z)) == 0) return 0;
	if(insideValid(vec3(min.x,max.y,min.z)) == 0) return 0;
	if(insideValid(vec3(max.x,max.y,min.z)) == 0) return 0;
	if(insideValid(vec3(min.x,min.y,max.z)) == 0) return 0;
	if(insideValid(vec3(max.x,min.y,max.z)) == 0) return 0;
	if(insideValid(vec3(min.x,max.y,max.z)) == 0) return 0;
	if(insideValid(vec3(max.x,max.y,max.z)) == 0) return 0;
	return 1;
}

/*
 */
INLINE int BoundBox::insideValid(const vec3 &point) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_cmpgt_ps(min.vec,point.vec);
		__m128 res_1 = _mm_cmplt_ps(max.vec,point.vec);
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,point.vec);
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,point.vec);
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		uint32x4_t res_0 = vcgtq_f32(min.vec,point.vec);
		uint32x4_t res_1 = vcltq_f32(max.vec,point.vec);
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		if(min.x > point.x || max.x < point.x) return 0;
		if(min.y > point.y || max.y < point.y) return 0;
		if(min.z > point.z || max.z < point.z) return 0;
		return 1;
	#endif
}

INLINE int BoundBox::insideValid(const vec3 &point,float radius) const {
	#ifdef USE_SSE
		__m128 r = vec3(radius).vec;
		__m128 res_0 = _mm_cmpgt_ps(min.vec,_mm_add_ps(point.vec,r));
		__m128 res_1 = _mm_cmplt_ps(max.vec,_mm_sub_ps(point.vec,r));
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_float4 r = vec_splats(radius);
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,vec_add(point.vec,r));
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,vec_sub(point.vec,r));
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		float32x4_t r = vdupq_n_f32(radius);
		uint32x4_t res_0 = vcgtq_f32(min.vec,vaddq_f32(point.vec,r));
		uint32x4_t res_1 = vcltq_f32(max.vec,vsubq_f32(point.vec,r));
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		if(min.x > point.x + radius || max.x < point.x - radius) return 0;
		if(min.y > point.y + radius || max.y < point.y - radius) return 0;
		if(min.z > point.z + radius || max.z < point.z - radius) return 0;
		return 1;
	#endif
}

INLINE int BoundBox::insideValid(const vec3 &min_,const vec3 &max_) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_cmpgt_ps(min.vec,max_.vec);
		__m128 res_1 = _mm_cmplt_ps(max.vec,min_.vec);
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,max_.vec);
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,min_.vec);
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		uint32x4_t res_0 = vcgtq_f32(min.vec,max_.vec);
		uint32x4_t res_1 = vcltq_f32(max.vec,min_.vec);
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		if(min.x > max_.x || max.x < min_.x) return 0;
		if(min.y > max_.y || max.y < min_.y) return 0;
		if(min.z > max_.z || max.z < min_.z) return 0;
		return 1;
	#endif
}

/*
 */
INLINE int BoundBox::insideValid(const BoundSphere &bs) const {
	const vec3 &center = bs.getCenter();
	#ifdef USE_SSE
		__m128 radius = _MM_SWIZZLE(center.vec,W,W,W,W);
		__m128 res_0 = _mm_cmpgt_ps(min.vec,_mm_add_ps(center.vec,radius));
		__m128 res_1 = _mm_cmplt_ps(max.vec,_mm_sub_ps(center.vec,radius));
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_float4 radius = VEC_SWIZZLE(center.vec,W,W,W,W);
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,vec_add(center.vec,radius));
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,vec_sub(center.vec,radius));
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		float32x4_t radius = vdupq_lane_f32(vget_high_f32(center.vec),1);
		uint32x4_t res_0 = vcgtq_f32(min.vec,vaddq_f32(center.vec,radius));
		uint32x4_t res_1 = vcltq_f32(max.vec,vsubq_f32(center.vec,radius));
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		float radius = bs.getRadius();
		if(min.x > center.x + radius || max.x < center.x - radius) return 0;
		if(min.y > center.y + radius || max.y < center.y - radius) return 0;
		if(min.z > center.z + radius || max.z < center.z - radius) return 0;
		return 1;
	#endif
}

INLINE int BoundBox::insideValid(const BoundBox &bb) const {
	return insideValid(bb.min,bb.max);
}

/*
 */
INLINE int BoundBox::insideAllValid(const BoundSphere &bs) const {
	const vec3 &center = bs.getCenter();
	#ifdef USE_SSE
		__m128 radius = _MM_SWIZZLE(center.vec,W,W,W,W);
		__m128 res_0 = _mm_cmpgt_ps(min.vec,_mm_sub_ps(center.vec,radius));
		__m128 res_1 = _mm_cmplt_ps(max.vec,_mm_add_ps(center.vec,radius));
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_float4 radius = VEC_SWIZZLE(center.vec,W,W,W,W);
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,vec_sub(center.vec,radius));
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,vec_add(center.vec,radius));
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		float32x4_t radius = vdupq_lane_f32(vget_high_f32(center.vec),1);
		uint32x4_t res_0 = vcgtq_f32(min.vec,vsubq_f32(center.vec,radius));
		uint32x4_t res_1 = vcltq_f32(max.vec,vaddq_f32(center.vec,radius));
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		float radius = bs.getRadius();
		if(min.x > center.x - radius || max.x < center.x + radius) return 0;
		if(min.y > center.y - radius || max.y < center.y + radius) return 0;
		if(min.z > center.z - radius || max.z < center.z + radius) return 0;
		return 1;
	#endif
	return 0;
}

INLINE int BoundBox::insideAllValid(const BoundBox &bb) const {
	#ifdef USE_SSE
		__m128 res_0 = _mm_cmpgt_ps(min.vec,bb.min.vec);
		__m128 res_1 = _mm_cmplt_ps(max.vec,bb.max.vec);
		return ((_mm_movemask_ps(_mm_or_ps(res_0,res_1)) & 0x07) == 0);
	#elif USE_ALTIVEC
		vec_uint4 res_0 = (vec_uint4)vec_cmpgt(min.vec,bb.min.vec);
		vec_uint4 res_1 = (vec_uint4)vec_cmplt(max.vec,bb.max.vec);
		return ((vec_perm(vec_or(res_0,res_1),res_1,VEC_PERM2(B0,B0,B0,B0))[0] & 0xffffff00) == 0);
	#elif USE_NEON
		uint32x4_t res_0 = vcgtq_f32(min.vec,bb.min.vec);
		uint32x4_t res_1 = vcltq_f32(max.vec,bb.max.vec);
		return ((vmaskq_u32(vorrq_u32(res_0,res_1)) & 0x07) == 0);
	#else
		if(min.x > bb.min.x || max.x < bb.max.x) return 0;
		if(min.y > bb.min.y || max.y < bb.max.y) return 0;
		if(min.z > bb.min.z || max.z < bb.max.z) return 0;
		return 1;
	#endif
	return 0;
}

/*
 */
INLINE int BoundBox::rayIntersectionValid(const vec3 &point,const vec3 &direction) const {
	return rayBoundBoxIntersection(point,direction,min,max);
}

INLINE int BoundBox::irayIntersectionValid(const vec3 &point,const vec3 &idirection) const {
	return irayBoundBoxIntersection(point,idirection,min,max);
}

INLINE int BoundBox::getIntersectionValid(const vec3 &p0,const vec3 &p1) const {
	return rayBoundBoxIntersection(p0,p1 - p0,min,max);
}

/*
 */
INLINE float BoundBox::distanceValid() const {
	#ifdef USE_SSE
		float ret;
		__m128 direction = _mm_min_ps(_mm_max_ps(vec3_zero.vec,min.vec),max.vec);
		direction = _mm_rcp_ss(_mm_rsqrt_ss(_mm_dot33_ps(direction,direction)));
		_mm_store_ss(&ret,direction);
		return ret;
	#elif USE_ALTIVEC
		vec4 ret;
		vec_float4 direction = vec_min(vec_max(vec3_zero.vec,min.vec),max.vec);
		direction = vec_re(vec_rsqrte(vec_dot33(direction,direction)));
		vec_ste(direction,0,&ret.x);
		return ret.x;
	#elif USE_NEON
		float32x4_t direction = vminq_f32(vmaxq_f32(vec3_zero.vec,min.vec),max.vec);
		direction = vrecpeq_f32(vrsqrteq_f32(vdot33q_f32(direction,direction)));
		return vgetq_lane_f32(direction,0);
	#else
		vec3 direction;
		if(min.x > 0.0f) direction.x = min.x;
		else if(max.x < 0.0f) direction.x = max.x;
		else direction.x = 0.0f;
		if(min.y > 0.0f) direction.y = min.y;
		else if(max.y < 0.0f) direction.y = max.y;
		else direction.y = 0.0f;
		if(min.z > 0.0f) direction.z = min.z;
		else if(max.z < 0.0f) direction.z = max.z;
		else direction.z = 0.0f;
		return Math::sqrtFast(direction.length2());
	#endif
}

INLINE float BoundBox::distanceValid(const vec3 &point) const {
	#ifdef USE_SSE
		float ret;
		__m128 direction = _mm_sub_ps(_mm_min_ps(_mm_max_ps(point.vec,min.vec),max.vec),point.vec);
		direction = _mm_rcp_ss(_mm_rsqrt_ss(_mm_dot33_ps(direction,direction)));
		_mm_store_ss(&ret,direction);
		return ret;
	#elif USE_ALTIVEC
		vec4 ret;
		vec_float4 direction = vec_sub(vec_min(vec_max(point.vec,min.vec),max.vec),point.vec);
		direction = vec_re(vec_rsqrte(vec_dot33(direction,direction)));
		vec_ste(direction,0,&ret.x);
		return ret.x;
	#elif USE_NEON
		float32x4_t direction = vsubq_f32(vminq_f32(vmaxq_f32(point.vec,min.vec),max.vec),point.vec);
		direction = vrecpeq_f32(vrsqrteq_f32(vdot33q_f32(direction,direction)));
		return vgetq_lane_f32(direction,0);
	#else
		vec3 direction;
		if(min.x > point.x) direction.x = min.x - point.x;
		else if(max.x < point.x) direction.x = max.x - point.x;
		else direction.x = 0.0f;
		if(min.y > point.y) direction.y = min.y - point.y;
		else if(max.y < point.y) direction.y = max.y - point.y;
		else direction.y = 0.0f;
		if(min.z > point.z) direction.z = min.z - point.z;
		else if(max.z < point.z) direction.z = max.z - point.z;
		else direction.z = 0.0f;
		return Math::sqrtFast(direction.length2());
	#endif
}

/******************************************************************************\
*
* BoundFrustum
*
\******************************************************************************/

/*
 */
class BoundFrustum {
		
	public:
		
		BoundFrustum();
		BoundFrustum(const mat4 &projection,const mat4 &modelview);
		BoundFrustum(const BoundFrustum &bf);
		BoundFrustum(const BoundFrustum &bf,const mat4 &itransform);
		~BoundFrustum();
		
		void clear();
		
		void set(const mat4 &projection,const mat4 &modelview);
		void set(const BoundFrustum &bf);
		void set(const BoundFrustum &bf,const mat4 &itransform);
		
		// transformation
		void setITransform(const mat4 &itransform);
		void setITransform(const dmat4 &itransform);
		
		// compare
		int compare(const BoundFrustum &bf) const;
		INLINE int operator==(const BoundFrustum &bf) const { return compare(bf); }
		INLINE int operator!=(const BoundFrustum &bf) const { return !compare(bf); }
		
		// portals
		int addPortal(const vec3 *points,int num_points,const mat4 &transform);
		void removePortal();
		
		int getNumPortals() const;
		const vec3 *getPortalPoints(int num) const;
		
		// expand
		void expand(float radius);
		
		// inside points
		int inside(const vec3 &point) const;
		int inside(const vec3 &point,float radius) const;
		int inside(const vec3 &min,const vec3 &max) const;
		int inside(const vec3 *points,int num) const;
		
		int insideFast(const vec3 &point) const;
		int insideFast(const vec3 &point,float radius) const;
		int insideFast(const vec3 &min,const vec3 &max) const;
		int insideFast(const vec3 *points,int num) const;
		
		// inside bounds
		int inside(const BoundSphere &bs) const;
		int inside(const BoundBox &bb) const;
		int inside(const BoundFrustum &bf) const;
		
		int insideValid(const BoundSphere &bs) const;
		int insideValid(const BoundBox &bb) const;
		int insideValid(const BoundFrustum &bf) const;
		
		int insideValidFast(const BoundSphere &bs) const;
		int insideValidFast(const BoundBox &bb) const;
		int insideValidFast(const BoundFrustum &bf) const;
		
		// inside all bounds
		int insideAll(const BoundSphere &bs) const;
		int insideAll(const BoundBox &bb) const;
		int insideAll(const BoundFrustum &bf) const;
		
		int insideAllValid(const BoundSphere &bs) const;
		int insideAllValid(const BoundBox &bb) const;
		int insideAllValid(const BoundFrustum &bf) const;
		
		int insideAllValidFast(const BoundSphere &bs) const;
		int insideAllValidFast(const BoundBox &bb) const;
		int insideAllValidFast(const BoundFrustum &bf) const;
		
		// inside planes bounds
		int insidePlanes(const BoundSphere &bs) const;
		int insidePlanes(const BoundBox &bb) const;
		int insidePlanes(const BoundFrustum &bf) const;
		
		int insidePlanesValid(const BoundSphere &bs) const;
		int insidePlanesValid(const BoundBox &bb) const;
		int insidePlanesValid(const BoundFrustum &bf) const;
		
		int insidePlanesValidFast(const BoundSphere &bs) const;
		int insidePlanesValidFast(const BoundBox &bb) const;
		int insidePlanesValidFast(const BoundFrustum &bf) const;
		
		// inside portals
		int insidePortals(const BoundSphere &bs) const;
		int insidePortals(const BoundBox &bb) const;
		int insidePortals(const BoundFrustum &bf) const;
		
		int insidePortalsValid(const BoundSphere &bs) const;
		int insidePortalsValid(const BoundBox &bb) const;
		int insidePortalsValid(const BoundFrustum &bf) const;
		
		// inside shadow
		int insideShadow(const vec3 &direction,const BoundSphere &object) const;
		int insideShadow(const BoundSphere &light,const vec3 &offset,const BoundSphere &object) const;
		
		// parameters
		INLINE int isValid() const { return valid; }
		INLINE const vec3 &getCamera() const { return camera; }
		INLINE const vec4 *getPlanes() const { return planes; }
		INLINE const vec3 *getPoints() const { return points; }
		
	private:
		
		struct Portal;
		
		enum {
			PLANE_L = 0,
			PLANE_R,
			PLANE_B,
			PLANE_T,
			PLANE_N,
			PLANE_F,
		};
		
		int inside_planes(const vec3 &point) const;
		int inside_planes(const vec3 &point,float radius) const;
		int inside_planes(const vec3 &min,const vec3 &max) const;
		int inside_planes(const vec3 *points,int num_points) const;
		
		int inside_planes_fast(const vec3 &point) const;
		int inside_planes_fast(const vec3 &point,float radius) const;
		int inside_planes_fast(const vec3 &min,const vec3 &max) const;
		int inside_planes_fast(const vec3 *points,int num_points) const;
		
		int inside_portal(const Portal &portal,const vec3 &point) const;
		int inside_portal(const Portal &portal,const vec3 &point,float radius) const;
		int inside_portal(const Portal &portal,const vec3 &min,const vec3 &max) const;
		int inside_portal(const Portal &portal,const vec3 *points,int num_points) const;
		
		int inside_portals(const vec3 &point) const;
		int inside_portals(const vec3 &point,float radius) const;
		int inside_portals(const vec3 &min,const vec3 &max) const;
		int inside_portals(const vec3 *points,int num_points) const;
		
		struct Portal {
			vec4 plane;					// portal plane
			vec4 planes[4];				// aos clipping planes
			vec4 tplanes[4];			// soa clipping planes
			vec3 points[4];				// bounding points
		};
		
		int valid;						// valid flag
		vec3 camera;					// camera position
		vec4 planes[6];					// aos clipping planes
		vec4 tplanes[8];				// soa clipping planes
		vec3 points[8];					// bounding points
		
		vec4 garbage[4];				// memory layout
		
		VectorStack<Portal,16> portals;	// portals
};

/*
 */
INLINE int BoundFrustum::inside(const vec3 &point) const {
	if(inside_portals(point) == 0) return 0;
	return inside_planes(point);
}

INLINE int BoundFrustum::inside(const vec3 &point,float radius) const {
	if(inside_portals(point,radius) == 0) return 0;
	return inside_planes(point,radius);
}

INLINE int BoundFrustum::inside(const vec3 &min,const vec3 &max) const {
	if(inside_portals(min,max) == 0) return 0;
	return inside_planes(min,max);
}

INLINE int BoundFrustum::inside(const vec3 *points,int num_points) const {
	if(inside_portals(points,num_points) == 0) return 0;
	return inside_planes(points,num_points);
}

/*
 */
INLINE int BoundFrustum::insideFast(const vec3 &point) const {
	if(inside_portals(point) == 0) return 0;
	return inside_planes_fast(point);
}

INLINE int BoundFrustum::insideFast(const vec3 &point,float radius) const {
	if(inside_portals(point,radius) == 0) return 0;
	return inside_planes_fast(point,radius);
}

INLINE int BoundFrustum::insideFast(const vec3 &min,const vec3 &max) const {
	if(inside_portals(min,max) == 0) return 0;
	return inside_planes_fast(min,max);
}

INLINE int BoundFrustum::insideFast(const vec3 *points,int num_points) const {
	if(inside_portals(points,num_points) == 0) return 0;
	return inside_planes_fast(points,num_points);
}

/*
 */
INLINE int BoundFrustum::insideValid(const BoundSphere &bs) const {
	return inside(bs.getCenter(),bs.getRadius());
}

INLINE int BoundFrustum::insideValid(const BoundBox &bb) const {
	return inside(bb.getMin(),bb.getMax());
}

INLINE int BoundFrustum::insideValid(const BoundFrustum &bf) const {
	return inside(bf.points,8);
}

/*
 */
INLINE int BoundFrustum::insideValidFast(const BoundSphere &bs) const {
	return insideFast(bs.getCenter(),bs.getRadius());
}

INLINE int BoundFrustum::insideValidFast(const BoundBox &bb) const {
	return insideFast(bb.getMin(),bb.getMax());
}

INLINE int BoundFrustum::insideValidFast(const BoundFrustum &bf) const {
	return insideFast(bf.points,8);
}

/*
 */
INLINE int BoundFrustum::insideAllValid(const BoundSphere &bs) const {
	return inside(bs.getCenter(),-bs.getRadius());
}

INLINE int BoundFrustum::insideAllValid(const BoundBox &bb) const {
	const vec3 &min = bb.getMin();
	const vec3 &max = bb.getMax();
	if(inside(vec3(min.x,min.y,min.z)) == 0) return 0;
	if(inside(vec3(max.x,min.y,min.z)) == 0) return 0;
	if(inside(vec3(min.x,max.y,min.z)) == 0) return 0;
	if(inside(vec3(max.x,max.y,min.z)) == 0) return 0;
	if(inside(vec3(min.x,min.y,max.z)) == 0) return 0;
	if(inside(vec3(max.x,min.y,max.z)) == 0) return 0;
	if(inside(vec3(min.x,max.y,max.z)) == 0) return 0;
	if(inside(vec3(max.x,max.y,max.z)) == 0) return 0;
	return 1;
}

INLINE int BoundFrustum::insideAllValid(const BoundFrustum &bf) const {
	for(int i = 0; i < 8; i++) {
		if(inside(bf.points[i]) == 0) return 0;
	}
	return 1;
}

/*
 */
INLINE int BoundFrustum::insideAllValidFast(const BoundSphere &bs) const {
	return insideFast(bs.getCenter(),-bs.getRadius());
}

INLINE int BoundFrustum::insideAllValidFast(const BoundBox &bb) const {
	const vec3 &min = bb.getMin();
	const vec3 &max = bb.getMax();
	if(insideFast(vec3(min.x,min.y,min.z)) == 0) return 0;
	if(insideFast(vec3(max.x,min.y,min.z)) == 0) return 0;
	if(insideFast(vec3(min.x,max.y,min.z)) == 0) return 0;
	if(insideFast(vec3(max.x,max.y,min.z)) == 0) return 0;
	if(insideFast(vec3(min.x,min.y,max.z)) == 0) return 0;
	if(insideFast(vec3(max.x,min.y,max.z)) == 0) return 0;
	if(insideFast(vec3(min.x,max.y,max.z)) == 0) return 0;
	if(insideFast(vec3(max.x,max.y,max.z)) == 0) return 0;
	return 1;
}

INLINE int BoundFrustum::insideAllValidFast(const BoundFrustum &bf) const {
	for(int i = 0; i < 8; i++) {
		if(insideFast(bf.points[i]) == 0) return 0;
	}
	return 1;
}

/*
 */
INLINE int BoundFrustum::insidePlanesValid(const BoundSphere &bs) const {
	return inside_planes(bs.getCenter(),bs.getRadius());
}

INLINE int BoundFrustum::insidePlanesValid(const BoundBox &bb) const {
	return inside_planes(bb.getMin(),bb.getMax());
}

INLINE int BoundFrustum::insidePlanesValid(const BoundFrustum &bf) const {
	return inside_planes(bf.points,8);
}

/*
 */
INLINE int BoundFrustum::insidePlanesValidFast(const BoundSphere &bs) const {
	return inside_planes_fast(bs.getCenter(),bs.getRadius());
}

INLINE int BoundFrustum::insidePlanesValidFast(const BoundBox &bb) const {
	return inside_planes_fast(bb.getMin(),bb.getMax());
}

INLINE int BoundFrustum::insidePlanesValidFast(const BoundFrustum &bf) const {
	return inside_planes_fast(bf.points,8);
}

/*
 */
INLINE int BoundFrustum::insidePortalsValid(const BoundSphere &bs) const {
	return inside_portals(bs.getCenter(),bs.getRadius());
}

INLINE int BoundFrustum::insidePortalsValid(const BoundBox &bb) const {
	return inside_portals(bb.getMin(),bb.getMax());
}

INLINE int BoundFrustum::insidePortalsValid(const BoundFrustum &bf) const {
	return inside_portals(bf.points,8);
}

#endif /* __BOUNDS_H__ */
