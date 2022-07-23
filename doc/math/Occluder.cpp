/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Occluder.cpp
 * Desc:    Software occluder
 * Version: 1.09
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

#include "Utils.h"
#include "Bounds.h"
#include "SimdLib.h"
#include "Occluder.h"

/*
 */
Occluder::Occluder() {
	
	has_occluders = 0;
	need_clear = 1;
	
	num_vertex = 0;
	vertex = NULL;
	
	data = new float[OCCLUDER_WIDTH * OCCLUDER_HEIGHT];
	
	clear();
}

Occluder::Occluder(const mat4 &projection,const mat4 &modelview) {
	
	has_occluders = 0;
	need_clear = 1;
	
	num_vertex = 0;
	vertex = NULL;
	
	data = new float[OCCLUDER_WIDTH * OCCLUDER_HEIGHT];
	
	set(projection,modelview);
}

Occluder::Occluder(const Occluder &occluder) {
	
	num_vertex = 0;
	vertex = NULL;
	
	data = new float[OCCLUDER_WIDTH * OCCLUDER_HEIGHT];
	
	for(int i = 0; i < NUM_PLANES; i++) {
		planes[i] = occluder.planes[i];
		tplanes[i] = occluder.tplanes[i];
	}
	
	modelviewprojection = occluder.modelviewprojection;
	tmodelviewprojection = occluder.tmodelviewprojection;
	
	has_occluders = occluder.has_occluders;
	need_clear = occluder.need_clear;
	
	if(has_occluders) Math::memcpy(data,occluder.data,sizeof(float) * OCCLUDER_WIDTH * OCCLUDER_HEIGHT);
}

Occluder::~Occluder() {
	
	delete [] vertex;
	
	delete [] data;
}

/*
 */
Occluder &Occluder::operator=(const Occluder &occluder) {
	
	if(this == &occluder) return *this;
	
	for(int i = 0; i < NUM_PLANES; i++) {
		planes[i] = occluder.planes[i];
		tplanes[i] = occluder.tplanes[i];
	}
	
	modelviewprojection = occluder.modelviewprojection;
	tmodelviewprojection = occluder.tmodelviewprojection;
	
	has_occluders = occluder.has_occluders;
	need_clear = occluder.need_clear;
	
	if(has_occluders) Math::memcpy(data,occluder.data,sizeof(float) * OCCLUDER_WIDTH * OCCLUDER_HEIGHT);
	
	return *this;
}

/*
 */
void Occluder::swap(Occluder &occluder) {
	
	if(this == &occluder) return;
	
	for(int i = 0; i < NUM_PLANES; i++) {
		::swap(planes[i],occluder.planes[i]);
		::swap(tplanes[i],occluder.tplanes[i]);
	}
	
	::swap(modelviewprojection,occluder.modelviewprojection);
	::swap(tmodelviewprojection,occluder.tmodelviewprojection);
	
	::swap(has_occluders,occluder.has_occluders);
	::swap(need_clear,occluder.need_clear);
	
	::swap(num_vertex,occluder.num_vertex);
	::swap(vertex,occluder.vertex);
	
	::swap(data,occluder.data);
}

/*
 */
void Occluder::clear() {
	
	for(int i = 0; i < NUM_PLANES; i++) {
		planes[i] = vec4_zero;
		tplanes[i] = vec4_zero;
	}
	modelviewprojection = mat4_identity;
	tmodelviewprojection = mat4_identity;
	
	has_occluders = 0;
	need_clear = 1;
}

/*
 */
void Occluder::set(const mat4 &projection,const mat4 &modelview) {
	
	clear();
	
	// clipping planes
	mat4 mvp = projection * modelview;
	planes[PLANE_L].set(mvp.m30 + mvp.m00,mvp.m31 + mvp.m01,mvp.m32 + mvp.m02,mvp.m33 + mvp.m03);
	planes[PLANE_R].set(mvp.m30 - mvp.m00,mvp.m31 - mvp.m01,mvp.m32 - mvp.m02,mvp.m33 - mvp.m03);
	planes[PLANE_T].set(mvp.m30 + mvp.m10,mvp.m31 + mvp.m11,mvp.m32 + mvp.m12,mvp.m33 + mvp.m13);
	planes[PLANE_B].set(mvp.m30 - mvp.m10,mvp.m31 - mvp.m11,mvp.m32 - mvp.m12,mvp.m33 - mvp.m13);
	planes[PLANE_N].set(mvp.m30 + mvp.m20,mvp.m31 + mvp.m21,mvp.m32 + mvp.m22,mvp.m33 + mvp.m23);
	for(int i = 0; i < NUM_PLANES; i++) {
		planes[i] /= length(vec3(planes[i]));
	}
	
	// clipping plane offsets
	float near_clip_plane = 1.0f;
	float far_clip_plane = 1.0f;
	decomposeProjection(projection,near_clip_plane,far_clip_plane);
	planes[PLANE_L].w += near_clip_plane;
	planes[PLANE_R].w += near_clip_plane;
	planes[PLANE_T].w += near_clip_plane;
	planes[PLANE_B].w += near_clip_plane;
	planes[PLANE_N].w -= near_clip_plane;
	
	// transformed clipping planes
	for(int i = 0; i < NUM_PLANES; i++) {
		tplanes[i] = planes[i];
	}
	
	// modelviewprojection matrix
	static const mat4 offset = scale(OCCLUDER_WIDTH / 2.0f,OCCLUDER_HEIGHT / 2.0f,-100.0f) * translate(1.0f,1.0f,-1.0f);
	modelviewprojection = offset * projection * modelview;
	tmodelviewprojection = modelviewprojection;
}

/*
 */
void Occluder::setITransform(const mat4 &itransform) {
	
	// near clipping plane
	for(int i = 0; i < NUM_PLANES; i++) {
		tplanes[i] = planes[i] * itransform;
	}
	
	// modelviewprojection matrix
	tmodelviewprojection = modelviewprojection * itransform;
}

/******************************************************************************\
*
* Render occluder
*
\******************************************************************************/

/*
 */
int Occluder::hasOccluders() const {
	return has_occluders;
}

/*
 */
void Occluder::addOccluder() {
	
	// clear occluder
	if(need_clear) {
		Math::memset(data,0,sizeof(float) * OCCLUDER_WIDTH * OCCLUDER_HEIGHT);
		need_clear = 0;
	}
	
	has_occluders = 1;
}

/*
 */
void Occluder::renderOccluder(const vec3 *vertex,int num,const unsigned short *indices,int num_indices) {
	
	// clear occluder
	if(need_clear) {
		Math::memset(data,0,sizeof(float) * OCCLUDER_WIDTH * OCCLUDER_HEIGHT);
		need_clear = 0;
	}
	
	// render triangles
	for(int i = 0; i < num_indices; i += 3) {
		const vec3 &v0 = vertex[indices[i + 0]];
		const vec3 &v1 = vertex[indices[i + 1]];
		const vec3 &v2 = vertex[indices[i + 2]];
		has_occluders |= clip_triangle(v0,v1,v2,0);
	}
}

/*
 */
int Occluder::clip_triangle(const vec3 &v0,const vec3 &v1,const vec3 &v2,int plane) const {
	
	// render triangle
	if(plane == NUM_PLANES) {
		vec3 p0,p1,p2;
		#ifdef USE_SSE
			__m128 col_0 = tmodelviewprojection.col0;
			__m128 col_1 = tmodelviewprojection.col1;
			__m128 col_2 = tmodelviewprojection.col2;
			__m128 col_3 = tmodelviewprojection.col3;
			#define PROJECT_POINT(RET,V) { \
				__m128 res_0 = _mm_mul_ps(col_0,_MM_SWIZZLE(V.vec,X,X,X,X)); \
				__m128 res_1 = _mm_mul_ps(col_1,_MM_SWIZZLE(V.vec,Y,Y,Y,Y)); \
				__m128 res_2 = _mm_mul_ps(col_2,_MM_SWIZZLE(V.vec,Z,Z,Z,Z)); \
				__m128 res_3 = _mm_add_ps(_mm_add_ps(res_0,res_1),_mm_add_ps(res_2,col_3)); \
				RET.vec = _mm_div_ps(res_3,_MM_SWIZZLE(res_3,W,W,W,W)); \
			}
			PROJECT_POINT(p0,v0);
			PROJECT_POINT(p1,v1);
			PROJECT_POINT(p2,v2);
			#undef PROJECT_POINT
		#elif USE_ALTIVEC
			vec_float4 col_0 = tmodelviewprojection.col0;
			vec_float4 col_1 = tmodelviewprojection.col1;
			vec_float4 col_2 = tmodelviewprojection.col2;
			vec_float4 col_3 = tmodelviewprojection.col3;
			#define PROJECT_POINT(RET,V) { \
				vec_float4 res_0 = vec_madd(col_0,VEC_SWIZZLE(V.vec,X,X,X,X),col_3); \
				vec_float4 res_1 = vec_madd(col_1,VEC_SWIZZLE(V.vec,Y,Y,Y,Y),res_0); \
				vec_float4 res_2 = vec_madd(col_2,VEC_SWIZZLE(V.vec,Z,Z,Z,Z),res_1); \
				RET.vec = vec_madd(res_2,vec_rcp_nr(VEC_SWIZZLE(res_2,W,W,W,W)),vec_splats(0.0f)); \
			}
			PROJECT_POINT(p0,v0);
			PROJECT_POINT(p1,v1);
			PROJECT_POINT(p2,v2);
			#undef PROJECT_POINT
		#elif USE_NEON
			float32x4_t col_0 = tmodelviewprojection.col0;
			float32x4_t col_1 = tmodelviewprojection.col1;
			float32x4_t col_2 = tmodelviewprojection.col2;
			float32x4_t col_3 = tmodelviewprojection.col3;
			#define PROJECT_POINT(RET,V) { \
				float32x2_t low = vget_low_f32(V.vec); \
				float32x2_t high = vget_high_f32(V.vec); \
				float32x4_t res_0 = vmlaq_lane_f32(col_3,col_0,low,0); \
				float32x4_t res_1 = vmlaq_lane_f32(res_0,col_1,low,1); \
				float32x4_t res_2 = vmlaq_lane_f32(res_1,col_2,high,0); \
				RET.vec = vmulq_lane_f32(res_2,vrcp_nr_f32(vget_high_f32(res_2)),1); \
			}
			PROJECT_POINT(p0,v0);
			PROJECT_POINT(p1,v1);
			PROJECT_POINT(p2,v2);
			#undef PROJECT_POINT
		#else
			proj(p0,tmodelviewprojection,v0);
			proj(p1,tmodelviewprojection,v1);
			proj(p2,tmodelviewprojection,v2);
		#endif
		return render_triangle(&p0,&p1,&p2);
	}
	
	// vertex distances
	const vec4 &tplane = tplanes[plane++];
	float d0 = dot(tplane,v0);
	float d1 = dot(tplane,v1);
	float d2 = dot(tplane,v2);
	
	// inside flags
	int i0 = (d0 > 0.0f);
	int i1 = (d1 > 0.0f);
	int i2 = (d2 > 0.0f);
	int num = i0 + i1 + i2;
	
	// three vertices are visible
	if(num == 3) {
		return clip_triangle(v0,v1,v2,plane);
	}
	
	// two vertices are visible
	if(num == 2) {
		
		if(i0 == 0) {
			vec3 v01,v02;
			mad(v01,v0 - v1,d0 / (d1 - d0),v0);
			mad(v02,v0 - v2,d0 / (d2 - d0),v0);
			return clip_triangle(v01,v1,v2,plane) | clip_triangle(v2,v02,v01,plane);
		}
		else if(i1 == 0) {
			vec3 v12,v10;
			mad(v12,v1 - v2,d1 / (d2 - d1),v1);
			mad(v10,v1 - v0,d1 / (d0 - d1),v1);
			return clip_triangle(v12,v2,v0,plane) | clip_triangle(v0,v10,v12,plane);
		}
		else {
			vec3 v20,v21;
			mad(v20,v2 - v0,d2 / (d0 - d2),v2);
			mad(v21,v2 - v1,d2 / (d1 - d2),v2);
			return clip_triangle(v20,v0,v1,plane) | clip_triangle(v1,v21,v20,plane);
		}
	}
	
	// one vertex is visible
	if(num == 1) {
		
		if(i0 != 0) {
			vec3 v01,v02;
			mad(v01,v0 - v1,d0 / (d1 - d0),v0);
			mad(v02,v0 - v2,d0 / (d2 - d0),v0);
			return clip_triangle(v0,v01,v02,plane);
		}
		else if(i1 != 0) {
			vec3 v12,v10;
			mad(v12,v1 - v2,d1 / (d2 - d1),v1);
			mad(v10,v1 - v0,d1 / (d0 - d1),v1);
			return clip_triangle(v1,v12,v10,plane);
		}
		else {
			vec3 v20,v21;
			mad(v20,v2 - v0,d2 / (d0 - d2),v2);
			mad(v21,v2 - v1,d2 / (d1 - d2),v2);
			return clip_triangle(v2,v21,v20,plane);
		}
	}
	
	return 0;
}

/*
 */
int Occluder::render_triangle(const vec3 *v0,const vec3 *v1,const vec3 *v2) const {
	
	if(v0->y > v1->y) ::swap(v0,v1);
	if(v0->y > v2->y) ::swap(v0,v2);
	if(v1->y > v2->y) ::swap(v1,v2);
	
	int v0_x = Math::ftoi(v0->x);
	int v1_x = Math::ftoi(v1->x);
	int v2_x = Math::ftoi(v2->x);
	int v0_y = Math::ftoi(v0->y);
	int v1_y = Math::ftoi(v1->y);
	int v2_y = Math::ftoi(v2->y);
	
	int v10_x = v1_x - v0_x;
	int v20_x = v2_x - v0_x;
	int v21_x = v2_x - v1_x;
	int v10_y = v1_y - v0_y;
	int v20_y = v2_y - v0_y;
	int v21_y = v2_y - v1_y;
	float v10_z = v1->z - v0->z;
	float v20_z = v2->z - v0->z;
	float iv10_y = Math::rcp(Math::itof(v10_y));
	float iv20_y = Math::rcp(Math::itof(v20_y));
	float iv21_y = Math::rcp(Math::itof(v21_y));
	
	float det = Math::itof(v20_x * v10_y - v20_y * v10_x);
	if(Math::abs(det) < EPSILON) return 0;
	float idet = Math::rcp(det);
	
	int y0 = clamp(v0_y,0,OCCLUDER_HEIGHT);
	int y1 = clamp(v1_y,0,OCCLUDER_HEIGHT);
	
	int ret = 0;
	
	for(int y = y0; y < y1; y++) {
		
		int yv0 = y - v0_y;
		int yv0_v20_x = yv0 * v20_x;
		int yv0_v10_x = yv0 * v10_x;
		int x0 = clamp(Math::ftoi(Math::ceil(v0_x + v20_x * yv0 * iv20_y)),0,OCCLUDER_WIDTH);
		int x1 = clamp(Math::ftoi(Math::ceil(v0_x + v10_x * yv0 * iv10_y)),0,OCCLUDER_WIDTH);
		
		if(x0 > x1) ::swap(x0,x1);
		
		for(int x = x0; x < x1; x++) {
			
			int xv0 = x - v0_x;
			float a = Math::itof(yv0_v20_x - xv0 * v20_y) * idet;
			float b = Math::itof(xv0 * v10_y - yv0_v10_x) * idet;
			float z = v0->z + v10_z * a + v20_z * b;
			
			int id = (y << OCCLUDER_SHIFT) + x;
			if(data[id] >= z) continue;
			
			data[id] = z;
			ret = 1;
		}
	}
	
	y0 = clamp(v1_y,0,OCCLUDER_HEIGHT);
	y1 = clamp(v2_y,0,OCCLUDER_HEIGHT);
	
	for(int y = y0; y < y1; y++) {
		
		int yv0 = y - v0_y;
		int yv1 = y - v1_y;
		int yv0_v20_x = yv0 * v20_x;
		int yv0_v10_x = yv0 * v10_x;
		int x0 = clamp(Math::ftoi(Math::ceil(v0_x + v20_x * yv0 * iv20_y)),0,OCCLUDER_WIDTH);
		int x1 = clamp(Math::ftoi(Math::ceil(v1_x + v21_x * yv1 * iv21_y)),0,OCCLUDER_WIDTH);
		
		if(x0 > x1) ::swap(x0,x1);
		
		for(int x = x0; x < x1; x++) {
			
			int xv0 = x - v0_x;
			float a = Math::itof(yv0_v20_x - xv0 * v20_y) * idet;
			float b = Math::itof(xv0 * v10_y - yv0_v10_x) * idet;
			float z = v0->z + v10_z * a + v20_z * b;
			
			int id = (y << OCCLUDER_SHIFT) + x;
			if(data[id] >= z) continue;
			
			data[id] = z;
			ret = 1;
		}
	}
	
	return ret;
}

/******************************************************************************\
*
* Inside occluder
*
\******************************************************************************/

/*
 */
int Occluder::inside(const vec3 &min,const vec3 &max) const {
	
	if(has_occluders == 0) return 1;
	
	float z;
	int x0,y0,x1,y1;
	if(project_bound_box(min,max,x0,y0,x1,y1,z) == 0) return 1;
	
	int num = x1 - x0 + 1;
	
	if(num & ~3) {
		for(int y = y0; y <= y1; y++) {
			const float *d = data + (y << OCCLUDER_SHIFT) + x0;
			int num = x1 - x0 + 1;
			for(int i = (num >> 2) - 1; i >= 0; i--) {
				if(d[0] < z) return 1;
				if(d[1] < z) return 1;
				if(d[2] < z) return 1;
				if(d[3] < z) return 1;
				d += 4;
			}
			for(int i = (num & 3) - 1; i >= 0; i--) {
				if(*d++ < z) return 1;
			}
		}
	}
	else {
		for(int y = y0; y <= y1; y++) {
			const float *d = data + (y << OCCLUDER_SHIFT) + x0;
			for(int i = num - 1; i >= 0; i--) {
				if(*d++ < z) return 1;
			}
		}
	}
	
	return 0;
}

/*
 */
int Occluder::insideAll(const vec3 &min,const vec3 &max) const {
	
	if(has_occluders == 0) return 1;
	
	float z;
	int x0,y0,x1,y1;
	if(project_bound_box(min,max,x0,y0,x1,y1,z) == 0) return 1;
	
	int num = x1 - x0 + 1;
	
	if(num & ~3) {
		for(int y = y0; y <= y1; y++) {
			const float *d = data + (y << OCCLUDER_SHIFT) + x0;
			int num = x1 - x0 + 1;
			for(int i = (num >> 2) - 1; i >= 0; i--) {
				if(d[0] > z) return 0;
				if(d[1] > z) return 0;
				if(d[2] > z) return 0;
				if(d[3] > z) return 0;
				d += 4;
			}
			for(int i = (num & 3) - 1; i >= 0; i--) {
				if(*d++ > z) return 0;
			}
		}
	}
	else {
		for(int y = y0; y <= y1; y++) {
			const float *d = data + (y << OCCLUDER_SHIFT) + x0;
			for(int i = num - 1; i >= 0; i--) {
				if(*d++ > z) return 0;
			}
		}
	}
	
	return 1;
}

/*
 */
int Occluder::inside(const BoundSphere &bs) const {
	if(bs.isValid()) {
		vec3 radius = vec3(bs.getRadius());
		return inside(bs.getCenter() - radius,bs.getCenter() + radius);
	}
	return 0;
}

int Occluder::inside(const BoundBox &bb) const {
	if(bb.isValid()) {
		return inside(bb.getMin(),bb.getMax());
	}
	return 0;
}

/*
 */
int Occluder::insideValid(const BoundSphere &bs) const {
	vec3 radius = vec3(bs.getRadius());
	return inside(bs.getCenter() - radius,bs.getCenter() + radius);
}

int Occluder::insideValid(const BoundBox &bb) const {
	return inside(bb.getMin(),bb.getMax());
}

/*
 */
int Occluder::insideAll(const BoundSphere &bs) const {
	if(bs.isValid()) {
		vec3 radius = vec3(bs.getRadius());
		return insideAll(bs.getCenter() - radius,bs.getCenter() + radius);
	}
	return 0;
}

int Occluder::insideAll(const BoundBox &bb) const {
	if(bb.isValid()) {
		return insideAll(bb.getMin(),bb.getMax());
	}
	return 0;
}

/*
 */
int Occluder::insideAllValid(const BoundSphere &bs) const {
	vec3 radius = vec3(bs.getRadius());
	return insideAll(bs.getCenter() - radius,bs.getCenter() + radius);
}

int Occluder::insideAllValid(const BoundBox &bb) const {
	return insideAll(bb.getMin(),bb.getMax());
}

/*
 */
INLINE int Occluder::project_bound_box(const vec3 &min,const vec3 &max,int &x0,int &y0,int &x1,int &y1,float &z) const {
	
	vec4 points_min = vec4_infinity;
	vec4 points_max = -vec4_infinity;
	
	#ifdef USE_SSE
		__m128 col_0 = tmodelviewprojection.col0;
		__m128 col_1 = tmodelviewprojection.col1;
		__m128 col_2 = tmodelviewprojection.col2;
		__m128 col_3 = tmodelviewprojection.col3;
		__m128 col_0_min = _mm_mul_ps(col_0,_MM_SWIZZLE(min.vec,X,X,X,X));
		__m128 col_1_min = _mm_mul_ps(col_1,_MM_SWIZZLE(min.vec,Y,Y,Y,Y));
		__m128 col_2_min = _mm_mul_ps(col_2,_MM_SWIZZLE(min.vec,Z,Z,Z,Z));
		__m128 col_0_max = _mm_mul_ps(col_0,_MM_SWIZZLE(max.vec,X,X,X,X));
		__m128 col_1_max = _mm_mul_ps(col_1,_MM_SWIZZLE(max.vec,Y,Y,Y,Y));
		__m128 col_2_max = _mm_mul_ps(col_2,_MM_SWIZZLE(max.vec,Z,Z,Z,Z));
		#define PROJECT_POINT(COL_0,COL_1,COL_2) { \
			__m128 point = _mm_add_ps(_mm_add_ps((COL_0),(COL_1)),_mm_add_ps((COL_2),col_3)); \
			if(_mm_movemask_ps(point) & 0x08) return 0; \
			point = _mm_div_ps(point,_MM_SWIZZLE(point,W,W,W,W)); \
			points_min.vec = _mm_min_ps(points_min.vec,point); \
			points_max.vec = _mm_max_ps(points_max.vec,point); \
		}
		PROJECT_POINT(col_0_min,col_1_min,col_2_min)
		PROJECT_POINT(col_0_min,col_1_min,col_2_max)
		PROJECT_POINT(col_0_min,col_1_max,col_2_min)
		PROJECT_POINT(col_0_min,col_1_max,col_2_max)
		PROJECT_POINT(col_0_max,col_1_min,col_2_min)
		PROJECT_POINT(col_0_max,col_1_min,col_2_max)
		PROJECT_POINT(col_0_max,col_1_max,col_2_min)
		PROJECT_POINT(col_0_max,col_1_max,col_2_max)
		#undef PROJECT_POINT
	#elif USE_NEON
		float32x4_t col_0 = tmodelviewprojection.col0;
		float32x4_t col_1 = tmodelviewprojection.col1;
		float32x4_t col_2 = tmodelviewprojection.col2;
		float32x4_t col_3 = tmodelviewprojection.col3;
		float32x4_t col_0_min = vmulq_n_f32(col_0,min.x);
		float32x4_t col_1_min = vmulq_n_f32(col_1,min.y);
		float32x4_t col_2_min = vmulq_n_f32(col_2,min.z);
		float32x4_t col_0_max = vmulq_n_f32(col_0,max.x);
		float32x4_t col_1_max = vmulq_n_f32(col_1,max.y);
		float32x4_t col_2_max = vmulq_n_f32(col_2,max.z);
		#define PROJECT_POINT(COL_0,COL_1,COL_2) { \
			float32x4_t point = vaddq_f32(vaddq_f32((COL_0),(COL_1)),vaddq_f32((COL_2),col_3)); \
			if(vmaskq_u32(vcltq_f32(point,vdupq_n_f32(0.0f))) & 0x08) return 0; \
			point = vmulq_lane_f32(point,vrcp_nr_f32(vget_high_f32(point)),1); \
			points_min.vec = vminq_f32(points_min.vec,point); \
			points_max.vec = vmaxq_f32(points_max.vec,point); \
		}
		PROJECT_POINT(col_0_min,col_1_min,col_2_min)
		PROJECT_POINT(col_0_min,col_1_min,col_2_max)
		PROJECT_POINT(col_0_min,col_1_max,col_2_min)
		PROJECT_POINT(col_0_min,col_1_max,col_2_max)
		PROJECT_POINT(col_0_max,col_1_min,col_2_min)
		PROJECT_POINT(col_0_max,col_1_min,col_2_max)
		PROJECT_POINT(col_0_max,col_1_max,col_2_min)
		PROJECT_POINT(col_0_max,col_1_max,col_2_max)
		#undef PROJECT_POINT
	#else
		vec4 col_0 = tmodelviewprojection.getColumn(0);
		vec4 col_1 = tmodelviewprojection.getColumn(1);
		vec4 col_2 = tmodelviewprojection.getColumn(2);
		vec4 col_3 = tmodelviewprojection.getColumn(3);
		vec4 col_0_min = col_0 * min.x;
		vec4 col_1_min = col_1 * min.y;
		vec4 col_2_min = col_2 * min.z;
		vec4 col_0_max = col_0 * max.x;
		vec4 col_1_max = col_1 * max.y;
		vec4 col_2_max = col_2 * max.z;
		#define PROJECT_POINT(COL_0,COL_1,COL_2) { \
			float w = (COL_0).w + (COL_1).w + (COL_2).w + col_3.w; \
			if(w < 0.0f) return 0; \
			float iw = Math::rcp(w); \
			float x = ((COL_0).x + (COL_1).x + (COL_2).x + col_3.x) * iw; \
			float y = ((COL_0).y + (COL_1).y + (COL_2).y + col_3.y) * iw; \
			float z = ((COL_0).z + (COL_1).z + (COL_2).z + col_3.z) * iw; \
			if(points_min.x > x) points_min.x = x; \
			if(points_min.y > y) points_min.y = y; \
			if(points_max.x < x) points_max.x = x; \
			if(points_max.y < y) points_max.y = y; \
			if(points_max.z < z) points_max.z = z; \
		}
		PROJECT_POINT(col_0_min,col_1_min,col_2_min)
		PROJECT_POINT(col_0_min,col_1_min,col_2_max)
		PROJECT_POINT(col_0_min,col_1_max,col_2_min)
		PROJECT_POINT(col_0_min,col_1_max,col_2_max)
		PROJECT_POINT(col_0_max,col_1_min,col_2_min)
		PROJECT_POINT(col_0_max,col_1_min,col_2_max)
		PROJECT_POINT(col_0_max,col_1_max,col_2_min)
		PROJECT_POINT(col_0_max,col_1_max,col_2_max)
		#undef PROJECT_POINT
	#endif
	
	x0 = clamp(Math::ftoi(Math::floor(points_min.x)),0,OCCLUDER_WIDTH - 1);
	y0 = clamp(Math::ftoi(Math::floor(points_min.y)),0,OCCLUDER_HEIGHT - 1);
	x1 = clamp(Math::ftoi(Math::ceil(points_max.x)),0,OCCLUDER_WIDTH - 1);
	y1 = clamp(Math::ftoi(Math::ceil(points_max.y)),0,OCCLUDER_HEIGHT - 1);
	
	z = points_max.z;
	
	return 1;
}
