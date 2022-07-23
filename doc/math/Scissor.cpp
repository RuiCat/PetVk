/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Scissor.cpp
 * Desc:    Scissor
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

#include "Log.h"
#include "Scissor.h"

/*
 */
Scissor::Scissor() {
	
	clear();
}

Scissor::~Scissor() {
	
}

/*
 */
void Scissor::clear() {
	
	left = 0.0f;
	right = 1.0f;
	bottom = 0.0f;
	top = 1.0f;
	znear = 0.0f;
	zfar = 1.0f;
}

void Scissor::set(float l,float r,float b,float t) {
	left = l;
	right = r;
	bottom = b;
	top = t;
}

void Scissor::set(const mat4 &p,const mat4 &m) {
	
	clear();
	
	// matrices
	projection = p;
	modelview = m;
	mul(modelviewprojection,projection,modelview);
	
	// projection parameters
	decomposeProjection(projection,near_clip_plane,far_clip_plane);
}

/*
 */
int Scissor::set(const vec3 &min,const vec3 &max) {
	
	// check visibility
	if(min.x > 1.0f || max.x < -1.0f) return 0;
	if(min.y > 1.0f || max.y < -1.0f) return 0;
	if(min.z > 1.0f || max.z < -1.0f) return 0;
	
	// scale and bias
	float l = min.x * 0.5f + 0.5f;
	float r = max.x * 0.5f + 0.5f;
	float b = min.y * 0.5f + 0.5f;
	float t = max.y * 0.5f + 0.5f;
	float n = min.z * 0.5f + 0.5f;
	float f = max.z * 0.5f + 0.5f;
	
	// clamp bounds
	if(l < 0.0f) l = 0.0f;
	if(r > 1.0f) r = 1.0f;
	if(b < 0.0f) b = 0.0f;
	if(t > 1.0f) t = 1.0f;
	if(n < 0.0f) n = 0.0f;
	if(f > 1.0f) f = 1.0f;
	
	// check scissors intersections
	if(l > right || r < left) return 0;
	if(b > top || t < bottom) return 0;
	if(n > zfar || f < znear) return 0;
	
	// intersect scissors
	if(left < l) left = l;
	if(right > r) right = r;
	if(bottom < b) bottom = b;
	if(top > t) top = t;
	if(znear < n) znear = n;
	if(zfar > f) zfar = f;
	
	return 1;
}

/*
 */
float Scissor::getX() const {
	return left;
}

float Scissor::getY() const {
	return bottom;
}

float Scissor::getWidth() const {
	return right - left;
}

float Scissor::getHeight() const {
	return top - bottom;
}

float Scissor::getZNear() const {
	return znear;
}

float Scissor::getZFar() const {
	return zfar;
}

/******************************************************************************\
*
* Portals
*
\******************************************************************************/

/*
 */
int Scissor::addPortal(const vec3 *points,int num_points,const mat4 &transform) {
	
	assert(num_points > 2 && "Scissor::addPortal(): bad number of points");
	
	// find bounds
	vec2 min = vec2_infinity;
	vec2 max = -vec2_infinity;
	mat4 m = modelviewprojection * transform;
	for(int i = 0; i < num_points; i++) {
		vec4 p = m * vec4(points[i],1.0f);
		min = ::min(min,vec2(p) / Math::abs(p.w));
		max = ::max(max,vec2(p) / Math::abs(p.w));
	}
	
	return set(vec3(min.x,min.y,-1.0f),vec3(max.x,max.y,1.0f));
}

/******************************************************************************\
*
* Lights
*
\******************************************************************************/

/*
 */
static void ellipsoid_projection(const vec3 *basis,const vec3 &radius,const vec3 &normal,vec3 &axis_a,vec3 &axis_b) {
	
	// projected ellipse basis
	float height[3];
	vec3 projection[3];
	vec3 radius2 = radius * radius;
	static const int next_1[3] = { 1, 2, 0 };
	static const int next_2[3] = { 2, 0, 1 };
	for(int i = 0; i < 3; i++) {
		int id_0 = i;
		int id_1 = next_1[i];
		int id_2 = next_2[i];
		vec3 dir = cross(basis[id_0],normal);
		if(dot(dir,dir) < EPSILON) {
			axis_a = basis[id_1] * radius[id_1];
			axis_b = basis[id_2] * radius[id_2];
			return;
		}
		projection[id_0] = normalize(cross(dir,normal));
		dir = normalize(cross(dir,basis[id_0]));
		float d1 = dot(dir,basis[id_1]);
		float d2 = dot(dir,basis[id_2]);
		float h2 = radius2[id_1] * radius2[id_2] / (radius2[id_2] * d1 * d1 + radius2[id_1] * d2 * d2);
		float c = dot(projection[id_0],dir);
		float c2 = c * c;
		height[id_0] = Math::sqrt(radius2[id_0] * (1.0f - c2) + h2 * c2);
	}
	
	if(Math::abs(dot(projection[0],projection[1])) > 1.0f - EPSILON) {
		axis_a = projection[0] * height[0];
		axis_b = projection[2] * height[2];
		return;
	}
	if(Math::abs(dot(projection[1],projection[2])) > 1.0f - EPSILON) {
		axis_a = projection[0] * height[0];
		axis_b = projection[1] * height[1];
		return;
	}
	if(Math::abs(dot(projection[0],projection[2])) > 1.0f - EPSILON) {
		axis_a = projection[1] * height[1];
		axis_b = projection[2] * height[2];
		return;
	}
	
	// 2d projected basis
	vec3 axis_x = projection[0];
	vec3 axis_y = normalize(cross(axis_x,normal));
	
	// ellipse equation solution
	mat3 matrix;
	for(int i = 0; i < 3; i++) {
		float x = dot(axis_x,projection[i]) * height[i];
		float y = dot(axis_y,projection[i]) * height[i];
		matrix[i + 0] = x * x;
		matrix[i + 4] = x * y;
		matrix[i + 8] = y * y;
	}
	vec3 root = inverse(matrix) * vec3_one;
	
	// orthonormalize basis
	mat3 rotate;
	float angle = (Math::abs(root.x - root.z) > EPSILON) ? Math::atan2(root.y,root.x - root.z) * RAD2DEG / 2.0f : 45.0f;
	if(normal.x == 1.0f) rotate.setRotateX(-angle);
	else if(normal.y == 1.0f) rotate.setRotateY(-angle);
	else rotate.setRotate(normal,-angle);
	axis_a = rotate * axis_x;
	axis_b = rotate * axis_y;
	axis_a = vec3(dot(axis_x,axis_a),dot(axis_y,axis_a),0.0f);
	axis_b = vec3(dot(axis_x,axis_b),dot(axis_y,axis_b),0.0f);
	axis_a *= Math::rsqrtFast(Math::abs(axis_a.x * axis_a.x * root.x + axis_a.x * axis_a.y * root.y + axis_a.y * axis_a.y * root.z));
	axis_b *= Math::rsqrtFast(Math::abs(axis_b.x * axis_b.x * root.x + axis_b.x * axis_b.y * root.y + axis_b.y * axis_b.y * root.z));
	axis_a = axis_x * axis_a.x + axis_y * axis_a.y;
	axis_b = axis_x * axis_b.x + axis_y * axis_b.y;
}

/*
 */
static INLINE float ellipse_left(const vec3 &axis_a,float dot_a,float dot_b,const vec3 &position,float x) {
	float z = Math::sqrt(1.0f - x * x);
	float c = axis_a.x * x + axis_a.z * z;
	float c2 = c * c;
	return position.x * x + position.z * z - Math::sqrt(dot_a * c2 + dot_b * (1.0f - c2));
}

static INLINE float ellipse_right(const vec3 &axis_a,float dot_a,float dot_b,const vec3 &position,float x) {
	float z = -Math::sqrt(1.0f - x * x);
	float c = axis_a.x * x + axis_a.z * z;
	float c2 = c * c;
	return position.x * x + position.z * z - Math::sqrt(dot_a * c2 + dot_b * (1.0f - c2));
}

static INLINE float ellipse_bottom(const vec3 &axis_a,float dot_a,float dot_b,const vec3 &position,float y) {
	float z = Math::sqrt(1.0f - y * y);
	float c = axis_a.y * y + axis_a.z * z;
	float c2 = c * c;
	return position.y * y + position.z * z - Math::sqrt(dot_a * c2 + dot_b * (1.0f - c2));
}

static INLINE float ellipse_top(const vec3 &axis_a,float dot_a,float dot_b,const vec3 &position,float y) {
	float z = -Math::sqrt(1.0f - y * y);
	float c = axis_a.y * y + axis_a.z * z;
	float c2 = c * c;
	return position.y * y + position.z * z - Math::sqrt(dot_a * c2 + dot_b * (1.0f - c2));
}

template <class Func>
static int ellipse_solver(const vec3 &axis_a,const vec3 &axis_b,const vec3 &position,Func func,float *roots) {
	
	int ret = 0;
	float epsilon = 1e-5f;
	float dot_a = dot(axis_a,axis_a);
	float dot_b = dot(axis_b,axis_b);
	vec3 axis = axis_a * Math::rsqrt(dot_a);
	
	// find intersection
	float x0 = -1.0f;
	float y0 = func(axis,dot_a,dot_b,position,x0);
	int y0_pos = (y0 > 0.0f);
	
	for(int i = 0; i < 10; i++) {
		
		float x1 = -1.0f + (i + 1) * 0.2f;
		float y1 = func(axis,dot_a,dot_b,position,x1);
		int y1_pos = (y1 > 0.0f);
		
		// more accuracy
		if(y0_pos != y1_pos) {
			
			float root = 0.0f;
			for(int i = 0; i < 10; i++) {
				
				float x_diff = x1 - x0;
				if(y1_pos) root = -x_diff * y0 / (y1 - y0) + x0;
				else root = x_diff * y0 / (y0 - y1) + x0;
				
				if(root - x0 < epsilon) break;
				if(x1 - root < epsilon) break;
				
				float y = func(axis,dot_a,dot_b,position,root);
				int y_pos = (y > 0.0f);
				
				if(y0_pos == y_pos) {
					x0 = root;
					y0 = y;
					y0_pos = y_pos;
				} else {
					x1 = root;
					y1 = y;
					y1_pos = y_pos;
				}
			}
			
			roots[ret++] = root;
		}
		
		x0 = x1;
		y0 = y1;
		y0_pos = y1_pos;
	}
	
	return ret;
}

/*
 */
int Scissor::addLight(float radius,const mat4 &transform) {
	
	// light bounds
	vec3 min = vec3(-1.0f,-1.0f,INFINITY);
	vec3 max = vec3(1.0f,1.0f,-INFINITY);
	
	// light transformation
	mat4 matrix = modelview * transform;
	vec3 position = matrix.getColumn3(3);
	
	// depth bounds
	vec4 depth[2];
	depth[0] = projection * vec4(position - vec3(0.0f,0.0f,radius),1.0f);
	depth[1] = projection * vec4(position + vec3(0.0f,0.0f,radius),1.0f);
	for(int i = 0; i < 2; i++) {
		depth[i].z /= ::max(Math::abs(depth[i].w),EPSILON);
		min.z = ::min(min.z,depth[i].z);
		max.z = ::max(max.z,depth[i].z);
	}
	if(min.z > 1.0f || max.z < -1.0f) return 0;
	
	/////////////////////////////////
	// project light on xz plane
	/////////////////////////////////
	
	// circle
	float a = position.x * position.x + position.z * position.z;
	float b = -radius * position.x;
	float c = radius * radius - position.z * position.z;
	float d = b * b - a * c;
	
	if(d > 0.0f && Math::abs(position.z) > radius * 1e-4f) {
		
		vec3 n0,n1;
		n0.x = (-b + Math::sqrt(d)) / a;
		n1.x = (-b - Math::sqrt(d)) / a;
		n0.z = (radius - n0.x * position.x) / position.z;
		n1.z = (radius - n1.x * position.x) / position.z;
		
		vec4 p0 = projection * vec4(n0.z * near_clip_plane / n0.x,0.0f,-near_clip_plane,1.0f);
		vec4 p1 = projection * vec4(n1.z * near_clip_plane / n1.x,0.0f,-near_clip_plane,1.0f);
		p0.x /= Math::abs(p0.w);
		p1.x /= Math::abs(p1.w);
		
		if(n0.x > 0.0f) min.x = ::max(min.x,p0.x);
		if(n0.x < 0.0f) max.x = ::min(max.x,p0.x);
		if(n1.x > 0.0f) min.x = ::max(min.x,p1.x);
		if(n1.x < 0.0f) max.x = ::min(max.x,p1.x);
	}
	
	/////////////////////////////////
	// project light on yz plane
	/////////////////////////////////
	
	// circle
	a = position.y * position.y + position.z * position.z;
	b = -radius * position.y;
	c = radius * radius - position.z * position.z;
	d = b * b - a * c;
	
	if(d > 0.0f && Math::abs(position.z) > radius * 1e-4f) {
		
		vec3 n0,n1;
		n0.y = (-b + Math::sqrt(d)) / a;
		n1.y = (-b - Math::sqrt(d)) / a;
		n0.z = (radius - n0.y * position.y) / position.z;
		n1.z = (radius - n1.y * position.y) / position.z;
		
		vec4 p0 = projection * vec4(0.0f,n0.z * near_clip_plane / n0.y,-near_clip_plane,1.0f);
		vec4 p1 = projection * vec4(0.0f,n1.z * near_clip_plane / n1.y,-near_clip_plane,1.0f);
		p0.y /= Math::abs(p0.w);
		p1.y /= Math::abs(p1.w);
		
		if(n0.y > 0.0f) min.y = ::max(min.y,p0.y);
		if(n0.y < 0.0f) max.y = ::min(max.y,p0.y);
		if(n1.y > 0.0f) min.y = ::max(min.y,p1.y);
		if(n1.y < 0.0f) max.y = ::min(max.y,p1.y);
	}
	
	return set(min,max);
}

/*
 */
int Scissor::addLight(const vec3 &radius,const mat4 &transform) {
	
	// light bounds
	vec3 min = vec3(-1.0f,-1.0f,INFINITY);
	vec3 max = vec3(1.0f,1.0f,-INFINITY);
	
	// light type
	int is_sphere = (compare(radius.x,radius.y) && compare(radius.x,radius.z));
	
	// light transformation
	mat4 matrix = modelview * transform;
	vec3 position = matrix.getColumn3(3);
	
	// light basis
	vec3 basis[3];
	if(is_sphere == 0) {
		basis[0] = normalize(matrix.getColumn3(0));
		basis[1] = normalize(matrix.getColumn3(1));
		basis[2] = normalize(matrix.getColumn3(2));
	}
	
	/////////////////////////////////
	// project light on xz plane
	/////////////////////////////////
	
	vec3 axis_a;
	vec3 axis_b;
	float dot_a = 0.0f;
	float dot_b = 0.0f;
	float height = 0.0f;
	
	// projection
	if(is_sphere) {
		height = radius.y;
	} else {
		ellipsoid_projection(basis,radius,vec3(0.0f,1.0f,0.0f),axis_a,axis_b);
		dot_a = dot(axis_a,axis_a);
		dot_b = dot(axis_b,axis_b);
		height = Math::sqrt(axis_a.z * axis_a.z + dot_b * (1.0f - axis_a.z * axis_a.z / dot_a));
	}
	
	// depth bounds
	vec4 depth[4];
	depth[0] = projection * vec4(position - vec3(0.0f,0.0f,height),1.0f);
	depth[1] = projection * vec4(position + vec3(0.0f,0.0f,height),1.0f);
	
	// circle
	if(is_sphere || compare(dot_a,dot_b,1e-2f)) {
		
		float r = (is_sphere) ? radius.y : Math::sqrt(::max(dot_a,dot_b));
		float a = position.x * position.x + position.z * position.z;
		float b = -r * position.x;
		float c = r * r - position.z * position.z;
		float d = b * b - a * c;
		
		if(d > 0.0f && Math::abs(position.z) > r * 1e-4f) {
			
			vec3 n0,n1;
			n0.x = (-b + Math::sqrt(d)) / a;
			n1.x = (-b - Math::sqrt(d)) / a;
			n0.z = (r - n0.x * position.x) / position.z;
			n1.z = (r - n1.x * position.x) / position.z;
			
			vec4 p0 = projection * vec4(n0.z * near_clip_plane / n0.x,0.0f,-near_clip_plane,1.0f);
			vec4 p1 = projection * vec4(n1.z * near_clip_plane / n1.x,0.0f,-near_clip_plane,1.0f);
			p0.x /= Math::abs(p0.w);
			p1.x /= Math::abs(p1.w);
			
			if(n0.x > 0.0f) min.x = ::max(min.x,p0.x);
			if(n0.x < 0.0f) max.x = ::min(max.x,p0.x);
			if(n1.x > 0.0f) min.x = ::max(min.x,p1.x);
			if(n1.x < 0.0f) max.x = ::min(max.x,p1.x);
		}
	}
	// ellipse
	else {
		
		vec3 normals[4];
		int num_normals = 0;
		
		// left plane
		float roots[4];
		int num_roots = ellipse_solver(axis_a,axis_b,position,ellipse_left,roots);
		for(int i = 0; i < num_roots; i++) {
			normals[num_normals++] = vec3(roots[i],0.0f,Math::sqrt(1.0f - roots[i] * roots[i]));
		}
		
		// right plane
		num_roots = ellipse_solver(axis_a,axis_b,position,ellipse_right,roots);
		for(int i = 0; i < num_roots; i++) {
			normals[num_normals++] = vec3(roots[i],0.0f,-Math::sqrt(1.0f - roots[i] * roots[i]));
		}
		
		if(num_normals == 2) {
			
			vec4 p0 = projection * vec4(normals[0].z * near_clip_plane / normals[0].x,0.0f,-near_clip_plane,1.0f);
			vec4 p1 = projection * vec4(normals[1].z * near_clip_plane / normals[1].x,0.0f,-near_clip_plane,1.0f);
			p0.x /= Math::abs(p0.w);
			p1.x /= Math::abs(p1.w);
			
			if(normals[0].x > 0.0f) min.x = ::max(min.x,p0.x);
			if(normals[0].x < 0.0f) max.x = ::min(max.x,p0.x);
			if(normals[1].x > 0.0f) min.x = ::max(min.x,p1.x);
			if(normals[1].x < 0.0f) max.x = ::min(max.x,p1.x);
		}
	}
	
	/////////////////////////////////
	// project light on yz plane
	/////////////////////////////////
	
	// projection
	if(is_sphere) {
		height = radius.x;
	} else {
		ellipsoid_projection(basis,radius,vec3(1.0f,0.0f,0.0f),axis_a,axis_b);
		dot_a = dot(axis_a,axis_a);
		dot_b = dot(axis_b,axis_b);
		height = Math::sqrt(axis_a.z * axis_a.z + dot_b * (1.0f - axis_a.z * axis_a.z / dot_a));
	}
	
	// depth bounds
	depth[2] = projection * vec4(position - vec3(0.0f,0.0f,height),1.0f);
	depth[3] = projection * vec4(position + vec3(0.0f,0.0f,height),1.0f);
	
	// calculate depth bounds
	for(int i = 0; i < 4; i++) {
		depth[i].z /= ::max(Math::abs(depth[i].w),EPSILON);
		min.z = ::min(min.z,depth[i].z);
		max.z = ::max(max.z,depth[i].z);
	}
	if(min.z > 1.0f || max.z < -1.0f) return 0;
	
	// circle
	if(is_sphere || compare(dot_a,dot_b,1e-2f)) {
		
		float r = (is_sphere) ? radius.x : Math::sqrt(::max(dot_a,dot_b));
		float a = position.y * position.y + position.z * position.z;
		float b = -r * position.y;
		float c = r * r - position.z * position.z;
		float d = b * b - a * c;
		
		if(d > 0.0f && Math::abs(position.z) > r * 1e-4f) {
			
			vec3 n0,n1;
			n0.y = (-b + Math::sqrt(d)) / a;
			n1.y = (-b - Math::sqrt(d)) / a;
			n0.z = (r - n0.y * position.y) / position.z;
			n1.z = (r - n1.y * position.y) / position.z;
			
			vec4 p0 = projection * vec4(0.0f,n0.z * near_clip_plane / n0.y,-near_clip_plane,1.0f);
			vec4 p1 = projection * vec4(0.0f,n1.z * near_clip_plane / n1.y,-near_clip_plane,1.0f);
			p0.y /= Math::abs(p0.w);
			p1.y /= Math::abs(p1.w);
			
			if(n0.y > 0.0f) min.y = ::max(min.y,p0.y);
			if(n0.y < 0.0f) max.y = ::min(max.y,p0.y);
			if(n1.y > 0.0f) min.y = ::max(min.y,p1.y);
			if(n1.y < 0.0f) max.y = ::min(max.y,p1.y);
		}
	}
	// ellipse
	else {
		
		vec3 normals[4];
		int num_normals = 0;
		
		// bottom plane
		float roots[4];
		int num_roots = ellipse_solver(axis_a,axis_b,position,ellipse_bottom,roots);
		for(int i = 0; i < num_roots; i++) {
			normals[num_normals++] = vec3(0.0f,roots[i],Math::sqrt(1.0f - roots[i] * roots[i]));
		}
		
		// top plane
		num_roots = ellipse_solver(axis_a,axis_b,position,ellipse_top,roots);
		for(int i = 0; i < num_roots; i++) {
			normals[num_normals++] = vec3(0.0f,roots[i],-Math::sqrt(1.0f - roots[i] * roots[i]));
		}
		
		if(num_normals == 2) {
			
			vec4 p0 = projection * vec4(0.0f,normals[0].z * near_clip_plane / normals[0].y,-near_clip_plane,1.0f);
			vec4 p1 = projection * vec4(0.0f,normals[1].z * near_clip_plane / normals[1].y,-near_clip_plane,1.0f);
			p0.y /= Math::abs(p0.w);
			p1.y /= Math::abs(p1.w);
			
			if(normals[0].y > 0.0f) min.y = ::max(min.y,p0.y);
			if(normals[0].y < 0.0f) max.y = ::min(max.y,p0.y);
			if(normals[1].y > 0.0f) min.y = ::max(min.y,p1.y);
			if(normals[1].y < 0.0f) max.y = ::min(max.y,p1.y);
		}
	}
	
	return set(min,max);
}

/*
 */
int Scissor::addLight(const mat4 &modelviewprojection) {
	
	// find light bounds
	vec3 min = vec3_infinity;
	vec3 max = -vec3_infinity;
	
	// unproject light
	vec4 points[8];
	mat4 imodelviewprojection = modelview * inverse(modelviewprojection);
	points[0] = imodelviewprojection * vec4(-1.0f,-1.0f,-1.0f,1.0f);
	points[1] = imodelviewprojection * vec4( 1.0f,-1.0f,-1.0f,1.0f);
	points[2] = imodelviewprojection * vec4( 1.0f, 1.0f,-1.0f,1.0f);
	points[3] = imodelviewprojection * vec4(-1.0f, 1.0f,-1.0f,1.0f);
	points[4] = imodelviewprojection * vec4(-1.0f,-1.0f, 1.0f,1.0f);
	points[5] = imodelviewprojection * vec4( 1.0f,-1.0f, 1.0f,1.0f);
	points[6] = imodelviewprojection * vec4( 1.0f, 1.0f, 1.0f,1.0f);
	points[7] = imodelviewprojection * vec4(-1.0f, 1.0f, 1.0f,1.0f);
	
	// edges
	static const int indices[12][2] = {
		{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
		{ 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
		{ 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
	};
	
	// normalize coordinates
	for(int i = 0; i < 8; i++) {
		points[i] /= Math::abs(points[i].w);
	}
	
	// visible points
	for(int i = 0; i < 8; i++) {
		if(points[i].z <= -near_clip_plane) {
			vec4 p = projection * points[i];
			min = ::min(min,vec3(p) / Math::abs(p.w));
			max = ::max(max,vec3(p) / Math::abs(p.w));
		}
	}
	
	// near clip plane intersections
	for(int i = 0; i < 12; i++) {
		const vec4 &p0 = points[indices[i][0]];
		const vec4 &p1 = points[indices[i][1]];
		if((p0.z > -near_clip_plane) != (p1.z > -near_clip_plane)) {
			vec4 p = projection * (p0 + (p1 - p0) * (-near_clip_plane - p0.z) / (p1.z - p0.z));
			min = ::min(min,vec3(p) / Math::abs(p.w));
			max = ::max(max,vec3(p) / Math::abs(p.w));
		}
	}
	
	return set(min,max);
}
