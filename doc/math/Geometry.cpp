/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Geometry.cpp
 * Desc:    Geometry utils
 * Version: 1.10
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

#include "Geometry.h"

/******************************************************************************\
*
* Orthonormal transformation
*
\******************************************************************************/

/*
 */
void orthoTransform(const vec3 &v,mat4 &transform) {
	if(Math::abs(v.z) > 0.7f) {
		float length2 = v.y * v.y + v.z * v.z;
		float ilength = Math::rsqrt(length2);
		transform.m00 = 0.0f;
		transform.m10 = -v.z * ilength;
		transform.m20 = v.y * ilength;
		transform.m01 = length2 * ilength;
		transform.m11 = -v.x * transform.m20;
		transform.m21 = v.x * transform.m10;
	} else {
		float length2 = v.x * v.x + v.y * v.y;
		float ilength = Math::rsqrt(length2);
		transform.m00 = -v.y * ilength;
		transform.m10 = v.x * ilength;
		transform.m20 = 0.0f;
		transform.m01 = -v.z * transform.m10;
		transform.m11 = v.z * transform.m00;
		transform.m21 = length2 * ilength;
	}
	transform.m02 = v.x;
	transform.m12 = v.y;
	transform.m22 = v.z;
}

void orthoTransform(const dvec3 &v,dmat4 &transform) {
	if(Math::abs(v.z) > 0.7) {
		double length2 = v.y * v.y + v.z * v.z;
		double ilength = Math::rsqrt(length2);
		transform.m00 = 0.0;
		transform.m10 = -v.z * ilength;
		transform.m20 = v.y * ilength;
		transform.m01 = length2 * ilength;
		transform.m11 = -v.x * transform.m20;
		transform.m21 = v.x * transform.m10;
	} else {
		double length2 = v.x * v.x + v.y * v.y;
		double ilength = Math::rsqrt(length2);
		transform.m00 = -v.y * ilength;
		transform.m10 = v.x * ilength;
		transform.m20 = 0.0;
		transform.m01 = -v.z * transform.m10;
		transform.m11 = v.z * transform.m00;
		transform.m21 = length2 * ilength;
	}
	transform.m02 = v.x;
	transform.m12 = v.y;
	transform.m22 = v.z;
}

/******************************************************************************\
*
* Triangle parameters
*
\******************************************************************************/

/*
 */
float triangleArea(const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	float a = length(v1 - v2);
	float b = length(v2 - v0);
	float c = length(v0 - v1);
	float s = (a + b + c) * 0.5f;
	return Math::sqrt(s * (s - a) * (s - b) * (s - c));
}

double triangleArea(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	double a = length(v1 - v2);
	double b = length(v2 - v0);
	double c = length(v0 - v1);
	double s = (a + b + c) * 0.5;
	return Math::sqrt(s * (s - a) * (s - b) * (s - c));
}

/*
 */
vec3 triangleNormal(const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	vec3 v10,v20,normal;
	cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
	normal.normalize();
	return normal;
}

dvec3 triangleNormal(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	dvec3 v10,v20,normal;
	cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
	normal.normalize();
	return normal;
}

/*
 */
vec4 trianglePlane(const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	vec3 v10,v20,normal;
	cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
	normal.normalize();
	return vec4(normal,-dot(normal,v0));
}

dvec4 trianglePlane(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	dvec3 v10,v20,normal;
	cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
	normal.normalize();
	return dvec4(normal,-dot(normal,v0));
}

/******************************************************************************\
*
* Polygon parameters
*
\******************************************************************************/

/*
 */
vec4 polygonPlane(const vec3 *vertex,const unsigned short *indices,int num_indices) {
	float area = 0.0f;
	vec3 v10,v20,normal;
	vec4 plane = vec4_zero;
	for(int i = 2; i < num_indices; i++) {
		const vec3 &v0 = vertex[indices[i - 2]];
		const vec3 &v1 = vertex[indices[i - 1]];
		const vec3 &v2 = vertex[indices[i]];
		cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
		float a = length(normal);
		if(area < a) {
			normal *= Math::rcp(a);
			plane.set(normal,-dot(normal,v0));
			area = a;
		}
	}
	return plane;
}

dvec4 polygonPlane(const dvec3 *vertex,const unsigned short *indices,int num_indices) {
	double area = 0.0;
	dvec3 v10,v20,normal;
	dvec4 plane = dvec4_zero;
	for(int i = 2; i < num_indices; i++) {
		const dvec3 &v0 = vertex[indices[i - 2]];
		const dvec3 &v1 = vertex[indices[i - 1]];
		const dvec3 &v2 = vertex[indices[i]];
		cross(normal,sub(v10,v1,v0),sub(v20,v2,v0));
		double a = length(normal);
		if(area < a) {
			normal *= Math::rcp(a);
			plane.set(normal,-dot(normal,v0));
			area = a;
		}
	}
	return plane;
}

/******************************************************************************\
*
* Point tringle intersections
*
\******************************************************************************/

/*
 */
int pointTriangleInside(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	vec3 axis,normal;
	vec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	cross(axis,normal,v20);
	float det = dot(v10,axis);
	if(det > 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s < 0.0f || s > det) return 0;
		cross(axis,pv0,v10);
		float t = dot(normal,axis);
		if(t < 0.0f || t + s > det) return 0;
		return 1;
	} else if(det < 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s > 0.0f || s < det) return 0;
		cross(axis,pv0,v10);
		float t = dot(normal,axis);
		if(t > 0.0f || t + s < det) return 0;
		return 1;
	}
	return 0;
}

int pointTriangleInside(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	dvec3 axis,normal;
	dvec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	cross(axis,normal,v20);
	double det = dot(v10,axis);
	if(det > 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s < 0.0 || s > det) return 0;
		cross(axis,pv0,v10);
		double t = dot(normal,axis);
		if(t < 0.0 || t + s > det) return 0;
		return 1;
	} else if(det < 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s > 0.0 || s < det) return 0;
		cross(axis,pv0,v10);
		double t = dot(normal,axis);
		if(t > 0.0 || t + s < det) return 0;
		return 1;
	}
	return 0;
}

/*
 */
int pointTriangleInside(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,const vec3 &normal) {
	vec3 axis;
	vec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,normal,v20);
	float det = dot(v10,axis);
	if(det > 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s < 0.0f || s > det) return 0;
		cross(axis,pv0,v10);
		float t = dot(normal,axis);
		if(t < 0.0f || t + s > det) return 0;
		return 1;
	} else if(det < 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s > 0.0f || s < det) return 0;
		cross(axis,pv0,v10);
		float t = dot(normal,axis);
		if(t > 0.0f || t + s < det) return 0;
		return 1;
	}
	return 0;
}

int pointTriangleInside(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,const dvec3 &normal) {
	dvec3 axis;
	dvec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,normal,v20);
	double det = dot(v10,axis);
	if(det > 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s < 0.0 || s > det) return 0;
		cross(axis,pv0,v10);
		double t = dot(normal,axis);
		if(t < 0.0 || t + s > det) return 0;
		return 1;
	} else if(det < 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s > 0.0 || s < det) return 0;
		cross(axis,pv0,v10);
		double t = dot(normal,axis);
		if(t > 0.0 || t + s < det) return 0;
		return 1;
	}
	return 0;
}

/*
 */
float pointTriangleDistance(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,const vec4 &plane) {
	vec3 axis,normal = vec3(plane);
	vec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,normal,v20);
	float det = dot(v10,axis);
	if(det > 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s >= 0.0f && s <= det) {
			cross(axis,pv0,v10);
			float t = dot(normal,axis);
			if(t >= 0.0f && t + s <= det) {
				return Math::abs(dot(plane,point));
			}
		}
	} else if(det < 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s <= 0.0f && s >= det) {
			cross(axis,pv0,v10);
			float t = dot(normal,axis);
			if(t <= 0.0f && t + s >= det) {
				return Math::abs(dot(plane,point));
			}
		}
	}
	vec3 v21;
	sub(v21,v2,v1);
	float k0 = saturate(dot(v10,sub(pv0,point,v0)) * Math::rcp(length2(v10)));
	float k1 = saturate(dot(v20,sub(pv0,point,v0)) * Math::rcp(length2(v20)));
	float k2 = saturate(dot(v21,sub(pv0,point,v1)) * Math::rcp(length2(v21)));
	k0 = length2(sub(pv0,point,mad(pv0,v10,k0,v0)));
	k1 = length2(sub(pv0,point,mad(pv0,v20,k1,v0)));
	k2 = length2(sub(pv0,point,mad(pv0,v21,k2,v1)));
	return Math::sqrt(min(k0,min(k1,k2)));
}

double pointTriangleDistance(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,const dvec4 &plane) {
	dvec3 axis,normal = dvec3(plane);
	dvec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,normal,v20);
	double det = dot(v10,axis);
	if(det > 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s >= 0.0 && s <= det) {
			cross(axis,pv0,v10);
			double t = dot(normal,axis);
			if(t >= 0.0 && t + s <= det) {
				return Math::abs(dot(plane,point));
			}
		}
	} else if(det < 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s <= 0.0 && s >= det) {
			cross(axis,pv0,v10);
			double t = dot(normal,axis);
			if(t <= 0.0 && t + s >= det) {
				return Math::abs(dot(plane,point));
			}
		}
	}
	dvec3 v21;
	sub(v21,v2,v1);
	double k0 = saturate(dot(v10,sub(pv0,point,v0)) * Math::rcp(length2(v10)));
	double k1 = saturate(dot(v20,sub(pv0,point,v0)) * Math::rcp(length2(v20)));
	double k2 = saturate(dot(v21,sub(pv0,point,v1)) * Math::rcp(length2(v21)));
	k0 = length2(sub(pv0,point,mad(pv0,v10,k0,v0)));
	k1 = length2(sub(pv0,point,mad(pv0,v20,k1,v0)));
	k2 = length2(sub(pv0,point,mad(pv0,v21,k2,v1)));
	return Math::sqrt(min(k0,min(k1,k2)));
}

/*
 */
void pointTriangleCoordinates(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,float &a,float &b) {
	vec3 area;
	vec3 v20,v10,v0p;
	sub(v20,v2,v0);
	sub(v10,v1,v0);
	float iarea = Math::rsqrt(length2(cross(area,v20,v10)));
	sub(v20,v2,point);
	sub(v10,v1,point);
	sub(v0p,v0,point);
	a = length(cross(area,v20,v0p)) * iarea;
	b = length(cross(area,v10,v0p)) * iarea;
}

void pointTriangleCoordinates(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,double &a,double &b) {
	dvec3 area;
	dvec3 v20,v10,v0p;
	sub(v20,v2,v0);
	sub(v10,v1,v0);
	double iarea = Math::rsqrt(length2(cross(area,v20,v10)));
	sub(v20,v2,point);
	sub(v10,v1,point);
	sub(v0p,v0,point);
	a = length(cross(area,v20,v0p)) * iarea;
	b = length(cross(area,v10,v0p)) * iarea;
}

/******************************************************************************\
*
* Point polygon intersections
*
\******************************************************************************/

/*
 */
int pointPolygonInside(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices) {
	vec3 clip,normal;
	vec3 v10,v20;
	const vec3 &v0 = vertex[indices[0]];
	const vec3 &v1 = vertex[indices[1]];
	const vec3 &v2 = vertex[indices[2]];
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const vec3 &v0 = vertex[indices[i]];
		const vec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		float d = dot(clip,v20);
		if(d > 0.0f) return 0;
	}
	return 1;
}

int pointPolygonInside(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices) {
	dvec3 clip,normal;
	dvec3 v10,v20;
	const dvec3 &v0 = vertex[indices[0]];
	const dvec3 &v1 = vertex[indices[1]];
	const dvec3 &v2 = vertex[indices[2]];
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const dvec3 &v0 = vertex[indices[i]];
		const dvec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		double d = dot(clip,v20);
		if(d > 0.0) return 0;
	}
	return 1;
}

/*
 */
int pointPolygonInside(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices,const vec3 &normal) {
	vec3 clip;
	vec3 v10,v20;
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const vec3 &v0 = vertex[indices[i]];
		const vec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		float d = dot(clip,v20);
		if(d > 0.0f) return 0;
	}
	return 1;
}

int pointPolygonInside(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices,const dvec3 &normal) {
	dvec3 clip;
	dvec3 v10,v20;
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const dvec3 &v0 = vertex[indices[i]];
		const dvec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		double d = dot(clip,v20);
		if(d > 0.0) return 0;
	}
	return 1;
}

/*
 */
float pointPolygonDistance(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices,const vec4 &plane) {
	vec3 clip,normal = vec3(plane);
	vec3 v10,v20;
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const vec3 &v0 = vertex[indices[i]];
		const vec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		float d = dot(clip,v20);
		if(d > 0.0f) {
			float distance = INFINITY;
			for(; i < num_indices; i++, j++) {
				if(j == num_indices) j = 0;
				const vec3 &v0 = vertex[indices[i]];
				const vec3 &v1 = vertex[indices[j]];
				sub(v10,v1,v0);
				sub(v20,point,v0);
				float k = saturate(dot(v10,v20) * Math::rcp(length2(v10)));
				float d = length2(sub(v20,point,mad(v20,v10,k,v0)));
				if(distance > d) distance = d;
			}
			return Math::sqrt(distance);
		}
	}
	return Math::abs(dot(plane,point));
}

double pointPolygonDistance(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices,const dvec4 &plane) {
	dvec3 clip,normal = dvec3(plane);
	dvec3 v10,v20;
	for(int i = 0, j = 1; i < num_indices; i++, j++) {
		if(j == num_indices) j = 0;
		const dvec3 &v0 = vertex[indices[i]];
		const dvec3 &v1 = vertex[indices[j]];
		sub(v10,v1,v0);
		sub(v20,point,v0);
		cross(clip,v10,normal);
		double d = dot(clip,v20);
		if(d > 0.0) {
			double distance = INFINITY;
			for(; i < num_indices; i++, j++) {
				if(j == num_indices) j = 0;
				const dvec3 &v0 = vertex[indices[i]];
				const dvec3 &v1 = vertex[indices[j]];
				sub(v10,v1,v0);
				sub(v20,point,v0);
				double k = saturate(dot(v10,v20) * Math::rcp(length2(v10)));
				double d = length2(sub(v20,point,mad(v20,v10,k,v0)));
				if(distance > d) distance = d;
			}
			return Math::sqrt(distance);
		}
	}
	return Math::abs(dot(plane,point));
}

/******************************************************************************\
*
* Bounding box intersections
*
\******************************************************************************/

/*
 */
int rayBoundBoxIntersection(const vec3 &point,const vec3 &direction,const vec3 &min,const vec3 &max) {
	float tmin,tmax;
	float idirectionx = Math::rcp(direction.x);
	float idirectiony = Math::rcp(direction.y);
	if(direction.x >= 0.0f) {
		tmin = (min.x - point.x) * idirectionx;
		tmax = (max.x - point.x) * idirectionx;
	} else {
		tmin = (max.x - point.x) * idirectionx;
		tmax = (min.x - point.x) * idirectionx;
	}
	float tymin,tymax;
	if(direction.y >= 0.0f) {
		tymin = (min.y - point.y) * idirectiony;
		tymax = (max.y - point.y) * idirectiony;
	} else {
		tymin = (max.y - point.y) * idirectiony;
		tymax = (min.y - point.y) * idirectiony;
	}
	if((tmin > tymax) || (tmax < tymin)) return 0;
	if(tmin < tymin) tmin = tymin;
	if(tmax > tymax) tmax = tymax;
	float tzmin,tzmax;
	float idirectionz = Math::rcp(direction.z);
	if(direction.z >= 0.0f) {
		tzmin = (min.z - point.z) * idirectionz;
		tzmax = (max.z - point.z) * idirectionz;
	} else {
		tzmin = (max.z - point.z) * idirectionz;
		tzmax = (min.z - point.z) * idirectionz;
	}
	if((tmin > tzmax) || (tmax < tzmin)) return 0;
	if(tmin < tzmin) tmin = tzmin;
	if(tmax > tzmax) tmax = tzmax;
	return (tmax > 0.0f && tmin < 1.0f);
}

int rayBoundBoxIntersection(const dvec3 &point,const dvec3 &direction,const dvec3 &min,const dvec3 &max) {
	double tmin,tmax;
	double idirectionx = Math::rcp(direction.x);
	double idirectiony = Math::rcp(direction.y);
	if(direction.x >= 0.0) {
		tmin = (min.x - point.x) * idirectionx;
		tmax = (max.x - point.x) * idirectionx;
	} else {
		tmin = (max.x - point.x) * idirectionx;
		tmax = (min.x - point.x) * idirectionx;
	}
	double tymin,tymax;
	if(direction.y >= 0.0) {
		tymin = (min.y - point.y) * idirectiony;
		tymax = (max.y - point.y) * idirectiony;
	} else {
		tymin = (max.y - point.y) * idirectiony;
		tymax = (min.y - point.y) * idirectiony;
	}
	if((tmin > tymax) || (tmax < tymin)) return 0;
	if(tmin < tymin) tmin = tymin;
	if(tmax > tymax) tmax = tymax;
	double tzmin,tzmax;
	double idirectionz = Math::rcp(direction.z);
	if(direction.z >= 0.0) {
		tzmin = (min.z - point.z) * idirectionz;
		tzmax = (max.z - point.z) * idirectionz;
	} else {
		tzmin = (max.z - point.z) * idirectionz;
		tzmax = (min.z - point.z) * idirectionz;
	}
	if((tmin > tzmax) || (tmax < tzmin)) return 0;
	if(tmin < tzmin) tmin = tzmin;
	if(tmax > tzmax) tmax = tzmax;
	return (tmax > 0.0 && tmin < 1.0);
}

/*
 */
int irayBoundBoxIntersection(const vec3 &point,const vec3 &idirection,const vec3 &min,const vec3 &max) {
	float tmin,tmax;
	if(idirection.x >= 0.0f) {
		tmin = (min.x - point.x) * idirection.x;
		tmax = (max.x - point.x) * idirection.x;
	} else {
		tmin = (max.x - point.x) * idirection.x;
		tmax = (min.x - point.x) * idirection.x;
	}
	float tymin,tymax;
	if(idirection.y >= 0.0f) {
		tymin = (min.y - point.y) * idirection.y;
		tymax = (max.y - point.y) * idirection.y;
	} else {
		tymin = (max.y - point.y) * idirection.y;
		tymax = (min.y - point.y) * idirection.y;
	}
	if((tmin > tymax) || (tmax < tymin)) return 0;
	if(tmin < tymin) tmin = tymin;
	if(tmax > tymax) tmax = tymax;
	float tzmin,tzmax;
	if(idirection.z >= 0.0f) {
		tzmin = (min.z - point.z) * idirection.z;
		tzmax = (max.z - point.z) * idirection.z;
	} else {
		tzmin = (max.z - point.z) * idirection.z;
		tzmax = (min.z - point.z) * idirection.z;
	}
	if((tmin > tzmax) || (tmax < tzmin)) return 0;
	if(tmin < tzmin) tmin = tzmin;
	if(tmax > tzmax) tmax = tzmax;
	return (tmax > 0.0f && tmin < 1.0f);
}

int irayBoundBoxIntersection(const dvec3 &point,const dvec3 &idirection,const dvec3 &min,const dvec3 &max) {
	double tmin,tmax;
	if(idirection.x >= 0.0) {
		tmin = (min.x - point.x) * idirection.x;
		tmax = (max.x - point.x) * idirection.x;
	} else {
		tmin = (max.x - point.x) * idirection.x;
		tmax = (min.x - point.x) * idirection.x;
	}
	double tymin,tymax;
	if(idirection.y >= 0.0) {
		tymin = (min.y - point.y) * idirection.y;
		tymax = (max.y - point.y) * idirection.y;
	} else {
		tymin = (max.y - point.y) * idirection.y;
		tymax = (min.y - point.y) * idirection.y;
	}
	if((tmin > tymax) || (tmax < tymin)) return 0;
	if(tmin < tymin) tmin = tymin;
	if(tmax > tymax) tmax = tymax;
	double tzmin,tzmax;
	if(idirection.z >= 0.0) {
		tzmin = (min.z - point.z) * idirection.z;
		tzmax = (max.z - point.z) * idirection.z;
	} else {
		tzmin = (max.z - point.z) * idirection.z;
		tzmax = (min.z - point.z) * idirection.z;
	}
	if((tmin > tzmax) || (tmax < tzmin)) return 0;
	if(tmin < tzmin) tmin = tzmin;
	if(tmax > tzmax) tmax = tzmax;
	return (tmax > 0.0 && tmin < 1.0);
}

/*
 */
int rayTriangleIntersection(const vec3 &point,const vec3 &direction,const vec3 &v0,const vec3 &v1,const vec3 &v2) {
	vec3 axis;
	vec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,direction,v20);
	float det = dot(v10,axis);
	if(det > 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s < 0.0f || s > det) return 0;
		cross(axis,pv0,v10);
		float t = dot(direction,axis);
		if(t < 0.0f || t + s > det) return 0;
		return 1;
	} else if(det < 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s > 0.0f || s < det) return 0;
		cross(axis,pv0,v10);
		float t = dot(direction,axis);
		if(t > 0.0f || t + s < det) return 0;
		return 1;
	}
	return 0;
}

int rayTriangleIntersection(const dvec3 &point,const dvec3 &direction,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2) {
	dvec3 axis;
	dvec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(axis,direction,v20);
	double det = dot(v10,axis);
	if(det > 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s < 0.0 || s > det) return 0;
		cross(axis,pv0,v10);
		double t = dot(direction,axis);
		if(t < 0.0 || t + s > det) return 0;
		return 1;
	} else if(det < 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s > 0.0 || s < det) return 0;
		cross(axis,pv0,v10);
		double t = dot(direction,axis);
		if(t > 0.0 || t + s < det) return 0;
		return 1;
	}
	return 0;
}

/******************************************************************************\
*
* Closest points
*
\******************************************************************************/

/*
 */
void getClosestPointOnLine(const vec3 &point,const vec3 &p0,const vec3 &p1,vec3 &ret) {
	vec3 direction;
	sub(direction,p1,p0);
	float d = length2(direction);
	if(d < EPSILON) {
		ret = p0;
		return;
	}
	vec3 v0;
	float k = dot(direction,sub(v0,point,p0)) * Math::rcp(d);
	mad(ret,direction,saturate(k),p0);
}

void getClosestPointOnLine(const dvec3 &point,const dvec3 &p0,const dvec3 &p1,dvec3 &ret) {
	dvec3 direction;
	sub(direction,p1,p0);
	double d = length2(direction);
	if(d < EPSILON) {
		ret = p0;
		return;
	}
	dvec3 v0;
	double k = dot(direction,sub(v0,point,p0)) * Math::rcp(d);
	mad(ret,direction,saturate(k),p0);
}

/*
 */
int getClosestPointOnTriangle(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,vec3 &ret) {
	vec3 axis,normal;
	vec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	cross(axis,normal,v20);
	float det = dot(v10,axis);
	if(det > 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s >= 0.0f && s <= det) {
			cross(axis,pv0,v10);
			float t = dot(normal,axis);
			if(t >= 0.0f && t + s <= det) {
				mad(ret,normal,dot(normal,v0) - dot(normal,point),point);
				return 1;
			}
		}
	} else if(det < 0.0f) {
		sub(pv0,point,v0);
		float s = dot(pv0,axis);
		if(s <= 0.0f && s >= det) {
			cross(axis,pv0,v10);
			float t = dot(normal,axis);
			if(t <= 0.0f && t + s >= det) {
				mad(ret,normal,dot(normal,v0) - dot(normal,point),point);
				return 1;
			}
		}
	}
	vec3 v21,p0,p1;
	sub(v21,v2,v1);
	float k0 = saturate(dot(v10,sub(pv0,point,v0)) * Math::rcp(length2(v10)));
	float k1 = saturate(dot(v20,sub(pv0,point,v0)) * Math::rcp(length2(v20)));
	float k2 = saturate(dot(v21,sub(pv0,point,v1)) * Math::rcp(length2(v21)));
	k0 = length2(sub(pv0,point,mad(p0,v10,k0,v0)));
	k1 = length2(sub(pv0,point,mad(p1,v20,k1,v0)));
	k2 = length2(sub(pv0,point,mad(ret,v21,k2,v1)));
	if(k0 < k1) { if(k0 < k2) ret = vec3(p0); }
	else { if(k1 < k2) ret = vec3(p1); }
	return 0;
}

int getClosestPointOnTriangle(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,dvec3 &ret) {
	dvec3 axis,normal;
	dvec3 v10,v20,pv0;
	sub(v10,v1,v0);
	sub(v20,v2,v0);
	cross(normal,v10,v20);
	cross(axis,normal,v20);
	double det = dot(v10,axis);
	if(det > 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s >= 0.0 && s <= det) {
			cross(axis,pv0,v10);
			double t = dot(normal,axis);
			if(t >= 0.0 && t + s <= det) {
				mad(ret,normal,dot(normal,v0) - dot(normal,point),point);
				return 1;
			}
		}
	} else if(det < 0.0) {
		sub(pv0,point,v0);
		double s = dot(pv0,axis);
		if(s <= 0.0 && s >= det) {
			cross(axis,pv0,v10);
			double t = dot(normal,axis);
			if(t <= 0.0 && t + s >= det) {
				mad(ret,normal,dot(normal,v0) - dot(normal,point),point);
				return 1;
			}
		}
	}
	dvec3 v21,p0,p1;
	sub(v21,v2,v1);
	double k0 = saturate(dot(v10,sub(pv0,point,v0)) * Math::rcp(length2(v10)));
	double k1 = saturate(dot(v20,sub(pv0,point,v0)) * Math::rcp(length2(v20)));
	double k2 = saturate(dot(v21,sub(pv0,point,v1)) * Math::rcp(length2(v21)));
	k0 = length2(sub(pv0,point,mad(p0,v10,k0,v0)));
	k1 = length2(sub(pv0,point,mad(p1,v20,k1,v0)));
	k2 = length2(sub(pv0,point,mad(ret,v21,k2,v1)));
	if(k0 < k1) { if(k0 < k2) ret = dvec3(p0); }
	else { if(k1 < k2) ret = dvec3(p1); }
	return 0;
}

/*
 */
int getClosestPointsOnLines(const vec3 &p00,const vec3 &p01,const vec3 &p10,const vec3 &p11,vec3 &ret_0,vec3 &ret_1) {
	if(p00 == p01) {
		ret_0 = p00;
		getClosestPointOnLine(p00,p10,p11,ret_1);
		return 0;
	}
	if(p10 == p11) {
		ret_1 = p10;
		getClosestPointOnLine(p10,p00,p01,ret_0);
		return 0;
	}
	vec3 v0,v1,v2;
	sub(v0,p01,p00);
	sub(v1,p11,p10);
	sub(v2,p10,p00);
	float d20 = dot(v2,v0);
	float d21 = dot(v2,v1);
	if(d20 < EPSILON && d21 > -EPSILON) {
		ret_0 = p00;
		ret_1 = p10;
		return 0;
	}
	vec3 v3;
	sub(v3,p11,p01);
	float d30 = dot(v3,v0);
	float d31 = dot(v3,v1);
	if(d30 > -EPSILON && d31 < EPSILON) {
		ret_0 = p01;
		ret_1 = p11;
		return 0;
	}
	vec3 v4;
	sub(v4,p11,p00);
	float d40 = dot(v4,v0);
	float d41 = dot(v4,v1);
	if(d40 < EPSILON && d41 < EPSILON) {
		ret_0 = p00;
		ret_1 = p11;
		return 0;
	}
	vec3 v5;
	sub(v5,p10,p01);
	float d50 = dot(v5,v0);
	float d51 = dot(v5,v1);
	if(d50 > -EPSILON && d51 > -EPSILON) {
		ret_0 = p01;
		ret_1 = p10;
		return 0;
	}
	vec3 v6;
	float d00 = dot(v0,v0);
	if(d20 > -EPSILON && d50 < EPSILON) {
		float k = d20 / d00;
		if(dot(mad(v6,v0,-k,v2),v1) > -EPSILON) {
			mad(ret_0,v0,k,p00);
			ret_1 = p10;
			return 0;
		}
	}
	if(d40 > -EPSILON && d30 < EPSILON) {
		float k = d40 / d00;
		if(dot(mad(v6,v0,-k,v4),v1) < EPSILON) {
			mad(ret_0,v0,k,p00);
			ret_1 = p11;
			return 0;
		}
	}
	float d11 = dot(v1,v1);
	if(d21 < EPSILON && d41 > -EPSILON) {
		float k = -d21 / d11;
		if(dot(mad(v6,v1,k,v2),v0) < EPSILON) {
			ret_0 = p00;
			mad(ret_1,v1,k,p10);
			return 0;
		}
	}
	if(d51 < EPSILON && d31 > -EPSILON) {
		float k = -d51 / d11;
		if(dot(mad(v6,v1,k,v5),v0) > -EPSILON) {
			ret_0 = p01;
			mad(ret_1,v1,k,p10);
			return 0;
		}
	}
	float d10 = dot(v1,v0);
	float det = d10 * d10 - d00 * d11;
	if(Math::abs(det) < EPSILON) {
		float distance = INFINITY;
		float d2 = length2(v2);
		if(distance > d2) {
			distance = d2;
			ret_0 = p00;
			ret_1 = p10;
		}
		float d3 = length2(v3);
		if(distance > d3) {
			distance = d3;
			ret_0 = p01;
			ret_1 = p11;
		}
		float d4 = length2(v4);
		if(distance > d4) {
			distance = d4;
			ret_0 = p00;
			ret_1 = p11;
		}
		float d5 = length2(v5);
		if(distance > d5) {
			distance = d5;
			ret_0 = p01;
			ret_1 = p10;
		}
		return 0;
	}
	float idet = Math::rcp(det);
	mad(ret_0,v0,saturate((d10 * d21 - d11 * d20) * idet),p00);
	mad(ret_1,v1,saturate((d00 * d21 - d10 * d20) * idet),p10);
	return 1;
}

int getClosestPointsOnLines(const dvec3 &p00,const dvec3 &p01,const dvec3 &p10,const dvec3 &p11,dvec3 &ret_0,dvec3 &ret_1) {
	if(p00 == p01) {
		ret_0 = p00;
		getClosestPointOnLine(p00,p10,p11,ret_1);
		return 0;
	}
	if(p10 == p11) {
		ret_1 = p10;
		getClosestPointOnLine(p10,p00,p01,ret_0);
		return 0;
	}
	dvec3 v0,v1,v2;
	sub(v0,p01,p00);
	sub(v1,p11,p10);
	sub(v2,p10,p00);
	double d20 = dot(v2,v0);
	double d21 = dot(v2,v1);
	if(d20 < EPSILON && d21 > -EPSILON) {
		ret_0 = p00;
		ret_1 = p10;
		return 0;
	}
	dvec3 v3;
	sub(v3,p11,p01);
	double d30 = dot(v3,v0);
	double d31 = dot(v3,v1);
	if(d30 > -EPSILON && d31 < EPSILON) {
		ret_0 = p01;
		ret_1 = p11;
		return 0;
	}
	dvec3 v4;
	sub(v4,p11,p00);
	double d40 = dot(v4,v0);
	double d41 = dot(v4,v1);
	if(d40 < EPSILON && d41 < EPSILON) {
		ret_0 = p00;
		ret_1 = p11;
		return 0;
	}
	dvec3 v5;
	sub(v5,p10,p01);
	double d50 = dot(v5,v0);
	double d51 = dot(v5,v1);
	if(d50 > -EPSILON && d51 > -EPSILON) {
		ret_0 = p01;
		ret_1 = p10;
		return 0;
	}
	dvec3 v6;
	double d00 = dot(v0,v0);
	if(d20 > -EPSILON && d50 < EPSILON) {
		double k = d20 / d00;
		if(dot(mad(v6,v0,-k,v2),v1) > -EPSILON) {
			mad(ret_0,v0,k,p00);
			ret_1 = p10;
			return 0;
		}
	}
	if(d40 > -EPSILON && d30 < EPSILON) {
		double k = d40 / d00;
		if(dot(mad(v6,v0,-k,v4),v1) < EPSILON) {
			mad(ret_0,v0,k,p00);
			ret_1 = p11;
			return 0;
		}
	}
	double d11 = dot(v1,v1);
	if(d21 < EPSILON && d41 > -EPSILON) {
		double k = -d21 / d11;
		if(dot(mad(v6,v1,k,v2),v0) < EPSILON) {
			ret_0 = p00;
			mad(ret_1,v1,k,p10);
			return 0;
		}
	}
	if(d51 < EPSILON && d31 > -EPSILON) {
		double k = -d51 / d11;
		if(dot(mad(v6,v1,k,v5),v0) > -EPSILON) {
			ret_0 = p01;
			mad(ret_1,v1,k,p10);
			return 0;
		}
	}
	double d10 = dot(v1,v0);
	double det = d10 * d10 - d00 * d11;
	if(Math::abs(det) < EPSILON) {
		double distance = INFINITY;
		double d2 = length2(v2);
		if(distance > d2) {
			distance = d2;
			ret_0 = p00;
			ret_1 = p10;
		}
		double d3 = length2(v3);
		if(distance > d3) {
			distance = d3;
			ret_0 = p01;
			ret_1 = p11;
		}
		double d4 = length2(v4);
		if(distance > d4) {
			distance = d4;
			ret_0 = p00;
			ret_1 = p11;
		}
		double d5 = length2(v5);
		if(distance > d5) {
			distance = d5;
			ret_0 = p01;
			ret_1 = p10;
		}
		return 0;
	}
	double idet = Math::rcp(det);
	mad(ret_0,v0,saturate((d10 * d21 - d11 * d20) * idet),p00);
	mad(ret_1,v1,saturate((d00 * d21 - d10 * d20) * idet),p10);
	return 1;
}
