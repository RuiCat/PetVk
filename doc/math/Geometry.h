/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Geometry.h
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

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "MathLib.h"

/******************************************************************************\
*
* Angle normalization
*
\******************************************************************************/

/*
 */
INLINE float normalizeAngle(float angle) {
	if(angle >= 180.0f) angle = Math::mod(angle + 180.0f,360.0f) - 180.0f;
	else if(angle <= -180.0f) angle = -Math::mod(-angle - 180.0f,360.0f) + 180.0f;
	return angle;
}

/******************************************************************************\
*
* Orthonormal basis
*
\******************************************************************************/

/*
 */
INLINE void orthoBasis(const vec3 &v,vec3 &tangent,vec3 &binormal) {
	if(Math::abs(v.z) > 0.7f) {
		float length2 = v.y * v.y + v.z * v.z;
		float ilength = Math::rsqrt(length2);
		tangent.x = 0.0f;
		tangent.y = -v.z * ilength;
		tangent.z = v.y * ilength;
		binormal.x = length2 * ilength;
		binormal.y = -v.x * tangent.z;
		binormal.z = v.x * tangent.y;
	} else {
		float length2 = v.x * v.x + v.y * v.y;
		float ilength = Math::rsqrt(length2);
		tangent.x = -v.y * ilength;
		tangent.y = v.x * ilength;
		tangent.z = 0.0f;
		binormal.x = -v.z * tangent.y;
		binormal.y = v.z * tangent.x;
		binormal.z = length2 * ilength;
	}
}

INLINE void orthoBasis(const dvec3 &v,dvec3 &tangent,dvec3 &binormal) {
	if(Math::abs(v.z) > 0.7) {
		double length2 = v.y * v.y + v.z * v.z;
		double ilength = Math::rsqrt(length2);
		tangent.x = 0.0;
		tangent.y = -v.z * ilength;
		tangent.z = v.y * ilength;
		binormal.x = length2 * ilength;
		binormal.y = -v.x * tangent.z;
		binormal.z = v.x * tangent.y;
	} else {
		double length2 = v.x * v.x + v.y * v.y;
		double ilength = Math::rsqrt(length2);
		tangent.x = -v.y * ilength;
		tangent.y = v.x * ilength;
		tangent.z = 0.0;
		binormal.x = -v.z * tangent.y;
		binormal.y = v.z * tangent.x;
		binormal.z = length2 * ilength;
	}
}

/*
 */
void orthoTransform(const vec3 &v,mat4 &transform);
void orthoTransform(const dvec3 &v,dmat4 &transform);

/******************************************************************************\
*
* Triangle parameters
*
\******************************************************************************/

/*
 */
float triangleArea(const vec3 &v0,const vec3 &v1,const vec3 &v2);
double triangleArea(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2);

/*
 */
vec3 triangleNormal(const vec3 &v0,const vec3 &v1,const vec3 &v2);
dvec3 triangleNormal(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2);

/*
 */
vec4 trianglePlane(const vec3 &v0,const vec3 &v1,const vec3 &v2);
dvec4 trianglePlane(const dvec3 &v0,const dvec3 &v1,const dvec3 &v2);

/******************************************************************************\
*
* Polygon parameters
*
\******************************************************************************/

/*
 */
vec4 polygonPlane(const vec3 *vertex,const unsigned short *indices,int num_indices);
dvec4 polygonPlane(const dvec3 *vertex,const unsigned short *indices,int num_indices);

/******************************************************************************\
*
* Point tringle intersections
*
\******************************************************************************/

/*
 */
int pointTriangleInside(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2);
int pointTriangleInside(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2);
int pointTriangleInside(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,const vec3 &normal);
int pointTriangleInside(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,const dvec3 &normal);

/*
 */
float pointTriangleDistance(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,const vec4 &plane);
double pointTriangleDistance(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,const dvec4 &plane);

/*
 */
void pointTriangleCoordinates(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,float &a,float &b);
void pointTriangleCoordinates(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,double &a,double &b);

/******************************************************************************\
*
* Point polygon intersections
*
\******************************************************************************/

/*
 */
int pointPolygonInside(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices);
int pointPolygonInside(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices);
int pointPolygonInside(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices,const vec3 &normal);
int pointPolygonInside(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices,const dvec3 &normal);

/*
 */
float pointPolygonDistance(const vec3 &point,const vec3 *vertex,const unsigned short *indices,int num_indices,const vec4 &plane);
double pointPolygonDistance(const dvec3 &point,const dvec3 *vertex,const unsigned short *indices,int num_indices,const dvec4 &plane);

/******************************************************************************\
*
* Bounding box intersections
*
\******************************************************************************/

/*
 */
int rayBoundBoxIntersection(const vec3 &point,const vec3 &direction,const vec3 &min,const vec3 &max);
int rayBoundBoxIntersection(const dvec3 &point,const dvec3 &direction,const dvec3 &min,const dvec3 &max);

/*
 */
int irayBoundBoxIntersection(const vec3 &point,const vec3 &idirection,const vec3 &min,const vec3 &max);
int irayBoundBoxIntersection(const dvec3 &point,const dvec3 &idirection,const dvec3 &min,const dvec3 &max);

/*
 */
int rayTriangleIntersection(const vec3 &point,const vec3 &direction,const vec3 &v0,const vec3 &v1,const vec3 &v2);
int rayTriangleIntersection(const dvec3 &point,const dvec3 &direction,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2);

/******************************************************************************\
*
* Closest points
*
\******************************************************************************/

/*
 */
void getClosestPointOnLine(const vec3 &point,const vec3 &p0,const vec3 &p1,vec3 &ret);
void getClosestPointOnLine(const dvec3 &point,const dvec3 &p0,const dvec3 &p1,dvec3 &ret);

/*
 */
int getClosestPointOnTriangle(const vec3 &point,const vec3 &v0,const vec3 &v1,const vec3 &v2,vec3 &ret);
int getClosestPointOnTriangle(const dvec3 &point,const dvec3 &v0,const dvec3 &v1,const dvec3 &v2,dvec3 &ret);

/*
 */
int getClosestPointsOnLines(const vec3 &p00,const vec3 &p01,const vec3 &p10,const vec3 &p11,vec3 &ret_0,vec3 &ret_1);
int getClosestPointsOnLines(const dvec3 &p00,const dvec3 &p01,const dvec3 &p10,const dvec3 &p11,dvec3 &ret_0,dvec3 &ret_1);

#endif /* __GEOMETRY_H__ */
