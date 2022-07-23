/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Concave.cpp
 * Desc:    Concave creator
 * Version: 1.05
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
#include "Utils.h"
#include "Convex.h"
#include "Bounds.h"
#include "SimdLib.h"
#include "Concave.h"

/*
 */
Concave::Concave() {
	
	clear();
}

Concave::Concave(const Concave &concave) {
	
	clear();
	
	error = concave.getError();
	threshold = concave.getThreshold();
	convexes.resize(concave.getNumConvexes());
	for(int i = 0; i < concave.getNumConvexes(); i++) {
		convexes[i] = new Convex();
		*convexes[i] = *concave.getConvex(i);
	}
}

Concave::~Concave() {
	
	clear();
}

/*
 */
Concave &Concave::operator=(const Concave &concave) {
	
	if(this == &concave) return *this;
	
	clear();
	
	error = concave.getError();
	threshold = concave.getThreshold();
	convexes.resize(concave.getNumConvexes());
	for(int i = 0; i < concave.getNumConvexes(); i++) {
		convexes[i] = new Convex(*concave.getConvex(i));
	}
	
	return *this;
}

/*
 */
void Concave::clear() {
	
	error = 0.0f;
	threshold = 0.0f;
	for(int i = 0; i < convexes.size(); i++) {
		delete convexes[i];
	}
	convexes.clear();
}

/******************************************************************************\
*
* Create concave
*
\******************************************************************************/

/*
 */
int Concave::create(const vec3 *vertex,int num_vertex,int depth,float e,float t) {
	
	clear();
	
	// parameters
	error = e;
	threshold = clamp(t,EPSILON,1.0f);
	
	// create convex
	Convex *convex = new Convex();
	if(convex->create(vertex,num_vertex,error) == 0) {
		Log::error("Concave::create(): can't create convex\n");
		delete convex;
		return 0;
	}
	
	// convex decomposition
	decompose(vertex,num_vertex,convex,depth);
	
	// convex merging
	while(merge());
	
	return 1;
}

/******************************************************************************\
*
* Convex decomposition
*
\******************************************************************************/

/*
 */
double Concave::get_concavity(const vec3 *vertex,int num_vertex,Convex *convex) const {
	
	double concavity = 0.0;
	
	// convex face planes
	VectorStack<dvec4> planes(convex->getNumFaces());
	for(int i = 0; i < convex->getNumFaces(); i++) {
		const dvec3 &v0 = convex->getVertex(convex->getFaceVertex(i,0));
		const dvec3 &v1 = convex->getVertex(convex->getFaceVertex(i,1));
		const dvec3 &v2 = convex->getVertex(convex->getFaceVertex(i,2));
		dvec3 normal = normalize(cross(v1 - v0,v2 - v0));
		planes[i] = dvec4(normal,-dot(normal,v0));
	}
	
	// nearest distance to convex face
	for(int i = 0; i < num_vertex; i += 3) {
		const vec3 &v0 = vertex[i + 0];
		const vec3 &v1 = vertex[i + 1];
		const vec3 &v2 = vertex[i + 2];
		vec3 v = (v0 + v1 + v2) / 3.0f;
		double volume = 0.0;
		double distance = -INFINITY;
		for(int j = 0; j < planes.size(); j++) {
			double d = dot(planes[j],dvec3(v));
			if(distance < d) {
				distance = d;
				volume = distance * length(cross(v1 - v0,v2 - v0)) / 2.0;
			}
		}
		concavity -= volume;
	}
	
	return concavity;
}

/*
 */
vec4 Concave::get_split_plane(const vec3 *vertex,int num_vertex,Convex *convex) const {
	
	// convex vertices
	VectorStack<vec3> src(convex->getNumVertex());
	VectorStack<vec3> dest(convex->getNumVertex());
	for(int i = 0; i < convex->getNumVertex(); i++) {
		src[i] = vec3(convex->getVertex(i));
	}
	
	// decompose convex inertia
	mat3 eigenvector;
	jacobi(convex->getInertia(),eigenvector);
	
	// calculate bounding box
	Simd::mulMat4Vec3(dest.get(),sizeof(vec3),mat4(eigenvector),src.get(),sizeof(vec3),src.size());
	BoundBox bound_box(dest.get(),dest.size());
	
	// choose splitting plane
	vec3 normal = vec3_zero;
	vec3 size = bound_box.getMax() - bound_box.getMin();
	if(size.x > size.y) {
		if(size.x > size.z) normal = eigenvector.getRow(0);
		else normal = eigenvector.getRow(2);
	} else {
		if(size.y > size.z) normal = eigenvector.getRow(1);
		else normal = eigenvector.getRow(2);
	}
	
	return vec4(normal,-dot(normal,vec3(convex->getCenter())));
}

/*
 */
void Concave::split(const vec3 *vertex,int num_vertex,const vec4 &plane,Vector<vec3> &left_vertex,Vector<vec3> &right_vertex) const {
	
	for(int i = 0; i < num_vertex; i += 3) {
		
		const vec3 &v0 = vertex[i + 0];
		const vec3 &v1 = vertex[i + 1];
		const vec3 &v2 = vertex[i + 2];
		
		float d0 = dot(plane,v0);
		float d1 = dot(plane,v1);
		float d2 = dot(plane,v2);
		
		// left triangle
		if(d0 <= 0.0f && d1 <= 0.0f && d2 <= 0.0f) {
			left_vertex.append(v0);
			left_vertex.append(v1);
			left_vertex.append(v2);
		}
		// right triangle
		else if(d0 >= 0.0f && d1 >= 0.0f && d2 >= 0.0f) {
			right_vertex.append(v0);
			right_vertex.append(v1);
			right_vertex.append(v2);
		}
		// splitted triangle
		else {
			VectorStack<vec3> left;
			VectorStack<vec3> right;
			const vec3 *v0 = &vertex[i];
			float d0 = dot(plane,*v0);
			for(int j = 0; j < 3; j++) {
				const vec3 *v1 = &vertex[i + (j + 1) % 3];
				float d1 = dot(plane,*v1);
				if(d0 <= 0.0f) left.append(*v0);
				if(d0 >= 0.0f) right.append(*v0);
				if((d0 > 0.0f) != (d1 > 0.0f)) {
					float k = -d0 / (d1 - d0);
					left.append(lerp(*v0,*v1,k));
					right.append(lerp(*v0,*v1,k));
				}
				v0 = v1;
				d0 = d1;
			}
			for(int j = 2; j < left.size(); j++) {
				left_vertex.append(left[0]);
				left_vertex.append(left[j - 1]);
				left_vertex.append(left[j]);
			}
			for(int j = 2; j < right.size(); j++) {
				right_vertex.append(right[0]);
				right_vertex.append(right[j - 1]);
				right_vertex.append(right[j]);
			}
		}
	}
}

/*
 */
int Concave::decompose(const vec3 *vertex,int num_vertex,Convex *convex,int depth) {
	
	// recursion depth
	if(depth <= 0) {
		convexes.append(convex);
		return 1;
	}
	
	// convex volume
	double volume = convex->getVolume();
	if(volume < 0.0) {
		Log::error("Concave::decompose(): bad convex volume %f\n",volume);
		return 0;
	}
	
	// mesh concavity
	double concavity = get_concavity(vertex,num_vertex,convex);
	if(concavity < volume * threshold) {
		convexes.append(convex);
		return 1;
	}
	
	// splitting plane
	vec4 plane = get_split_plane(vertex,num_vertex,convex);
	
	// split mesh
	VectorStack<vec3> left_vertex;
	VectorStack<vec3> right_vertex;
	Convex *left_convex = new Convex();
	Convex *right_convex = new Convex();
	split(vertex,num_vertex,plane,left_vertex,right_vertex);
	
	// check split plane
	if(left_vertex.size() < 4 || right_vertex.size() < 4) {
		Log::error("Concave::decompose(): can't split mesh\n");
		convexes.append(convex);
		return 1;
	}
	
	// create convexes
	if(left_convex->create(left_vertex.get(),left_vertex.size(),error) == 0 ||
		right_convex->create(right_vertex.get(),right_vertex.size(),error) == 0) {
		Log::error("Concave::decompose(): can't create convex\n");
		convexes.append(convex);
		return 1;
	}
	
	// recursion
	decompose(left_vertex.get(),left_vertex.size(),left_convex,depth - 1);
	decompose(right_vertex.get(),right_vertex.size(),right_convex,depth - 1);
	
	// delete input convex
	delete convex;
	
	return 1;
}

/*
 */
struct ConcaveConvexCompare {
	INLINE int operator()(const Convex *c0,const Convex *c1) {
		return (c0->getVolume() > c1->getVolume());
	}
};

/*
 */
int Concave::merge() {
	
	// sort convexes by volume
	ConcaveConvexCompare convex_compare;
	quickSort(convexes.get(),convexes.size(),convex_compare);
	
	// cache convex parameters
	dvec3 min,max;
	VectorStack<double> volumes(convexes.size());
	VectorStack<BoundBox> bounds(convexes.size());
	for(int i = 0; i < convexes.size(); i++) {
		volumes[i] = convexes[i]->getVolume();
		convexes[i]->getBoundBox(min,max);
		bounds[i].set(vec3(min),vec3(max));
	}
	
	// find merge candidats
	int ret = 0;
	Convex *convex = NULL;
	VectorStack<dvec3> vertex;
	for(int i = 0; i < convexes.size(); i++) {
		for(int j = i + 1; j < convexes.size(); j++) {
			if(i == j) continue;
			if(bounds[i].inside(bounds[j]) == 0) continue;
			vertex.clear();
			vertex.append(convexes[i]->getVertex(),convexes[i]->getNumVertex());
			vertex.append(convexes[j]->getVertex(),convexes[j]->getNumVertex());
			if(convex == NULL) convex = new Convex();
			if(convex->create(vertex.get(),vertex.size(),error) == 0) continue;
			double volume = convex->getVolume();
			if(Math::abs(volume - volumes[i] - volumes[j]) > volume * threshold) continue;
			convexes.remove(j);
			convexes.remove(i);
			volumes.remove(j);
			volumes.remove(i);
			bounds.remove(j);
			bounds.remove(i);
			convexes.append(convex);
			volumes.append(convex->getVolume());
			convex->getBoundBox(min,max);
			bounds.append().set(vec3(min),vec3(max));
			convex = NULL;
			ret = 1;
		}
	}
	
	delete convex;
	
	return ret;
}
