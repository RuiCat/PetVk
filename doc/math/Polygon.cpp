/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Polygon.cpp
 * Desc:    Polygon creator
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

#include "Geometry.h"
#include "Polygon.h"

/*
 */
Polygon::Polygon() {
	
}

Polygon::~Polygon() {
	
}

/*
 */
void Polygon::clear() {
	vertex.clear();
	indices.clear();
}

/******************************************************************************\
*
* Create convex
*
\******************************************************************************/

/*
 */
int Polygon::createConvex(const vec3 *v,int num_vertex,const vec3 &normal) {
	
	clear();
	
	if(num_vertex < 3) return 0;
	
	// vertex buffer
	VectorStack<vec3,NUM_VERTEX> buffer(num_vertex);
	for(int i = 0; i < num_vertex; i++) {
		buffer[i] = v[i];
	}
	
	// first vertex
	int index = 0;
	for(int i = 1; i < buffer.size(); i++) {
		if(buffer[index].x > buffer[i].x) index = i;
	}
	vertex.append(buffer[index]);
	buffer.removeFast(index);
	
	// create convex polygon
	for(int i = 0; i < buffer.size(); i++) {
		for(int j = 0; j < buffer.size(); j++) {
			vec3 left = cross(buffer[j] - vertex[vertex.size() - 1],normal);
			vec4 plane = vec4(left,-dot(left,buffer[j]));
			int is_convex = 1;
			for(int k = 0; k < buffer.size(); k++) {
				if(j != k && dot(plane,buffer[k]) < 0.0f) {
					is_convex = 0;
					break;
				}
			}
			if(is_convex) {
				vertex.append(buffer[j]);
				buffer.removeFast(j--);
				i--;
			}
		}
	}
	
	// create polygon indices
	for(int i = 2; i < vertex.size(); i++) {
		indices.append(0);
		indices.append(i - 1);
		indices.append(i + 0);
	}
	
	return 1;
}

/******************************************************************************\
*
* Create concave
*
\******************************************************************************/

/*
 */
int Polygon::createConcave(const vec3 *v,int num_vertex,const vec3 &normal) {
	
	clear();
	
	if(num_vertex < 3) return 0;
	
	// copy vertices
	vertex.resize(num_vertex);
	for(int i = 0; i < num_vertex; i++) {
		vertex[i] = v[i];
	}
	
	// polygon area
	vec3 area = cross(vertex[0],vertex[num_vertex - 1]);
	for(int i = 1; i < num_vertex; i++) {
		area += cross(vertex[i],vertex[i - 1]);
	}
	
	// indices buffer
	VectorStack<short,NUM_INDICES> buffer(num_vertex);
	if(dot(area,normal) < 0.0f) {
		for(int i = 0; i < num_vertex; i++) {
			buffer[i] = num_vertex - i - 1;
		}
	} else {
		for(int i = 0; i < num_vertex; i++) {
			buffer[i] = i;
		}
	}
	
	// create indices
	int i1 = buffer.size() - 1;
	int iterations = buffer.size() * 2;
	while(iterations-- > 0 && buffer.size() > 2) {
		
		int i0,i2;
		i0 = i1; if(i0 >= buffer.size()) i0 = 0;
		i1 = i0 + 1; if(i1 >= buffer.size()) i1 = 0;
		i2 = i1 + 1; if(i2 >= buffer.size()) i2 = 0;
		
		const vec3 &v0 = vertex[buffer[i0]];
		const vec3 &v1 = vertex[buffer[i1]];
		const vec3 &v2 = vertex[buffer[i2]];
		
		if(dot(cross(v2 - v0,v1 - v0),normal) < 0.0f) continue;
		
		int inside = 0;
		for(int j = 0; j < buffer.size(); j++) {
			if(j == i0 || j == i1 || j == i2) continue;
			if(pointTriangleInside(vertex[buffer[j]],v0,v1,v2,normal)) {
				inside = 1;
				break;
			}
		}
		if(inside) continue;
	  	
		indices.append(buffer[i0]);
		indices.append(buffer[i1]);
		indices.append(buffer[i2]);
		
		buffer.remove(i1);
		
		iterations = buffer.size() * 2;
	}
	
	return (iterations > 0);
}
