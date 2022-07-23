/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Polygon.h
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

#ifndef __POLYGON_H__
#define __POLYGON_H__

#include "Vector.h"
#include "MathLib.h"

/*
 */
class Polygon {
		
	public:
		
		Polygon();
		~Polygon();
		
		// clear
		void clear();
		
		// create polygon
		int createConvex(const vec3 *vertex,int num_vertex,const vec3 &normal);
		int createConcave(const vec3 *vertex,int num_vertex,const vec3 &normal);
		
		// vertices
		INLINE int getNumVertex() const { return vertex.size(); }
		INLINE const vec3 *getVertex() const { return vertex.get(); }
		INLINE const vec3 &getVertex(int num) const { return vertex[num]; }
		
		// indices
		INLINE int getNumIndices() const { return indices.size(); }
		INLINE const short *getIndices() const { return indices.get(); }
		INLINE short getIndex(int num) const { return indices[num]; }
		
	private:
		
		enum {
			NUM_VERTEX = 64,
			NUM_INDICES = 128,
		};
		
		VectorStack<vec3,NUM_VERTEX> vertex;
		VectorStack<short,NUM_INDICES> indices;
};

#endif /* __POLYGON_H__ */
