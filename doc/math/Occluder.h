/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Occluder.h
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

#ifndef __OCCLUDER_H__
#define __OCCLUDER_H__

#include "MathLib.h"

/*
 */
class BoundSphere;
class BoundBox;

/*
 */
class Occluder {
		
	public:
		
		Occluder();
		Occluder(const mat4 &projection,const mat4 &modelview);
		Occluder(const Occluder &occluder);
		~Occluder();
		
		Occluder &operator=(const Occluder &occluder);
		
		// swap occluders
		void swap(Occluder &occluder);
		
		// clear occluder
		void clear();
		
		// occluder matrices
		void set(const mat4 &projection,const mat4 &modelview);
		
		// occluder transformation
		void setITransform(const mat4 &itransform);
		
		// check occluder
		int hasOccluders() const;
		
		// add occluder
		void addOccluder();
		
		// render occluder
		void renderOccluder(const vec3 *vertex,int num_vertex,const unsigned short *indices,int num_indices);
		
		// inside points
		int inside(const vec3 &min,const vec3 &max) const;
		
		int insideAll(const vec3 &min,const vec3 &max) const;
		
		// inside bounds
		int inside(const BoundSphere &bs) const;
		int inside(const BoundBox &bb) const;
		
		int insideValid(const BoundSphere &bs) const;
		int insideValid(const BoundBox &bb) const;
		
		int insideAll(const BoundSphere &bs) const;
		int insideAll(const BoundBox &bb) const;
		
		int insideAllValid(const BoundSphere &bs) const;
		int insideAllValid(const BoundBox &bb) const;
		
		// parameters
		enum {
			OCCLUDER_WIDTH = 128,
			OCCLUDER_HEIGHT = 64,
			OCCLUDER_SHIFT = 7,
		};
		
		INLINE float *getData() const { return data; }
		INLINE const mat4 &getModelviewProjection() const { return tmodelviewprojection; }
		
	private:
		
		enum {
			PLANE_L = 0,
			PLANE_R,
			PLANE_B,
			PLANE_T,
			PLANE_N,
			NUM_PLANES,
		};
		
		// project bounding box
		int project_bound_box(const vec3 &min,const vec3 &max,int &x0,int &y0,int &x1,int &y1,float &z) const;
		
		// clip triangle
		int clip_triangle(const vec3 &v0,const vec3 &v1,const vec3 &v2,int plane) const;
		
		// render triangle
		int render_triangle(const vec3 *v0,const vec3 *v1,const vec3 *v2) const;
		
		vec4 planes[NUM_PLANES];		// clipping planes
		vec4 tplanes[NUM_PLANES];		// transformed clipping planes
		
		mat4 modelviewprojection;		// modelviewprojection matrix
		mat4 tmodelviewprojection;		// transformed modelviewprojection matrix
		
		int has_occluders;				// has occluders flag
		int need_clear;					// need clear flag
		
		int num_vertex;					// number of vertices
		vec4 *vertex;					// vertices buffer
		
		float *data;					// depth buffer
};

#endif /* __OCCLUDER_H__ */
