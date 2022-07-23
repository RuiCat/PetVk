/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Concave.h
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

#ifndef __CONCAVE_H__
#define __CONCAVE_H__

#include "Vector.h"
#include "MathLib.h"

/*
 */
class Convex;

/*
 */
class Concave {
		
	public:
		
		Concave();
		Concave(const Concave &concave);
		~Concave();
		
		Concave &operator=(const Concave &concave);
		
		// clear
		void clear();
		
		// create concave
		int create(const vec3 *vertex,int num_vertex,int depth = 4,float error = 0.01f,float threshold = 0.01f);
		
		// parameters
		INLINE float getError() const { return error; }
		INLINE float getThreshold() const { return threshold; }
		
		// convexes
		INLINE int getNumConvexes() const { return convexes.size(); }
		INLINE const Convex *getConvex(int num) const { return convexes[num]; }
		
	private:
		
		double get_concavity(const vec3 *vertex,int num_vertex,Convex *convex) const;
		
		vec4 get_split_plane(const vec3 *vertex,int num_vertex,Convex *convex) const;
		
		void split(const vec3 *vertex,int num_vertex,const vec4 &plane,Vector<vec3> &left_vertex,Vector<vec3> &right_vertex) const;
		
		int decompose(const vec3 *vertex,int num_vertex,Convex *convex,int depth);
		
		int merge();
		
		float error;						// convex error
		float threshold;					// merge threshold
		Vector<Convex*> convexes;			// convex hulls
};

#endif /* __CONCAVE_H__ */
