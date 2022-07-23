/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Scissor.h
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

#ifndef __SCISSOR_H__
#define __SCISSOR_H__

#include "MathLib.h"

/*
 */
class Scissor {
		
	public:
		
		Scissor();
		~Scissor();
		
		void clear();
		
		void set(float left,float right,float bottom,float top);
		void set(const mat4 &projection,const mat4 &modelview);
		
		// portals
		int addPortal(const vec3 *points,int num_points,const mat4 &transform);
		
		// lights
		int addLight(float radius,const mat4 &transform);
		int addLight(const vec3 &radius,const mat4 &transform);
		int addLight(const mat4 &modelviewprojection);
		
		// bounds
		float getX() const;
		float getY() const;
		float getWidth() const;
		float getHeight() const;
		float getZNear() const;
		float getZFar() const;
		
	private:
		
		// create scissor from bounds
		int set(const vec3 &min,const vec3 &max);
		
		mat4 projection;			// matrices
		mat4 modelview;
		mat4 modelviewprojection;
		
		float near_clip_plane;		// distance to near clipping plane
		float far_clip_plane;		// distance to far clipping plane
		
		float left;					// scissor
		float right;
		float bottom;
		float top;
		float znear;
		float zfar;
};

#endif /* __SCISSOR_H__ */
