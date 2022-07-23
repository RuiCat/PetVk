/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Convex.h
 * Desc:    Convex creator
 * Version: 1.13
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

#ifndef __CONVEX_H__
#define __CONVEX_H__

#include "Vector.h"
#include "MathLib.h"

/*
 */
class Convex {
		
	public:
		
		Convex();
		Convex(const Convex &convex);
		~Convex();
		
		Convex &operator=(const Convex &convex);
		
		// clear
		void clear();
		
		// create convex
		int create(const vec3 *vertex,int num_vertex,float error = 0.01f);
		int create(const dvec3 *vertex,int num_vertex,float error = 0.01f);
		
		// vertices
		INLINE int getNumVertex() const { return vertex.size(); }
		INLINE const dvec3 *getVertex() const { return vertex.get(); }
		INLINE const dvec3 &getVertex(int num) const { return vertex[num]; }
		
		// faces
		INLINE int getNumFaces() const { return faces.size(); }
		INLINE int getNumFaceVertex(int face) const { return faces[face].vertex.size(); }
		INLINE int getFaceVertex(int face,int num) const { return faces[face].vertex[num]; }
		
		// parameters
		INLINE double getVolume() const { return volume; }
		INLINE double getThreshold() const { return threshold; }
		dvec3 getCenter() const;
		mat3 getInertia() const;
		
		// closest point
		dvec3 getClosestPoint(const dvec3 &point) const;
		
		// bounding box
		void getBoundBox(dvec3 &min,dvec3 &max) const;
		
	private:
		
		struct Triangle;
		
		static int vertex_compare(const dvec3 &v0,const dvec3 &v1);
		static int triangle_compare(const Triangle *t0,const Triangle *t1);
		
		void add_vertex(Triangle *t,int v);
		Triangle *add_triangle(int v0,int v1,int v2);
		
		int create_edges(const Vector<Triangle*> &triangles);
		int create_convex();
		int create_vertex();
		int create_faces();
		
		struct Vertex {
			double distance;				// vertex distance
			int id;							// vertex number
		};
		
		struct Edge {
			int v[2];						// edge vertices
			int counter;					// edge counter
		};
		
		struct Triangle {
			int v[3];						// triangle vertices
			int e[3];						// triangle edges
			dvec4 plane;					// triangle plane
			VectorStack<Vertex,16> vertex;	// separate vertices
		};
		
		struct Face {
			VectorStack<int,16> vertex;		// face vertices
		};
		
		double volume;						// convex volume
		double threshold;					// distance threshold
		Vector<dvec3> vertex;				// convex vertices
		Vector<Edge> edges;					// convex edges
		Vector<Triangle*> triangles;		// convex triangles
		Vector<Face> faces;					// convex faces
};

#endif /* __CONVEX_H__ */
