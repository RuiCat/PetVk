/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Geometry.spu.h
 * Desc:    Geometry spu utils
 * Version: 1.01
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

#ifndef __GEOMETRY_SPU_H__
#define __GEOMETRY_SPU_H__

#include "MathLib.spu.h"

/******************************************************************************\
*
* Orthonormal basis
*
\******************************************************************************/

/*
 */
INLINE void orthoBasis(vec_float4 v,vec_float4 &tangent,vec_float4 &binormal) {
	vec_uchar16 yzxw = SPU_PERM2(Y,Z,X,W);
	vec_uchar16 zxyw = SPU_PERM2(Z,X,Y,W);
	vec_uint4 mask = spu_cmpabsgt(v,spu_splats(0.7f));
	mask = spu_splats(spu_extract(mask,2));
	vec_float4 v_yzxw = spu_shuffle(v,v,yzxw);
	vec_float4 v_zxyw = spu_shuffle(v,v,zxyw);
	vec_float4 b_yzxw = spu_sel(SPU_FLOAT4(0.0f,1.0f,0.0f,0.0f),SPU_FLOAT4(0.0f,0.0f,1.0f,0.0f),mask);
	vec_float4 b_zxyw = spu_sel(SPU_FLOAT4(1.0f,0.0f,0.0f,0.0f),SPU_FLOAT4(0.0f,1.0f,0.0f,0.0f),mask);
	tangent = spu_nmsub(b_zxyw,v_yzxw,spu_mul(b_yzxw,v_zxyw));
	vec_float4 tangent_yzxw = spu_shuffle(tangent,tangent,yzxw);
	vec_float4 tangent_zxyw = spu_shuffle(tangent,tangent,zxyw);
	binormal = spu_nmsub(v_zxyw,tangent_yzxw,spu_mul(v_yzxw,tangent_zxyw));
	tangent = spu_normalize3(tangent);
	binormal = spu_normalize3(binormal);
}

#endif /* __GEOMETRY_SPU_H__ */
