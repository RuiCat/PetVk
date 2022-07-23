/* Copyright (C) 2005-2013, Unigine Corp. All rights reserved.
 *
 * File:    Numerical.h
 * Desc:    Numerical library
 * Version: 1.11
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

#ifndef __NUMERICAL_H__
#define __NUMERICAL_H__

#include "MathLib.h"
#include "SimdLib.h"

/*
 */
template <int Size> class vecX;
template <int Rows,int Columns> class matX;

/******************************************************************************\
*
* vecX
*
\******************************************************************************/

/*
 */
template <int Size> class vecX {
		
	public:
		
		vecX() {
			
		}
		explicit vecX(float v) {
			set(v);
		}
		explicit vecX(const float *v) {
			set(v);
		}
		vecX(const vecX<Size> &v) {
			set(v.get());
		}
		~vecX() {
			
		}
		
		vecX<Size> &operator=(const vecX<Size> &v) {
			if(this == &v) return *this;
			set(v.get());
			return *this;
		}
		
		vecX<Size> operator-() const {
			vecX<Size> ret;
			for(int i = 0; i < Size; i++) {
				ret.data[i] = -data[i];
			}
			return ret;
		}
		vecX<Size> &operator*=(float v) {
			return mul(*this,*this,v);
		}
		vecX<Size> &operator*=(const vecX<Size> &v) {
			return mul(*this,this,v);
		}
		vecX<Size> &operator+=(const vecX<Size> &v) {
			return add(*this,*this,v);
		}
		vecX<Size> &operator-=(const vecX<Size> &v) {
			return sub(*this,*this,v);
		}
		
		INLINE float &operator[](int index) {
			assert((unsigned int)index < (unsigned int)Size && "vecX::operator[](): bad index");
			return data[index];
		}
		INLINE float operator[](int index) const {
			assert((unsigned int)index < (unsigned int)Size && "vecX::operator[](): bad index");
			return data[index];
		}
		
		void set(float v) {
			for(int i = 0; i < Size; i++) {
				data[i] = v;
			}
		}
		void set(const float *v) {
			Math::memcpy(data,v,sizeof(float) * Size);
		}
		void get(float *v) const {
			Math::memcpy(v,data,sizeof(float) * Size);
		}
		
		INLINE float *get() { return data; }
		INLINE const float *get() const { return data; }
		
		INLINE int getSize() const { return Size; }
		
	private:
		
		ATTRIBUTE_ALIGNED16(float data[Size]);
};

/*
 */
template <int Size> int operator==(const vecX<Size> &v0,const vecX<Size> &v1) {
	return compare(v0,v1);
}

template <int Size> int operator!=(const vecX<Size> &v0,const vecX<Size> &v1) {
	return !compare(v0,v1);
}

template <int Size> vecX<Size> operator*(const vecX<Size> &v0,float v1) {
	vecX<Size> ret;
	return mul(ret,v0,v1);
}

template <int Size> vecX<Size> operator*(const vecX<Size> &v0,const vecX<Size> &v1) {
	vecX<Size> ret;
	return mul(ret,v0,v1);
}

template <int Size> vecX<Size> operator+(const vecX<Size> &v0,const vecX<Size> &v1) {
	vecX<Size> ret;
	return add(ret,v0,v1);
}

template <int Size> vecX<Size> operator-(const vecX<Size> &v0,const vecX<Size> &v1) {
	vecX<Size> ret;
	return sub(ret,v0,v1);
}

/*
 */
template <int Size> int compare(const vecX<Size> &v0,const vecX<Size> &v1) {
	for(int i = 0; i < Size; i++) {
		if(compare(v0[i],v1[i]) == 0) return 0;
	}
	return 0;
}

template <int Size> int compare(const vecX<Size> &v0,const vecX<Size> &v1,float epsilon) {
	for(int i = 0; i < Size; i++) {
		if(compare(v0[i],v1[i],epsilon) == 0) return 0;
	}
	return 0;
}

template <int Size> float dot(const vecX<Size> &v0,const vecX<Size> &v1) {
	float ret = 0.0;
	for(int i = 0; i < Size; i++) {
		ret += v0[i] * v1[i];
	}
	return ret;
}

template <int Size> vecX<Size> &mul(vecX<Size> &ret,const vecX<Size> &v0,float v1) {
	Simd::mul(ret.get(),v0.get(),v1,Size);
	return ret;
}

template <int Size> vecX<Size> &mul(vecX<Size> &ret,const vecX<Size> &v0,const vecX<Size> &v1) {
	Simd::mul(ret.get(),v0.get(),v1.get(),Size);
	return ret;
}

template <int Size> vecX<Size> &mad(vecX<Size> &ret,const vecX<Size> &v0,float v1,const vecX<Size> &v2) {
	Simd::mad(ret.get(),v0.get(),v1,v2.get(),Size);
	return ret;
}

template <int Size> vecX<Size> &add(vecX<Size> &ret,const vecX<Size> &v0,const vecX<Size> &v1) {
	Simd::add(ret.get(),v0.get(),v1.get(),Size);
	return ret;
}

template <int Size> vecX<Size> &sub(vecX<Size> &ret,const vecX<Size> &v0,const vecX<Size> &v1) {
	Simd::sub(ret.get(),v0.get(),v1.get(),Size);
	return ret;
}

template <int Size> vecX<Size> &clamp(vecX<Size> &ret,const vecX<Size> &v,float v0,float v1) {
	for(int i = 0; i < Size; i++) {
		ret[i] = clamp(v[i],v0,v1);
	}
	return ret;
}

template <int Size> vecX<Size> &clamp(vecX<Size> &ret,const vecX<Size> &v,const vecX<Size> &v0,const vecX<Size> &v1) {
	for(int i = 0; i < Size; i++) {
		ret[i] = clamp(v[i],v0[i],v1[i]);
	}
	return ret;
}

template <int Size> vecX<Size> &lerp(vecX<Size> &ret,const vecX<Size> &v0,const vecX<Size> &v1,float k) {
	for(int i = 0; i < Size; i++) {
		ret[i] = lerp(v0[i],v1[i],k);
	}
	return ret;
}

/*
 */
template <int Size> vecX<Size> clamp(const vecX<Size> &v,float v0,float v1) {
	vecX<Size> ret;
	return clamp(ret,v,v0,v1);
}

template <int Size> vecX<Size> clamp(const vecX<Size> &v,const vecX<Size> &v0,const vecX<Size> &v1) {
	vecX<Size> ret;
	return clamp(ret,v,v0,v1);
}

template <int Size> vecX<Size> lerp(const vecX<Size> &v0,const vecX<Size> &v1,float k) {
	vecX<Size> ret;
	return lerp(ret,v0,v1,k);
}

/*
 */
int biSolve(float *ret,const float *c,int size,float iepsilon = 1e6f);

/*
 */
template <int Size> int biSolve(vecX<Size> &ret,const vecX<Size> &c,float iepsilon = 1e6f) {
	return biSolve(ret.get(),c.get(),Size,iepsilon);
}

/******************************************************************************\
*
* matX
*
\******************************************************************************/

/*
 */
template <int Rows,int Columns> class matX {
		
	public:
		
		matX() {
			
		}
		explicit matX(float v) {
			set(v);
		}
		explicit matX(const float *m) {
			set(m);
		}
		matX(const matX<Rows,Columns> &m) {
			set(m.get());
		}
		~matX() {
			
		}
		
		matX<Rows,Columns> &operator=(const matX<Rows,Columns> &m) {
			if(this == &m) return *this;
			set(m.get());
			return *this;
		}
		
		matX<Rows,Columns> operator-() const {
			matX<Rows,Columns> ret;
			for(int i = 0; i < Size; i++) {
				ret.data[i] = -data[i];
			}
			return ret;
		}
		matX<Rows,Columns> &operator*=(float v) {
			return mul(*this,*this,v);
		}
		matX<Rows,Columns> &operator*=(const matX<Rows,Columns> &m) {
			return mul(*this,matX<Rows,Columns>(*this),m);
		}
		matX<Rows,Columns> &operator+=(const matX<Rows,Columns> &m) {
			return add(*this,*this,m);
		}
		matX<Rows,Columns> &operator-=(const matX<Rows,Columns> &m) {
			return sub(*this,*this,m);
		}
		
		INLINE float &operator[](int index) {
			assert((unsigned int)index < (unsigned int)Size && "matX::operator[](): bad index");
			return data[index];
		}
		INLINE float operator[](int index) const {
			assert((unsigned int)index < (unsigned int)Size && "matX::operator[](): bad index");
			return data[index];
		}
		
		INLINE void set(int row,int column,float v) {
			assert((unsigned int)row < (unsigned int)Rows && "matX::set(): bad row");
			assert((unsigned int)column < (unsigned int)Columns && "matX::set(): bad column");
			data[Rows * column + row] = v;
		}
		INLINE float &get(int row,int column) {
			assert((unsigned int)row < (unsigned int)Rows && "matX::get(): bad row");
			assert((unsigned int)column < (unsigned int)Columns && "matX::get(): bad column");
			return data[Rows * column + row];
		}
		INLINE float get(int row,int column) const {
			assert((unsigned int)row < (unsigned int)Rows && "matX::get(): bad row");
			assert((unsigned int)column < (unsigned int)Columns && "matX::get(): bad column");
			return data[Rows * column + row];
		}
		
		void set(float v) {
			for(int i = 0; i < Size; i++) {
				data[i] = v;
			}
		}
		void set(const float *m) {
			Math::memcpy(data,m,sizeof(float) * Size);
		}
		void get(float *m) const {
			Math::memcpy(m,data,sizeof(float) * Size);
		}
		
		INLINE float *get() { return data; }
		INLINE const float *get() const { return data; }
		
		INLINE int getSize() const { return Size; }
		INLINE int getRows() const { return Rows; }
		INLINE int getColumns() const { return Columns; }
		
		void setRow(int row,const vecX<Columns> &v) {
			assert((unsigned int)row < (unsigned int)Rows && "matX::setRow(): bad row");
			float *d = data + row;
			for(int i = 0; i < Columns; i++) {
				*d = v[i];
				d += Rows;
			}
		}
		vecX<Columns> &getRow(vecX<Columns> &ret,int row) const {
			assert((unsigned int)row < (unsigned int)Rows && "matX::getRow(): bad row");
			const float *s = data + row;
			for(int i = 0; i < Columns; i++) {
				ret[i] = *s;
				s += Rows;
			}
			return ret;
		}
		vecX<Columns> getRow(int row) const {
			vecX<Columns> ret;
			return getRow(ret,row);
		}
		
		void setColumn(int column,const vecX<Rows> &v) {
			assert((unsigned int)column < (unsigned int)Columns && "matX::setColumn(): bad column");
			Math::memcpy(data + Rows * column,v.get(),sizeof(float) * Rows);
		}
		vecX<Rows> &getColumn(vecX<Rows> &ret,int column) const {
			assert((unsigned int)column < (unsigned int)Columns && "matX::getColumn(): bad column");
			Math::memcpy(ret.get(),data + Rows * column,sizeof(float) * Rows);
			return ret;
		}
		vecX<Rows> getColumn(int column) const {
			vecX<Rows> ret;
			return getColumn(ret,column);
		}
		
		void swapRows(int row_0,int row_1) {
			assert((unsigned int)row_0 < (unsigned int)Rows && "matX::swapRows(): bad row");
			assert((unsigned int)row_1 < (unsigned int)Rows && "matX::swapRows(): bad row");
			float *d0 = data + row_0;
			float *d1 = data + row_1;
			for(int i = 0; i < Columns; i++) {
				float d = *d0;
				*d0 = *d1;
				*d1 = d;
				d0 += Rows;
				d1 += Rows;
			}
		}
		void swapColumns(int column_0,int column_1) {
			assert((unsigned int)column_0 < (unsigned int)Columns && "matX::swapColumns(): bad column");
			assert((unsigned int)column_1 < (unsigned int)Columns && "matX::swapColumns(): bad column");
			float *d0 = data + Rows * column_0;
			float *d1 = data + Rows * column_1;
			for(int i = 0; i < Rows; i++) {
				float d = *d0;
				*d0 = *d1;
				*d1 = d;
				d0++;
				d1++;
			}
		}
		
		void setBlock(int row,int column,const mat3 &m) {
			assert(row >= 0 && row + 2 < Rows && "matX::setBlock(): bad row");
			assert(column >= 0 && column + 2 < Columns && "matX::setBlock(): bad column");
			float *d = data + Rows * column + row;
			const float *s = m.get();
			for(int i = 0; i < 3; i++) {
				*d++ = *s++;
				*d++ = *s++;
				*d++ = *s++;
				d += Rows - 3;
				s++;
			}
		}
		mat3 &getBlock(mat3 &ret,int row,int column) const {
			assert(row >= 0 && row + 2 < Rows && "matX::getBlock(): bad row");
			assert(column >= 0 && column + 2 < Columns && "matX::getBlock(): bad column");
			float *d = ret.get();
			const float *s = data + Rows * column + row;
			for(int i = 0; i < 3; i++) {
				*d++ = *s++;
				*d++ = *s++;
				*d++ = *s++;
				s += Rows - 3;
				d++;
			}
			return ret;
		}
		
		void setBlock(int row,int column,const mat3 &m,int rows,int columns) {
			assert(row >= 0 && row + rows - 1 < Rows && "matX::setBlock(): bad row");
			assert(column >= 0 && column + columns - 1 < Columns && "matX::setBlock(): bad column");
			assert(rows > 0 && rows < 4 && "matX::setBlock(): bad rows");
			assert(columns > 0 && columns < 4 && "matX::setBlock(): bad columns");
			float *d = data + Rows * column + row;
			const float *s = m.get();
			for(int i = 0; i < columns; i++) {
				for(int j = 0; j < rows; j++) {
					*d++ = *s++;
				}
				d += Rows - rows;
				s += 4 - rows;
			}
		}
		mat3 &getBlock(mat3 &ret,int row,int column,int rows,int columns) const {
			assert(row >= 0 && row + rows - 1 < Rows && "matX::getBlock(): bad row");
			assert(column >= 0 && column + columns - 1 < Columns && "matX::getBlock(): bad column");
			assert(rows > 0 && rows < 4 && "matX::getBlock(): bad rows");
			assert(columns > 0 && columns < 4 && "matX::getBlock(): bad columns");
			float *d = ret.get();
			const float *s = data + Rows * column + row;
			for(int i = 0; i < columns; i++) {
				for(int j = 0; j < rows; j++) {
					*d++ = *s++;
				}
				d += 4 - rows;
				s += Rows - rows;
			}
			return ret;
		}
		
		void setZero() {
			Math::memset(data,0,sizeof(float) * Size);
		}
		void setIdentity() {
			assert(Rows == Columns && "matX::setIdentity(): bad matrix size");
			Math::memset(data,0,sizeof(float) * Size);
			float *d = data;
			for(int i = 0; i < Rows; i++) {
				*d = 1.0f;
				d += Rows + 1;
			}
		}
		
	private:
		
		enum {
			Size = Rows * Columns,
		};
		
		ATTRIBUTE_ALIGNED16(float data[Size]);
};

/*
 */
template <int Rows,int Columns> int operator==(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	return compare(m0,m1);
}

template <int Rows,int Columns> int operator!=(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	return !compare(m0,m1);
}

template <int Rows,int Columns> matX<Rows,Columns> operator*(const matX<Rows,Columns> &m,float v) {
	matX<Rows,Columns> ret;
	return mul(ret,m,v);
}

template <int Rows,int Columns> vecX<Rows> operator*(const matX<Rows,Columns> &m,const vecX<Columns> &v) {
	vecX<Rows> ret;
	return mul(ret,m,v);
}

template <int Rows,int Columns> vecX<Columns> operator*(const vecX<Rows> &v,const matX<Rows,Columns> &m) {
	vecX<Columns> ret;
	return mul(ret,v,m);
}

template <int Rows,int Columns> matX<Rows,Columns> operator*(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	matX<Rows,Columns> ret;
	return mul(ret,m0,m1);
}

template <int Rows,int Columns> matX<Rows,Columns> operator+(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	matX<Rows,Columns> ret;
	return add(ret,m0,m1);
}

template <int Rows,int Columns> matX<Rows,Columns> operator-(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	matX<Rows,Columns> ret;
	return sub(ret,m0,m1);
}

/*
 */
template <int Rows,int Columns> int compare(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	for(int i = 0; i < m0.getSize(); i++) {
		if(compare(m0[i],m1[i]) == 0) return 0;
	}
	return 1;
}

template <int Rows,int Columns> int compare(const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1,float epsilon) {
	for(int i = 0; i < m0.getSize(); i++) {
		if(compare(m0[i],m1[i],epsilon) == 0) return 0;
	}
	return 1;
}

template <int Size> float trace(const matX<Size,Size> &m) {
	float trace = m[0];
	for(int i = 1; i < Size; i++) {
		trace += m[Size * i + i];
	}
	return trace;
}

template <int Size> float determinant(const matX<Size,Size> &m) {
	int index[Size];
	float det = 0.0f;
	matX<Size,Size> temp = m;
	int ret = luDecompose(temp,index);
	if(ret) det = luDeterminant(temp,ret);
	return det;
}

template <int Rows,int Columns> matX<Rows,Columns> &mul(matX<Rows,Columns> &ret,const matX<Rows,Columns> &m,float v) {
	Simd::mul(ret.get(),m.get(),v,ret.getSize());
	return ret;
}

template <int Rows,int Columns> vecX<Rows> &mul(vecX<Rows> &ret,const matX<Rows,Columns> &m,const vecX<Columns> &v) {
	float *d = ret.get();
	const float *s = m.get();
	Simd::mul(d,s,v[0],Rows);
	for(int i = 1; i < Columns; i++) {
		s += Rows;
		Simd::mad(d,s,v[i],d,Rows);
	}
	return ret;
}

template <int Rows,int Columns> vecX<Columns> &mul(vecX<Columns> &ret,const vecX<Rows> &v,const matX<Rows,Columns> &m) {
	const float *s = m.get();
	for(int i = 0; i < Columns; i++) {
		Simd::dot(ret[i],s,v.get(),Rows);
		s += Rows;
	}
	return ret;
}

template <int Rows,int Columns,int Size> matX<Rows,Columns> &mul(matX<Rows,Columns> &ret,const matX<Rows,Size> &m0,const matX<Size,Columns> &m1) {
	float *d = ret.get();
	const float *s1 = m1.get();
	for(int j = 0; j < Columns; j++) {
		for(int i = 0; i < Rows; i++) {
			const float *s0 = m0.get() + i;
			float sum = *s0 * *s1;
			for(int k = 1; k < Size; k++) {
				s0 += Rows;
				sum += *s0 * s1[k];
			}
			*d++ = sum;
		}
		s1 += Size;
	}
	return ret;
}

template <int Rows,int Columns> matX<Rows,Columns> &mad(matX<Rows,Columns> &ret,const matX<Rows,Columns> &m0,float v,const matX<Rows,Columns> &m1) {
	Simd::mad(ret.get(),m0.get(),v,m1.get(),ret.getSize());
	return ret;
}

template <int Rows,int Columns> matX<Rows,Columns> &add(matX<Rows,Columns> &ret,const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	Simd::add(ret.get(),m0.get(),m1.get(),ret.getSize());
	return ret;
}

template <int Rows,int Columns> matX<Rows,Columns> &sub(matX<Rows,Columns> &ret,const matX<Rows,Columns> &m0,const matX<Rows,Columns> &m1) {
	Simd::sub(ret.get(),m0.get(),m1.get(),ret.getSize());
	return ret;
}

template <int Rows,int Columns> matX<Rows,Columns> &transpose(matX<Rows,Columns> &ret,const matX<Columns,Rows> &m) {
	float *d = ret.get();
	for(int j = 0; j < Columns; j++) {
		const float *s = m.get() + j;
		*d++ = *s;
		for(int i = 1; i < Rows; i++) {
			s += Columns;
			*d++ = *s;
		}
	}
	return ret;
}

template <int Size> matX<Size,Size> &inverse(matX<Size,Size> &ret,const matX<Size,Size> &m) {
	int index[Size];
	matX<Size,Size> temp = m;
	if(luDecompose(temp,index)) luInverse(ret,temp,index);
	else ret.set(0.0f);
	return ret;
}

/*
 */
template <int Rows,int Columns> matX<Rows,Columns> transpose(const matX<Columns,Rows> &m) {
	matX<Rows,Columns> ret;
	return transpose(ret,m);
}

template <int Size> matX<Size,Size> inverse(const matX<Size,Size> &m) {
	matX<Size,Size> ret;
	return inverse(ret,m);
}

/*
 */
int luDecompose(float *m,int size,int *index);
float luDeterminant(const float *m,int size,int ret);
void luSolve3(float *ret,const float *m,const float *b,const int *index);
void luSolve4(float *ret,const float *m,const float *b,const int *index);
void luSolve5(float *ret,const float *m,const float *b,const int *index);
void luSolve6(float *ret,const float *m,const float *b,const int *index);
void luSolve(float *ret,const float *m,const float *b,int size,const int *index);

/*
 */
template <int Size> int luDecompose(matX<Size,Size> &m,int *index) {
	return luDecompose(m.get(),Size,index);
}

template <int Size> float luDeterminant(const matX<Size,Size> &m,int ret) {
	return luDeterminant(m.get(),Size,ret);
}

INLINE vecX<3> &luSolve(vecX<3> &ret,const matX<3,3> &m,const vecX<3> &b,const int *index) {
	luSolve3(ret.get(),m.get(),b.get(),index);
	return ret;
}

INLINE vecX<4> &luSolve(vecX<4> &ret,const matX<4,4> &m,const vecX<4> &b,const int *index) {
	luSolve4(ret.get(),m.get(),b.get(),index);
	return ret;
}

INLINE vecX<5> &luSolve(vecX<5> &ret,const matX<5,5> &m,const vecX<5> &b,const int *index) {
	luSolve5(ret.get(),m.get(),b.get(),index);
	return ret;
}

INLINE vecX<6> &luSolve(vecX<6> &ret,const matX<6,6> &m,const vecX<6> &b,const int *index) {
	luSolve6(ret.get(),m.get(),b.get(),index);
	return ret;
}

template <int Size> vecX<Size> &luSolve(vecX<Size> &ret,const matX<Size,Size> &m,const vecX<Size> &b,const int *index) {
	luSolve(ret.get(),m.get(),b.get(),Size,index);
	return ret;
}

template <int Size> matX<Size,Size> &luInverse(matX<Size,Size> &ret,const matX<Size,Size> &m,const int *index) {
	vecX<Size> x;
	vecX<Size> b(0.0f);
	for(int i = 0; i < Size; i++) {
		b[i] = 1.0f;
		luSolve(x.get(),m.get(),b.get(),Size,index);
		ret.setColumn(i,x);
		b[i] = 0.0f;
	}
	return ret;
}

#endif /* __NUMERICAL_H__ */
