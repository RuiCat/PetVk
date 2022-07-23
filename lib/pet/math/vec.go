package math

import (
	"reflect"
	"unsafe"
)

// Vec 向量
type Vec[T NumNber] interface {
	MatrixFace[T]
	Tran() Vec[T]               // 转置
	Cross(a, b Vec[T]) Vec[T]   // 叉积
	Reflect(a, b Vec[T]) Vec[T] // 反射
}

type vec[T NumNber] struct{ *matrix[T] }

func (v vec[T]) Tran() Vec[T] {
	// TODO:  向量直接转置
	if v.Column == 1 || v.Row == 1 {
		v.Row, v.Column = v.Column, v.Row
	}
	return v
}
func (v vec[T]) Cross(a, b Vec[T]) Vec[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	mV := v.GetMatrix()
	if len(mV.Data) == len(mA.Data) && len(mV.Data) == len(mB.Data) && len(mB.Data) >= 3 {
		mV.Data[0] = mA.Data[1]*mB.Data[2] - mA.Data[2]*mB.Data[1]
		mV.Data[1] = mA.Data[2]*mB.Data[0] - mA.Data[0]*mB.Data[2]
		mV.Data[2] = mA.Data[0]*mB.Data[1] - mA.Data[1]*mB.Data[0]
		for i, n := 3, len(v.Data); i < n; i++ {
			mV.Data[i] = 1
		}
	}
	return v
}
func (v vec[T]) Reflect(a, b Vec[T]) Vec[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	mV := v.GetMatrix()
	if len(mV.Data) == len(mA.Data) && len(mV.Data) == len(mB.Data) {
		p := 2 * a.Dot(b)
		for i := range v.Data {
			v.Data[i] = mA.Data[i] - p*mB.Data[i]
		}
	}
	return v
}
func VecN[T NumNber](v ...T) Vec[T] {
	return &vec[T]{matrix: &matrix[T]{
		Data:   v,
		Row:    len(v),
		Column: 1,
	}}
}
func Vec2[T NumNber]() Vec[T] {
	return VecN[T](0, 0)
}
func Vec3[T NumNber]() Vec[T] {
	return VecN[T](0, 0, 0)
}
func Vec4[T NumNber]() Vec[T] {
	return VecN[T](0, 0, 0, 0)
}
func VecPtr[T NumNber](ptr unsafe.Pointer, len int) Vec[T] {
	return &vec[T]{matrix: &matrix[T]{
		Data: *(*[]T)((unsafe.Pointer)(&reflect.SliceHeader{
			Data: uintptr(ptr),
			Len:  len,
			Cap:  len},
		)),
		Row:    len,
		Column: 1,
	}}
}
