package math

import (
	"bytes"
	"fmt"
)

// MatrixFace 矩阵接口定义
type MatrixFace[T NumNber] interface {
	Identity() MatrixFace[T]                 // 单位矩阵
	Min(a, b MatrixFace[T]) MatrixFace[T]    // 小于
	Max(a, b MatrixFace[T]) MatrixFace[T]    // 大于
	Add(a, b MatrixFace[T]) MatrixFace[T]    // 加法
	Sub(a, b MatrixFace[T]) MatrixFace[T]    // 减法
	Mult(a, b MatrixFace[T]) MatrixFace[T]   // 乘法
	Len() (v T)                              // 范数
	Norm() MatrixFace[T]                     // 化一
	Scale(num T) MatrixFace[T]               // 数乘
	Dot(a MatrixFace[T]) (v T)               // 点积
	Distance(a MatrixFace[T]) (v T)          // 距离
	Transpose(a MatrixFace[T]) MatrixFace[T] // 转置

	String() string
	GetMatrix() *matrix[T]
}

// matrix 矩阵
type matrix[T NumNber] struct {
	Data   []T
	Row    int
	Column int
}

func (m *matrix[T]) GetMatrix() *matrix[T] {
	return m
}

func Matrix[T NumNber](Row, Column int, a ...T) MatrixFace[T] {
	if len(a) == 0 {
		return &matrix[T]{
			Data:   make([]T, Row*Column),
			Row:    Row,
			Column: Column,
		}
	}
	// TODO: 拷贝数据,防止出错
	return &matrix[T]{
		Data:   append(make([]T, 0, Row*Column), a[:]...),
		Row:    Row,
		Column: Column,
	}
}
func (m *matrix[T]) Identity() MatrixFace[T] {
	for i, n := 0, 0; i < m.Row; i++ {
		for j := 0; j < m.Column; j++ {
			if i == j {
				m.Data[n] = 1
			} else {
				m.Data[n] = 0
			}
			n++
		}
	}
	return m
}
func (m *matrix[T]) Min(a, b MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	if len(m.Data) == len(mA.Data) && len(m.Data) == len(mB.Data) {
		for i := range m.Data {
			if mA.Data[i] < mB.Data[i] {
				m.Data[i] = mA.Data[i]
			} else {
				m.Data[i] = mB.Data[i]
			}
		}
	}
	return m
}
func (m *matrix[T]) Max(a, b MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	if len(m.Data) == len(mA.Data) && len(m.Data) == len(mB.Data) {
		for i := range m.Data {
			if mA.Data[i] > mB.Data[i] {
				m.Data[i] = mA.Data[i]
			} else {
				m.Data[i] = mB.Data[i]
			}
		}
	}
	return m
}
func (m *matrix[T]) Add(a, b MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	if len(m.Data) == len(mA.Data) && len(m.Data) == len(mB.Data) {
		for i := range m.Data {
			m.Data[i] = mA.Data[i] + mB.Data[i]
		}
	}
	return m
}
func (m *matrix[T]) Sub(a, b MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	if len(m.Data) == len(mA.Data) && len(m.Data) == len(mB.Data) {
		for i := range m.Data {
			m.Data[i] = mA.Data[i] - mB.Data[i]
		}
	}
	return m
}
func (m *matrix[T]) Mult(a, b MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	mB := b.GetMatrix()
	if len(m.Data) == mA.Row*mB.Column && mA.Column == mB.Row {
		for i, n := 0, 0; i < m.Row; i++ {
			for j := 0; j < m.Column; j++ {
				m.Data[n] = 0
				for o := 0; o < mA.Column; o++ {
					m.Data[n] += mA.Data[i*m.Row+o] * mB.Data[o*m.Column+j]
				}
				n++
			}
		}
	}
	return m
}
func (m *matrix[T]) Scale(num T) MatrixFace[T] {
	for i := range m.Data {
		m.Data[i] = m.Data[i] * num
	}
	return m
}
func (m *matrix[T]) Dot(a MatrixFace[T]) (v T) {
	mA := a.GetMatrix()
	if len(m.Data) == len(mA.Data) {
		for i, k := range mA.Data {
			v += m.Data[i] * k
		}
	}
	return v
}
func (m *matrix[T]) Distance(a MatrixFace[T]) (v T) {
	mA := a.GetMatrix()
	if len(m.Data) == len(mA.Data) {
		for i, o := range mA.Data {
			q := (m.Data[i] - o)
			v += q * q
		}
	}

	return Sqrt[T](v)
}
func (m *matrix[T]) Transpose(a MatrixFace[T]) MatrixFace[T] {
	mA := a.GetMatrix()
	// TODO: 转置矩阵限制条件
	le := len(mA.Data)
	if (mA.Row == m.Row || m.Row == mA.Column) && le == len(m.Data) {
		var n int
		for j := 0; j < mA.Row; j++ {
			for i := j; i < le; i += mA.Row {
				m.Data[n] = mA.Data[i]
				n++
			}
		}
	}
	return m
}
func (m *matrix[T]) Len() (v T) {
	return Sqrt[T](m.Dot(m))
}
func (m *matrix[T]) Norm() MatrixFace[T] {
	return m.Scale(1.0 / m.Len())
}
func (m *matrix[T]) String() string {
	str := bytes.NewBuffer(make([]byte, len(m.Data)))
	k := len(m.Data)
	// TODO: 第一行
	if m.Column == 1 {
		str.WriteString("[")
	} else {
		str.WriteString("╭")
	}
	var n int
	for o := m.Row - 1; ; n++ {
		str.WriteString(fmt.Sprintf("%4v", m.Data[n]))
		if n == o {
			n++
			break
		}
		str.WriteString(",")
	}
	if m.Column == 1 {
		str.WriteString("]")
		return str.String()
	}
	if m.Column > 2 {
		str.WriteString(" ╮\n│")
		// TODO: 中间行
		i := k - m.Row
		for o := 0; n < i; n++ {
			if m.Row == o {
				str.WriteString(" │\n│")
				o = 0
			} else if o != 0 {
				str.WriteString(",")
			}
			str.WriteString(fmt.Sprintf("%4v", m.Data[n]))
			o++
		}
		str.WriteString(" │\n")
	} else {
		str.WriteString(" ╮\n")
	}
	// TODO: 最后行
	str.WriteString("╰")
	for k := k - 1; ; n++ {
		str.WriteString(fmt.Sprintf("%4v", m.Data[n]))
		if n == k {
			break
		}
		str.WriteString(",")
	}
	str.WriteString(" ╯")
	return str.String()
}
