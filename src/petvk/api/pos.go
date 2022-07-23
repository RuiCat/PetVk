package api

import (
	"pet/math"
)

// Pos3 位置信息
type Pos3 struct {
	math.Vec[float64]
	X, Y, Z float64
}

// Pos4 位置信息
type Pos4 struct {
	Pos3
	W float64
}
