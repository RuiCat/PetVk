package math

import "math"

// NumNber 泛型数值定义
type NumNber interface {
	float32 | float64 | ~int
}

func Sqrt[T NumNber](v T) T {
	return T(math.Sqrt(float64(v)))
}
func Sin[T NumNber](v T) T {
	return T(math.Sin(float64(v)))
}
func Cos[T NumNber](v T) T {
	return T(math.Cos(float64(v)))
}
func Tan[T NumNber](v T) T {
	return T(math.Tan(float64(v)))
}
func DegreesToRadians[T NumNber](angleDegrees T) T {
	return T(float64(angleDegrees) * float64(math.Pi) / 180.0)
}
func RadiansToDegrees[T NumNber](angleRadians T) T {
	return T(float64(angleRadians) * 180.0 / float64(math.Pi))
}
