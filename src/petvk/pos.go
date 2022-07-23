package petvk

import (
	"pet/math"
	"petvk/api"
	"unsafe"
)

func Pos3() *api.Pos3 {
	pos := &api.Pos3{}
	pos.Vec = math.VecPtr[float64](unsafe.Pointer(&pos.X), 3)
	return pos
}
func Pos4() *api.Pos4 {
	pos := &api.Pos4{}
	pos.Vec = math.VecPtr[float64](unsafe.Pointer(&pos.X), 4)
	return pos
}
