package petvk

import (
	"pet/get"
	"unsafe"
)

// offsetGid 偏移
var offsetGid uintptr

func init() {
	if v, ok := get.GetSymbol("*runtime.g"); ok {
		if sf, ok := v.Elem().FieldByName("goid"); ok {
			offsetGid = uintptr(sf.Offset)
			return
		}
	}
	panic("Symbol: *runtime.g is nil")
}

// Gid 携程 ID 用于非锁数据安全
func Gid() int64 {
	return *(*int64)(unsafe.Pointer((uintptr)(get.GetG()) + uintptr(offsetGid)))
}
