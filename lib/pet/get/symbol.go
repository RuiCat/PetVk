package get

import (
	"reflect"
	"unsafe"
)

//go:linkname activeModules runtime.activeModules
func activeModules() []*uintptr

type typeOff int32

type Types uintptr

// FIXME: 非版本安全代码
var offsetTypes uintptr = 280
var offsetLinks uintptr = 336
var offsetTypemap uintptr = 528

//go:linkname (*Types).Name runtime.(*_type).name
func (*Types) Name() string

//go:linkname (*Types).String runtime.(*_type).string
func (*Types) String() string

//go:linkname (*Types).Pkgpath runtime.(*_type).pkgpath
func (*Types) Pkgpath() string

//go:linkname ToType reflect.toType
func ToType(t *Types) reflect.Type

func getSymbol() map[string]reflect.Type {
	m := make(map[string]reflect.Type)
	ptr := uintptr(unsafe.Pointer(activeModules()[0]))
	// 结构体元素偏移
	types := *(*uintptr)(unsafe.Pointer(ptr + offsetTypes))
	links := *(*[]int32)(unsafe.Pointer(ptr + offsetLinks))
	typemap := *(*map[typeOff]*Types)(unsafe.Pointer(ptr + offsetTypemap))
	// 输出结果
	for _, tl := range links {
		var t *Types
		if typemap == nil {
			t = (*Types)(unsafe.Pointer(types + uintptr(tl)))
		} else {
			t = typemap[typeOff(tl)]
		}
		m[t.String()] = ToType(t)
	}
	return m
}

// SymbolList 符号列表
var SymbolList = getSymbol()

// UpdateSymbol 更新符号列表
func UpdateSymbol() {
	SymbolList = getSymbol()
}

// GetSymbol 获取符号
func GetSymbol(symbol string) (reflect.Type, bool) {
	if sym, ok := SymbolList[symbol]; ok {
		return sym, true
	}
	return nil, false
}

// NewSymbol 通过符号值创建符号
func NewSymbol(symbol string, ptr unsafe.Pointer) (reflect.Value, bool) {
	if sym, ok := SymbolList[symbol]; ok {
		if ptr != nil {
			return reflect.NewAt(sym, ptr), true
		}
		return reflect.New(sym), true
	}
	return reflect.Value{}, false
}
