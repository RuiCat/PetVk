package peterror

import (
	"fmt"
	"reflect"
	"unsafe"

	// 处理错误
	. "pet/error"
)

// Member 结构体成员信息定义
type Member interface {
	Any
	Name() string // 获取当前参数指向的成员名称
	Value() Any   // 返回当前值
}

// MemberHook 结构成员初始化拦截
type MemberHook interface {
	Hook(t any, v Member) bool
}

type member struct {
	name  *string
	value *Any
}

func (v member) Name() string { return *v.name }
func (v member) Value() Any   { return *v.value }

// NewMember 创建通用结构体成员定义
func NewMember(name string) func(v Any) Member {
	return func(v Any) Member {
		return &member{name: &name, value: &v}
	}
}

// Pattern 构建结构体
func Pattern[Out any](v ...Member) (out *Out, err error) {
	// 创建对象
	out = new(Out)
	// 反射
	of := reflect.ValueOf(out).Elem()
	switch of.Kind() {
	case reflect.Struct: // 处理结构体
		for _, vl := range v {
			// 得到成员
			tv := of.FieldByName(vl.Name())
			if tv.Kind() == reflect.Invalid {
				return out, fmt.Errorf(" Name: %s.%s == nil", of.Type().Name(), vl.Name())
			}
			// 从地址创建成员类型
			tv = reflect.NewAt(tv.Type(), unsafe.Pointer(tv.UnsafeAddr())).Elem()
			// 处理输入类型
			nv := reflect.ValueOf(vl.Value())
			if tv.Kind() != nv.Kind() {
				return out, fmt.Errorf(" In: %v != Set: %v", nv.Kind(), tv.Kind())
			}
			// 设置成员数据
			tv.Set(nv)
		}
	default:
		ty := of.Type()
		return nil, fmt.Errorf(" Type: %s != struct", ty.String())
	}
	return
}
