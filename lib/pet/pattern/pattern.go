package peterror

import (
	"fmt"
	"reflect"
	"unsafe"

	// 处理错误
	err "pet/errors"
)

// Member 结构体成员信息定义
type Member interface {
	err.Any
	Name() string   // 获取当前参数指向的成员名称
	Value() err.Any // 返回当前值
}

// MemberHook 结构成员初始化拦截
type MemberHook interface {
	Hook(in any, v Member) bool // 拦截设置值,对值处理后通知是否继续设置值
}

type member struct {
	name  *string
	value *err.Any
}

func (v member) Name() string   { return *v.name }
func (v member) Value() err.Any { return *v.value }

// NewMember 创建构建过程
func NewMember(name string) func(v err.Any) Member {
	return func(v err.Any) Member {
		return &member{name: &name, value: &v}
	}
}

type memberHook struct {
	member
	MemberHook
}

// NewMemberHook 创建带拦截构建过程
func NewMemberHook(name string, hook MemberHook) func(v err.Any) Member {
	return func(v err.Any) Member {
		return &memberHook{MemberHook: hook, member: member{name: &name, value: &v}}
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
			// 获取值
			val := vl.Value()
			// TODO: 断言是否可以拦截
			if hook, ok := vl.(MemberHook); ok {
				if !hook.Hook(val, vl) {
					continue
				}
			}
			// 得到成员
			tv := of.FieldByName(vl.Name())
			if tv.Kind() == reflect.Invalid {
				return out, fmt.Errorf(" Name: %s.%s == nil", of.Type().Name(), vl.Name())
			}
			// 从地址创建成员类型
			tv = reflect.NewAt(tv.Type(), unsafe.Pointer(tv.UnsafeAddr())).Elem()
			// 处理输入类型
			nv := reflect.ValueOf(val)
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
