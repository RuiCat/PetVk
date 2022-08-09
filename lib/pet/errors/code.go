package errors

import "fmt"

// Code 错误代码
type Code uint64

// Map 错误信息
type Map map[Code]func(v any) error

// Set 设置错误
func (c *Map) Set(code Code, fu func(v any) error) {
	(*c)[code] = fu
}

// Call 处理错误
func (c *Map) Call(code Code, v any) error {
	if _, ok := (*c)[code]; !ok {
		return fmt.Errorf("code %d not found", code)
	}
	return (*c)[code](v)
}

// Is 错误代码是否存在
func (c *Map) Is(code Code) bool {
	_, ok := (*c)[code]
	return ok
}
