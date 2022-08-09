package errors

// Any 基础类型
type Any interface {
	any | *any
}

// Is 错误相关处理
func Is[Out Any](fn func() (Out, error)) Out {
	out, err := fn()
	IsE(err)
	return out
}
func IsE(err error) {
	if err != nil {
		panic(err)
	}
}
func IsV[In Any](in In, err error) In {
	IsE(err)
	return in
}
func IsP[In Any](in *In, err error) *In {
	IsE(err)
	return in
}
func Is1[In0 Any, Out Any](fn func(in0 In0) (Out, error), in0 In0) Out {
	out, err := fn(in0)
	IsE(err)
	return out
}
func Is2[In0 Any, In1 Any, Out Any](fn func(in0 In0, in1 In1) (Out, error), in0 In0, in1 In1) Out {
	out, err := fn(in0, in1)
	IsE(err)
	return out
}
func Is3[In0 Any, In1 Any, In2 Any, Out Any](fn func(in0 In0, in1 In1, in2 In2) (Out, error), in0 In0, in1 In1, in2 In2) Out {
	out, err := fn(in0, in1, in2)
	IsE(err)
	return out
}

// Transmit 错误传递
type Transmit interface {
	Be() bool         // 存在错误
	Error() string    // 错误接口
	Recover()         // 拦截错误
	BindMap(list Map) // 绑定错误代码信息
	// 设置错误,当拦截被触发则输出当前设置信息与错误信息
	Set(code Code, a ...any)
	// 添加错误
	Add(code Code)
	AddE(code Code, err error)
	AddS(code Code, str string)
	AddF(code Code, format string, a ...any)
}
