package decode

import "pet/errors"

// Type 解码过程
type Type uint64

// Call 解码回调
type Call[V any] func(this Decode[V], out V) bool

// Decode 解码定义
type Decode[V any] interface {
	errors.Transmit                         // 共用: 错误传递
	Bind(t Type, fn func(t Type, v V) bool) // 共用: 绑定一个解析过程
	Process(t Type, call Call[V]) bool      // 层级: 定义一个解析
}

// Coding 变长编码树
type Coding[Data any, Value any] interface {
	Seek(data Data) (Value, bool)         // 查找数据
	Increase(data Data, value Value)      // 以当前节点添加数据构建新树
	ForEach(func(data Data, varlu Value)) // 递归可达数据
}
