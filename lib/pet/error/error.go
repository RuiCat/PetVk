package error

// Any 基础类型
type Any interface {
	any | *any
}

type Is[Out Any] func(fn func() (Out, error)) Out
type Is1[In0 Any, Out Any] func(fn func(in0 In0) (Out, error), in0 In0) Out
type Is2[In0 Any, In1 Any, Out Any] func(fn func(in0 In0, in1 In1) (Out, error), in0 In0, in1 In1) Out
type Is3[In0 Any, In1 Any, In2 Any, Out Any] func(fn func(in0 In0, in1 In1, in2 In2) (Out, error), in0 In0, in1 In1, in2 In2) Out

type IsE func(err error)
type IsV[In Any] func(in In, err error) In
type IsP[In Any] func(in *In, err error) *In

// *err, _ = recover().(error)
