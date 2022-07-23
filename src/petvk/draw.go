package petvk

/*
	Darw 实现模块功能
		1. vulkan 相关初始化
		2. 数据间 携程通信
*/
import (
	"glfw"
	"petvk/api"
)

// NewDraw 创建显示处理
func NewDraw() api.ShowDraw {
	return &draw{}
}

type draw struct{}

func (d *draw) Init(w *glfw.Window) {}
func (d *draw) Draw()               {}
func (d *draw) Destroy()            {}
