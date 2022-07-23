package api

import "glfw"

// Show 构建显示窗口
type Show interface {
	Init() error // 窗口初始化
	Run()
}

// ShowDraw 窗口绘画
type ShowDraw interface {
	Init(*glfw.Window)
	Draw()
	Destroy()
}
