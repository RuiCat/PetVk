package petvk

import (
	"fmt"
	"glfw"
	"petvk/api"
	"runtime"
	"time"
	vk "vulkan"
)

// NewShow 创建窗口
func NewShow(width, height int, title string, Darw api.ShowDraw) api.Show {
	return &show{
		Width:  width,
		Height: height,
		Title:  title,
		Darw:   Darw,
		run:    make(chan struct{}, 1),
	}
}

type show struct {
	Window *glfw.Window
	Width  int
	Height int
	Title  string
	Darw   api.ShowDraw
	run    chan struct{}
	isinit bool
}

func (s *show) Run() {
	<-s.run // 等待软件退出
}
func (s *show) Init() error {
	if s.Darw == nil {
		panic("Show.Darw is nil")
	}
	if s.isinit {
		s.isinit = true
		return nil
	}
	// 初始化
	re := make(chan error, 1)
	go func(s *show) {
		// 错误处理
		defer func() {
			if err := recover(); err != nil {
				re <- fmt.Errorf("%s", err)
			} else {
				close(re)
			}
			// 重新初始化
			if s.isinit {
				s.run = make(chan struct{}, 1)
				s.isinit = false
			}
		}()
		orError := func(err any) {
			switch v := err.(type) {
			case error:
				if v != nil {
					panic(err)
				}
			case vk.Result:
				if err := vk.Error(v); err != nil {
					panic(err)
				}
			}
		}
		// 绑定线程
		runtime.LockOSThread()
		// 初始化
		orError(glfw.Init())
		orError(vk.Init())
		// 窗口初始化
		glfw.WindowHint(glfw.ClientAPI, glfw.NoAPI)
		var err error
		s.Window, err = glfw.CreateWindow(s.Width, s.Height, s.Title, nil, nil)
		orError(err)
		// 初始化渲染
		s.Darw.Init(s.Window)
		// 携程逃逸并绑定线程,进入窗口事件处理主循环
		go func() {
			fpsDelay := time.Second / 60
			fpsTicker := time.NewTicker(fpsDelay)
			for range fpsTicker.C {
				if s.Window.ShouldClose() {
					break
				}
				s.Darw.Draw()
				glfw.PollEvents()
			}
			// 销毁
			fpsTicker.Stop()
			s.Darw.Destroy()
			s.Window.Destroy()
			glfw.Terminate()
			// 退出
			close(s.run)
		}()
	}(s)
	return <-re
}
