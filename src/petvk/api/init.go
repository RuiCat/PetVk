package api

import (
	"glfw"
	vk "vulkan"
)

func init() {
	// 绑定地址
	procAddr := glfw.GetVulkanGetInstanceProcAddress()
	if procAddr == nil {
		panic("GetInstanceProcAddress is nil")
	}
	vk.SetGetInstanceProcAddr(procAddr)
}
