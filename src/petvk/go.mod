module petvk

go 1.18

require (
	vulkan v0.0.0
	glfw v0.0.0
	pet/math v0.0.0
	pet/get v0.0.0
)

replace (
	vulkan => ../../lib/vulkan
	glfw => ../../lib/glfw
	pet/math => ../../lib/pet/math
	pet/get => ../../lib/pet/get
)
