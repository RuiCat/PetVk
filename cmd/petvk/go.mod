module main

go 1.18

require (
	petvk v0.0.0
	vulkan v0.0.0
	glfw v0.0.0
    pet/math v0.0.0
    pet/get v0.0.0
)
replace (
	petvk => ../../src/petvk
	vulkan => ../../lib/vulkan
	glfw => ../../lib/glfw
	pet/math => ../../lib/pet/math
    pet/get => ../../lib/pet/get
)