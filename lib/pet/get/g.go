// Copyright 2018 Huan Du. All rights reserved.
// Licensed under the MIT license that can be found in the LICENSE file.

// Package g exposes goroutine struct g to user space.
package get

import (
	"unsafe"
)

func getg() unsafe.Pointer

// GetG returns current g (the goroutine struct) to user space.
func GetG() unsafe.Pointer {
	return getg()
}
