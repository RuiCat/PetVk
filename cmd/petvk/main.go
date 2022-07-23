package main

import (
	"fmt"
	"petvk"
)

func main() {
	show := petvk.NewShow(300, 200, "测试窗口", petvk.NewDraw())
	if err := show.Init(); err != nil {
		fmt.Println(err)
		return
	}
	show.Run()
}
