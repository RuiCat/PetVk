package scanin

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

type Fn func(code ...string)

type Scanner interface {
	Add(code string, fn Fn) Scanner
	Fprint(a ...any) (n int, err error)
	Fprintf(format string, a ...any) (n int, err error)
}

type wr struct {
	io.Writer
	*bufio.Scanner
}

func (w *wr) Fprint(a ...any) (n int, err error) {
	return fmt.Fprint(w.Writer, a...)
}
func (w *wr) Fprintf(format string, a ...any) (n int, err error) {
	return fmt.Fprintf(w.Writer, format, a...)
}
func (w *wr) Scan() (scan []string, ok bool) {
	ok = w.Scanner.Scan()
	if ok {
		scan = strings.Split(w.Text(), " ")
	}
	return scan, ok
}

type scanner struct {
	*wr
	fn   Fn
	list map[string]*scanner
}

// Add 添加命令
func (s *scanner) Add(code string, fn Fn) Scanner {
	sc := &scanner{
		wr:   s.wr,
		fn:   fn,
		list: make(map[string]*scanner),
	}
	s.list[code] = sc
	return sc
}

// Fn 添加命令
func (s *scanner) Fn(code ...string) {
	(*s).fn(code...)
	if len(code) > 0 {
		if f, ok := s.list[code[0]]; ok {
			f.Fn(code[1:]...)
		}
	}
}

// NewScanner 创建命令解析器
func NewScanner(r io.Reader, w io.Writer) Scanner {
	s := &scanner{
		wr: &wr{
			Writer:  w,
			Scanner: bufio.NewScanner(r),
		},
		list: make(map[string]*scanner),
	}
	go func() {
		s.Fprint("> ")
		scan, ok := s.Scan()
		for ; ok; scan, ok = s.Scan() {
			if fn, ok := s.list[scan[0]]; ok {
				(*fn).Fn(scan[1:]...)
			} else if scan[0] != "" {
				s.Fprintf("command not found: %s\n", scan[0])
			}
			s.Fprint("> ")
		}
	}()
	return s
}
