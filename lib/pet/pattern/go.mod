module pattern

go 1.18

require (
    pet/error v0.0.0
)
replace (
	pet/error => ../error
)