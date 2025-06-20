pkgname <- "MultiDiscreteRNG"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MultiDiscreteRNG')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("validation.Bparameters")
### * validation.Bparameters

flush(stderr()); flush(stdout())

### Name: validation.Bparameters
### Title: Validate if the input Binomial parameters are within feasible
###   range
### Aliases: validation.Bparameters

### ** Examples

validation.Bparameters(n.vec = c(10, 15), p.vec = c(0.4, 0.2))



cleanEx()
nameEx("validation.GPDparameters")
### * validation.GPDparameters

flush(stderr()); flush(stdout())

### Name: validation.GPDparameters
### Title: Validate if the input GPD parameters are within feasible range
### Aliases: validation.GPDparameters

### ** Examples

validation.GPDparameters(theta.vec = c(3, 2), lambda.vec = c(0.4, 0.2))



cleanEx()
nameEx("validation.NBparameters")
### * validation.NBparameters

flush(stderr()); flush(stdout())

### Name: validation.NBparameters
### Title: Validate if the input NB parameters are within feasible range
### Aliases: validation.NBparameters

### ** Examples

validation.NBparameters(r.vec = c(10, 15), prob.vec = c(0.7, 0.5))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
