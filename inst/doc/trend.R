### R code from vignette source 'trend.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: trend.Rnw:83-87
###################################################
require(trend)
data(maxau)
Q <- maxau[,"Q"]
mk.test(Q)


###################################################
### code chunk number 2: trend.Rnw:107-109
###################################################
require(trend)
smk.test(nottem)


###################################################
### code chunk number 3: trend.Rnw:117-119
###################################################
require(trend)
csmk.test(nottem)


###################################################
### code chunk number 4: trend.Rnw:129-132
###################################################
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
cor.test(s,Q, meth="spearman")


###################################################
### code chunk number 5: trend.Rnw:137-141
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.mk.test(s,Q)


###################################################
### code chunk number 6: trend.Rnw:150-154
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.cor.trend.test(s,Q, "spearman")


###################################################
### code chunk number 7: trend.Rnw:185-188
###################################################
require(trend)
s <- maxau[,"s"]
sens.slope(s)


###################################################
### code chunk number 8: trend.Rnw:203-205
###################################################
require(trend)
sea.sens.slope(nottem)


###################################################
### code chunk number 9: trend.Rnw:232-235
###################################################
require(trend)
data(PagesData)
pettitt.test(PagesData)


