### R code from vignette source 'trend.Rnw'

###################################################
### code chunk number 1: trend.Rnw:86-90
###################################################
require(trend)
data(maxau)
Q <- maxau[,"Q"]
mk.test(Q)


###################################################
### code chunk number 2: trend.Rnw:110-112
###################################################
require(trend)
smk.test(nottem)


###################################################
### code chunk number 3: trend.Rnw:120-122
###################################################
require(trend)
csmk.test(nottem)


###################################################
### code chunk number 4: trend.Rnw:133-136
###################################################
require(trend)
data(hcb)
plot(hcb)


###################################################
### code chunk number 5: trend.Rnw:139-144
###################################################
## Single site trends
site <- c("we", "ka", "mz", "ko", "bh", "bi")
for (i in 1:6) {print(site[i]) ; print(mk.test(hcb[,site[i]], continuity = TRUE))}
## Regional trend (all stations including covariance between stations
mult.mk.test(hcb)


###################################################
### code chunk number 6: trend.Rnw:153-156
###################################################
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
cor.test(s,Q, meth="spearman")


###################################################
### code chunk number 7: trend.Rnw:161-165
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.mk.test(s,Q)


###################################################
### code chunk number 8: trend.Rnw:174-178
###################################################
require(trend)
data(maxau)
s <- maxau[,"s"]; Q <- maxau[,"Q"]
partial.cor.trend.test(s,Q, "spearman")


###################################################
### code chunk number 9: trend.Rnw:187-197
###################################################
## Example from Schoenwiese (1992, p. 114)
## Number of frost days in April at Munich from 1957 to 1968
## z = -0.5, Accept H0
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
cs.test(frost)
     
## Example from Sachs (1997, p. 486-487)
## z ~ 2.1, Reject H0 on a level of p = 0.0357
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
cs.test(x)


###################################################
### code chunk number 10: trend.Rnw:226-229
###################################################
require(trend)
s <- maxau[,"s"]
sens.slope(s)


###################################################
### code chunk number 11: trend.Rnw:244-246
###################################################
require(trend)
sea.sens.slope(nottem)


###################################################
### code chunk number 12: trend.Rnw:273-276
###################################################
require(trend)
data(PagesData)
pettitt.test(PagesData)


###################################################
### code chunk number 13: trend.Rnw:312-314
###################################################
require(trend)
(res <- br.test(Nile))


###################################################
### code chunk number 14: trend.Rnw:317-319
###################################################
par(mfrow=c(2,1))
plot(Nile); plot(res)


###################################################
### code chunk number 15: trend.Rnw:340-342
###################################################
require(trend)
(res <- bu.test(Nile))


###################################################
### code chunk number 16: trend.Rnw:345-347
###################################################
par(mfrow=c(2,1))
plot(Nile); plot(res)


###################################################
### code chunk number 17: trend.Rnw:375-377
###################################################
require(trend)
(res <- snh.test(Nile))


###################################################
### code chunk number 18: trend.Rnw:380-382
###################################################
par(mfrow=c(2,1))
plot(Nile); plot(res)


###################################################
### code chunk number 19: trend.Rnw:391-401
###################################################
## Example from Schoenwiese (1992, p. 113)
## Number of frost days in April at Munich from 1957 to 1968
## z = -0.124, Accept H0
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
wm.test(frost)
     
## Example from Sachs (1997, p. 486)
## z = 2.56, Reject H0 on a level of p < 0.05
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
wm.test(x)


###################################################
### code chunk number 20: trend.Rnw:408-421
###################################################
## Example from Schoenwiese (1992, p. 113)
## Number of frost days in April at Munich from 1957 to 1968
## 
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
bartels.test(frost)

## Example from Sachs (1997, p. 486)
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
bartels.test(x)
     
## Example from Bartels (1982, p. 43)
x <- c(4, 7, 16, 14, 12, 3, 9, 13, 15, 10, 6, 5, 8, 2, 1, 11, 18, 17)
bartels.test(x)


###################################################
### code chunk number 21: trend.Rnw:427-440
###################################################
## Example from Schoenwiese (1992, p. 113)
## Number of frost days in April at Munich from 1957 to 1968
## 
frost <- ts(data=c(9,12,4,3,0,4,2,1,4,2,9,7), start=1957)
ww.test(frost)

## Example from Sachs (1997, p. 486)
x <- c(5,6,2,3,5,6,4,3,7,8,9,7,5,3,4,7,3,5,6,7,8,9)
ww.test(x)
     
## Example from Bartels (1982, p. 43)
x <- c(4, 7, 16, 14, 12, 3, 9, 13, 15, 10, 6, 5, 8, 2, 1, 11, 18, 17)
ww.test(x)


