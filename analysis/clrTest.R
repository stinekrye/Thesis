library(compositions)
m <- matrix(c(1,2,3,4,5,6,7,8,0,0,0,0,0,1,2,3,4,5,6,7,8,0,0,0), nrow = 12, ncol = 2)
clrM <- clr(as.data.frame(m[,1]))
m[,1] <- clr(m[,1])
m[,2] <- clr(m[,2])
n <- apply(m, 2, clr)

rM <- t(t(m)/colSums(m))
rM[,1] <- clr(rM[,1])
rM[,2] <- clr(rM[,2])

n <- apply(m, 2, clr)


# Test for first column
GM <- exp(mean(log(m[,1])))
test <- log(m[,1]/GM)
