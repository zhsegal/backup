



#########################################################################################################

#	0.  READ DATA and ARRANGE in Tables

#########################################################################################################


case.df		<- read.table("~/Dropbox/Zvi Segal/case.csv",sep = ",",header = T)
control.df	<- read.table("~/Dropbox/Zvi Segal/control.csv",sep = ",",header = T)

case.tbl				<- array(dim = c(100,252,14))
dimnames(case.tbl)	<- list(paste("Patient",1:100),paste("month",1:252),dimnames(case.df)[[2]][3:16])
for(i in 1:14)	case.tbl[,,i]	<- t(array(case.df[[i+2]],dim = c(252,100)))

control.tbl				<- array(dim = c(100,252,14))
dimnames(control.tbl)	<- list(paste("Patient",1:100),paste("month",1:252),dimnames(control.df)[[2]][3:16])
for(i in 1:14)	control.tbl[,,i]	<- t(array(control.df[[i+2]],dim = c(252,100)))


apply(is.na(case.tbl),3,mean)
apply(is.na(control.tbl),3,mean)

apply(case.tbl,3,mean,na.rm = T)
apply(control.tbl,3,mean,na.rm = T)


#	1. Function that performs imputation to a blood-work-outcome sequence 


imp.bwo.seq	 <- function(bwo.seq,lowess.f = 2/3,plt = F)
{
	mnth.ind <- 1:252
	
	ind.y 	<- which(!is.na(bwo.seq))
	pred.x	<- mnth.ind[ind.y]
	pred.y	<- bwo.seq[ind.y]
	
	if(is.na(bwo.seq[1]))
	{
		pred.x <- c(1,pred.x)
		pred.y <- c(median(pred.y[1:3]),pred.y)
	}
	
	if(is.na(bwo.seq[252]))
	{
		pred.x <- c(pred.x,252)
		pred.y <- c(median(pred.y[(length(pred.y)-2):length(pred.y)]),pred.y)
	}
	
	low.pred	<- approx(lowess(pred.x,pred.y,f = lowess.f),xout = mnth.ind)	

	res.sd	<- sd(bwo.seq - low.pred$y,na.rm = T)
	imp.seq	<- low.pred$y + rnorm(252,mean = 0, sd = res.sd)	
	imp.seq[ind.y] <- bwo.seq[ind.y]
	names(imp.seq)	<- paste("month",1:252)
	
	if(plt)
	{
			par(mfcol=c(1,1))
			plot(mnth.ind[is.na(bwo.seq)],imp.seq[is.na(bwo.seq)],pch = 4, col=2, xlim = c(1,252))
			lines(low.pred,col=4)
			points(mnth.ind[ind.y],imp.seq[ind.y],pch=1,col=3)
	}
	names(imp.seq)	<- 
	return(imp.seq)
}
	

	
#	Example: Imputation plot of hb profile for Case-Patient 1, 

aa <- imp.bwo.seq(case.tbl[1,,1],plt = T)


	
#	1. Agregate and impute t3 profiles for the cases and the controls

t3.imp	<- t(apply(case.tbl[,,14],1,imp.bwo.seq))
t3.imp	<- rbind(t3.imp, t(apply(control.tbl[,,14],1,imp.bwo.seq)))

#	Normalize t3 profile for each patient

t3.imp	<- 	t3.imp - apply(t3.imp,1,mean)
t3.imp	<- 	t3.imp / apply(t3.imp,1,sd)
	

#	Perform eigen value decomposition for estimated covariance matrix of t3 profiles

aa <- eigen(t(t3.imp) %*% t3.imp)	
	
par(mfcol=c(1,1))
plot(aa$values)	
sum(aa$values)	

par(mfcol=c(2,2))
for(i in 1:4)
{
	plot(1:252,aa$vectors[,i],main = paste("Eigen vector",i),xlab = "Month")
	lines(lowess(1:252,aa$vectors[,i],f = 1/3),lwd = 2,col = 2)
}

#	Save  -- SMOOTHED -- top 4 eigen vectors 

#t3.PC.matrix <- aa[[2]][,1:4]
t(t3.PC.matrix) %*% t3.PC.matrix

t3.PC.matrix <- array(dim = c(252,4))
for(i in 1:4)	t3.PC.matrix[,i] <- lowess(1:252,aa$vectors[,i],f = 1/3)$y
t(t3.PC.matrix) %*% t3.PC.matrix

#	Express t3 profile for cases and controls in terms of top 3 smoothed eigen vector regression coefficients 

summary(lm(case.tbl[1,,14] ~ t3.PC.matrix[,1:3]))

control.coef <- array(dim = c(100,4))
dimnames(control.coef) <- list(paste("patient",1:100),c("Int","PC 1","PC 2","PC 3"))

case.coef <- array(dim = c(100,4))
dimnames(control.coef) <- list(paste("patient",1:100),c("Int","PC 1","PC 2","PC 3"))


for(i in 1:100)
{
	case.coef[i,] <- lm(case.tbl[i,,14] ~ t3.PC.matrix[,1:3])$coef
	control.coef[i,] <- lm(control.tbl[i,,14] ~ t3.PC.matrix[,1:3])$coef
}

par(mfcol = c(1,2))
boxplot(split(c(control.coef),rep(c("Int","PC 1","PC 2","PC 3"),each = 100)))
boxplot(split(c(case.coef),rep(c("Int","PC 1","PC 2","PC 3"),each = 100)))

par(mfcol = c(1,1))
pairs(case.coef,main = "Cases")

par(mfcol = c(1,1))
pairs(control.coef,,main = "Controls")


which(control.coef[,3] < -2)
mean(is.na(control.tbl[74,,14]))




