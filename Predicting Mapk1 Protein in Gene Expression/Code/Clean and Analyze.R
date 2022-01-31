#Data Mining and Machine Learning (STSCI 4740)
#Authors: Binglin Wang (bw474) and Nivedita Vatsa (nkv4)

#######################
#1. Organizing the Data
#######################

#import data
dir = #UPDATE DIRECTORY NAME
Kidney <- read.csv(paste0(dir, "Predicting Mapk1 Protein in Gene Expression/Data/Kidney_2.csv"))

seed = 55
set.seed(seed)

#transposing the matrix
Kidney <- t(Kidney)
colnames(Kidney) = Kidney[1, ] # the first row will be the header
Kidney = Kidney[-1, ]          # removing the first row.
row.names(Kidney) <- NULL

#destringing the matrix
Kidney.new = matrix(rep(0,24*40),ncol=24,nrow=40)
for(i in 1:24){
  Kidney.new[,i] = as.numeric(as.character(Kidney[,i]))}
colnames(Kidney.new) = colnames(Kidney)

#creating a matrix of all the predictors
x.matrix = matrix(Kidney.new[,-5],ncol=23,nrow=40)
colnames(x.matrix) = colnames(Kidney.new[,-5])

#create data.frame
dat = as.data.frame(cbind(Kidney.new[,5],x.matrix))
colnames(dat)[1] ="Mapk1"

#splitting the data - training and test
train = sample(c(rep(TRUE,32),rep(FALSE,8)),replace=FALSE)
test = (!train)

#easily usable names
dat.train = dat[train,]
dat.test = dat[test,]

x.matrix.train = x.matrix[train,]
x.matrix.test = x.matrix[test,]

y.train = Kidney.new[train,5]
y.test = Kidney.new[test,5]

#######################
#2. Exploring the Data
#######################

#means
round(colMeans(Kidney.new),2)

#sd
round(apply(Kidney.new, 2, sd),2)

#correlations

#3. Best Subset Selection

library(leaps)
regfit.full=regsubsets(x=x.matrix.train,y=y.train,nvmax=23)
#summary(regfit.full)
reg.summary=summary(regfit.full)

# a) Plot the RSS, adj.R2, Cp,  BIC for the selected models
#   (these are approximate to the test MSE)
par(mfrow=c(2,3))

#(i) RSS
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="l") # With more variables, RSS always decreases

#(ii) Adj-R^2
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="l")
which.max(reg.summary$adjr2) # % to see which is the best. max
points(which.max(reg.summary$adjr2),reg.summary$adjr2[which.max(reg.summary$adjr2)], col="red",cex=2,pch=20) # % mark the points max in the curve. using BIC to choose, fewer variables.

#(iii) Mallow Cp
plot(reg.summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(reg.summary$cp)
points(which.min(reg.summary$cp),reg.summary$cp[which.min(reg.summary$cp)],col="red",cex=2,pch=20)

#(iv) BIC is measured on full training data
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
which.min(reg.summary$bic)
points(which.min(reg.summary$bic),reg.summary$bic[which.min(reg.summary$bic)],col="red",cex=2,pch=20)

#(v) AIC is measured on full training data
AIC=rep(0,23)
for (i in 1:23) {
  AIC[i]=glm(Mapk1~.,data=dat[,reg.summary$which[i,]])$aic}

plot(AIC,xlab="Number of Variables",ylab="AIC",type='l')
which.min(AIC)
points(which.min(AIC),AIC[which.min(AIC)],col="red",cex=2,pch=20)

#vi) 10-fold CV (done manually)
# define the function for prediction for subset selection, because there is no predict() method for regsubsets()
predict.regsubsets =function(object ,newdata ,id ,...){
  form=as.formula(object$call[[2]])
  mat=model.matrix(form,newdata)
  coefi=coef(object ,id=id)
  xvars=names(coefi)
  return(mat[,xvars]%*%coefi)
}
# create the folds and shuffle the data rowise, so that for each seed the folds and the corresponding results
# are different.
k=10
folds = c(rep(1:8,each=3),rep(9:10,each=4))
table(folds)
dat.train.ramdomized=dat.train[sample(nrow(dat.train)),]

cv.errors=matrix(NA,k,23, dimnames=list(NULL, paste(1:23)))

for(j in 1:k){
  best.fit = regsubsets(Mapk1 ~., data=dat.train.ramdomized[folds != j,], nvmax = 23)
  
  for (i in 1:23){
    pred = predict.regsubsets(best.fit, dat.train.ramdomized[folds == j, ], id = i)
    cv.errors[j, i] = mean((dat.train.ramdomized$Mapk1[folds == j] - pred)^2)
  }
}
mean.cv.errors = apply(cv.errors,2,mean)
mean.cv.errors

plot(mean.cv.errors,xlab="Number of Variables",ylab="cv errors",type='l')
which.min(mean.cv.errors)
points(which.min(mean.cv.errors),mean.cv.errors[which.min(mean.cv.errors)],col="red",cex=2,pch=20)

#(c)test prediction error and best nump
best.nump = as.numeric(which.min(mean.cv.errors))
glm.fit.best=glm(Mapk1~.,data=dat.train[,reg.summary$which[best.nump,]])
tmse.glm=mean(y.test-predict(glm.fit.best,dat.test,type="response")^2)


summary(glm.fit.best)$coef
par(mfrow=c(1,1))

################
#4. Lasso Model
################
set.seed(seed)
#(a) basic model
library(glmnet) #by default performs 10-fold CV
lasso.model = cv.glmnet(x=x.matrix.train,y=y.train,alpha=1)
round(coef(lasso.model),4)

#just checking how many predictors are non-zero (new stuff)
lasso.coef.extra.nonzero=rep(0,100)
for(i in 1:100){
  set.seed(i)
  lasso.model.extra = cv.glmnet(x=x.matrix.train,y=y.train,alpha=1)
  lasso.coef.extra = coef(lasso.model.extra)[1:24]
  lasso.coef.extra.nonzero[i] = sum(lasso.coef.extra!=0)
}

hist(lasso.coef.extra.nonzero)
sum(lasso.coef.extra.nonzero==1)
sum(lasso.coef.extra.nonzero==2)
sum(lasso.coef.extra.nonzero==3)
sum(lasso.coef.extra.nonzero==4)
sum(lasso.coef.extra.nonzero==5)
sum(lasso.coef.extra.nonzero>5)


set.seed(seed) #new stuff end

#i. basic plot
plot(lasso.model)

#ii. all lambdas
lambda.lasso = lasso.model$lambda

#iii. best lambda
bestlambda.lasso = lasso.model$lambda.min

#iv. highest lambda within 1 SE of best lambda
se1.lambda.lasso = lasso.model$lambda.1se

#v. the 10-fold CV error for each lambda
cverror.lasso = lasso.model$cvm

#plotting CV error v. lambda - variable based on single seed
#the error looks negative in the plot, but it's just very small!
plot(cverror.lasso~lambda.lasso, main="Lasso CV Error v. Lambda",xlab="Lambda",ylab="10-fold CV Error")
points(bestlambda.lasso,min(cverror.lasso),col="green",pch=19)
points(se1.lambda.lasso,cverror.lasso[lambda.lasso==se1.lambda.lasso],col="dark green",pch=19)
legend("topright",col=c("green","dark green"),lty=c(2,2),legend=c("Best Lambda","Best Lambda within 1-SE"))
#(b) sanity check -- running lasso model multiple times to see what the best lambda is
#the randomization process changes for the 10-fold CV error, which selects the best lambda
#normally, this would not be practical

bestlambda.lasso.vec = rep(0,100)
se1.lambda.lasso.vec = rep(0,100)
cv.err.lasso.vec = rep(0,100)

for(i in 1:100){ #start loop
set.seed(i)

lasso.model.loop = cv.glmnet(x=x.matrix.train,y=y.train,alpha=1)
bestlambda.lasso.loop = lasso.model.loop$lambda.min
se1.lambda.lasso.loop = lasso.model.loop$lambda.1se
cv.err.lasso.loop = lasso.model.loop$cvm

bestlambda.lasso.vec[i] = bestlambda.lasso.loop
se1.lambda.lasso.vec[i] = se1.lambda.lasso.loop
cv.err.lasso.vec[i] = min(cv.err.lasso.loop)} #end loop

set.seed(seed) #resetting the seed

#plots
hist(bestlambda.lasso.vec, main="Histogram of Best Lasso Lambdas", xlab="Lambda",col="dark green")
abline(v=bestlambda.lasso,lwd=2,col="yellow")

#hist(se1.lambda.lasso.vec, main="Histogram of the Highest Lasso Lambda within 1 SE of Best Lambda", xlab="Lambda",col="light green")

boxplot(cv.err.lasso.vec, names= c("Lasso"),ylab="10-fold CV Error",main="Lasso CV Errors",sub="for 100 randomization processes",col=c("green"))
 
#(c) test prediction error

pred.lasso = predict(lasso.model,s=bestlambda.lasso,newx=x.matrix.test)
tmse.lasso = mean((pred.lasso-y.test)^2)

####################
#5. Ridge Regression
####################

#(a) basic model
library(glmnet) #by default performs 10-fold CV
ridge.model = cv.glmnet(x=x.matrix.train,y=y.train,alpha=0)
coef(ridge.model)
#i. basic plot
plot(ridge.model)

#ii. all lambdas
lambda.ridge= ridge.model$lambda

#iii. best lambda
bestlambda.ridge = ridge.model$lambda.min

#iv. highest lambda within 1 SE of best lambda
se1.lambda.ridge = ridge.model$lambda.1se

#v. 10-fold CV error for each lambda
cverror.ridge = ridge.model$cvm

#plotting CV error v. lambda - variable based on seed
plot(cverror.ridge~lambda.ridge, main="Ridge CV Error v. Lambda",xlab="Lambda",ylab="10-fold CV Error")
points(bestlambda.ridge,min(cverror.ridge),col="purple",pch=19)
points(se1.lambda.ridge,cverror.ridge[lambda.ridge==se1.lambda.ridge],col="pink",pch=19)
legend("bottomright",col=c("purple","pink"),lty=c(2,2),legend=c("Best Lambda","Best Lambda within 1-SE"))

#b) sanity check - running ridge model multiple times to see what the best lambda is
#the randomization process changes for the 10-fold CV error, which selects the best lambda
bestlambda.ridge.vec = rep(0,100)
se1.lambda.ridge.vec = rep(0,100)
cv.err.ridge.vec = rep(0,100)

for(i in 1:100){ #start loop
  set.seed(i)
  
  ridge.model.loop = cv.glmnet(x=x.matrix.train,y=y.train,alpha=0)
  bestlambda.ridge.loop= ridge.model.loop$lambda.min
  se1.lambda.ridge.loop = ridge.model.loop$lambda.1se
  cv.err.ridge.loop = ridge.model.loop$cvm
  
  bestlambda.ridge.vec[i] = bestlambda.ridge.loop
  se1.lambda.ridge.vec[i] = se1.lambda.ridge.loop
  cv.err.ridge.vec[i] = min(cv.err.ridge.loop)} #end loop

#reset seet
set.seed(seed)

#plots
hist(bestlambda.ridge.vec, main="Histgram of Best Ridge Lambdas", xlab="Lambda",
     col="purple",breaks=3)
abline(v=bestlambda.ridge,lwd=2,col="yellow")
#hist(se1.lambda.ridge.vec, main="Histgram of the Highest Ridge Lambda within 1 SE of Best Lambda", 
#     xlab="Lambda",col="lavender")

boxplot(cv.err.lasso.vec,cv.err.ridge.vec,names=c("Lasso","Ridge"),ylab="10-fold CV Error",
        main="Comparing CV Errors",sub="for 100 randomization processes",col=c("green","purple"))

#(c) test prediction error

pred.ridge = predict(ridge.model,s=bestlambda.ridge,newx=x.matrix.test)
tmse.ridge = mean((pred.ridge-y.test)^2)

########
#6. PCR 
########

library(pls)
set.seed(seed)

#(a) PCR on training data
pcr.model = pcr(Mapk1~.,data=dat.train,scale=TRUE,validation="CV")
MSEP = MSEP(pcr.model)

plot(seq(0:23),MSEP$val[1,1,],col="dark blue",type="l",xlab="Number of Components",
     ylab="Training CV MSE",main="Training CV Mean Squared Error for PCR")
lines(seq(0:23),MSEP$val[2,1,],col="grey") #adjusted CV

best.M = MSEP$comps[MSEP$val[1,1,]==min(MSEP$val[1,1,])]

#(b) sanity check for M
M.vec = rep(0,200)
cv.err.pcr.vec = rep(0,200)
for(i in 1:200){
  set.seed(i)
  pcr.model.loop = pcr(Mapk1~.,data=dat.train,scale=TRUE,validation="CV")
  MSEP.loop = MSEP(pcr.model.loop)
  best.M.loop = MSEP.loop$comps[MSEP.loop$val[1,1,]==min(MSEP.loop$val[1,1,])]
  M.vec[i] = best.M.loop
  cv.err.pcr.vec[i] = min(MSEP.loop$val[1,1,])
}

#reset seed
set.seed(seed)

#plots
hist(M.vec,main="Histogram of Best Number of Components",sub="for 200 randomization processes",xlab="Number of Components",
     col="orange",breaks=3)
abline(v=best.M,lwd=3,col="brown")

boxplot(cv.err.lasso.vec,cv.err.ridge.vec,cv.err.pcr.vec,names=c("Lasso","Ridge","PCR"),ylab="10-fold CV Error",
        main="Comparing CV Errors",sub="for 100 randomization processes for lasso/ridge and 200 for PCR",col=c("green","purple","orange"))



#(c) test prediction error

pred.pcr = predict(pcr.model,x.matrix.test,ncomp=best.M)
tmse.pcr = mean((pred.pcr-y.test)^2)
                   
#comparing all the tmses
barplot(c(tmse.glm,tmse.lasso,tmse.ridge,tmse.pcr),names.arg=c("Subset Selection","Lasso","Ridge","PCR"),
        main="Comparing Test MSE for the Methods",ylab="Test MSE",col=c("red","green","purple","orange"))

#comparing all the tmses
barplot(c(tmse.lasso,tmse.ridge,tmse.pcr),names.arg=c("Lasso","Ridge","PCR"),
        main="Comparing Test MSE for the Methods",ylab="Test MSE",col=c("green","purple","orange"))

##########################################################################################
#8. Comparing test errors of all the models, with our selected parameters (lambda, M etc.).
#We generate multiple training and test subsets to compile the test MSE
##########################################################################################

iter2 = 100
tmse.lasso.loop = rep(0,iter2)
tmse.ridge.loop =  rep(0,iter2)
tmse.pcr.loop =  rep(0,iter2)
tmse.glm.loop = rep(0,iter2)
for(i in 1:iter2){
  
  set.seed(i)
  #splitting the data - training and test
  train.loop = sample(c(rep(TRUE,32),rep(FALSE,8)),replace=FALSE)
  test.loop = (!train.loop)
  
  #easily usable names
  dat.train.loop = dat[train.loop,]
  dat.test.loop = dat[test.loop,]
  
  x.matrix.train.loop = x.matrix[train.loop,]
  x.matrix.test.loop = x.matrix[test.loop,]
  
  y.train.loop = Kidney.new[train.loop,5]
  y.test.loop = Kidney.new[test.loop,5]
  
  #(a) Best subset
  # fit a linear model
  reg.fit.loop = regsubsets(x=x.matrix.train.loop,y=y.train.loop,nvmax=23)
  reg.summary.loop = summary(reg.fit.loop)
  glm.fit.loop2=glm(Mapk1~.,data=dat.train.loop[,reg.summary.loop$which[best.nump,]])
  tmse.glm.loop[i]=mean(y.test.loop-predict(glm.fit.loop2,dat.test.loop,type="response")^2)
  
  #(b) Lasso
  lasso.model.loop2 = cv.glmnet(x=x.matrix.train.loop,y=y.train.loop,alpha=1)
  pred.lasso.loop = predict(lasso.model.loop2,s=bestlambda.lasso,newx=x.matrix.test.loop)
  tmse.lasso.loop[i] = mean((pred.lasso.loop-y.test.loop)^2)
  
  #(c) Ridge
  ridge.model.loop2 = cv.glmnet(x=x.matrix.train.loop,y=y.train.loop,alpha=0)
  pred.ridge.loop = predict(ridge.model.loop2,s=bestlambda.ridge,newx=x.matrix.test.loop)
  tmse.ridge.loop[i] = mean((pred.ridge.loop-y.test.loop)^2)

  #(d) PCR
  pcr.model.loop2 = pcr(Mapk1~.,data=dat.train.loop,scale=TRUE,validation="CV")
  pred.pcr.loop = predict(pcr.model.loop2,x.matrix.test.loop,ncomp=best.M)
  tmse.pcr.loop[i] = mean((pred.pcr.loop-y.test.loop)^2)
}

#comparing TMSE for all 4 methods
boxplot(tmse.glm.loop,tmse.lasso.loop,tmse.ridge.loop,tmse.pcr.loop,names=c("Subset","Lasso","Ridge","PCR"),
        col=c("red","green","purple","orange"),ylab="Test MSE",main="Comparing Test MSE for selected Parameter (number of parameters,lambda, or M)",
              sub="Over 100 Randomization Processes")

#comparing TMSE for all methods except best subset
boxplot(tmse.lasso.loop,tmse.ridge.loop,tmse.pcr.loop,names=c("Lasso","Ridge","PCR"),
        col=c("green","purple","orange"),ylab="Test MSE",main="Comparing Test MSE for selected Parameter (number of parameters,lambda, or M)",
        sub="Over 100 Randomization Processes")
