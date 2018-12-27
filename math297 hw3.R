library(readr)
dat <- read_csv("a3-20181118.csv")
x = unlist(dat[,1])
y = unlist(dat[,2])

## linear splines
require(splines)


splinesboot = function(x,y,k) {
  B = 100
  n = length(x)
  cov = matrix(data = 0, nrow = length(y), ncol = 1)

  dat = data.frame(x,y)
  n = length(x)
  fit <- lmspline(x, y, nknots = 3, na.rm = FALSE)
  predict = pred.spline(x,fit)
  err = sum((predict - dat$y)^2) # calculate err
  fit2 = lmspline(x, y, nknots = k, na.rm = FALSE)
  predict2 = pred.spline(x,fit2)
  err2 = sum((predict2 - dat$y)^2)
  
  y.star = NULL
  mu.star = NULL
  y.star.mean = NULL
  for (b in 1:B) {
    y.boot = predict + rnorm(n, 0, sqrt(err/(n-1)))
    fit = lmspline(x,y.boot, nknots = k, na.rm = FALSE)
    predict = pred.spline(x,fit)
    y.star = cbind(y.star, y.boot)
    mu.star = cbind(mu.star, predict)
  }
  
  for (i in 1:length(x)){ 
    y.star.mean[i] = mean(y.star[i,])
  }
  for (i in 1:length(x)){
    cov[i] = sum(mu.star[i,] * (y.star[i,]-y.star.mean[i]))
  }
  cov = cov/(B-1)
  err.hat = err2 + 2*sum(cov)
  return(err.hat)
}


err.out=NULL

for (k in seq(2,50)){
  err.out = c(err.out,splinesboot(x,y,k))
}
best.k = seq(2,50)[which.min(err.out)]

plot(seq(2,50),err.out)


k = seq(min(x),max(x),len = best.k+2)
k=k[seq(2,length(k)-1)]
k = quantile(x, k, type=1)

fit1 <- lmspline(x, y, nknots = best,k, na.rm = FALSE)

fit1 = lm(y ~ bs(x, degree=1, knots=k))


summary(fit1)
#plot(fit1)
  
plot(x, y, pch = 20, main="Linear splines fit")
pts = seq(min(x), max(x), len=length(x))
val = predict(fit1, data.frame(x = pts))
lines(pts, val, col="red", lwd = 3)

# df of obs


k = seq(min(x),max(x),len = 10)
k = quantile(x, k, type=1)
fit.best = lm(y ~ bs(x, degree=1, knots=k[-1]))
mu = fit.best$fitted.values

df = c()
for (i in 1:length(x)) {
  df = c(df, cov(mu[1:i],y[1:i]))
}
plot(df, pch = 20, col = 'blue')



# cross-validation
cv = function(X,Y){
  n = length(x)
  N_cv = 1
  k = 7
  cv_lab = sample(n,n,replace=F) %% k
  ## randomly split all the indices into k numbers
  h_seq = seq(from=2,to=50, by=1)
  CV_err_h = rep(0,length(h_seq))
  for(i_tmp in 1:N_cv){
    CV_err_h_tmp = rep(0, length(h_seq))
    cv_lab = sample(n,n,replace=F) %% k
    for(i in 1:length(h_seq)){
      h0 = h_seq[i]
      CV_err =0
      for(i_cv in 1:k){
        w_val = which(cv_lab==(i_cv-1))
        X_tr = X[-w_val]
        Y_tr = Y[-w_val]
        X_val = X[w_val]
        Y_val = Y[w_val]
        fit = lmspline(x = X_tr, y = Y_tr, nknots = h0,na.rm = FALSE)
        pred = pred.spline(X_val,fit)
        #kernel_reg = ksmooth(x = X_tr,y=Y_tr,kernel = "normal",bandwidth=h0,x.points=X_val)
        CV_err = CV_err+mean((Y_val-pred)^2,na.rm=TRUE)
        # na.rm = T: remove the case of 'NA'
      }
      CV_err_h_tmp[i] = CV_err/k
    }
    CV_err_h = CV_err_h+CV_err_h_tmp
  }
  CV_err_h = CV_err_h/N_cv
  print(h_seq[which.min(CV_err_h)])
  plot(h_seq,CV_err_h, type="b", lwd=4, col="blue", xlab="knots", ylab="7-CV Error")
  
  return(CV_err_h)
}


cv(x,y)

## df
# e's method
B = 100
n = length(x)
cov = matrix(data = 0, nrow = length(y), ncol = 1)


dat = data.frame(x,y)
#dat = dat[order(dat$x),]
n = length(x)
K = quantile(x, c(0.5), type=1) #knots
#fit <- lm(y ~ bs(x,degree = 1, knots = K))
fit <- lmspline(x, y, nknots = 3, na.rm = FALSE)
predict = pred.spline(x,fit)
#predict = fit$fitted.values
#mu.hat = mean(predict)
err = sum((predict - dat$y)^2) # calculate err
#k = seq(min(x),max(x),len = k+2)
#k=k[seq(2,length(k)-1)]
#k = quantile(x, k, type=1)
#fit2 <- lm(y ~ bs(x,degree = 1, knots = k))
#predict2 = fit2$fitted.values
fit2 = lmspline(x, y, nknots = 10, na.rm = FALSE)
predict2 = pred.spline(x,fit2)
err2 = sum((predict2 - dat$y)^2)

y.star = NULL
mu.star = NULL
y.star.mean = NULL
#prob = dnorm(y, mean = mu.hat,sd = sqrt(err/(n-1)))
for (b in 1:B) {
  y.boot = predict + rnorm(n, 0, sqrt(err/(n-1)))
  #boot.sample = sample(c(1:n),size = n, replace = TRUE,prob = prob) # follow the normal distribution
  #boot.sample = cbind(x, y)[sample(1:n, size = n, replace = TRUE,prob = prob), ]
  #sample = dat[boot.sample,]
  #sample = sample[order(sample$x),]
  #fit = lm(sample[,2] ~ bs(sample[,1], degree = 1, knots = k))
  #fit = lm(y.boot ~ bs(x, degree = 1, knots = k))
  fit = lmspline(x,y.boot, nknots = 10, na.rm = FALSE)
  predict = pred.spline(x,fit)
  #predict = fit$fitted.values
  y.star = cbind(y.star, y.boot)
  mu.star = cbind(mu.star, predict)
  # err2 = sum((predict - dat$y)^2)
}

for (i in 1:length(x)){ 
  y.star.mean[i] = mean(y.star[i,])
}
for (i in 1:length(x)){
  cov[i] = sum(mu.star[i,] * (y.star[i,]-y.star.mean[i]))
}
cov = cov/(B-1)

df=NULL  
for (i in 1:n) {
  df[i] = sum(cov[1:i])
}

plot(df, pch = 20, col = 'red')

## cv



### part 2 ###
require(locfit)
## algotithm for part 2

n = length(x)
fit <- locfit(y~x, alpha = c(0,0.01), deg = 1)
predict = predict(fit, newdata = x)
err = sum((predict - dat$y)^2) # calculate err

localboot = function(x,y,h) {
  B = 100
  n = length(x)
  cov = matrix(data = 0, nrow = length(y), ncol = 1)
  
  dat = data.frame(x,y)
  #dat = dat[order(dat$x),]

  
  fit2 <- locfit(y~x, alpha = c(0,h), deg = 1)
  predict2 = predict(fit2, newdata = x)
  err2 = sum((predict2 - dat$y)^2)
  
  y.star = NULL
  mu.star = NULL
  y.star.mean = NULL
  
  for (b in 1:B) {
    y.boot = predict + rnorm(n, 0, sqrt(err/(n-1)))
    #boot.sample = sample(c(1:n),size = n, replace = TRUE,prob = prob) # follow the normal distribution
    #boot.sample = cbind(x, y)[sample(1:n, size = n, replace = TRUE,prob = prob), ]
    #sample = dat[boot.sample,]
    #sample = sample[order(sample$x),]
    #fit = lm(sample[,2] ~ bs(sample[,1], degree = 1, knots = k))
    fit = locfit(y.boot~x, alpha = c(0,h), deg = 1)
    predict = predict(fit, newdata = x)
    y.star = cbind(y.star, y.boot)
    mu.star = cbind(mu.star, predict)
    # err2 = sum((predict - dat$y)^2)
  }
  
  for (i in 1:length(x)){ 
    y.star.mean[i] = mean(y.star[i,])
  }
  for (i in 1:length(x)){
    cov[i] = sum(mu.star[i,] * (y.star[i,]-y.star.mean[i]))
  }
  cov = cov/(B-1)
  err.hat = err2 + 2*sum(cov)
  return(err.hat)
}


err.out=NULL

for (h in seq(0.01,0.5, 0.01)){
  err.out = c(err.out,localboot(x,y,h))
}
best.h = seq(0.01,0.5, 0.01)[which.min(err.out)]

plot(err.out)



### 7-fold cv

# cross-validation
cv = function(X,Y){
  n = length(x)
  N_cv = 1
  k = 7
  cv_lab = sample(n,n,replace=F) %% k
  ## randomly split all the indices into k numbers
  h_seq = seq(from=0.01,to=0.5, by=0.01)
  CV_err_h = rep(0,length(h_seq))
  for(i_tmp in 1:N_cv){
    CV_err_h_tmp = rep(0, length(h_seq))
    cv_lab = sample(n,n,replace=F) %% k
    for(i in 1:length(h_seq)){
      h0 = h_seq[i]
      CV_err =0
      for(i_cv in 1:k){
        w_val = which(cv_lab==(i_cv-1))
        X_tr = X[-w_val]
        Y_tr = Y[-w_val]
        X_val = X[w_val]
        Y_val = Y[w_val]
        fit = locfit(Y_tr~X_tr, alpha = c(0,h0), deg = 1)
        pred = predict(fit, newdata = X_val)
        #kernel_reg = ksmooth(x = X_tr,y=Y_tr,kernel = "normal",bandwidth=h0,x.points=X_val)
        CV_err = CV_err+mean((Y_val-pred)^2,na.rm=TRUE)
        # na.rm = T: remove the case of 'NA'
      }
      CV_err_h_tmp[i] = CV_err/k
    }
    CV_err_h = CV_err_h+CV_err_h_tmp
  }
  CV_err_h = CV_err_h/N_cv
  print(h_seq[which.min(CV_err_h)])
  plot(h_seq,CV_err_h, type="b", lwd=4, col="blue", xlab="knots", ylab="7-CV Error")
  
  return(CV_err_h)
}


cv(x,y)

## plot degree of freedom

n = length(x)
fit <- locfit(y~x, alpha = c(0,0.01), deg = 1)
predict = predict(fit, newdata = x)
err = sum((predict - dat$y)^2) # calculate err


  B = 100
  n = length(x)
  cov = matrix(data = 0, nrow = length(y), ncol = 1)
  
  dat = data.frame(x,y)
  #dat = dat[order(dat$x),]
  
  
  fit2 <- locfit(y~x, alpha = c(0,0.04), deg = 1)
  predict2 = predict(fit2, newdata = x)
  err2 = sum((predict2 - dat$y)^2)
  
  y.star = NULL
  mu.star = NULL
  y.star.mean = NULL
  
  for (b in 1:B) {
    y.boot = predict + rnorm(n, 0, sqrt(err/(n-1)))
    #boot.sample = sample(c(1:n),size = n, replace = TRUE,prob = prob) # follow the normal distribution
    #boot.sample = cbind(x, y)[sample(1:n, size = n, replace = TRUE,prob = prob), ]
    #sample = dat[boot.sample,]
    #sample = sample[order(sample$x),]
    #fit = lm(sample[,2] ~ bs(sample[,1], degree = 1, knots = k))
    fit = locfit(y.boot~x, alpha = c(0,0.04), deg = 1)
    predict = predict(fit, newdata = x)
    y.star = cbind(y.star, y.boot)
    mu.star = cbind(mu.star, predict)
    # err2 = sum((predict - dat$y)^2)
  }
  
  for (i in 1:length(x)){ 
    y.star.mean[i] = mean(y.star[i,])
  }
  for (i in 1:length(x)){
    cov[i] = sum(mu.star[i,] * (y.star[i,]-y.star.mean[i]))
  }
  cov = cov/(B-1)
  err.hat = err2 + 2*sum(cov)

df=NULL  
for (i in 1:n) {
  df[i] = sum(cov[1:i])
}

plot(df, pch = 20, col = 'red')





### part 3
# split point(0.47, 0.78, and 0.90)
splinesboot3 = function(x,y,k) {
  
  B = 100
  n = length(x)
  cov = matrix(data = 0, nrow = length(y), ncol = 1)
  
  
  dat = data.frame(x,y)
  #dat = dat[order(dat$x),]
  n = length(x)
  K = quantile(x, c(0.5), type=1) #knots
  fit <- lm(y ~ bs(x,degree = 1, knots = K))
  predict = fit$fitted.values
  mu.hat = mean(predict)
  err = sum((predict - dat$y)^2) # calculate err
  k = quantile(x, k, type=1)
  fit2 <- lm(y ~ bs(x,degree = 1, knots = k))
  predict2 = fit2$fitted.values
  err2 = sum((predict2 - dat$y)^2)
  
  y.star = NULL
  mu.star = NULL
  y.star.mean = NULL
  #prob = dnorm(y, mean = mu.hat,sd = sqrt(err/(n-1)))
  for (b in 1:B) {
    y.boot = predict + rnorm(n, 0, sqrt(err/(n-1)))
    fit = lm(y.boot ~ bs(x, degree = 1, knots = k))
    #predict = pred.spline(x,fit)
    predict = fit$fitted.values
    y.star = cbind(y.star, y.boot)
    mu.star = cbind(mu.star, predict)
    # err2 = sum((predict - dat$y)^2)
  }
  
  for (i in 1:length(x)){ 
    y.star.mean[i] = mean(y.star[i,])
  }
  for (i in 1:length(x)){
    cov[i] = sum(mu.star[i,] * (y.star[i,]-y.star.mean[i]))
  }
  cov = cov/(B-1)
  err.hat = err2 + 2*sum(cov)
  return(err.hat)
}


a = splinesboot3(x,y,0.47)
b = splinesboot3(x,y,0.78)
c = splinesboot3(x,y,0.9)

# plot(c(a,b,c))

ggplot(data=data.frame(c(0.47,0.8,0.9),c(a,b,c)), aes(x='split points', y='risk')) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()


cv3 = function(X,Y,t){
  n = length(x)
  N_cv = 1
  k = 7
  cv_lab = sample(n,n,replace=F) %% k
  ## randomly split all the indices into k numbers
  h_seq = t
  CV_err_h = rep(0,length(h_seq))
  for(i_tmp in 1:N_cv){
    CV_err_h_tmp = rep(0, length(h_seq))
    cv_lab = sample(n,n,replace=F) %% k
    for(i in 1:length(h_seq)){
      h0 = h_seq[i]
      CV_err =0
      for(i_cv in 1:k){
        w_val = which(cv_lab==(i_cv-1))
        X_tr = X[-w_val]
        Y_tr = Y[-w_val]
        X_val = X[w_val]
        Y_val = Y[w_val]
        fit = locfit(Y_tr~X_tr, alpha = c(0,h0), deg = 1)
        pred = predict(fit, newdata = X_val)
        #kernel_reg = ksmooth(x = X_tr,y=Y_tr,kernel = "normal",bandwidth=h0,x.points=X_val)
        CV_err = CV_err+mean((Y_val-pred)^2,na.rm=TRUE)
        # na.rm = T: remove the case of 'NA'
      }
      CV_err_h_tmp[i] = CV_err/k
    }
    CV_err_h = CV_err_h+CV_err_h_tmp
  }
  CV_err_h = CV_err_h/N_cv
  print(h_seq[which.min(CV_err_h)])
  #plot(h_seq,CV_err_h, type="b", lwd=4, col="blue", xlab="knots", ylab="7-CV Error")
  
  return(CV_err_h)
}


aa = cv3(x,y,0.47)
bb = cv3(x,y,0.78)
cc = cv3(x,y,0.9)
plot(c(0.47,0.78,0.9),c(aa,bb,cc), col = 'red', pch = 16)
barplot(c(aa,bb,cc))


# autometic splitting 
# refer to http://statweb.lsu.edu/faculty/marx/SKiP.pdf