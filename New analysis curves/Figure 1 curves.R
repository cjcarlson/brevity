
library(codependent)



rob <- curve.df(rob1929,iter=100)

rob <- log(rob)

m1 <- nls(n.par ~ a + b1 * n.host, 
          start=list(a=0,b1=1), data = rob)
m2 <- nls(n.par ~ a + b1 * n.host + b2 * n.host^2, 
          start=list(a=0,b1=1,b2=0), data = rob)
m3 <- nls(n.par ~ a + b1 * exp(b2*n.host), 
          start=list(a=0,b1=1,b2=1), data = rob)
m4 <- nls(n.par ~ a + b1 * n.host + b2 * exp (b3 * n.host), 
          start=list(a=0,b1=1,b2=1,b3=0), data = rob)

# m3 and m4 don't converge 

AIC(m1,m2)

plot(rob)
lines(pr1 <- data.frame(host = rob$n.host,pred = predict(m1)),col='red')
lines(pr2 <- data.frame(host = rob$n.host,pred = predict(m2)),col='blue')

plot(exp(rob))
lines(exp(pr1),col='red')
lines(exp(pr2),col='blue')

# why's rd line so bad

rob <- exp(rob)

m1 <- nls(n.par ~ a * n.host^b, 
          start=list(a=1,b=1), data = rob)
m2 <- nls(n.par ~ a * n.host^(b1+b2*log(n.host)),
          start=list(a=1,b1=1,b2=0), data = rob)

plot(rob)
lines(pr1 <- data.frame(host = rob$n.host,pred = predict(m1)),col='red')
lines(pr2 <- data.frame(host = rob$n.host,pred = predict(m2)),col='blue')

AIC(m1,m2)

coef(m1)[[1]] * (5000)^(coef(m1)[[2]])
coef(m2)[[1]] * (5000)^(coef(m2)[[2]] + coef(m2)[[3]]*log(5000))

# glm approach didnt work

rob <- log(rob)

m1 <- glm(n.par ~ n.host, data = rob)
m2 <- glm(n.par ~ n.host + (n.host^2), data = rob)
m3 <- glm(n.par ~ n.host + exp(n.host), data = rob)

AIC(m1,m2,m3)


plot(rob)
lines(pr1 <- data.frame(host = rob$n.host,pred = predict(m1)),col='red')
lines(pr2 <- data.frame(host = rob$n.host,pred = predict(m2)),col='blue')
lines(pr3 <- data.frame(host = rob$n.host,pred = predict(m3)),col='green')

plot(exp(rob))
lines(exp(pr1),col='red')
lines(exp(pr2),col='blue')
lines(exp(pr3),col='green')
