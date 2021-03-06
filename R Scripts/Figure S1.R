library(cowplot)
library(reshape)

# 1. Robertson 1929: plant-pollinator interactions

rob1929_raw <- read.csv("~/Github/brevity/Figure 1 Demo Networks/robertson1929.csv")
rob1929 <- rob1929_raw[,c('plant','poll')]
colnames(rob1929) <- c("Plant","Pollinator")

# 2. Schleuning 2010: seed-disperser interactions

sch2010_raw <- read.csv("~/Github/brevity/Figure 1 Demo Networks/schleuning2010.csv")
sch2010_raw[sch2010_raw>1] <- 1
sch2010 <- melt(sch2010_raw, id=c("Plant.species"))
sch2010 <- sch2010[sch2010$value==1,]
sch2010 <- sch2010[,c(2,1)]
colnames(sch2010) <- c("Plant","Disperser")

# 3. Toju 2018: plant-arbuscular mycorrhizae

tojuraw.1 <- read.csv('~/Github/brevity/Figure 1 Demo Networks/toju descriptors.csv')
tojuraw.2 <- read.csv('~/Github/brevity/Figure 1 Demo Networks/toju.csv')
otu.list <- unique(tojuraw.1[tojuraw.1$Category.In.This.Study=='Arbuscular_Mycorrhizal',]$OTU.code)

tojuraw.2$Plant <- gsub('_TES', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_TOMA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YKS', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_SGD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YSD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YAKU', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YONA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_IRI', '', tojuraw.2$Plant)

toju <- aggregate(. ~ Plant, data=tojuraw.2, FUN=sum)
toju <- melt(toju, id=c("Plant"))
toju <- toju[toju$value>1,]
colnames(toju)[1:2] <- c("Plant","Microbe")
toju$Microbe <- gsub('X', '', toju$Microbe)
toju$Microbe <- gsub('\\.', ':', toju$Microbe)
toju <- toju[toju$Microbe %in% as.character(otu.list),]
toju <- toju[,c(1:2)]

# 4. Host-helminth relationships

helminths.raw <- read.csv('~/Github/brevity/Figure 1 Demo Networks/helminths.csv')
helminths <- helminths.raw[helminths.raw$group=='Nematoda',]
helminths <- helminths[helminths$hostgroup=='Mammalia',]
helminths <- na.omit(helminths)
helminths <- helminths[,c(2,3)]
helminths <- unique(helminths)
##########################################################

# Analyses for each

library(codependent)

df.poll <- curve.df(rob1929, 100)
df.disp <- curve.df(sch2010, 100)
df.myco <- curve.df(toju, 100)
df.helm <- curve.df(helminths, 100)


########

# Read in data 
virus <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')
# Exclude humans
virus <- virus[!(virus$Host %in% c('Homo_sapiens')),]
virusmeta <- read.csv("~/Github/brevity/Olival Nature 2017 Raw Data/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)

# Compile sub-association lists
rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

df.rna <- curve.df(rnavirus,100)
df.dna <- curve.df(dnavirus,100)

model.rna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.rna)
predicted.rna <- data.frame(pred = predict(model.rna), host = df.rna$n.host)

model.dna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.dna)
predicted.dna <- data.frame(pred = predict(model.dna), host = df.dna$n.host)

###################

ggplotter <- function(data, xlab, ylab, colorname) {
  
  model <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=data)
  model2 <- lm(n.par~n.host,data=data)
  model3 <- lm(n.par~log(n.host),data=data)
  pred <- data.frame(pred = predict(model), host = data$n.host)
  pred2 <- data.frame(pred = predict(model2), host = data$n.host)
  pred3 <- data.frame(pred = predict(model3), host = data$n.host)
  g <- ggplot(data, aes(n.host, n.par)) + xlab(xlab) + ylab(ylab) + 
    geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c(colorname)) + theme_bw() +
    geom_line(color='black',lwd=1,data = pred, aes(x=host, y=pred))+
    geom_line(color='black',lwd=1,lty=2,data = pred2, aes(x=host, y=pred))+
    geom_line(color='black',lwd=1,lty=3,data = pred3, aes(x=host, y=pred))
  g
  
}

g1 <- ggplotter(df.poll,'Plants','Pollinators','steelblue1')
g2 <- ggplotter(df.disp,'Plants','Seed dispersers','red')
g3 <- ggplotter(df.myco,'Plants','Microbe OTUs','orange')
g4 <- ggplotter(df.helm,'Hosts','Helminths','seagreen1')
g5 <- ggplotter(df.rna,'Hosts','RNA viruses','mediumorchid1')
g6 <- ggplotter(df.dna,'Hosts','DNA viruses','slateblue1')

plot2by2 <- plot_grid(g2, g3, g1, g4, g5, g6,
                      labels=c("A", "B", "C", "D", "E", "F"), nrow=3, ncol = 2)
plot2by2

aic.doer <- function(data) {
  
  model <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=data)
  model2 <- lm(n.par~n.host,data=data)
  model3 <- lm(n.par~log(n.host),data=data)
  AIC(model, model2, model3)
  
}

t(data.frame(aic.doer(df.poll)$AIC,
aic.doer(df.disp)$AIC,
aic.doer(df.myco)$AIC,
aic.doer(df.helm)$AIC,
aic.doer(df.rna)$AIC,
aic.doer(df.dna)$AIC))


z.doer <- function(data) {
  
  model <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=data)
  
  q <- stats::coef(model)
  p.cis <- nlstools::confint2(model)
  
  return(c(q[2],p.cis[2,1:2]))
  
}


round(t(data.frame(z.doer(df.poll),
             z.doer(df.disp),
             z.doer(df.myco),
             z.doer(df.helm),
             z.doer(df.rna),
             z.doer(df.dna))),4)





###################
load(file.choose())

model.poll <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.poll)
predicted.poll <- data.frame(pred = predict(model.poll), host = df.poll$n.host)
g1 <- ggplot(log(df.poll), aes(n.host, n.par)) + xlab('Plants') + ylab('Pollinators') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('steelblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.poll), aes(x=host, y=pred))

model.disp <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.disp)
predicted.disp <- data.frame(pred = predict(model.disp), host = df.disp$n.host)
g2 <- ggplot(log(df.disp), aes(n.host, n.par)) + xlab('Plants') + ylab('Seed dispersers') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('red')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.disp), aes(x=host, y=pred))

model.myco <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.myco)
predicted.myco <- data.frame(pred = predict(model.myco), host = df.myco$n.host)
g3 <- ggplot(log(df.myco), aes(n.host, n.par)) + xlab('Plants') + ylab('Microbe OTUs') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('orange')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.myco), aes(x=host, y=pred))

model.helm <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.helm)
predicted.helm <- data.frame(pred = predict(model.helm), host = df.helm$n.host)
g4 <- ggplot(log(df.helm), aes(n.host, n.par)) + xlab('Hosts') + ylab('Helminths') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('seagreen1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.helm), aes(x=host, y=pred))


model.rna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.rna)
predicted.rna <- data.frame(pred = predict(model.rna), host = df.rna$n.host)
g5 <- ggplot(log(df.rna), aes(n.host, n.par)) + xlab('Hosts') + ylab('RNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('mediumorchid1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.rna), aes(x=host, y=pred))

model.dna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.dna)
predicted.dna <- data.frame(pred = predict(model.dna), host = df.dna$n.host)
g6 <- ggplot(log(df.dna), aes(n.host, n.par)) + xlab('Hosts') + ylab('DNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .20, color = c('slateblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = log(predicted.dna), aes(x=host, y=pred))


plot2by2 <- plot_grid(g2, g3, g1, g4, g5, g6,
                      labels=c("A", "B", "C", "D", "E", "F"), nrow=3, ncol = 2)
plot2by2



############# FIGURE S2




resid.poll <- data.frame(n.host = df.poll$n.host,n.par = residuals(model.poll))
ss.p <- data.frame(n.host = smooth.spline(resid.poll)$x, n.par = smooth.spline(resid.poll)$y)
flat <- data.frame(n.host = df.poll$n.host, n.par=0)
g1 <- ggplot(resid.poll, aes(n.host, n.par)) + xlab('Plants') + ylab('Pollinators') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .10, color = c('steelblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))

resid.disp <- data.frame(n.host = df.disp$n.host,n.par = residuals(model.disp))
ss.p <- data.frame(n.host = smooth.spline(resid.disp)$x, n.par = smooth.spline(resid.disp)$y)
flat <- data.frame(n.host = df.disp$n.host, n.par=0)
g2 <- ggplot(resid.disp, aes(n.host, n.par)) + xlab('Plants') + ylab('Seed dispersers') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('red')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))

resid.myco <- data.frame(n.host = df.myco$n.host,n.par = residuals(model.myco))
ss.p <- data.frame(n.host = smooth.spline(resid.myco)$x, n.par = smooth.spline(resid.myco)$y)
flat <- data.frame(n.host = df.myco$n.host, n.par=0)
g3 <- ggplot(resid.myco, aes(n.host, n.par)) + xlab('Plants') + ylab('Microbe OTUs') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('orange')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))

resid.helm <- data.frame(n.host = df.helm$n.host,n.par = residuals(model.helm))
ss.p <- data.frame(n.host = smooth.spline(resid.helm)$x, n.par = smooth.spline(resid.helm)$y)
flat <- data.frame(n.host = df.helm$n.host, n.par=0)
g4 <- ggplot(resid.helm, aes(n.host, n.par)) + xlab('Hosts') + ylab('Helminths') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('seagreen1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))

resid.rna <- data.frame(n.host = df.rna$n.host,n.par = residuals(model.rna))
ss.p <- data.frame(n.host = smooth.spline(resid.rna)$x, n.par = smooth.spline(resid.rna)$y)
flat <- data.frame(n.host = df.rna$n.host, n.par=0)
g5 <- ggplot(resid.rna, aes(n.host, n.par)) + xlab('Hosts') + ylab('RNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('mediumorchid')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))

resid.dna <- data.frame(n.host = df.dna$n.host,n.par = residuals(model.dna))
ss.p <- data.frame(n.host = smooth.spline(resid.dna)$x, n.par = smooth.spline(resid.dna)$y)
flat <- data.frame(n.host = df.dna$n.host, n.par=0)
g6 <- ggplot(resid.dna, aes(n.host, n.par)) + xlab('Hosts') + ylab('DNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('slateblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = ss.p, aes(x=n.host, y=n.par)) + 
  geom_line(color='black',lwd=0.5,lty=2,data = flat, aes(x=n.host, y=n.par))


plot2by2 <- plot_grid(g2, g3, g1, g4, g5, g6,
                      labels=c("A", "B", "C", "D", "E", "F"), nrow=3, ncol = 2)
plot2by2



##################

power.fit <- function(df.input,pred=FALSE,ext=0) {
  
  df <- df.input; names(df) <- c('h','p')
  
  lin <- lm(p ~ h, data=df)
  pow <- nls(p ~ b0*h^b1, start=list(b0=1,b1=1), control=nls.control(maxiter=500), data = df)
  pow.q <- nls(p ~ exp(b0+b1*log(h)+b2*log(h)^2), start=list(b0=1,b1=1,b2=1), control=nls.control(maxiter=1000), data=df)
  epm1 <- nls(p ~ a*h^(b*h^(-d)), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  epm2 <- nls(p ~ a*h^(b-(d/h)), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  p1 <- nls(p ~ a*h^(b)*exp(-d*h), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  p2 <- nls(p ~ a*h^(b)*exp(-d/h), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  
  results.df <- data.frame(AIC(lin,pow,pow.q,epm1,epm2,p1,p2))

  results.df$min <- as.numeric(results.df$AIC == min(results.df$AIC))
  if(pred==TRUE) {
    
    pred.1 <- predict(lin,list(h=ext))
    pred.2 <- predict(pow,list(h=ext))
    pred.3 <- predict(pow.q,list(h=ext))
    pred.4 <- predict(epm1,list(h=ext))
    pred.5 <- predict(epm2,list(h=ext))
    pred.6 <- predict(p1,list(h=ext))
    pred.7 <- predict(p2,list(h=ext))
    
    results.df$pred <- c(pred.1,
                         pred.2,
                         pred.3,
                         pred.4,
                         pred.5,
                         pred.6,
                         pred.7)
    
    return(results.df)
    
  } else {
    
    return(results.df)
    
  }
  
}
  
power.fit.noepm <- function(df.input,pred=FALSE,ext=0) {
  
  df <- df.input; names(df) <- c('h','p')
  
  lin <- lm(p ~ h, data=df)
  pow <- nls(p ~ b0*h^b1, start=list(b0=1,b1=1), control=nls.control(maxiter=500), data = df)
  pow.q <- nls(p ~ exp(b0+b1*log(h)+b2*log(h)^2), start=list(b0=1,b1=1,b2=1), control=nls.control(maxiter=1000), data=df)
  epm2 <- nls(p ~ a*h^(b-(d/h)), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  p1 <- nls(p ~ a*h^(b)*exp(-d*h), start=list(a=1,b=1,d=0), control=nls.control(maxiter=1000),  data = df)
  
  results.df <- data.frame(AIC(lin,pow,pow.q,epm2,p1))
  
  results.df$min <- as.numeric(results.df$AIC == min(results.df$AIC))
  if(pred==TRUE) {
    
    pred.1 <- predict(lin,list(h=ext))
    pred.2 <- predict(pow,list(h=ext))
    pred.3 <- predict(pow.q,list(h=ext))
    pred.5 <- predict(epm2,list(h=ext))
    pred.6 <- predict(p1,list(h=ext))
    
    results.df$pred <- c(pred.1,
                         pred.2,
                         pred.3,
                         pred.5,
                         pred.6)
    
    return(results.df)
    
  } else {
    
    return(results.df)
    
  }
  
}

power.fit(df.disp)
power.fit.noepm(df.myco)
power.fit.noepm(df.poll)
power.fit(df.helm)
power.fit(df.dna)
power.fit(df.rna)

power.fit(df.dna,pred=TRUE,ext=5291)
power.fit(df.rna,pred=TRUE,ext=5291)

