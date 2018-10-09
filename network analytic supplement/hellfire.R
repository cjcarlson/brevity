library(codependent)
library(bipartite)

data("rob1929")
sch <- as.matrix(table(rob1929))
sch <- sch[,-1]
sch <- sch[-nrow(sch),]
null <- null.distr(N=1, sch, distr="negbin")[[1]]

colnames(null) <- c(1:ncol(null))
rownames(null) <- c(1:nrow(null))
df <- data.frame(null)
df <- df[df$Freq>0,][,1:2]
codependent::binera(df, 10, plots=TRUE)
codependent::binera(rob1929, 10, plots=TRUE)

g <- genweb(N1 = nrow(sch), N2 = ncol(sch), dens = 1)
df2 <- data.frame(as.table(g))
df2 <- df2[df2$Freq>0,][,1:2]
codependent::binera(df2, 10, plots=TRUE)
