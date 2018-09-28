library(ggplot2)
library(cowplot)

virus <- read.csv("~/Github/brevity/olival nature 2017/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')
virus <- virus[!(virus$Host %in% c('Homo_sapiens')),]

copredict(5291,virus,200,1)

length(unique(virus$Host))
length(unique(virus$Parasite))

mean(data.frame(table(virus$Host))$Freq)

14.93607*1431
14.93607*1435
14.93607*1439

virusmeta <- read.csv("~/Github/brevity/olival nature 2017/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)

rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

copredict(5291,rnavirus,200,1)
copredict(5291,dnavirus,200,1)

14.93607*892
14.93607*894
14.93607*897


14.93607*1612.7
14.93607*1591.6
14.93607*1633.8


table(virusmeta$vDNAoRNA,virusmeta$IsZoonotic)
24087*29/(29+176)
13353*159/(159+222)




#### prop zoonotic
virus$zoonotic=0
for (i in 1:nrow(virus)) {
virus$zoonotic[i] <- virusmeta[virusmeta$vVirusNameCorrected==virus$vVirusNameCorrected[i],]$IsZoonotic
}

hostmeta <- read.csv("~/Github/brevity/olival nature 2017/hosts.csv")
virus$hostgroup = ''

for (i in 1:nrow(virus)) {
  virus$hostgroup[i] <- as.character(hostmeta[hostmeta$hHostNameFinal==virus$hHostNameFinal[i],]$hOrder)
}

virus <- virus[,c(1,7,8)]
virus <- unique(virus)
df <- table(virus$hostgroup,virus$zoonotic)


#####


copredict.ci(5291,virus,200,1)

virusmeta <- read.csv("~/Github/brevity/olival nature 2017/viruses.csv")

dnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='DNA',]$vVirusNameCorrected)
rnalist <- unique(virusmeta[virusmeta$vDNAoRNA=='RNA',]$vVirusNameCorrected)

rnavirus <- virus[virus$Parasite %in% rnalist,]
dnavirus <- virus[virus$Parasite %in% dnalist,]

copredict.ci(5291,rnavirus,200,1)

# [1] "True number of species in entire network is 345"
# [1] "Estimated number of species in entire network is 393.590434663421"
# [1] "The lower 95% CI is 389.705103384448"
# [1] "The upper 95% CI is 397.514502410087"
# [1] "Extrapolated estimated number of species is 1156.54367524801"

14.92537*1157

17269*159/(159+222)

# [1] "The lower 95% CI is 1142.10898172487"

14.92537*1142 

17045*159/(159+222)
# [1] "The upper 95% CI is 1171.16080353039"

14.92537*1171

17478*159/(159+222)

copredict.ci(5291,dnavirus,200,1)


# [1] "Extrapolated estimated number of species is 2092.91700248616"

14.92537*2093

31239*159/(159+222)
# [1] "The lower 95% CI is 2035.41012410777"

14.92537*2035

30373*159/(159+222)

# [1] "The upper 95% CI is 2152.0486350218"

14.92537*2152

32119*159/(159+222)



df.virus <- curve.df(virus,100)
df.rna <- curve.df(rnavirus,100)
df.dna <- curve.df(dnavirus,100)

model.all <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.virus)
predicted.all <- data.frame(pred = predict(model.all), host = df.virus$n.host)

g1 <- ggplot(df.virus, aes(n.host, n.par)) + xlim(0,775) + ylim(0,650) + xlab('Mammals') + ylab('All viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('mediumpurple1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.all, aes(x=host, y=pred))
g1
model.all

model.rna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.rna)
predicted.rna <- data.frame(pred = predict(model.rna), host = df.rna$n.host)
g2 <- ggplot(df.rna, aes(n.host, n.par)) + xlim(0,740) + ylim(0,370) + xlab('Mammals') + ylab('RNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('mediumorchid1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.rna, aes(x=host, y=pred))
g2
model.rna

model.dna <- nls(n.par~b*n.host^z,start = list(b = 1, z = 0.5),data=df.dna)
predicted.dna <- data.frame(pred = predict(model.dna), host = df.dna$n.host)
g3 <- ggplot(df.dna, aes(n.host, n.par)) + xlim(0,180) + ylim(0,170) + xlab('Mammals') + ylab('DNA viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('slateblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.dna, aes(x=host, y=pred))
g3
model.dna

plot1 <- plot_grid(g1, g2, g3,
                      labels=c("A", "B", "C"), nrow=1, ncol = 3)

plot2 <- plot_grid(NULL, NULL,
                      labels=c("D", "E"), nrow=1, ncol = 2)

plot3 <- plot_grid(blank, blank, nrow=1, ncol = 2)

plot_grid(plot1,plot2,plot2,nrow=3)




#### supp figure


virus <- read.csv("~/Github/brevity/olival nature 2017/associations.csv")
virus <- virus[,c(2:1)]
names(virus) <- c('Host','Parasite')

df.virus <- curve.df(virus,100)
