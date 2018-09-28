
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
g1 <- ggplot(df.virus, aes(n.host, n.par)) + xlim(0,500) + ylim(0,1500) + xlab('Mammals') + ylab('All viruses') + 
  geom_point(shape = 16, size = 2.5, show.legend = FALSE, alpha = .05, color = c('steelblue1')) + theme_bw() +
  geom_line(color='black',lwd=1,data = predicted.all, aes(x=host, y=pred))



