
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
